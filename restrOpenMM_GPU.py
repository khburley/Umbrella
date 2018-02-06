#!/usr/bin/env python


import os, sys, glob
import numpy as np
import math
import re 

import openeye.oechem as oechem
import openeye.oeomega as oeomega

import parmed
from parmed.openmm import topsystem
from parmed import unit as u

import simtk.openmm as mm
from simtk.openmm import Platform
from simtk.openmm import app
from simtk.unit import *


def writeUpdatedMol(Mol, fname):
    """

    Parameters
    ----------
    Mol: an OEChem molecule
    fname: str - name of the output mol2 file

    Returns
    -------
    True if the function completed successfully

    """

    ofs = oechem.oemolostream()
    if not ofs.open(fname):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % fname)
    oechem.OEWriteConstMolecule(ofs, Mol)
    ofs.close()
    return True

def prepGAFFx(parm, solvate):
    """

    Creates topology, system, and coordinates for AMBER
       Prmtop and Inpcrd input files. This function should work
       with either GAFF or GAFF2 files.

    Parameters
    ----------
    parm: parmed.structure.Structure instance from an OpenMM topology

    Returns
    -------
    topology:  OpenMM topology for this mol's Prmtop and Inpcrd files
    system:    OpenMM system for this mol's Prmtop and Inpcrd files
    positions: OpenMM positions for this mol's Prmtop and Inpcrd files

    """

    topology = parm.topology
    positions = parm.positions
    if solvate==True:
        system = parm.createSystem(nonbondedMethod=app.PME,constraints=app.HBonds,rigidWater=False,removeCMMotion=True)
    else:
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,constraints=app.HBonds,rigidWater=False,removeCMMotion=True)

    return topology, system, positions

def runOpenMM(parm, topology, system, positions, rad, K, Indices, solvate):
    """

    Minimize molecule with specified topology, system, and positions
       using OpenMM. Return the positions of the optimized moelcule.

    Parameters
    ----------
    topology:  OpenMM topology for this mol's Prmtop and Inpcrd files
    system:    OpenMM system for this mol's Prmtop and Inpcrd files
    positions: OpenMM positions for this mol's Prmtop and Inpcrd files
    rad:       float of reference angle around which to restrain
    K:         float of force constant for harmonic restraint
    Indices:   list of integers of zero-based atom indices for restraint

    Returns
    -------
    topology.positions: OpenMM topology positions for minimized mol

    """

    def newIntegrator():
        integrator = mm.LangevinIntegrator(
                300.0 * u.kelvin,
                10.0 / u.picosecond,
                1.0 * u.femtosecond)
        return integrator


    # harmonically restrain dihedral angle
    # see units, http://docs.openmm.org/6.3.0/userguide/theory.html
    pi = np.pi
    harmonic = mm.CustomTorsionForce("k*min(dtheta, 2*pi-dtheta)^2 + pi; dtheta = abs(theta-theta0);");
    harmonic.addPerTorsionParameter("theta0");
    harmonic.addPerTorsionParameter("k");
    harmonic.addGlobalParameter("pi",pi);
    #harmonic.setParameter("pi",pi);
    system.addForce(harmonic)
    harmonic.addTorsion(Indices[0], Indices[1], Indices[2], Indices[3], (rad, K))

    #Restrain backbone atoms
    force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", 5.0*kilocalories_per_mole/angstroms**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for i, atom_crd in enumerate(parm.positions):
            if parm.atoms[i].name in ('CA', 'C', 'N'):
                        force.addParticle(i, atom_crd.value_in_unit(u.nanometers))
    system.addForce(force)


    # build simulaion
    #platform = mm.Platform.getPlatformByName('CPU')
    platform = mm.Platform.getPlatformByName('CUDA')
    integ1 = newIntegrator()
    simulation = app.Simulation(topology, system, integ1)
    simulation.context.setPositions(positions)

    # perform minimization
    print('Minimizing...')
    simulation.minimizeEnergy()

    # NVT equilibration
    simulation.context.setVelocitiesToTemperature(300*u.kelvin)
    simulation.reporters.append(app.DCDReporter('nvt01.dcd', 1000))   # write every 1000 steps
    simulation.reporters.append(app.StateDataReporter('data01.csv', 1000, step=True, potentialEnergy=True, volume=True,temperature=True, separator='\t'))
    print('Equilibrating at NVT...')
    simulation.step(10000) # 10 ps

    if solvate==True:
        positionsNVT = simulation.context.getState(getPositions=True).getPositions()
        velocitiesNVT = simulation.context.getState(getVelocities=True).getVelocities()

        # NPT equilibration
        barostat = mm.MonteCarloBarostat(1.0*u.bar, 300.0*u.kelvin)
        system.addForce(barostat)
        ### bc barostat, need new simulation and associated properties
        integ2 = newIntegrator()
        simulation = app.Simulation(topology, system, integ2) 
        simulation.context.setPositions(positionsNVT)
        simulation.context.setVelocities(velocitiesNVT)
        simulation.reporters.append(app.DCDReporter('npt01.dcd', 1000))
        simulation.reporters.append(app.StateDataReporter('data02.csv', 1000, step=True, potentialEnergy=True, volume=True,temperature=True, separator='\t'))
        print('Equilibrating at NPT...')
        simulation.step(10000) # 10 ps
    
        # NPT production
        print('Production run at NPT...')
        simulation.step(3000000) # 100 ps

    else:
        print('Production run at NVT...')
        simulation.step(3000000) # 100 ps


    topology.positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    return topology.positions

# --------------------------- Main Function ------------------------- #


def load_and_minimize(infiles,dogaff,dogaff2,gaffdir,atomlist,solvate):
    molfile = glob.glob(os.path.join(infiles, '0.mol2'))
    print('found mol2 file')
    print(molfile)
    rad = 0
    f = molfile[0]
        ### Load molecule from .mol2 file.
    print(f)
    ifs = oechem.oemolistream()
    if not ifs.open(f):
         oechem.OEThrow.Warning("Unable to open %s for reading" % f)
    try:
        mol = next(ifs.GetOEMols())
    except StopIteration:
            print('No mol loaded for %s (StopIteration exception)' % mol.GetTitle())
    ifs.close()

    # identify the input files 
    prmFile = 'vacDivaline.prmtop'
    inpFile = 'vacDivaline.inpcrd'

    # set umbrella windows/angle restraints
    for i in range(-180,180,10):
        ### Set output file name and make a copy of molecule for opt.
        startdir = os.getcwd()
        print('This is the starting directory')
        print(startdir)
        print('This is the mol')
        print(mol)
        if i < 0:
            #convert i to generate 0 - 365 intead of -180 to 180 folders
            adj_i = i + 360
            moldir = os.path.join(startdir, opt.fftype.upper()+'/',str(adj_i))
        else:
            moldir = os.path.join(startdir, opt.fftype.upper()+'/',str(i))
        print('This is the mol directory')
        print(moldir)
        if not os.path.exists(moldir):
            os.makedirs(moldir)
            print('mol directory created')
        tmpmol = oechem.OEGraphMol( mol)
        if dogaff: print('==> Starting on GAFF optimization for %s' % tmpmol.GetTitle())
        if dogaff2: print('==> Starting on GAFF2 optimization for %s' % tmpmol.GetTitle())

        ### Get relevant files and file names.
        ffield = opt.fftype.upper()

        ### Check: gaff dependencies present, gaff files have content
        if not (os.path.exists(prmFile) and os.path.exists(inpFile)):
            print('%s.inpcrd or %s.prmtop files do not exist' % (mol.GetTitle(), mol.GetTitle()))
            continue
        if not (os.path.getsize(prmFile) and os.path.getsize(inpFile) > 30):
            print('%s.inpcrd or %s.prmtop files are empty' % (mol.GetTitle(), mol.GetTitle()))
            continue
        parm = parmed.load_file(prmFile, inpFile)
        print('Parm object generated')
        ### Generate topology, system, and position.
        top, syst, pos = prepGAFFx(parm, solvate)
        print('print system prepped (gaffx)')
        
        # restrain according to angle specified in range
        rad_real = math.radians(i)
        print(rad_real)
        os.chdir(moldir)
        rad = math.radians(i)
        # Set the force constants for each angle
        # there's probably a much more efficent way to do this (ie - dictionary)
        smallk = []
        try:
            if i in smallk:
                parm.positions = runOpenMM(parm, top, syst, pos, rad, 50., atomlist, solvate)
                print('using K = 50')
                print('Restraining around degree: %f' % (math.degrees(rad)))
            else:
                print('Trying to execute runOpenMM line')
                parm.positions = runOpenMM(parm, top, syst, pos, rad, 200., atomlist, solvate)
                print('using K = 200')
                print('Restraining around degree: %f' % (math.degrees(rad)))

        except Exception as e:
            if str(e) == 'Particle coordinate is nan':
                print("Simulation failed with ValueError for conformation %s" % f)
                os.chdir(startdir)
                continue
            else:
                print("Simulation failed for conformation %s" % f)
                print(e)
                os.chdir(startdir)
                continue
        os.chdir(startdir)
        


        ### Format OpenMM positions to a form readable to oemol
        concat_coords = []
        for atomi in parm.positions:
            concat_coords += [float(i) for i in atomi._value]
        tmpmol.SetCoords(oechem.OEFloatArray(concat_coords))

        #for debugging to make sure angle was restrained during minimization
        dihedcrds = []
        for ind in atomlist: 
            xyz = oechem.OEFloatArray(3)
            for atom in tmpmol.GetAtoms():
                #print atom.GetIdx(), atom.GetName()
                if atom.GetIdx() == ind: 
                    tmpmol.GetCoords( atom, xyz)
                    dihedcrds.append([xyz[0],xyz[1],xyz[2]])
        ang = parmed.geometry.dihedral(dihedcrds[0],dihedcrds[1],dihedcrds[2],dihedcrds[3])
        rad2 = math.radians(ang)
        print(ang,rad2)

        #fulln = moldir+"/%s.mol2" % (mol.GetTitle())
        #writeUpdatedMol(tmpmol, fulln)

        #rad = rad + math.radians(15)

# ------------------------- Parse Inputs ----------------------- #



if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage = " ... update this from description at top of python script ... ")

    parser.add_option('-i', '--inmols',
            help = "REQUIRED! Path to directory containing all .mol2 files\
 to be minimized.",
            type = "string",
            dest = 'inmols')

    parser.add_option('-f','--fftype',
            help = "REQUIRED! Force field type. Supported options:\
 'mmff', 'gaff', 'gaff2', 'smirff'.",
            type = "string",
            dest = 'fftype')

    parser.add_option('-g', '--gaffdir',
            help = "Directory containing GAFFx parameter/topology files\
 (*.prmtop) and coordinates files (*.inpcrd). Required if and only if --fftype\
 is 'gaff' or 'gaff2'.",
            type = "string",
            default = None,
            dest = 'gaffdir')

    parser.add_option('--solvate',
            help = "Do input systems have explicit water molecules?",
            action="store_true",
            default = False,
            dest = 'solvate')

    ### Check required fields.
    (opt, args) = parser.parse_args()
    if opt.fftype == None:
        parser.error("ERROR: No force field was specified.")
    if opt.inmols == None:
        parser.error("ERROR: No directory for input mol2 files provided.")
    if opt.fftype == 'gaff' or opt.fftype == 'gaff2':
        if opt.gaffdir == None:
            parser.error("ERROR: No GAFF* files' directory provided for GAFF* force field.")
    dogaff = opt.fftype == 'gaff'      # is True if opt.type == 'gaff'
    dogaff2 = opt.fftype == 'gaff2'    # is True if opt.type == 'gaff2'

    atomlist=[0,4,6,8]
    print('starting load/min')
    load_and_minimize(opt.inmols, dogaff, dogaff2, opt.gaffdir, atomlist, opt.solvate==False)

