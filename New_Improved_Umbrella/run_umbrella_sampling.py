## This is a very quick skeleton driver script for umbrella sampling, created by David Mobley
## It is currently NOT FUNCTIONAL but attempts to outline what is probably a simpler
## way to do the calculations than how they are done in `restrOpenMM_GPU.py`


import numpy as np
import os
import parmed
from parmed.openmm import topsystem
from parmed import unit as u

import simtk.openmm as mm
from simtk.openmm import Platform
from simtk.openmm import app
from simtk.unit import *
# There need to be more imports here relating to OpenMM


## Load umbrella centers and force constants
umbrella_data = np.loadtxt('centers.dat')
centers = umbrella_data[:,0]*np.pi/180. #Convert to radians
force_constants = umbrella_data[:,1]

n_centers = len(centers)

# Let's assume we already have input input files
prmFile = 'watVA.prmtop'
inpFile = 'watVA.inpcrd'

if not os.path.isfile(prmFile):
    raise IOError('Missing input file %s' %prmFile)
if not os.path.isfile(inpFile):
    raise IOError('Missing input file %s' %inpFile)


def runOpenMM(parm, system, rad, K, Indices, solvate, out_dcd, out_csv):
    """

    Minimize molecule with specified topology, system, and positions
       using OpenMM. Return the positions of the optimized moelcule.

    Parameters
    ----------
    parm:      ParmEd object for system (parameters/corodinates)
    system:    OpenMM system for this mol's Prmtop and Inpcrd files
    rad:       float of reference angle around which to restrain
    K:         float of force constant for harmonic restraint
    Indices:   list of integers of zero-based atom indices for restraint
    out_dcd:   filename for output dcd file for constant pressure equilibration/production
    out_csv:   filename for output csv file for constant pressure equilibration/production

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
    harmonic = mm.CustomTorsionForce("k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = %.5f" % pi);
    harmonic.addPerTorsionParameter("theta0");
    harmonic.addPerTorsionParameter("k");
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
    simulation = app.Simulation(parm.topology, system, integ1)
    simulation.context.setPositions(parm.positions)

    # perform minimization
    print('Minimizing...')
    simulation.minimizeEnergy()

    # NVT equilibration
    simulation.context.setVelocitiesToTemperature(300*u.kelvin)
    simulation.reporters.append(app.StateDataReporter('data01.csv', 1000, step=True, potentialEnergy=True, volume=True,temperature=True, separator='\t'))
    print('Equilibrating at NVT...')
    simulation.step(10000) # 10 ps

    if solvate==True:
        positionsNVT = simulation.context.getState(getPositions=True).getPositions()
        velocitiesNVT = simulation.context.getState(getVelocities=True).getVelocities()

        # NPT equilibration/production
        barostat = mm.MonteCarloBarostat(1.0*u.bar, 300.0*u.kelvin)
        system.addForce(barostat)
        ### bc barostat, need new simulation and associated properties
        integ2 = newIntegrator()
        simulation = app.Simulation(parm.topology, system, integ2)
        simulation.context.setPositions(positionsNVT)
        simulation.context.setVelocities(velocitiesNVT)
        simulation.reporters.append(app.DCDReporter(out_dcd, 1000))
        simulation.reporters.append(app.StateDataReporter(out_csv, 1000, step=True, potentialEnergy=True, volume=True,temperature=True, separator='\t'))
        print('Constant pressure equilibration/production...')
        simulation.step(3000000) # 10 ps


    else:
        print('Production run at NVT...')
        simulation.step(3000000) # 100 ps


    parm.positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    return parm.positions




# Which atoms to umbrella sample
atomlist=[0,4,6,8]

# Note that if you want to make this run faster, you can trivially parallelize by
# just making this script run one particular center and having a separate
# shell script run each umbrella
for i,umbr_num in enumerate(range(n_centers)):
    #Retrieve force constants, file names, etc.
    fc = force_constants[i]
    center = centers[i]
    prefix = 'umbrella%s' % umbr_num
    dcd_file = prefix+'.dcd'
    csv_file = prefix+'.csv'

    # Load inputs; I put this inside the loop just to be safe (in case anything gets changed)
    parm = parmed.load_file(prmFile, inpFile)

    # Create OpenMM system
    system = parm.createSystem(nonbondedMethod=app.PME,constraints=app.HBonds,rigidWater=True,removeCMMotion=True)

    # I would probably add a check here that the atoms actually have the names you expect them to have.

    # Run simulation
    runOpenMM(parm, system, center, fc, atomlist, True, dcd_file, csv_file)



# I THINK something roughly like that should be about all you need.
# IMPORTANT NOTES:
# - CHECK UNITS: I have not yet checked the units of the force constants to make sure they are OpenMM units or get converted.
# - I am certainly missing some OpenMM-related imports
# - I changed the arguments to your runOpenMM function considerably, stripping out redundant information
# - One thing I'm worried about in your umbrella sampling script is that it had two topologies
# that it was carting around -- one embedded in parm, and one created by the GAFF thing.
# therewas no checking these were equivalent. Likewise two sets of positions.


# AFTER THIS: I would just modify the umbrella-sampling1.py file to DIRECTLY load the dcd files with mdtraj and do the anlaysis on them
# At the same time I would re-check units, etc.
