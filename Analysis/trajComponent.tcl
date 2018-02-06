

# Usage: vmdt -e trajComponent.tcl -args startFrame file.parm file.dcd
# Measure COOH dihed component over all frames of trajectory.

# ====================================================================================


# open data output file for making plots
set outFile [open "valineAngle_first3000.dat" w]
puts $outFile "# frame\tvaline_ang "

set startFrame [lindex $argv 0] 
set parm [lindex $argv 1]
set dcd [lindex $argv 2]

mol new $parm type {parm7} first 0 last -1 step 1 waitfor -1
mol addfile $dcd type {dcd} first $startFrame last -1 step 1 waitfor -1 0


set nframes [molinfo 0 get numframes]
for {set i 0} {$i < $nframes} {incr i} {
    set dang [measure dihed "0 4 6 8" frame $i]
    puts $outFile "$i\t$dang"
}

close $outFile
exit
