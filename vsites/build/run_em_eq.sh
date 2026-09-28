#!/bin/bash
set -e
# Build-time record: this ran inside the original working tree.  Set G to your gmx.
G=${GMX:-gmx}
MD="$G mdrun -ntmpi 1 -ntomp 4 -nb cpu -pin off"
$G grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
$MD -deffnm em -resetstep 2000
$G grompp -f eq.mdp -c em.gro -p topol.top -o eq.tpr
$MD -deffnm eq -resetstep 2000
echo "BUILD_EQ_DONE"
