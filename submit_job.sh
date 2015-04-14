#!/bin/bash
# The number of threads for each blastn job

num_procs=$1
shift
Nlev=$1
shift
Nsub=$1
shift
node=$1
shift
tasks_per_node=$1
shift

echo "num_procs = $num_procs"
echo "Nlev = $Nlev"
echo "Nsub = $Nsub"
echo "node = $node"
echo "tasks_per_node = $tasks_per_node"


if [ $num_procs -ne $(expr $node \* $tasks_per_node) ]
then
	echo "Error: wrong argument. $num_procs != $node * $tasks_per_node"
	exit
fi
wall_clock_limit=72:00:00
llscript_name=simulate
queuetime=`date +%G%m%d_%H%M%S`
outputdir=${basedir}/out-${queuetime}
email=katharina.dormanns@pg.canterbury.ac.nz
lltemplate=simulate.ll.tpl #LoadLeveler script template
 
 
# generates batcher_input_command1....X
 
# generates LL scripts and submit them
 
llfile=${llscript_name}_${num_procs}_${Nlev}_${Nsub}.ll
#replace variables in the template file with actual value
cat ${lltemplate} |sed -e  "s/\${num_procs}/$num_procs/" -e "s/\${Nlev}/$Nlev/"   -e "s/\${Nsub}/$Nsub/"   -e "s/\${tasks_per_node}/$tasks_per_node/"  -e "s/\${node}/$node/"   -e "s/\${wall_clock_limit}/$wall_clock_limit/" -e "s/\${email}/$email/" > ${llfile}
#        llsubmit ${llfile}
echo "LL script: ${llfile} has been generated"
cat ${llfile}
llsubmit ${llfile}
