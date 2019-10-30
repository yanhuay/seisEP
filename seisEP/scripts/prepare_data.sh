#!/bin/bash

# parameters
source parameter
data_tag='DATA_obs'
velocity_dir=$1
SAVE_FORWARD=false

# local id (from 0 to $ntasks-1)
if [ $system == 'slurm' ]; then
    iproc=$SLURM_PROCID  
elif [ $system == 'pbs' ]; then
    iproc=$PBS_VNODENUM
fi
#iproc=${OMPI_COMM_WORLD_RANK}

# allocate tasks over all sources
# ntasks in parallel and nsrc in total
# nsrc_per_task=$(( $NSRC / $ntasks ))
# take ceiling 
nsrc_per_task_ceiling=$(echo $(echo "$NSRC $ntasks" | awk '{ print $1/$2 }') | awk '{printf("%d\n",$0+=$0<0?0:0.999)}')
ntasks_ceiling=$(echo $(echo "$NSRC $ntasks" | awk '{print $1%$2}')) 
# take floor 
nsrc_per_task_floor=$(echo $(echo "$NSRC $ntasks" | awk '{ print int($1/$2) }'))

# allocate nsrc for each task
if [ $iproc -lt $ntasks_ceiling ]; then
    nsrc_this_task=$nsrc_per_task_ceiling
    isource_start=$(echo $(echo "$iproc $nsrc_per_task_ceiling" | awk '{ print $1*$2 }'))
else
    nsrc_this_task=$nsrc_per_task_floor
    isource_start=$(echo $(echo "$iproc $nsrc_per_task_floor $ntasks_ceiling $nsrc_per_task_ceiling" | awk '{ print ($1-$3)*$2+$3*$4 }'))
fi

if [ $iproc -eq  0 ]; then
    echo "allocate $NSRC sources over $ntasks tasks"
    echo "iproc 0-$(($ntasks_ceiling-1)): nsrc_per_task=$nsrc_per_task_ceiling"
    echo "iproc $ntasks_ceiling-$(($ntasks-1)): nsrc_per_task=$nsrc_per_task_floor"
fi
echo "iproc = $iproc, isource_start = $isource_start, nsrc_this_task = ${nsrc_this_task}, source: $(($isource_start+1))-$(($isource_start+$nsrc_this_task)) "

# source for this task
for ((isrc_this_task=1; isrc_this_task<=${nsrc_this_task}; isrc_this_task++));
do 
    let isource=`echo $isource_start, $isrc_this_task |awk '{print $1 + $2 }'` 
    echo "iproc = $iproc, isource = $isource"

    STARTTIME=$(date +%s)
    if  $ExistDATA && [ -d "$DATA_DIR" ]; then
        sh $SCRIPTS_DIR/copy_data.sh $isource $data_tag $data_list $WORKING_DIR $DATA_DIR 2>./job_info/error_copy_data
    else
        sh $SCRIPTS_DIR/Forward_${solver}.sh $isource $NPROC_SPECFEM $data_tag $data_list \
            $velocity_dir $SAVE_FORWARD $WORKING_DIR $SUBMIT_RESULT $job 2>./job_info/error_Forward_simulation
        if [ $isource -eq 1 ] && [ -d "$m_target" ]; then
            cp -r $velocity_dir $SUBMIT_RESULT/m_target
        fi
    fi
 if [ $isource -eq 1 ] ; then
     ENDTIME=$(date +%s)
     Ttaken=$(($ENDTIME - $STARTTIME))
     echo "Data preparation took $Ttaken seconds for one source"
 fi
done

