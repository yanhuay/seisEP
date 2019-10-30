#!/bin/bash

echo 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo
echo " source parameter file ..." 
source parameter

echo 
echo " create new job_info file ..."
rm -rf job_info
mkdir job_info

echo 
echo " create result file ..."
mkdir -p RESULTS

echo
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo
echo

echo " select workflow ..."
workflow_DIR="$package_path/workflow"

var1="modeling"
var2="kernel"
var3="inversion"
var4="fwi"
var5="misfit"
if [ "${job,,}" == "${var1,,}"  ]
then
    echo " ########################################################"
    echo " Forward modeling ..." 
    echo " ########################################################"
    cp $workflow_DIR/Modeling.sh $Job_title.sh

elif [ "${job,,}" == "${var2,,}"  ]
then
    echo " ########################################################"
    echo " Kernel Construction ..." 
    echo " ########################################################"
    cp $workflow_DIR/Kernel.sh $Job_title.sh  
elif [ "${job,,}" == "${var3,,}"  ] || [ "${job,,}" == "${var4,,}"  ]
then
    echo " ########################################################"
    echo " Adjoint Inversion ..." 
    echo " ########################################################"
    cp $workflow_DIR/AdjointInversion.sh $Job_title.sh
elif [ "${job,,}" == "${var5,,}"  ]
then
    echo " ########################################################"
    echo " Misfit Evaluation ..."
    echo " ########################################################"
    cp $workflow_DIR/misfit.sh $Job_title.sh
else
    echo "Wrong job: $job"
    exit
fi
echo 
#read -p "Is the correct workflow selected (y/n)?" yn
yn='y'
if [ $yn == 'n' ]; then
    exit
fi
# check compiler
echo 
echo "check the compiler for mpif90"
echo "----------------------------"
mpif90 -v
echo "----------------------------"
#read -p "the compiler is $compiler (y/n)?" yn
if [ $yn == 'n' ]; then
    exit
fi

## not modeling 
if [ "${job,,}" != "${var1,,}"  ]; then
    mkdir -p bin
    echo
    echo " renew parameter file ..."
    cp $package_path/SRC/seismo_parameters.f90 ./bin/
    cp $package_path/lib/src/constants.f90 ./bin/
    cp $package_path/scripts/renew_parameter.sh ./
    ./renew_parameter.sh

    echo 
    #read -p "Do you wish to compile lib codes (y/n)?" yn
    if [ $yn == 'y' ]; then 
        rm -rf *.mod make_*
        cp $package_path/lib/make_lib ./make_lib
        FILE="make_lib"
        sed -e "s#^FC=.*#FC=${compiler}#g"  $FILE > temp;  mv temp $FILE
        sed -e "s#^SRC_DIR=.*#SRC_DIR=$package_path/lib/src#g"  $FILE > temp;  mv temp $FILE
        sed -e "s#^MOD_DIR=.*#MOD_DIR=./bin#g"  $FILE > temp;  mv temp $FILE
        sed -e "s#^LIB_preprocess=.*#LIB_preprocess=./bin/seismo.a#g"  $FILE > temp;  mv temp $FILE
        make -f make_lib clean
        make -f make_lib
    else
        echo 
        #read -p "Do you wish to compile source codes (y/n)?" yn
    fi
    if [ $yn == 'y' ]; then
        cp $package_path/make/make_file ./make_file
        FILE="make_file"
        sed -e "s#^FC=.*#FC=${compiler}#g"  $FILE > temp;  mv temp $FILE
        sed -e "s#^SRC_DIR=.*#SRC_DIR=$package_path/SRC#g"  $FILE > temp;  mv temp $FILE
        sed -e "s#^LIB_seismo=.*#LIB_seismo=./bin/seismo.a#g"  $FILE > temp;  mv temp $FILE
        make -f make_file clean
        make -f make_file
    fi 
fi ## not modeling 

nproc=$NPROC_SPECFEM
ntaskspernode=$(echo "$max_nproc_per_node $nproc" | awk '{ print $1/$2 }')
nodes=$(echo $(echo "$ntasks $nproc $max_nproc_per_node" | awk '{ print $1*$2/$3 }') | awk '{printf("%d\n",$0+=$0<0?0:0.999)}')
echo
echo "Request $nodes nodes, $ntasks tasks, $ntaskspernode tasks per node, $nproc cpus per task "
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#read -p "ready to submit job (y/n)?" yn
if [ $yn == 'y' ]; then
    if [ $system == 'slurm' ]; then
        echo "slurm system ..."
        echo "sbatch -p $queue -N $nodes -n $ntasks --cpus-per-task=$nproc -t $WallTime -e job_info/error -o job_info/output $Job_title.sh"
        sbatch -p $queue -N $nodes -n $ntasks --cpus-per-task=$nproc -t $WallTime -e job_info/error -o job_info/output $Job_title.sh

    elif [ $system == 'pbs' ]; then
        # Benchmarked by Chao Zhang (dkzhangchao@gmail.com) 
        echo "pbs system ..."
        echo "qsub -q $queue -l nodes=$nodes:ppn=$max_nproc_per_node -l --walltime=$WallTime -e job_info/error -o job_info/output  $Job_title.sh"
        qsub -q $queue -l nodes=$nodes:ppn=$max_nproc_per_node -l --walltime=$WallTime -e job_info/error -o job_info/output  $Job_title.sh
    fi
fi
echo
