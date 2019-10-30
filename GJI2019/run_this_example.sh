#!/bin/bash

source parameter
currentdir=`pwd`

echo "Configure and compile specfem2D ..."
cd $specfem_path
#make clean
if [ $NPROC_SPECFEM == 1 ]; then
    ./configure FC=$compiler 
else
    ./configure FC=$compiler --with-mpi
fi
#make all

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
read -rsp $'Press any key to run this example ...\n' -n1 key
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

cd $currentdir
if [ -z "$submit_dir" ]; then
    export submit_dir="./submit_job"
fi

rm -rf $submit_dir
mkdir $submit_dir

cp -r $specfem_path/bin $submit_dir/
cp -r submit.sh $submit_dir/
cp -r DATA $submit_dir/
cp -r parameter $submit_dir/

echo "ready to submit job in the directory: $submit_dir ..."
echo "cd $submit_dir"
echo './submit.sh'
cd $currentdir

