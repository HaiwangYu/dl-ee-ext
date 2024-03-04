#!/usr/bin/env bash

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v08_00_00_59 -q e17:prof
unsetup mrb;
setup mrb -o
#source /uboone/data/users/wgu/larsoft_j11/localProducts_larsoft_v08_05_00_17_e17_prof/setup
#mrbslp

# DL stuff ref: uboonecode v08_00_00_29e_dl
setup numpy v1_14_3 -q e17:openblas:p2714b:prof
setup libtorch v1_0_1 -q e17:prof
setup SparseConvNet 8422a6f -q e17:prof

#WCP_FQ_DIR=/uboone/app/users/$USER/products/wire-cell/
#path-prepend $WCP_FQ_DIR/lib/ LD_LIBRARY_PATH
#path-prepend $WCP_FQ_DIR/lib64/ LD_LIBRARY_PATH
#path-prepend $WCP_FQ_DIR/bin/ PATH
#path-prepend $WCP_FQ_DIR/python PYTHONPATH

find-fhicl(){
    fhicl_file=$1
    for path in `echo $FHICL_FILE_PATH  | sed -e 's/:/\n/g'`;do find $path -name "$fhicl_file"  2>/dev/null;done
}

locality () {
    f=$1;
    dir=`samweb locate-file $f | egrep 'enstore:|dcache:' | cut -d: -f2 | cut -d\( -f1`;
    pushd $dir;
    cat ".(get)($f)(locality)";
    popd;
}

