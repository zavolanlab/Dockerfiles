#!/bin/bash

export CURRENT=`pwd`
export LIB=$CURRENT/libs

#installing gsl
export GSLLIB=$CURRENT/libs/gsl-1.15
#cd $LIB
#tar xzvf gsl-1.15.tar.gz
#mkdir gsl
#cd $GSLLIB
#./configure --prefix=$LIB/gsl
#make
#make install
#rm -rf $GSLLIB

#installing hdf
export HDFLIB=$CURRENT/libs/hdf5-1.8.8
#cd $LIB
#tar xzvf hdf5-1.8.8.tar.gz
#mkdir hdf
#cd $HDFLIB
#./configure --prefix=$LIB/hdf
#make
#make install
#rm -rf $HDFLIB

#installing o2scl
export SCLLIB=$CURRENT/libs/o2scl-0.906
#cd $LIB
#tar xzvf o2scl-0.906.tar.gz
#mkdir o2scl
#cd $SCLLIB
#export LDFLAGS="-L$CURRENT/libs/gsl/lib -L$CURRENT/libs/hdf/lib"
#export CPPFLAGS="-O3 -I$CURRENT/libs/gsl/include -I$CURRENT/libs/hdf/include"
#./configure --prefix=$LIB/o2scl
#make 
#make install
#rm -rf $SCLLIB

#installing mirza
#export CPPFLAGS="-O3 -I$CURRENT/libs/gsl/include -I$CURRENT/libs/hdf/include -I$CURRENT/libs/o2scl/include"
cd $CURRENT

make all

rm $CURRENT/*.o
mv MIRZA ./bin/
rm $CURRENT/src/*.o
./launch_MIRZA

