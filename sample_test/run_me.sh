#!/bin/bash


## Please first dowload the demo data from https://sourceforge.net/projects/tiglon/files/DemoData/

# !!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!

export LD_LIBRARY_PATH=/your/boost/dir/lib:$LD_LIBRARY_PATH

#export LD_LIBRARY_PATH=/yuting/local/boost/lib:$LD_LIBRARY_PATH

../Tiglon -B bamlist -s first -o tiglon_oudir 
