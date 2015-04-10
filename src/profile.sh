#!/bin/bash
DEFAULT_NUM_PLOTS=3
NUM_PLOTS_REQUESTED=${1:-$DEFAULT_NUM_PLOTS}

./run_N1000_gene645
cd ../CPU/profile/src/
make
./run_N1000
NUM_SUBNETS="$(cat n_subnets.out)"
NUM_PLOTS=$(($NUM_SUBNETS<$NUM_PLOTS_REQUESTED?$NUM_SUBNETS:$NUM_PLOTS_REQUESTED))
# extract 2 column of data
#sed -n '13,24p' profile.out > ../out/N1000/profile.data
cd ../out/N1000
# plot using gnuplot and output to a png
gnuplot <<EOF
set terminal png
set output "profile.png"
plot 'actual_profile' using 1:2 w points ps 3, for [i=1 : $NUM_PLOTS] 'run'.i.'.out' using 1:2 w lp title 'run '.i
EOF
# this should popup a window with the plot
display profile.png

#plot 'actual_profile' using 1:2 w points, 'profile.data' using 1:2 w lp


