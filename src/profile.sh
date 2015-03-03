#!/bin/bash
./run_N1000_gene645
cd ../CPU/profile/src/
./run_N1000 > profile.out
# extract 2 column of data
sed -n '13,24p' profile.out > ../out/N1000/profile.data
cd ../out/N1000
# plot using gnuplot and output to a png
gnuplot <<EOF
set terminal png
set output "profile.png"
plot 'actual_profile' using 1:2 w points, 'profile.data' using 1:2 w points
EOF
# this should popup a window with the plot
display profile.png


