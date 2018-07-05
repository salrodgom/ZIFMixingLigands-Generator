!e_min=$( tail -n1 log.txt | awk '{print $2}' ) ; echo "e_min=${e_min}" > c
load 'c'
set logscale xy
plot 'log.txt' u 1:(abs($2/e_min - 1 )) w lp
