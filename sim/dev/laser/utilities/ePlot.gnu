#set datafile separator "\t"

plot "meanE.dat" using 1:2 with points title 'Mean Energy',\
    "meanE.dat" using 1:3 with points title 'Temperature'
set logscale xy
pause 1
reread
