#! /software/gnuplot/5.2.0b2/bin/gnuplot

set pm3d map
splot '2DES_T0.txt' u 1:2:3
set terminal png crop enhanced
set output '2DES_T0_real_g08.png'

set palette defined (-6 "blue", 0 "white", 3 "red")
set size square
set contour base
set cntrparam bspline
set cntrparam level incremental -6, 0.9, 3
set style line  1 lc rgb "black" lw 1
set style line  2 lc rgb "black" lw 1
set style line  3 lc rgb "black" lw 1
set style line  4 lc rgb "black" lw 1
set style line  5 lc rgb "black" lw 1
set style line  6 lc rgb "black" lw 1
set style line  7 lc rgb "black" lw 1
set style line  8 lc rgb "black" lw 1
set style line  9 lc rgb "black" lw 1
set style line  10 lc rgb "black" lw 1
set style line  11 lc rgb "black" lw 1
set style line  12 lc rgb "black" lw 1
set style line  13 lc rgb "black" lw 1
set style line  14 lc rgb "black" lw 1
set style line  15 lc rgb "black" lw 1
set style line  16 lc rgb "black" lw 1
set style line  17 lc rgb "black" lw 1
set style line  18 lc rgb "black" lw 1
set style line  19 lc rgb "black" lw 1
set style line  20 lc rgb "black" lw 1
set style line  21 lc rgb "black" lw 1
set style line  22 lc rgb "black" lw 1
set style line  23 lc rgb "black" lw 1
set style line  24 lc rgb "black" lw 1
set style line  25 lc rgb "black" lw 1
set style line  26 lc rgb "black" lw 1
set style line  27 lc rgb "black" lw 1
set style line  28 lc rgb "black" lw 1
set style line  29 lc rgb "black" lw 1
set style line  30 lc rgb "black" lw 1
set style line  31 lc rgb "black" lw 1
set style line  32 lc rgb "black" lw 1
set style line  33 lc rgb "black" lw 1
set style line  34 lc rgb "black" lw 1
set style line  35 lc rgb "black" lw 1
set style line  36 lc rgb "black" lw 1
set style line  37 lc rgb "black" lw 1
set style line  38 lc rgb "black" lw 1
set style line  39 lc rgb "black" lw 1
set style line  40 lc rgb "black" lw 1
set style increment user
set xrange [1.8 to 2.55]
set yrange [1.8 to 2.55]
set cbrange[-6 to 3]
set xlabel 'excitation frequency (eV)'
set ylabel 'detection frequency (eV)'
#set isosample 250, 250
set pm3d interpolate 0,0
set arrow from 1.8,1.8 to 2.55,2.55 nohead front lt 0 lw 2.5 lc rgb "black"
unset key
replot
