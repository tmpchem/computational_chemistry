#!/usr/bin/gnuplot
set terminal postscript eps enhanced 'Helvetica' 10.5 solid color
set output "/home/tmpchem/coding/compchem/graph/ch2o_12.eps"
set size 0.67,0.67

unset clip points

# Point Type (pt)
# 1 cross : 2 x : 3 asterick
# 4 open box : 5 filled box
# 6 open circle : 7 filled circle
# 8 open triangle : 9 filled triangle
# 10 open upside down triangle
# 11 filled upside down triangle
# 12 open diamond : 13 filled diamond

set style line 1  lt 1 linecolor rgb "#000000" lw 1.2 pt 9  ps 1.0
set style line 2  lt 3 linecolor rgb "#05A000" lw 0.5 pt 9  ps 1.0
set style line 3  lt 3 linecolor rgb "#FF0000" lw 0.5 pt 13 ps 1.0
set style line 4  lt 1 linecolor rgb "#0000FF" lw 0.7 pt 13 ps 1.0
set style line 5  lt 1 linecolor rgb "#FF6800" lw 0.7 pt 13 ps 1.0
set style line 6  lt 1 linecolor rgb "#800080" lw 0.7 pt 13 ps 1.0
set style line 7  lt 1 linecolor rgb "#FF00B6" lw 0.5 pt 13 ps 1.0
set style line 8  lt 3 linecolor rgb "#DDDD00" lw 0.5 pt 9  ps 1.0
set style line 9  lt 3 linecolor rgb "#FF8888" lw 0.9 pt 9  ps 1.0
set style line 10 lt 3 linecolor rgb "#88FF88" lw 0.9 pt 13 ps 1.0
set style line 11 lt 3 linecolor rgb "#8888FF" lw 0.9 pt 13 ps 1.0
set style line 12 lt 3 linecolor rgb "#FF8EF9" lw 0.9 pt 13 ps 1.0

set title "Molecular Dynamics Simulation of 12 Formaldehyde Molecules"
#set title offset -9.5, -2.6
set ytics 100
set mytics 4
set yrange[-90:620]
set xtics 10
set mxtics 2
set xzeroaxis
set grid x
set xrange [0.000:100.000]
set format x "%.0f"
set format y "%.0f"
set xlabel "Time (ps)"
set ylabel "Energy Terms (kcal mol^{-1})"
set key below spacing 1.0

plot \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:13 notitle w l ls 12, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:12 notitle w l ls 11, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:11 notitle w l ls 10, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:10 notitle w l ls 9, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:9  notitle w l ls 8, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:8  notitle w l ls 7, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:7  notitle w l ls 6, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:6  notitle w l ls 5, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:5  notitle w l ls 4, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:4  notitle w l ls 3, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:3  notitle w l ls 2, \
  "/home/tmpchem/coding/compchem/data/ch2o_12.dat" using 1:2  notitle w l ls 1, \
  10000 title "Total"       w l ls 1  lw 3.0, \
  10000 title "Kinetic"     w l ls 2  lw 3.0, \
  10000 title "Potential"   w l ls 3  lw 3.0, \
  10000 title "Non-bonded"  w l ls 4  lw 3.0, \
  10000 title "Bonded"      w l ls 5  lw 3.0, \
  10000 title "Boundary"    w l ls 6  lw 3.0, \
  10000 title "Vdw"         w l ls 7  lw 3.0, \
  10000 title "Elst"        w l ls 8  lw 3.0, \
  10000 title "Bonds"       w l ls 9  lw 3.0, \
  10000 title "Angles"      w l ls 10 lw 3.0, \
  10000 title "Torsions"    w l ls 11 lw 3.0, \
  10000 title "Outofplanes" w l ls 12 lw 3.0

