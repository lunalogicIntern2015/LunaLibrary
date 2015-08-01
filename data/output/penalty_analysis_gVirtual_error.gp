### gnuplot file
### grid xtics , 0, 0.001, 0.01, 0.1, 1
### grid ytics , 0, 1e-006, 0.0001, 0.01, 1set xlabel "xLabel"
set ylabel "yLabel"
set title "Graph title" 
set style data lines
set hidden3d
set pm3d ; set palette
set cntrparam levels 10
set contour base
splot "penalty_analysis_gVirtual_error.dat"
pause -1 " Hit return to save image" 
set terminal postscript eps
set output "penalty_analysis_gVirtual_error.eps"
replot
