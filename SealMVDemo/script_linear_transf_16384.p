# Set the output terminal
set terminal canvas
set output "canvas_linear_transf_16384.html"
set title "Linear Transformation Benchmark 16384"
set xlabel 'Dimension'
set ylabel 'Time (microseconds)'
set logscale
set ytics nomirror
set xtics nomirror
set grid
set key outside

# Set the styling 
set style line 1\
linecolor rgb '#0060ad'\
linetype 1 linewidth 2\
pointtype 7 pointsize 1.5

set style line 2\
linecolor rgb '#dd181f'\
linetype 1 linewidth 2\
pointtype 5 pointsize 1.5


plot 'linear_transf_16384.dat' index 0 title "C_Vec * P_Mat" with linespoints ls 1, \
'' index 1 title "C_Vec * C_Mat"  with linespoints ls 2