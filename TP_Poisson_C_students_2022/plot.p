#Type of output and name of it
set terminal pdf
set output "Convergence_Gauss_Siedel.pdf"

set logscale y 10
set xlabel "Iterations"
set ylabel "log_1_0(Residual)"

plot 'bin/RESVEC.dat' w l title 'Residual Value', \
     
pause -1 "Hit any key to continue"