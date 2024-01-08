#Type of output and name of it
set terminal pdf
set output "DGBSV_Results.pdf"

set xlabel "Size of the matrix (NxN)"
set ylabel "Relative Forward Error"

plot 'DGBSV_Solve/result.dat' w l title 'Results Precision', \
     
pause -1 "Hit any key to continue"