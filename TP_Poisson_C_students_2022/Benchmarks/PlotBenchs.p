#Type of output and name of it
set terminal pdf
set output "DGBSV_Solve.pdf"

set xlabel "Size of the matrix (NxN)"
set ylabel "Time in Nanoseconds"

plot 'DGBSV_Solve/DGBSV.dat' w l title 'DGBSV Performances', \

#'DGBTRFTRI_DGBTRS_Solve/DGBTRS.dat' w l title 'DGBTRS Performances', \
     
pause -1 "Hit any key to continue"