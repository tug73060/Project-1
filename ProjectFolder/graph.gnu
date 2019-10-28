set terminal png size 1000,800
set output 'out.png'

set xlabel 'Matrix Size'
set ylabel 'Time (sec)'
set title 'Matrix Multiplcation Speed on Wolfgand Cluster'
plot "test_regular.txt" using 1:2 title 'SIMD w/o -O3' with linespoint, "test_simd.txt" using 1:2 title 'SIMD /w -O3' with linespoint, "test_omp.txt" using 1:2 title 'OpenMP' with linespoint, "test_mpi.txt" using 1:2 title "MPI(10)" with linespoint
