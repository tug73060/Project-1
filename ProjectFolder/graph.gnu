set terminal png size 800,600
set output 'out.png'

set xlabel 'Matrix Size'
set ylabel 'Time (sec)'
set title 'Matrix Multiplcation Speed on Wolfgand Cluster'
plot "test_regular.txt" using 1:2 title 'No -O3' with linespoint, "test_simd.txt" using 1:2 title 'With -O3' with linespoint, "test_omp.txt" using 1:2 title 'OpenMP' with linespoint
