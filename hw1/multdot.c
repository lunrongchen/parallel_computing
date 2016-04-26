/*
  This program computes a matrix-vector multiplication using
  the dot-product method.
  
  Run as "time ./multdot".
  
  Output should be something like:

12.253u 11.784s 0:08.11 296.3% 0+0k 1+5io 1pf+0w <-- (output from time command)
  ^       ^      ^
  |       |      (elapsed real runtime used by our program)
  |       (CPU-seconds used by system on behalf of user program)
  (CPU-seconds used directly by user program)
*/

#include <omp.h>

/* We wish to compute c = A * b */
double a[15000][15000], b[15000], c[15000];

main()
{
  int row, i;
  double row_sum;

#pragma omp parallel private(row,i,row_sum) num_threads(8)
  {
    /* We want to compute each row solution in parallel */
#pragma omp for schedule(static)
    for (row = 0; row < 15000; row++) {
      row_sum = 0;

      /* This next loop performs the dot-product */
      for (i = 0; i < 15000; i++) {
	row_sum += a[row][i] * b[i];
      }
      
      c[row] = row_sum;
    }
  }
}

