/*
  This program computes the value of "pi".
  
  Run as "time ./pi1".
  
  Output should be something like:

Pi is 3.141593                        <-- (output from pi1.c)
29.7u 0.0s 0:29 99% 0+0k 0+0io 0pf+0w <-- (output from time command)
  ^    ^    ^
  |    |    (elapsed real runtime used by our program)
  |    (CPU-seconds used by system on behalf of user program)
  (CPU-seconds used directly by user program)
*/

/* This version uses parallelization directives */

#include <stdio.h>
#include <omp.h>

const long num_intervals = 1000000000;
                           
double f( double a )
{
  return (4.0 / (1.0 + a*a));
}

main()
{
  double h, sum, partial_sum, x, pi;
  long i;

  h   = 1.0 / (double)num_intervals;
  sum = 0.0;

#pragma omp parallel private(i,x,partial_sum) num_threads(4)
  {
    partial_sum = 0.0;
    
#pragma omp for schedule(static)
    for (i = 0; i < num_intervals; i++) {
      x = h * ((double)i + 0.5);
      partial_sum += f(x);
    }

#pragma omp critical
    {
      sum += partial_sum;
    }
  }

  pi = h * sum;

  printf( "Pi is %f\n", pi );
}

