/*
  This program computes the value of "pi".
  
  Run as "time ./pi2".
  
  Output should be something like:

Pi is 3.141593                        <-- (output from pi2.c)
29.7u 0.0s 0:29 99% 0+0k 0+0io 0pf+0w <-- (output from time command)
  ^    ^    ^
  |    |    (elapsed real runtime used by our program)
  |    (CPU-seconds used by system on behalf of user program)
  (CPU-seconds used directly by user program)
*/

/* This version uses explicit parallel programming */

#include <stdio.h>
#include <pthread.h>

#define CHUNK_SIZE 1000

const long num_intervals = 1000000000;
double h;
double sum;
long global_i;
pthread_mutex_t global_i_lock, sum_lock;
int NTHREADS=4;
void *compute_pi();

double f( double a )
{
  return (4.0 / (1.0 + a*a));
}

void *compute_pi(void *threadid)
{
  long i, imax;
  double x, partial_sum;
  unsigned int tryout;
  pthread_t self;
  long tid;
  tid = (long)threadid;

  fprintf (stdout, "Thread %ld has started\n", tid);

  i = 0;
  partial_sum = 0.0;

  /* We wish to reserve a chunk of "i" values every time we
     lock the global_i variable */
  while (i < num_intervals) {

    pthread_mutex_lock(&global_i_lock);
    //m_lock();
    i = global_i;
    global_i += CHUNK_SIZE;
    //m_unlock();
    pthread_mutex_unlock(&global_i_lock);

    imax = i + CHUNK_SIZE - 1;
    for ( ; i <= imax; i++) {
     	x = h * ((double)i + 0.5);
    	partial_sum += f(x);
    }
  }

  pthread_mutex_lock(&sum_lock);
  //m_lock();
  sum += partial_sum;
  //m_unlock();
  pthread_mutex_unlock(&sum_lock);
}

main()
{
  double pi;
  int i=0; 
  pthread_t threads[NTHREADS];
  h = 1.0 / (double)num_intervals;
  
  global_i = 0;
  sum      = 0.0;
  
  pthread_mutex_init(&global_i_lock,NULL);
  pthread_mutex_init(&sum_lock,NULL);

  //m_set_procs( 4 );
  for (i=0; i < NTHREADS; i++)
  {
     pthread_create(&threads[i],NULL,&compute_pi,(void*)i);
  }
  //m_fork( compute_pi ); /* perform parallel computation here */
  //m_sync(); /* wait until all processors are done */
  for(i=0; i < NTHREADS; i++)
  {
     pthread_join( threads[i], NULL); 
  } 

  pi = h * sum;

  printf( "Pi is %f\n", pi );
  pthread_mutex_destroy(&global_i_lock);
  pthread_mutex_destroy(&sum_lock);
  pthread_exit(NULL);
}
