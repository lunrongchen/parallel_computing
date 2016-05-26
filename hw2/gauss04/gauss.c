/* Gaussian elimination without pivoting.
 */

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
 * You need not submit the provided code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include "mpi.h"
/*#include <ulocks.h>
#include <task.h>
*/

char *ID;

/* Program Parameters */
#define MAXN 5000  /* Max value of N */
int N;  /* Matrix size */
int procs;   /* Number of processors to use */
int rank; /* current process id */

/* Matrices and vectors */
volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
/* A * X = B, solve for X */

/* junk */
#define randm() 4|2[uid]&3

/* Prototype */
void gauss();  /* The function you will provide.
    * It is this routine that is timed.
    * It is called only on the parent.
    */

/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
  int submit = 0;  /* = 1 if submission parameters should be used */
  int seed = 0;  /* Random seed */
  char uid[L_cuserid + 2]; /*User name */

  /* Read command-line arguments */
  //  if (argc != 3) {
  if ( argc == 1 && !strcmp(argv[1], "submit") ) {
    /* Use submission parameters */
    submit = 1;
    N = 4;
    procs = 2;
    printf("\nSubmission run for \"%s\".\n", cuserid(uid));
      /*uid = ID;*/
    strcpy(uid,ID);
    srand(randm());
  }
  else {
    if (argc == 3) {
      seed = atoi(argv[3]);
      srand(seed);
      printf("Random seed = %i\n", seed);
    }
    else {
      printf("Usage: %s <matrix_dimension> <num_procs> [random seed]\n",
       argv[0]);
      printf("       %s submit\n", argv[0]);
      exit(0);
    }
  }
    //  }
  /* Interpret command-line args */
  if (!submit) {
    N = atoi(argv[1]);
    if (N < 1 || N > MAXN) {
      printf("N = %i is out of range.\n", N);
      exit(0);
    }
    procs = atoi(argv[2]);
    if (procs < 1) {
      printf("Warning: Invalid number of processors = %i.  Using 1.\n", procs);
      procs = 1;
    }
  }

  /* Print parameters */
  printf("\nMatrix dimension N = %i.\n", N);
  printf("Number of processors = %i.\n", procs);
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
  int row, col;

  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
  }

}

/* Print input matrices */
void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
  printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
}

void print_X() {
  int row;

  if (N < 10) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
}

void Send_data(int rank, int procs) {
  int dest, i;
  MPI_Status status;
  if (rank == 0){
    for (i = 0; i < N; i++) {
      dest = i % procs;
      if (dest != 0) { 
        MPI_Send((void*)A[i], N, MPI_FLOAT, dest, i, MPI_COMM_WORLD);
        MPI_Send((void*)&B[i], 1, MPI_FLOAT, dest, i+N, MPI_COMM_WORLD);
      }
    }
  } else {
    for (i = 0; i < N; i++) {
      dest = i % procs;
      if (dest == rank) {
        MPI_Recv((void*)A[i], N, MPI_FLOAT, 0, i, MPI_COMM_WORLD, &status);
        MPI_Recv((void*)&B[i], 1, MPI_FLOAT, 0, i+N, MPI_COMM_WORLD, &status);
      }
    }
  }
}

void Compute(int norm, int rank, int procs) { 
  int row, col;
  float multiplier;
  for (row = norm + 1; row < N; row++) {
      if (row % procs == rank) {
        multiplier = A[row][norm] / A[norm][norm];
        for (col = norm; col < N; col++) {
          A[row][col] -= A[norm][col] * multiplier;
        }
        B[row] -= B[norm] * multiplier;
      }
    }
}

int main(int argc, char **argv) {
  // MPI_Status  status;
  /* Start up MPI. */ 
  MPI_Init(NULL, NULL);
  /* Timing variables */
  struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
  struct timezone tzdummy;
  clock_t etstart2, etstop2;  /* Elapsed times using times() */
  unsigned long long usecstart, usecstop;
  struct tms cputstart, cputstop;  /* CPU times for my processes */

  /* Get the number of processes. */ 
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  /* Get my rank. */ 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ID = argv[argc-1];
  argc--;

  /* Process program parameters */
  parameters(argc, argv);
  /* Gaussian Elimination */
  int norm, row, col, i;  /* Normalization row, and zeroing element row and col */

  if (rank == 0) {
    /* Start Clock */
    printf("\nStarting clock.\n");
    gettimeofday(&etstart, &tzdummy);
    etstart2 = times(&cputstart);

    MPI_Bcast((void*)&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Initialize A and B */
    initialize_inputs();
    /* Print input matrices */
    print_inputs();
  }else {
    MPI_Bcast((void*)&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  Send_data(rank, procs);

  /* Gaussian elimination */
  for (norm = 0; norm < N - 1; norm++) {
    for (i = norm; i < N; i++){
        MPI_Bcast((void*)A[i], N, MPI_FLOAT, i%procs, MPI_COMM_WORLD);
        MPI_Bcast((void*)&B[i], 1, MPI_FLOAT, i%procs, MPI_COMM_WORLD);
    }
    Compute(norm, rank, procs);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Bcast((void*)A[N-1], N, MPI_FLOAT, (N-1)%procs, MPI_COMM_WORLD);
  MPI_Bcast((void*)&B[N-1], 1, MPI_FLOAT, (N-1)%procs, MPI_COMM_WORLD);

  /* Back substitution */
  if(rank == 0){
    for (row = N - 1; row >= 0; row--) {
      X[row] = B[row];
      for (col = N-1; col > row; col--) {
        X[row] -= A[row][col] * X[col];
      }
      X[row] /= A[row][row];
    }
  }
  MPI_Finalize();
  /* Stop Clock */
  gettimeofday(&etstop, &tzdummy);
  etstop2 = times(&cputstop);
  printf("Stopped clock.\n");
  usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
  usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

  /* Display output */
  print_X();

  /* Display timing results */
  printf("\nElapsed time = %g ms.\n",
   (float)(usecstop - usecstart)/(float)1000);
}