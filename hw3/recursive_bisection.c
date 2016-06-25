#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM_POINTS 524288
unsigned int X_axis[NUM_POINTS];
unsigned int Y_axis[NUM_POINTS];

int i = 0, j = 0, k = 0, z = 0;
//global_cost computed at process 0
double global_cost = 0;

int numprocs; /* Number of processors to use */
int myid;

int num_quadrants;

//find the kth smallest element in array v
unsigned int find_kth(unsigned int * v, int n, int k, unsigned int * y);

// swap elemnt at index i and j in array
void swap(unsigned int array[], int i, int j) {
    int temp = array[i];
    array[i] = array[j];
    array[j] = temp;
}

//initialize quadrant
void initialize(int * min_x, int * max_x, int * min_y, int * max_y) {
    for (i = 1; i < NUM_POINTS; i++) {
        if (X_axis[i] < ( * min_x)) {
            ( * min_x) = X_axis[i];
        }
        if (X_axis[i] > ( * max_x)) {
            ( * max_x) = X_axis[i];
        }
        if (Y_axis[i] < ( * min_y)) {
            ( * min_y) = Y_axis[i];
        }
        if (Y_axis[i] > ( * max_y)) {
            ( * max_y) = Y_axis[i];
        }
    }
}

//bisection on X 
void bisection_X(int * x_cut, int quadrants, int points, int * pivot_arr) {
    for (i = 0; i < quadrants; i++) {
        //find  median  X_axis[]
        int x_pivot = find_kth(X_axis + i * points, points, points / 2 - 1, Y_axis + i * points);
        k = i * points;
        j = i * points + points - 1;
        while (k < j && k < k + points / 2) {
            if (X_axis[k] > x_pivot) {
                while (X_axis[j] > x_pivot) {
                    j--;
                }
                if (k < j) {
                    swap(X_axis, k, j);
                    swap(Y_axis, k, j);
                }
            }
            k++;
        }
        //save the dividing point
        pivot_arr[i + quadrants - 1] = x_pivot;
    }
    //bisection_Y's turn
    ( * x_cut) = 1;
}

//bisection on Y
void bisection_Y(int * x_cut, int quadrants, int points, int * pivot_arr) {
    for (i = 0; i < quadrants; i++) {
        int y_pivot = find_kth(Y_axis + i * points, points, points / 2 - 1, X_axis + i * points);
        k = i * points;
        j = i * points + points - 1;
        while (k < j && k < k + points / 2) {
            if (Y_axis[k] > y_pivot) {
                while (Y_axis[j] > y_pivot) {
                    j--;
                }
                if (k < j) {
                    swap(X_axis, k, j);
                    swap(Y_axis, k, j);
                }
            }
            k++;
        }
        //save the dividing point
        pivot_arr[i + quadrants - 1] = y_pivot;
    }
    //bisection_X's turn
    ( * x_cut) = 0;
}

//divide the quadrant by X
void divide_x(int * top, int * bot, int * lef, int * rht, int * pivot_arr, int quadrants, int * temp, int * x_cut) {
    for (j = 0; j < quadrants * 2; j += 2) {
        top[j] = top[quadrants + j / 2];
        bot[j] = bot[quadrants + j / 2];
        lef[j] = lef[quadrants + j / 2];
        rht[j] = pivot_arr[( * temp)];

        top[j + 1] = top[quadrants + j / 2];
        bot[j + 1] = bot[quadrants + j / 2];
        lef[j + 1] = pivot_arr[( * temp)];
        rht[j + 1] = rht[quadrants + j / 2];
        ( * temp) ++;
    }
    ( * x_cut) = 1;
}

//divide the quadrant by Y
void divide_Y(int * top, int * bot, int * lef, int * rht, int * pivot_arr, int quadrants, int * temp, int * x_cut) {
    for (j = 0; j < quadrants * 2; j += 2) {
        top[j] = top[quadrants + j / 2];
        bot[j] = pivot_arr[( * temp)];
        lef[j] = lef[quadrants + j / 2];
        rht[j] = rht[quadrants + j / 2];

        top[j + 1] = pivot_arr[( * temp)];
        bot[j + 1] = bot[quadrants + j / 2];
        lef[j + 1] = lef[quadrants + j / 2];
        rht[j + 1] = rht[quadrants + j / 2];
        ( * temp) ++;
    }
    ( * x_cut) = 0;
}

//find the kth smallest element in array v
unsigned int find_kth(unsigned int * v, int n, int k, unsigned int * y) {
    int j0 = 0;
    int i1 = 0;
    int j1 = 0;
    if (n == 1 && k == 0) return v[0];
    int m = (n + 4) / 5;
    unsigned int * medians = (unsigned int * ) malloc(m * sizeof(int));
    for (i1 = 0; i1 < m; i1++) {
        if (5 * i1 + 4 < n) {
            unsigned int * w = v + 5 * i1;
            unsigned int * w1 = y + 5 * i1;
            for (j0 = 0; j0 < 3; j0++) {
                int jmin = j0;
                for (j1 = j0 + 1; j1 < 5; j1++) {
                    if (w[j1] < w[jmin]) jmin = j1;
                }
                swap(w, j0, jmin);
                swap(w1, j0, jmin);
            }
            medians[i1] = w[2];
        } else {
            medians[i1] = v[5 * i1];
        }
    }
    int pivot = find_kth(medians, m, m / 2, medians);
    free(medians);
    for (i1 = 0; i1 < n; i1++) {
        if (v[i1] == pivot) {
            swap(v, i1, n - 1);
            swap(y, i1, n - 1);
            break;
        }
    }
    int store = 0;
    for (i1 = 0; i1 < n - 1; i1++) {
        if (v[i1] < pivot) {
            swap(v, i1, store);
            swap(y, i1, store);
            store++;
        }
    }
    swap(v, store, n - 1);
    swap(y, store, n - 1);
    if (store == k) {
        return pivot;
    } else if (store > k) {
        return find_kth(v, store, k, y);
    } else {
        return find_kth(v + store + 1, n - store - 1, k - store - 1, y + store + 1);
    }
}

void find_quadrants(num_quadrants) {
    if (myid == 0) {
        int x_cut = 0;
        int quadrants = 1;
        //coordinates of the quadrants, 
        int top[num_quadrants]; //the smallest y axis
        int bot[num_quadrants]; //the largest  y axis
        int lef[num_quadrants]; //the smallest x axis
        int rht[num_quadrants]; //the largest  x axis
        //record the bisection coordinates
        int pivot_arr[num_quadrants];
        //Step1
        //ends when there are enough quadrants
        while (num_quadrants > quadrants) {
            int points = NUM_POINTS / quadrants;
            if (!x_cut) {
                bisection_X( & x_cut, quadrants, points, pivot_arr);
            } else {
                //bisection on Y  
                bisection_Y( & x_cut, quadrants, points, pivot_arr);
            }
            quadrants *= 2;
        }

        int min_x = X_axis[0], max_x = X_axis[0], min_y = Y_axis[0], max_y = Y_axis[0];
        initialize( & min_x, & max_x, & min_y, & max_y);

        //Step2
        //update the quadrants' coordinates
        i = 0;
        x_cut = 0;
        quadrants = 1;
        top[0] = min_y;
        bot[0] = max_y;
        lef[0] = min_x;
        rht[0] = max_x;
        while (i < num_quadrants - 1) {
            int temp = i;
            for (j = 0; j < quadrants; j++) {
                top[j + quadrants] = top[j];
                bot[j + quadrants] = bot[j];
                lef[j + quadrants] = lef[j];
                rht[j + quadrants] = rht[j];
            }
            if (!x_cut) {
                divide_x(top, bot, lef, rht, pivot_arr, quadrants, & temp, & x_cut);
            } else {
                divide_Y(top, bot, lef, rht, pivot_arr, quadrants, & temp, & x_cut);
            }
            i += quadrants;
            quadrants *= 2;
        }
        //print out the four coordinates
        printf("\nQuadrants coordinates:\nQuarant number :  (top-left) | (top-right) | (bottom-left) | (bottom-righ)\n");
        for (z = 0; z < num_quadrants; z++) {
            printf("Number %d : ", z);
            printf(" (%d, %d) |", lef[z], top[z]);
            printf(" (%d, %d) |", rht[z], top[z]);
            printf(" (%d, %d) |", lef[z], bot[z]);
            printf(" (%d, %d)\n", rht[z], bot[z]);
        }
    }

    //Step3
    //broadcast the X and Y to other proesses
    MPI_Bcast( & X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( & Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);

    //calculate the local_cost on each process
    double local_cost = 0;
    for (i = myid; i < num_quadrants; i += numprocs) {
        int points = NUM_POINTS / num_quadrants;
        for (j = 0; j < points - 1; j++) {
            for (k = j + 1; k < points; k++) {
                int x1 = points * i + j;
                int x2 = points * i + k;
                int y1 = points * i + j;
                int y2 = points * i + k;

                double diff_x = abs(X_axis[x1] - X_axis[x2]);
                double diff_y = abs(Y_axis[y1] - Y_axis[y2]);
                local_cost += sqrt((double) diff_x * diff_x + diff_y * diff_y);
            }
        }
    }
    //Step4
    //reduce and calculate the global_cost on process 0
    MPI_Reduce( & local_cost, & global_cost, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

int main(argc, argv)
int argc;
char * argv[]; {
    int num_quadrants;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    /*Time Variables*/
    double startwtime = 0.0, endwtime;

    MPI_Init( & argc, & argv);
    MPI_Comm_size(MPI_COMM_WORLD, & numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, & myid);
    MPI_Get_processor_name(processor_name, & namelen);

    if (argc != 2) {
        fprintf(stderr, "Usage: recursive_bisection <#of quadrants>\n");
        MPI_Finalize();
        exit(0);
    }

    fprintf(stderr, "Process %d on %s\n", myid, processor_name);

    num_quadrants = atoi(argv[1]);

    if (myid == 0) {
        fprintf(stdout, "Extracting %d quadrants with %d processors \n", num_quadrants, numprocs);
    }
    if (myid == 0) {
        int i;
        srand(10000);
        for (i = 0; i < NUM_POINTS; i++)
            X_axis[i] = (unsigned int) rand();
        for (i = 0; i < NUM_POINTS; i++)
            Y_axis[i] = (unsigned int) rand();
        //start timer at process 0
        printf("\nComputing Parallely Using MPI.\n");
        startwtime = MPI_Wtime();
    }

    MPI_Bcast( & X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( & Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);

    find_quadrants(num_quadrants);
    //barrier to guarantee the correctness of timing
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid == 0) {
        //end timer at process 0
        endwtime = MPI_Wtime();
        printf("\nelapsed time = %f\n", endwtime - startwtime);
        printf("\nTotal cost:  %lf \n", global_cost);

    }
    MPI_Finalize();
    return 0;
}