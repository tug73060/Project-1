#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#define min(x, y) ((x)<(y)?(x):(y))

#include "mmult.c"

/********************
This program tests matrix multiplication with MPI. 
It creates 2 output files for comparison.
One is the result from MPI. The other is from the simple solution defined in "mmult.c"
********************/
int main(int argc, char** argv)
{
  int nrows, ncols;
  double *A, *B, *C;
  double *CC;
  double *buffer;
  double *ans;
  double answer;
  double total_times;
  int run_index;
  int stripe;
  double starttime, endtime, timeTaken;
  int i, j, q, k, numsent, sender;
  int anstype, row;
  int myid, master, numprocs;
  MPI_Status status;
  srand(time(0));

  // Initialize the master and slaves
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  char me[255];
  gethostname(me, 254);
  printf("Hello from %s I am process %d of %d\n", me, myid, numprocs);

  FILE* outf1 = fopen("output_MPI.txt", "w+");
  FILE* outf2 = fopen("output_simple.txt", "w+");

  // Need one argument for the matrix dimension
  if (argc > 1) {
    
    nrows = atoi(argv[1]);
    ncols = nrows;
    A = (double*)malloc(sizeof(double) * nrows * ncols);
    B = (double*)malloc(sizeof(double) * nrows * ncols);
    C = (double*)malloc(sizeof(double) * nrows * ncols);

    // This is the output of the simple solution
    CC = (double*)malloc(sizeof(double) * nrows *ncols);

    // Split matrix A into stripes to send to the slaves
    // The slaves are all processes except the master
    // Below is the number of rows per stripe
    if (numprocs-1 < nrows) {
      stripe = nrows / (numprocs-1);
    } else {
      stripe = 1;
    }

    buffer = (double*)malloc(sizeof(double) * stripe * ncols);
    ans = (double*)malloc(sizeof(double) * stripe * ncols);    
    master = 0;

    if (myid == master) {
      // Master Code goes here

      // Initialize matrix A and B
      for (i = 0; i < nrows; i++) {
	      for (j = 0; j < ncols; j++) {
	        A[i*ncols + j] = (double)rand()/RAND_MAX;
          B[i*ncols + j] = (double)rand()/RAND_MAX;
          C[i*ncols + j] = 0;
	      }
      }
      
      numsent = 0;
      run_index = 0;
      starttime = MPI_Wtime();
            
      // Broadcast matrix B to all processes
      MPI_Bcast(B, ncols * nrows, MPI_DOUBLE, master, MPI_COMM_WORLD);

      printf("Hello from %s, broadcasting B successfully\n", me);
      printf("Stripe length = %d\n", stripe);

      // Send each stripe i to process i+1 
      for (i = 0; i < min(numprocs-1, nrows); i++) {
	      for (q = 0; q < stripe; q++) {
          for (j = 0; j < ncols; j++) {
            buffer[q*ncols + j] = A[i*stripe*ncols + q*ncols + j];
          }
        }    
	      MPI_Send(buffer, stripe * ncols, MPI_DOUBLE, i+1, i+1, MPI_COMM_WORLD);
	      numsent+=stripe;
        run_index = i;
        printf("%s: sent stripe %d. Numrows sent: %d\n", me, i, numsent);
      }

      // Assemble the answer matrix
      for (i = 0; i < (nrows/stripe); i++) {
        MPI_Recv(ans, stripe * ncols, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, 
		    MPI_COMM_WORLD, &status);
	      sender = status.MPI_SOURCE;
	      anstype = status.MPI_TAG - 1;
        for (q = 0; q < stripe; q++) {
          for (j = 0; j < ncols; j++) {
            C[anstype*stripe*ncols + q*ncols + j] = ans[q*ncols + j];
          }
        }
        printf("Finished with stripe %d\n", anstype);

        // If there are enough rows to fill the buffer, send them to the free slave
        // Otherwise, signal the slave to terminate
        if ((nrows - numsent) >= stripe) {
          run_index+=1;
          for (int q = 0; q < stripe; q++) {
            for (int j = 0; j < ncols; j++) {
              buffer[q*ncols + j] = A[run_index*stripe*ncols + q*ncols + j];
            }
          }
          MPI_Send(buffer, stripe * ncols, MPI_DOUBLE, sender, run_index+1,MPI_COMM_WORLD);
          numsent += stripe;
          printf("%s: sent stripe %d. Numrows sent = %d\n", me, run_index, numsent);
        } else {
	        MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD);
	      }
      } 

      // If there are rows left (less than the buffer size), the master takes care of them
      if (numsent < nrows) {
        int numleft = nrows - numsent;
        run_index++;
        buffer = (double*)realloc(buffer, sizeof(double)* ncols * numleft);
        ans = (double*)realloc(ans, sizeof(double)* ncols * numleft);
        for (int q = 0; q < numleft; q++) {
          for (int j = 0; j < ncols; j++) {
            buffer[q*ncols + j] = A[run_index*stripe*ncols + q*ncols + j];
          }
        }
        mmult(ans, buffer, numleft, ncols, B, nrows, ncols);
        for (q = 0; q < numleft; q++) {
          for (int j = 0; j < ncols; j++) {
            C[run_index*stripe*ncols + q*ncols + j] = ans[q*ncols + j];
          }
        }
        printf("Master: Finished with stripe %d\n", run_index);
      }

      // Print the total time 
      endtime = MPI_Wtime();
      timeTaken = endtime-starttime;
      printf("*** MPI runtime ***: %f\n",timeTaken);

      // The simple solution
      mmult(CC, A, nrows, ncols, B, nrows, ncols);
      fprintf(outf1,"[Matrix Result by MPI]\n");
      fprintf(outf2,"[Matrix Result: Simple Solution\n");
      for (k = 0; k < nrows; k++) {
        for (j = 0; j < ncols; j++) {
          fprintf(outf1, "%.2f\t", C[k*ncols + j]);
          fprintf(outf2, "%.2f\t", CC[k*ncols + j]);
        }
        fprintf(outf1, "\n");
        fprintf(outf2, "\n");
      }
      
      free(A);
      free(B);
      free(C); free(CC); free(buffer); free(ans);
      
    }

  // Slave Code
    else {
      // Slave Code goes here
      MPI_Bcast(B, ncols * nrows, MPI_DOUBLE, master, MPI_COMM_WORLD);
      if (myid <= numprocs) {
	      while(1) {
	        MPI_Recv(buffer, stripe * ncols, MPI_DOUBLE, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	        if (status.MPI_TAG == 0){
	          break;
          }
	        row = status.MPI_TAG;
        
          printf("-- Slave process %d: received stripe %d successfully\n", myid, row-1);
          /// Multiply the stripe by matrix B 
          mmult(ans, buffer, stripe, ncols, B, nrows, ncols);
                    
	        MPI_Send(ans, stripe * ncols, MPI_DOUBLE, master, row, MPI_COMM_WORLD);
	      }
      }
    }    

  } else {
    fprintf(stderr, "Usage matrix_times_vector <size>\n");
  }

  MPI_Finalize();
  //fclose(outf1);
  //fclose(outf2);
  return 0;
}

