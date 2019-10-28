#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#define min(x, y) ((x)<(y)?(x):(y))

#include "mmult.c"

typedef struct
{
	int startLen;
	int curLen;
	int numIte;
	double timeTaken;
	struct timespec start;
	struct timespec end; 
}TestData;


int nrows, ncols;
double *A, *B, *C;
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


////////////////////////////////////
/// call log to print to file
/// > prints the MATRIX SIZE
/// > prints the TIME IT TOOK
void logResult(FILE* f, TestData* data)
{
	fprintf(f, "%d\t%f\n", data->curLen, data->timeTaken);
}


/////////////////////////////////////////////////////
/// deltaTime() equation from mmul_omp_timing.c file
double deltaTime(struct timespec* start, struct timespec* end)
{
	double delta = (end->tv_sec - start->tv_sec) + (end->tv_nsec - start->tv_nsec)/1e9;
	return delta;
}


////////////////////////////////////////////////////////
/// testing single thread without any modifications
/// 1. Start with matrix of size Testdata.startLen
/// 2. Multiply the matrices
/// 3. Update TestData.curLen and call logResult()
/// 4. Increment the size by 1 and multiple the matrices
/// 5. Repeat 2-4 until the size reach TestData.numIte
void TestLinear(TestData* data, FILE* f, int numprocs)
{
    srand(time(0));
	int k;
	for(k = 0; k < data->numIte; k++)
	{
        nrows = data->startLen + k;
        ncols = nrows;
        data->curLen = nrows;
        A = (double*)malloc(sizeof(double) * nrows * ncols);
        B = (double*)malloc(sizeof(double) * nrows * ncols);
        C = (double*)malloc(sizeof(double) * nrows * ncols);

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
            
            clock_gettime(CLOCK_REALTIME, &(data->start));            
            // Broadcast matrix B to all processes
            MPI_Bcast(B, ncols * nrows, MPI_DOUBLE, master, MPI_COMM_WORLD);

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
                //printf("%s: sent stripe %d to sender %d successfully. Numrows sent: %d\n", me, i, i+1, numsent);
            }

            // Assemble the answer matrix
            for (i = 0; i < (nrows/stripe); i++) {
                MPI_Recv(ans, stripe * ncols, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	            sender = status.MPI_SOURCE;
	            anstype = status.MPI_TAG - 1;
                for (q = 0; q < stripe; q++) {
                    for (j = 0; j < ncols; j++) {
                        C[anstype*stripe*ncols + q*ncols + j] = ans[q*ncols + j];
                    }
                }

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
                    //printf("%s: sent stripe %d to sender %d successfully. Numrows sent = %d\n", me, run_index, sender, numsent);
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
            //printf("Finished outputing for stripe %d\n", run_index);
            }

            // Print the total time 
            clock_gettime(CLOCK_REALTIME, &(data->end));
            data->timeTaken = deltaTime(&data->start, &data->end);
            logResult(f, data);
      
            free(A);
            free(B);
            free(C); free(buffer); free(ans);
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
                /// Multiply the stripe by matrix B 
                mmult(ans, buffer, stripe, ncols, B, nrows, ncols);
                    
	            MPI_Send(ans, stripe * ncols, MPI_DOUBLE, master, row, MPI_COMM_WORLD);
	            }
            }
        }        
    }

}

int main(int argc, char** argv)
{
    /// Initialize the master and slaves
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	/// checking for correct arguments
	if(argc > 3)
	{
		/// setup output file and test data
		FILE* outf1 = fopen(argv[3], "w+");
        TestData data = {0};
		data.startLen = atol(argv[1]);
		data.numIte = atol(argv[2]);
		data.curLen = data.startLen;

		fprintf(outf1, "\n\n[MPI TIME TEST]\n");
        TestLinear(&data, outf1, numprocs);

		fclose(outf1);
		//fclose(outf2);
	}
	else
	{	
		/// usage of arguments

		printf("USAGE:\t[START LEN] [# OF TRIALS] [OUTPUT FILE NAME]\n");
		printf("W/ CLUSTER: mpiexec -f ~/hosts -n x [program] [args]\n");

		return 0;
	}
    MPI_Finalize();
}

