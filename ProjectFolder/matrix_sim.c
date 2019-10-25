#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#include "mmult.c"


///////////////////////////////////////////////////////////////
/// dynamically changing struct depending on testing variables
/// Update curLen and timeTaken before calling logResult()
typedef struct
{
	int startLen;
	int curLen;
	int numIte;
	double timeTaken;
	struct timespec start;
	struct timespec end; 
}TestData;


/////////////////////////////////////////////////////
/// deltaTime() equation from mmul_omp_timing.c file
double deltaTime(struct timespec* start, struct timespec* end)
{
	double delta = (end->tv_sec - start->tv_sec) + (end->tv_nsec - start->tv_nsec)/1e9;
	return delta;
}

////////////////////////////////////
/// call log to print to file
/// > prints the MATRIX SIZE
/// > prints the TIME IT TOOK
void logResult(FILE* f, TestData* data)
{
	fprintf(f, "%d\t%f\n", data->curLen, data->timeTaken);
}

////////////////////////////////////////////////////////
/// testing single thread without any modifications
/// 1. Start with matrix of size Testdata.startLen
/// 2. Multiply the matrices
/// 3. Update TestData.curLen and call logResult()
/// 4. Increment the size by 1 and multiple the matrices
/// 5. Repeat 2-4 until the size reach TestData.numIte
void TestLinear(TestData* data, FILE* f)
{
	double* M = NULL;
	double* A = NULL;
	double* C = NULL;	
	int templen = data->startLen;

	int i;
	for(i = 0; i < data->numIte; i++)
	{
		templen = data->startLen + i;
		printf("Matrix Len: %d\nNum Iterations: %d\n", templen, data->numIte);	
		M = gen_matrix(templen, templen);
		A = gen_matrix(templen, templen);
		C = gen_matrix(templen, templen);

		
		clock_gettime(CLOCK_REALTIME, &(data->start));
		mmult(C, M, templen, templen, A, templen, templen);
		clock_gettime(CLOCK_REALTIME, &(data->end));
		
			
		data->timeTaken = deltaTime(&data->start, &data->end);
		
		///updating the current matrix size before printing
		data->curLen = templen;
		logResult(f, data);
		
		free(M);
		free(A);

	}
	free(C);
}

///////////////////////////////////////////////////
/// Generate array of doubles using a file
/// file content example:
///
/// 2
/// 2.32    2423.5
/// 984.02  1.11
double* gen_matrix_file(FILE* f, int* size)
{
	char buff[255];
	if(fscanf(f, "%s", buff) == 1)
	{
		if((*size = atoi(buff)) != 0)
		{
			printf("%d\n", *size);
		}
		else{ exit(-1); }
	}
	else { exit(-1); }

	double* matrix = (double*)malloc((*size) * (*size) * sizeof(double));

	int dimension = (*size) * (*size);
	int i = 0;
	while(fscanf(f, "%s", buff) == 1 && (i < dimension))
	{
		matrix[i] = atof(buff);
		printf("Loaded %f\n", matrix[i]);
		i++;
	}
}


int main(int argc, char** argv)
{
	char me[255];
	gethostname(me, 254);

	printf("\nHost: %s\n", me);


	/// checking for correct arguments
	if(argc > 3)
	{
		/// setup output file and test data
		FILE* outf1 = fopen(argv[3], "w+");
		//FILE* outf2 = fopen("test_simd.txt", "w+");
		TestData data = {0};
		data.startLen = atol(argv[1]);
		data.numIte = atol(argv[2]);
		data.timeTaken = 0.0;
		data.curLen = data.startLen;
		
		/// Start linear test

		fprintf(outf1, "\n\n[LINEAR TIME TEST]\n");
		TestLinear(&data, outf1);

		/*
		fprintf(outf, "\n\n[LINEAR TIME TEST]\n");

		TestLinear(&data, outf);
		printf("\n");
		
		*/
		

		/// Start other tests here ?


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

}
