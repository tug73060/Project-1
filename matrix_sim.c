
#include <stdio.h>
#include <time.h>
#include <sys/times.h>

#include "mmult.c"

/****
* Generate random square matrices and mutliply them
* Time the time it takes for each iteration
* Output data as 
* 
* Current implementation: Single threaded
****/


int main(int argc, char** argv)
{
	if(argc > 3)
	{	
		int mlen = atol(argv[1]);
		int numIte = atol(argv[2]);
		FILE* outf = fopen(argv[3], "w+");

		double* M = NULL;
		double* A = NULL;
		double* C = gen_matrix(mlen, mlen);

		printf("Matrix Len: %d\nNum Iterations: %d\n", mlen, numIte);
		int i;
		for(i = 0; i < numIte; i++)
		{
			M = gen_matrix(mlen, mlen);
			A = gen_matrix(mlen, mlen);

			mmult(C, M, mlen, mlen, A, mlen, mlen);


			free(M);
			free(A);
			
			printf("Done %d\n", i);
		}	
		free(C);
		fclose(outf);
	}
	else{
		printf("ARGUMENTS: [Matrix Len]  [# of Trials]  [Output File Name]\n");
		return 0; 
	}
}
