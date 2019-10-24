int nrows, ncols;
	double *aa;	/* the A matrix */
	double *bb;	/* the B matrix */
	double *cc1;	/* A x B computed using the omp-mpi code you write */
	double *cc2;	/* A x B computed using the conventional algorithm */
	int myid, numprocs;
	double starttime, endtime;

	if(argc > 2)
	{
		FILE* f1 = fopen(argv[1], "r+");
		FILE* f2 = fopen(argv[2], "r+");
		FILE* f3 = fopen("out_C.txt", "w+");
		if(f1 == NULL || f2 == NULL) exit(-1);

		int size = 0;
		printf("%s\n%s\n", argv[1], argv[2]);		
	
		
		gen_matrix_file(f1, &size);
		gen_matrix_file(f2, &size);
		
	
			
		
		if (myid == 0) {
      			// Master Code goes here
			aa = gen_matrix_file(f1, nrows);
			bb = gen_matrix_file(f2, nrows);
			ncols = nrows;
			cc1 = malloc(sizeof(double) * nrows * nrows); 
			
			starttime = MPI_Wtime();
			/* Insert your master code here to store the product into cc1 */
			mmult(cc1, aa, nrows, ncols, bb, nrows, ncols);
			endtime = MPI_Wtime();
			printf("%f\n",(endtime - starttime));
			cc2  = malloc(sizeof(double) * nrows * nrows);
			mmult(cc2, aa, nrows, ncols, bb, ncols, nrows);
			compare_matrices(cc2, cc1, nrows, nrows);
      	
			/// Print resultant matrix to file ///
			fprintf(f3, "%d\n", nrows);
			int y = 0, x = 0;            
			
			for(y = 0; y < nrows; y++)
			{
				for(x = 0; x < ncols; x++)
				{
					fprintf(f3, "%f ", cc1[y*nrows + x]);
				}
				fprintf(f3, "\n");
			}

			free(aa);
			free(bb);
			free(cc);
		} else {
			// Slave Code goes here
    		}

		fclose(f1);
		fclose(f2);
		fclose(f3);

	} else {
		fprintf(stderr, "Usage matrix_times_vector <size>\n");
	}
