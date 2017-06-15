#include "solver.h" // contains other header files

//TODO: more output (e.g. % done, ETA)

//        printf("%s, %s, %d\n", __FUNCTION__, __FILE__, __LINE__);

// Main simulation program 
int main(int argc, char **argv)
{

	// Initialisation step. Construct and initialise the ODE workspace
    ode_workspace *odews;

    MPI_Init(&argc, &argv);			// initialise MPI
    odews = malloc(sizeof *odews); 	// allocate memory

    int verbose = 1;

    // Problem parameters
    odews->dt  	  = 1e-5; 			// time step  1e-5
    odews->t0     = 0;   			// initial time 0
    odews->tf     = 10;  			// final time  10
    odews->ftol   = 1e-3; 			// function evaluation tolerance for Newton convergence 1e-3
    odews->ytol   = 1e-3; 			// relative error tolerance for Newton convergence 1e-3
    odews->nconv  = 5;    			// Newton iteration threshold for Jacobian reevaluation 5
    odews->maxits = 100;   			// Maximum number of Newton iterations 100
    odews->dtwrite = 0.1; 			// Time step for writing to file (and screen)

    // If optional command line argument is passed change the default final time.
    if (argc > 3)
    {
    	odews->tf = atoi(argv[3]);
    }

    // Initialise the solver with all the bits and pieces
    solver_init(odews, argc, argv);

    // Print out the adjacency matrix A containing the structure of the H tree
//    if (odews->W->rank == 0)
//    {
//		sparseprint(odews->W->A); // print adjacency matrix
//		printf("\n ------------------ \n");
//    }

    // Begin computation. t0 and tf just measure elapsed simulation time
    double t0 = MPI_Wtime();
    back_euler(odews); // All of the simulation occurs in here
    double tf = MPI_Wtime();

    odews->W->ntimestamps = (odews->tf-odews->t0)/odews->dt;

    // Display diagnostics
    if (odews->W->rank == 0)
    {
        if (verbose)
        {
            printf("Levels: %d Subtree size: %d N procs: %d\n", odews->W->N, odews->W->Nsub, odews->W->n_procs);
            printf("Solution time:                %g seconds\n", tf - t0);
            printf("    # fevals:                 %d\n", odews->W->fevals);
            printf("    # Jacobians:              %d\n", odews->W->jacupdates);
            printf("    # feval time:             %g seconds\n", odews->W->tfeval);
            printf("    # Jacobian update time:   %g seconds\n", odews->W->tjacupdate);
            printf("    # Jacobian symbolic time: %g seconds\n", odews->W->tjacfactorize);
        }
        else
        {
            printf("%4d%4d%4d%12.4e%4d%4d%12.4e%12.4e%12.4e\n", odews->W->N, odews->W->Nsub, odews->W->n_procs, tf - t0, odews->W->fevals, odews->W->jacupdates, odews->W->tfeval, odews->W->tjacupdate, odews->W->tjacfactorize);
        }
    }

    // And clean up
    close_io(odews->W);
    free_var(odews);

    MPI_Finalize();

    return 0;
}
