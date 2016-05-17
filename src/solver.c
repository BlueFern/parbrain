#include "solver.h"

// Fixed step Backward Euler ODE solver
void back_euler(ode_workspace *odews)
{
    // Declare and initialise additional workspace variables
    double *beta, *w, *x; // arrays
    workspace *W;
    W = odews->W;
    int ny = W->nu;
    beta = zerosv(ny);
    w    = zerosv(ny);
    x    = zerosv(ny);

    // Time stuff used later to output progress of the program at each iteration
	time_t start_t = 0;
	time_t end_t = 0;
	double total_t = 0;
	// Record the time at the start of program
	time(&start_t);

    double t = odews->t0;
    double tnext;

    // Initial newton_matrix computation
    newton_matrix(odews);
    int jac_needed = 0; // flag to say that Jacobian is current

    int converged = 0;
    write_data(W, t, odews->y); // Write initial data to file
    write_flow(W, t, W->q, W->q0);
    write_pressure(W, t, W->p, W->p0);

    MPI_Barrier(MPI_COMM_WORLD);

    // Timestep loop
    for (int i = 0; t < odews->tf; i++)
    {
        // Perform a Jacobian update if necessary. This is in sync across
        // all processors
        if (jac_needed)
        {
            jacupdate(W, t, odews->y);
            jac_needed = 0;
            newton_matrix(odews);
        }

        // Copy values from previous completed timestep
        dcopy(ny, odews->y, beta);
        dcopy(ny, odews->y, w);
        tnext = t + odews->dt;

        // Indicate that we haven't converged yet
        W->flag[W->rank] = 0;
        converged = 0;

        // Newton loop
        for (int k = 0; k < odews->maxits; k++)
        {
            evaluate(W, tnext, w, odews->f); // f = g(w)

            // evaluate also exchanges convergence information. If everyone
            // has converged, then we can stop. Everyone likewise gets the
            // request for Jacobian update
            if (all(W->flag, W->n_procs))
            {
                converged = 1;

                if (k > odews->nconv)
                {
                    jac_needed = 1;
                }

                break; // w contains the correct value
            }

            // Form x = w - beta - dt g (our fcn value for Newton)
            dcopy(ny, w, x);
            daxpy(ny, -1, beta, x);
            daxpy(ny, -odews->dt, odews->f, x);

            // TEST x[0] = radius etc. - w is the state var value, x is passed on to Newton
//            for (int la = 0; la < 27; la++)
//            {
//            	printf("iteration %d, state variable %2d - x: %e w: %e\n", k, la, x[la], w[la] );
//            }

            W->flag[W->rank] = sizecheck(x, ny, odews->ftol); // function value size check
            lusoln(odews, x);  // solve (x is now increment)
            W->flag[W->rank] |= sizecheck(x, ny, odews->ytol); // increment size check
            daxpy(ny, -1, x, w); // update w with new value
        }

        if (!converged)
        {
            printf("Newton iteration failed to converge\n");
            exit(1);
        }

        t = tnext;
        dcopy(ny, w, odews->y); // update y values

        if (fmod(t, odews->dtwrite) < odews->dt)
        {
            write_data(W, t, odews->y);
            write_flow(W, t, W->q, W->q0);
            write_pressure(W, t, W->p, W->p0);

        	// Time taken since beginning
        	time(&end_t);

        	// Update total time since beginning and reset start_t
        	total_t += difftime(end_t, start_t);
        	start_t = end_t;

        	// Output progress in terms of percentage remaining, time elapsed and estimated time remaining
        	if (W->rank == 0)
        	{
        		printf("Time: %4.0f  | %4d min %2d sec elapsed\n", t, (int)(total_t / 60), ((int)total_t % 60));
        	}

        }
    }

    free(w);
    free(x);
    free(beta);
}







void solver_init(ode_workspace *odews, int argc, char **argv)
{
    workspace *W;

    // Initialise the workspace
    odews->W = workspace_init(argc, argv);
    W = odews->W;
    odews->mdeclared = 0;

    // Put initial conditions in to y
    odews->y = zerosv(W->nu);
    odews->f = zerosv(W->nu);
    set_initial_conditions(W, odews->y);

    // Initial Jacobian computation
    double t0 = MPI_Wtime();
    evaluate(W, odews->t0, odews->y, odews->f);

    double tf = MPI_Wtime();
    odews->W->tfeval = tf - t0;
    t0 = MPI_Wtime();
    jacupdate(W, odews->t0, odews->y);

    double ta = MPI_Wtime();
    odews->S = newton_sparsity(W->J);

    double tb = MPI_Wtime();
    newton_matrix(odews);
    tf = MPI_Wtime();

    odews->W->tjacupdate = (tf - t0) - (tb - ta);
    odews->W->tjacfactorize = (tb - ta);
}

int sizecheck(double *x, int n, double tol) { // n - # of equations total (nblocks*nequs)
    int smallenough = 1;

    double x0[27] =
    {		1,   	// 0
		 	1e-7,	// 1
			1e-4,	// 2
			1e-3,	// 3
			1e-4,	// 4
			1e-4,	// 5
			1e-3,	// 6
			1e-5, 	// 7
			1e-4,	// 8
			1e+3,	// 9
			1e-4,	// 10
			1e-1,	// 11 *
			1.0,	// 12
			1e-1,	// 13
			1e-1, 	// 14
			1e-1,	// 15 **
			1e+5,	// 16
			1.0,	// 17
			1e-1,	// 18
			1e1,	// 19 *
			1.0,	// 20 **
			1e-1,	// 21
			1e-1,	// 22
			1e-1,	// 23
			1,
			1,
			1
	};

    for (int i = 0; i < n; i++)
    {
// 	    for (int la = 0; la < 27; la++)
// 	    {
//            printf("***** tolerance check: var = %d: %e %e  %e \n", la, x[la], x0[la % 27], fabs(x[la] / x0[la % 27])); // TEST
//        }

        smallenough &= (fabs(x[i] / x0[i % 27]) < tol);  // W->neq hardcoded
        //smallenough &= (fabs(x[i]) < tol);
        //printf("%f \n", x[i]);

        if (!smallenough)
        {
            break;
        }
    }

    return smallenough;
}

css * newton_sparsity(cs *J)
{
    // Perform symbolic analysis of the Jacobian for subsequent LU
    // factorisation
    css *S;
    int order = 0;           // 0 = natural, 1 = min deg order of A + A'
    S = cs_sqr(order, J, 0); // 0 means we're doing LU and not QR
    return S;
}

void newton_matrix(ode_workspace *odews)
{
    // Create a Newton matrix from the given step gamma and Jacobian in W
    cs *M, *eye;
    if (odews->mdeclared)
    {
        cs_nfree(odews->N);
    }
    else
    {
        odews->mdeclared = 1;
    }

    eye = speye(odews->W->J->m);
    M = cs_add(eye, odews->W->J, 1, -odews->dt);
    cs_spfree(eye);

    odews->N = cs_lu(M, odews->S, 1);
    cs_spfree(M);
}

int lusoln(ode_workspace *odews, double *b)
{
    // Can only be called if newton_matrix has been called already
    double *x;
    int n = odews->W->J->n;
    x = cs_malloc (n, sizeof (*x));
    int ok = odews->S && odews->N && x;

    if (ok)
    {
        cs_ipvec(odews->N->pinv, b, x, n);
        cs_lsolve(odews->N->L, x);
        cs_usolve(odews->N->U, x);
        cs_ipvec(odews->S->q, x, b, n);
    }

    cs_free (x);

    return ok;
}

void free_var(ode_workspace *odews)
{
    if (odews->W->buf != NULL) free(odews->W->buf);
    if (odews->W->flag != NULL) free(odews->W->flag);
    if (odews->W->level != NULL) free(odews->W->level);
    if (odews->W->l != NULL) free(odews->W->l);
    if (odews->W->x != NULL) free(odews->W->x);
    if (odews->W->y != NULL) free(odews->W->y);
    if (odews->W->b != NULL) free(odews->W->b);
    if (odews->W->u != NULL) free(odews->W->u);
    if (odews->W->v != NULL) free(odews->W->v);
    if (odews->W->p != NULL) free(odews->W->p);
    if (odews->W->q != NULL) free(odews->W->q);
    if (odews->W->w != NULL) free(odews->W->w);
    if (odews->W->ucomm != NULL) free(odews->W->ucomm);
    if (odews->W->vcomm != NULL) free(odews->W->vcomm);
    if (odews->W->gcomm != NULL) free(odews->W->gcomm);
    if (odews->W->xn != NULL) free(odews->W->xn);
    if (odews->W->xm != NULL) free(odews->W->xm);
    if (odews->W->level0 != NULL) free(odews->W->level0);
    if (odews->W->b0 != NULL) free(odews->W->b0);
    if (odews->W->p0 != NULL) free(odews->W->p0);
    if (odews->W->q0 != NULL) free(odews->W->q0);
    if (odews->W->xn0 != NULL) free(odews->W->xn0);
    if (odews->W->J != NULL) cs_spfree(odews->W->J);
    if (odews->W->dfdx->r != NULL) free(odews->W->dfdx->r);
    if (odews->W->dfdx->g != NULL) free(odews->W->dfdx->g);
    if (odews->W->dfdx->A != NULL) cs_spfree(odews->W->dfdx->A);
    if (odews->W->dfdx != NULL) free(odews->W->dfdx);
    if (odews->W->A != NULL) cs_spfree(odews->W->A);
    if (odews->W->At != NULL) cs_spfree(odews->W->At);
    if (odews->W->A_A != NULL) cs_spfree(odews->W->A_A);
    if (odews->W->A_At != NULL) cs_spfree(odews->W->A_At);
    if (odews->W->dpdgneg != NULL) cs_spfree(odews->W->dpdgneg);
    if (odews->W->dfdp->r != NULL) free(odews->W->dfdp->r);
    if (odews->W->dfdp->g != NULL) free(odews->W->dfdp->g);
    if (odews->W->dfdp->A != NULL) cs_spfree(odews->W->dfdp->A);
    if (odews->W->dfdp != NULL) free(odews->W->dfdp);
    if (odews->W->G != NULL) cs_spfree(odews->W->G);
    if (odews->W->dgdx != NULL) cs_spfree(odews->W->dgdx);
    if (odews->W->Proj != NULL) cs_spfree(odews->W->Proj);
    if (odews->W->symbchol != NULL) cs_sfree(odews->W->symbchol);
    if (odews->W->symbchol_reduced != NULL) cs_sfree(odews->W->symbchol_reduced);
    if (odews->W->nvu_w->dfdx_pattern != NULL) cs_spfree(odews->W->nvu_w->dfdx_pattern);
    if (odews->W->nvu_w->dfdp_pattern != NULL) cs_spfree(odews->W->nvu_w->dfdp_pattern);
    if (odews->W->nvu_w != NULL) free(odews->W->nvu_w);
    if (odews->W->symbchol0 != NULL) cs_sfree(odews->W->symbchol0);
    if (odews->W->G0 != NULL) cs_spfree(odews->W->G0);
    if (odews->W->A0 != NULL) cs_spfree(odews->W->A0);
    if (odews->W->A0t != NULL) cs_spfree(odews->W->A0t);

    if (odews->W != NULL) free(odews->W);
    if (odews->N != NULL) cs_nfree(odews->N);
    if (odews->y != NULL) free(odews->y);
    if (odews->f != NULL) free(odews->f);
}



