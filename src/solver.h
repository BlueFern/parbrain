#ifndef SOLVER_H_
#define SOLVER_H_

#include "brain.h"

typedef struct ode_workspace {
    workspace *W;
    csn *N; // Newton matrix numeric factorisation
    css *S; // Newton matrix sybolic factorisation
    double *y; // Workspace variable
    double *p; //
    double *q; //
    double *f; // Workspace variable
    double dt;
    double t0;
    double tf;
    double ftol;
    double ytol;
    int    maxits;
    int    nconv;
    int mdeclared;
    double dtwrite;
} ode_workspace;

// prototypes
void newton_matrix(ode_workspace *odews);
int lusoln(ode_workspace *odews, double *b);
css * newton_sparsity(cs *J);
int sizecheck(double *x, int n, double tol);
void back_euler(ode_workspace *odews);
void solver_init(ode_workspace *odews, int argc, char **argv);
void free_var(ode_workspace *odews); //free remaining variables

#endif
