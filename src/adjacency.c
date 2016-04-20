#include <stdio.h>
#include <cs.h>
#include <math.h>
#include "matops.h"
#include "brain.h"

// Structure of the H tree?
cs * adjacency(int Np) {
	/* Size of A matrix dependent on N (i.e. Np which is dependent on the number of cores and number of levels)
	Np = # levels - log2 (# cores)
	 */

    // Initialise the sparse matrix for filling
    cs *A, *T;
    int *Ti, *Tj;
    double *Tx;
    int m, n;
    m = POW_OF_2(Np-1) - 1; 	// number of rows
    n = POW_OF_2(Np) - 1;		// number of cols
    T = cs_spalloc(m, n, 3*m, 1, 1);
    Ti = T->i; Tj = T->p; Tx = T->x;

    // Set the size of the lowest level grid
    int ncols = POW_OF_2((Np-1)/2); 				// equivalent to 2^((Np-1)/2)
    int nrows = POW_OF_2((Np-1)/2 + (Np-1)%2);
    
    int a, k = 0;
    int k1, k2, row = 0, col = POW_OF_2(Np-1);
    int xbranch = 0;
    for (int L = Np - 1; L > 0; L--) {
        a = POW_OF_2(Np) - POW_OF_2(L+1);
        //b = (1 << N) - (1 << (L  ));
        //c = (1 << N) - (1 << (L-1));

        if (xbranch) {
            for (int j = 0; j < ncols; j+=2) {
                for (int i = 0; i < nrows; i++) {
                    k1 = a + i + j*nrows;
                    k2 = a + i + (j+1)*nrows;
                    Ti[k] = row; Tj[k] = k1; Tx[k++] = 1;
                    Ti[k] = row; Tj[k] = k2; Tx[k++] = 1;
                    Ti[k] = row++; Tj[k] = col++; Tx[k++] = -1;
                }
            }
            ncols /= 2;
        } 
        else {
            for (int j = 0; j < ncols; j++) {
                for (int i = 0; i < nrows; i+=2) {
                    k1 = a + i + j*nrows;
                    k2 = k1 + 1;
                    Ti[k] = row; Tj[k] = k1; Tx[k++] = 1;
                    Ti[k] = row; Tj[k] = k2; Tx[k++] = 1;
                    Ti[k] = row++; Tj[k] = col++; Tx[k++] = -1;
                }
            }
            nrows /= 2;
        }
        xbranch = !xbranch; // switch xbranch 0 <--> 1
    } // L loop: from bottom level up to the top of the tree (internal nodes)
    T->nz = k;
    A = cs_compress(T);
    cs_spfree(T);
    //sparseprint(A); // good way to print sparse matrix!
    return A;
}








