#include <stdio.h>
#include <cs.h>
#include <math.h>
#include "matops.h"
#include "brain.h"

// Structure of the H tree defined by a matrix A consisting of 0, 1 and -1's
cs * adjacency(int Np)
{
	/* Size of A matrix dependent on Np: the number of levels of the local subtree
	 * which is dependent on the number of cores and number of levels)
	Np = # levels - log2 (# cores), also Np = n_blocks(local) - 1
	 */

    // Initialise the sparse matrix for filling
    cs *A, *T;
    int *Ti, *Tj;
    double *Tx;
    int m, n;
    m = POW_OF_2(Np-1) - 1; 	// number of rows of T
    n = POW_OF_2(Np) - 1;		// number of cols of T
    T = cs_spalloc(m, n, 3*m, 1, 1); // create T with size m x n
    Ti = T->i; 					// array of ints: row indices (size nzmax = max number of entries)
    Tj = T->p; 					// array of ints: column indices (size nzmax = max number of entries)
    Tx = T->x; 					// array of doubles: numerical values (size nzmax = max number of entries)

    // Set the size of the lowest level grid
    int ncols = POW_OF_2((Np-1)/2); 				// equivalent to 2^((Np-1)/2)
    int nrows = POW_OF_2((Np-1)/2 + (Np-1)%2);
    
    int a, k = 0;
    int k1, k2 = 0;
	int row = 0;
	int col = POW_OF_2(Np-1);
    int xbranch = 0;

    // L loop: from bottom level up to the top of the tree (internal nodes)
    for (int L = Np - 1; L > 0; L--)
    {
        a = POW_OF_2(Np) - POW_OF_2(L+1);
        //b = (1 << N) - (1 << (L  ));
        //c = (1 << N) - (1 << (L-1));

        if (xbranch)
        {
            for (int j = 0; j < ncols; j+=2)
            {
                for (int i = 0; i < nrows; i++)
                {
                    k1 = a + i + j*nrows;
                    k2 = a + i + (j+1)*nrows;
                    Ti[k] = row; Tj[k] = k1; Tx[k++] = 1;
                    Ti[k] = row; Tj[k] = k2; Tx[k++] = 1;
                    Ti[k] = row++; Tj[k] = col++; Tx[k++] = -1;
                }
            }
            ncols /= 2;
        } 
        else
        {
            for (int j = 0; j < ncols; j++)
            {
                for (int i = 0; i < nrows; i+=2)
                {
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
    }

    T->nz = k;
    A = cs_compress(T);
    cs_spfree(T);
    //sparseprint(A); // good way to print sparse matrix!
    return A;
}








