#include <stdio.h>

double b0[4] = {0.000001, 0.100001, 0.200001, 0.300001};
double b1[4] = {1.000001, 1.100001, 1.200001, 1.300001};
double b2[4] = {2.000001, 2.100001, 2.200001, 2.300001};
double b3[4] = {3.000001, 3.100001, 3.200001, 3.300001};



int main(int argc, char **argv) {

int ml = atoi(argv[1]);
int nl = atoi(argv[2]);

double N[ml * nl];
double E[ml * nl];
double S[ml * nl];
double W[ml * nl];

int i = 0;
int j = 0;

	for (j = 0; j < nl; j++) {
		for (i = 0; i < ml; i++) {

			if ((i + ml * j) < ml) {
				W[i + ml * j] = b0[i];
			} else {
			    W[i + ml * j] = i + ml * j - ml;
			}

			if (((i + ml * j + 1) % ml) == 0) {
				N[i + ml * j] = b1[j];
			} else {
			    N[i + ml * j] = i + ml * j + 1; 
			}

			if ((i + ml * j) >= (ml * (nl - 1))) {
				E[i + ml * j] = b2[i];
			} else {
			    E[i + ml * j] = i + ml * j + ml; 
			}

			if (((i + ml * j) % ml) == 0) {
				S[i + ml * j] = b3[j];
			} else {
			    S[i + ml * j] = i + ml * j - 1;
			}
		}
	
	}

	int k = 0;

	for (k = 0; k < (ml * nl); k++) {
		printf("block: %d \t W: %f \t N: %f \t E: %f \t S: %f \n", k, W[k], N[k], E[k], S[k]);
	} 

}

