#include <stdio.h>


int main(int argc, char **argv) {

	int mlocal = atoi(argv[1]);
	int nlocal = atoi(argv[2]);

	int neighbours[4*nlocal * mlocal];

	int b0[mlocal];
	int b1[nlocal];
	int b2[mlocal];
	int b3[nlocal];

	for (int j = 0; j < mlocal; j++) {
		b0[j] = -1 * (j + 1);
		b2[j] = -1 * (mlocal + nlocal + j + 1);
	}

	for (int k = 0; k < nlocal; k++) {
		b1[k] = -1 * (mlocal + k + 1);
		b3[k] = -1 * (mlocal + 2*nlocal + k + 1);
	}

	//int direction_offset_m = mlocal; //
	//int direction_offset_n = nlocal;

	int block_offset = 4; // neighbours per block

	int i = 0;
	int j = 0;

	for (j = 0; j < nlocal; j++) {
		for (i = 0; i < mlocal; i++) {

			if ((i + mlocal * j) < mlocal) {
				neighbours[block_offset * (i + mlocal * j)] = b0[i]; //W
			} else {
			    neighbours[block_offset * (i + mlocal * j)] = i + mlocal * j - mlocal; //W
			}

			if (((i + mlocal * j + 1) % mlocal) == 0) {
				neighbours[1 + block_offset * (i + mlocal * j)] = b1[j]; //N
			} else {
			    neighbours[1 + block_offset * (i + mlocal * j)] = i + mlocal * j + 1; //N
			}

			if ((i + mlocal * j) >= (mlocal * (nlocal - 1))) {
				neighbours[2 + block_offset * (i + mlocal * j)] = b2[i]; //E
			} else {
			    neighbours[2 + block_offset * (i + mlocal * j)] = i + mlocal * j + mlocal; //E
			}

			if (((i + mlocal * j) % mlocal) == 0) {
				neighbours[3 + block_offset * (i + mlocal * j)] = b3[j]; //S
			} else {
			    neighbours[3 + block_offset * (i + mlocal * j)] = i + mlocal * j - 1; //S
			}
		}
	}
	
	int l = 0;

	for (l = 0; l < (mlocal * nlocal); l++) {
		printf("block: %d \t W: %d \t N: %d \t E: %d \t S: %d \n", l, neighbours[block_offset * l], neighbours[1+block_offset * l], neighbours[2+block_offset * l],neighbours[3+block_offset * l]);
	} 
	
	int m = 0;
	printf("neighbours array: ");
	for (m = 0; m < (4 * mlocal * nlocal); m++) {
		printf("%d \t", neighbours[m]);
	} 
	
}

