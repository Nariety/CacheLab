/*
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"

void trans(int M, int N, int A[N][M], int B[M][N]);
void trans_1(int M, int N, int A[N][M], int B[M][N]);
void trans_test(int M, int N, int A[N][M], int B[M][N]);
int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N]) {

		int diagonal = 0;
		int temp = 0;
		int blockSize1 = 8, blockSize3 = 16; // For the case N ==64, blocksize =4 is the best
		// variables used to optimize the case N == 64
		// a1,a2,a3 for diagonals, a4,a5,a6 exploits spatial locality
		int a1,a2,a3,a4,a5,a6;//,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16;


		if (N == 32) {
			for (int row = 0; row < N; row += blockSize1) {
				for (int column = 0; column < N; column += blockSize1) {
					for (int i = row; i < row + blockSize1; i++) {
						for (int j = column; j < column + blockSize1; j++) {
							if (i != j) {
								B[j][i] = A[i][j];
							}
							else {
								temp = A[i][j];
								diagonal = i;
							}
						}
						if (row == column) {
							B[diagonal][diagonal] = temp;
						}
					}
				}
			}
		}
		// In order to increase the efficiency, we only used two for loops for this case
		// and used local variables to manually store and load.
		else if (N == 64) {
			for (int row = 0; row < N; row += 4) {
				for (int column = 0; column < N; column += 4) {
			/*
					a16 = A[row+3][column];
					a15 = A[row+3][column+1];
					a14 = A[row+3][column+2];
					a13 = A[row+3][column+3];

					a12 = A[row+2][column];
					a11 = A[row+2][column+1];
					a10 = A[row+2][column+2];
					a9 = A[row+2][column+3];

					a8 = A[row+1][column];
					a7 = A[row+1][column+1];
					a6 = A[row+1][column+2];
					a5 = A[row+1][column+3];

					a4 = A[row][column];
					a3 = A[row][column+1];
					a2 = A[row][column+2];
					a1 = A[row][column+3];


					B[column][row] = a4;
					B[column][row+1] = a8;
					B[column][row+2] = a12;
					B[column][row+3] = a16;

					B[column+1][row] = a3;
					B[column+1][row+1] = a7;
					B[column+1][row+2] = a11;
					B[column+1][row+3] = a15;

					B[column+2][row] = a2;
					B[column+2][row+1] = a6;
					B[column+2][row+2] = a10;
					B[column+2][row+3] = a14;

					B[column+3][row] = a1;
					B[column+3][row+1] = a5;
					B[column+3][row+2] = a9;
					B[column+3][row+3] = a13;
*/ // result in 1603 misses

					a6 = A[row+3][column];
					a5 = A[row+3][column+1];
					a4 = A[row+3][column+2];
					a3 = A[row+3][column+3];

					temp = A[row+2][column];
					B[column][row+3] = a6;
					a6 = A[row+2][column+1];// a6,a5,a4,a3 in use
					a2 = A[row+2][column+2];
					B[column+1][row+3] = a5;
					a5 = A[row+2][column+3];

					diagonal = A[row+1][column];
					a1 = A[row+1][column+1]; //a6,a5,a4,a3,a2,a1 in use
					blockSize3 = A[row+1][column+2];

					blockSize1 = A[row+1][column+3];

					B[column][row] = A[row][column];
					B[column][row+1] = diagonal;
					B[column][row+2] = temp;
					B[column+1][row] = A[row][column+1];
					B[column+1][row+1] = a1;
					B[column+1][row+2] = a6;
					B[column+2][row] = A[row][column+2];
					B[column+2][row+1] = blockSize3;
					B[column+2][row+2] = a2;
					B[column+2][row+3] = a4;
					B[column+3][row] = A[row][column+3];
					B[column+3][row+1] = blockSize1;
					B[column+3][row+2] = a5;
					B[column+3][row+3] = a3;

				}
			}
		}
		// for any case with non-symmetric matrix
		else {
			for (int row = 0; row < N; row += blockSize3) {
				for (int column = 0; column < N; column += blockSize3) {
					for (int i = row; (i < row + blockSize3) && (i < N); i++) {
						for (int j = column; (j < column + blockSize3) && (j < M); j++) {
							if (i != j) {
								B[j][i] = A[i][j];
							} else {
								temp = A[i][j];
								diagonal = i;
							}
						}
						if (row == column) {
							B[diagonal][diagonal] = temp;
						}
					}
		 		}
			}
	}
}

/*
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started.
 */

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N]) {
	int i, j, tmp;

	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}
}

/*
 * trans_1 -
 */
char trans_1_desc[] = "Simple transpose";
void trans_1(int M, int N, int A[N][M], int B[M][N]) {
	int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			B[j][i] = A[i][j];
		}
	}
}

/*
 * trans_test -
 */
char trans_test_desc[] = "Testing transpose";
void trans_test(int M, int N, int A[N][M], int B[M][N]) {

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions() {
	/* Register your solution function */
	registerTransFunction(transpose_submit, transpose_submit_desc);

	/* Register any additional transpose functions */
	registerTransFunction(trans, trans_desc);

	//registerTransFunction(trans_1,trans_1_desc);

}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N]) {
	int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < M; ++j) {
			if (A[i][j] != B[j][i]) {
				return 0;
			}
		}
	}
	return 1;
}


