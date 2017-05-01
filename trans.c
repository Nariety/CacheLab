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

	// when N == 32
	if (N == 32) {
		// for each row
		for (int row = 0; row < N; row += blockSize1) {
			// for each column
			for (int column = 0; column < N; column += blockSize1) {
				// for each row of the block
				for (int i = row; i < row + blockSize1; i++) {
					// for each column of the block
					for (int j = column; j < column + blockSize1; j++) {
						// when i != j, stores A[i][j] into B[j][i]
						if (i != j) {
							B[j][i] = A[i][j];
						}
						// else stores A[i][j] into temp and i into diagonal for later use
						else {
							temp = A[i][j];
							diagonal = i;
						}
					}
					// when row == column
					// stores temp into B[diagonal[diagonal] to avoid some misses
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
		// local variables used to optimize the case N == 64
		// a1,a2,a3,a4,a5,a6 to exploits spatial locality
		int a1, a2, a3, a4, a5, a6;
		for (int row = 0; row < N; row += 4) {
			for (int column = 0; column < N; column += 4) {
				//for each 4x4 block in matrix A

				//stores the first row of A into a3, a4, a5, a6
				a6 = A[row][column];
				a5 = A[row][column + 1];
				a4 = A[row][column + 2];
				a3 = A[row][column + 3];

				//stores the second row of A into a2, a1, temp, diagonal
				a2 = A[row + 1][column];
				a1 = A[row + 1][column + 1];
				temp = A[row + 1][column + 2];
				diagonal = A[row + 1][column + 3];

				//stores the first two elements in the third row of A into blockSize1, blockSize2
				blockSize1 = A[row + 2][column];
				blockSize3 = A[row + 2][column + 1];

				//storing the third row of B
				B[column + 2][row + 2] = A[row + 2][column + 2];
				B[column + 2][row] = a4;
				B[column + 2][row + 1] = temp;
				B[column + 2][row + 3] = A[row + 3][column + 2];

				//stores the first two elements of A into a4, temp
				temp = A[row + 3][column];
				a4 = A[row + 3][column + 1];

				//storing the last row of B
				B[column + 3][row + 3] = A[row + 3][column + 3];
				B[column + 3][row + 2] = A[row + 2][column + 3];
				B[column + 3][row + 1] = diagonal;
				B[column + 3][row] = a3;

				//storing the second row of B
				B[column + 1][row + 3] = a4;
				B[column + 1][row + 2] = blockSize3;
				B[column + 1][row + 1] = a1;
				B[column + 1][row] = a5;

				//storing the first row of B
				B[column][row + 3] = temp;
				B[column][row + 1] = a2;
				B[column][row + 2] = blockSize1;
				B[column][row] = a6;
			}
		}
	}
	// for any case with non-symmetric matrix
	else {
		//for each row
		for (int row = 0; row < N; row += blockSize3) {
			//for each column
			for (int column = 0; column < N; column += blockSize3) {
				//for each row of the block
				for (int i = row; (i < row + blockSize3) && (i < N); i++) {
					//for each column of the block
					for (int j = column; (j < column + blockSize3) && (j < M);
							j++) {
						// when i != j, stores A[i][j] into B[j][i]
						if (i != j) {
							B[j][i] = A[i][j];
						}
						// else stores A[i][j] into temp and i into diagonal for later use
						else {
							temp = A[i][j];
							diagonal = i;
						}
					}
					// when row == column
					// stores temp into B[diagonal[diagonal] to avoid some misses
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

