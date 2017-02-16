#include "matrix.hpp"

// g++ matrix.cpp matrix_det.cpp --std=gnu++11 -Wall

int main() {
	double arr1[] = {2, 4, 1, 1, 
					 1, 2, 3, 4, 
					 2, 1, 1, 3, 
					 4, 0, 2, 3};

	Matrix m1(4, 4, arr1);
	double det1 = m1.determinant();
	printf("determinant for my matrix:\n");
	printf("%g\n", det1);

	int N = 10;
	double *arr2 = new double[N * N];
	for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++) 
          arr2[i*N + j] = 1/(double)(i + j + 1);

	Matrix m2(N, N, arr2);
	double det = m2.determinant();
	printf("\n1/determinant G %d*%d :\n", N, N);
	printf("%g\n", 1/det);

	delete[] arr2;

	return 0;
}