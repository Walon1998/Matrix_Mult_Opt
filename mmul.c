#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>

unsigned long long readTSC() {
    // _mm_lfence();  // optionally wait for earlier insns to retire before reading the clock
    return __rdtsc();
    // _mm_lfence();  // optionally block later instructions until rdtsc retires
}

static int *randmatrix(int size) {
    int *matrix = malloc(size * size * sizeof(int));
    assert(matrix != NULL);

    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            matrix[row * size + col] = random() / (RAND_MAX / 100);
        }
    }

    return matrix;
}

static void printmatrix(int size, int *matrix, char name) {
//    printf("Matrix %c:\n", name);
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
//      printf("%10d", matrix[row * size + col]);
        }
//    printf("\n");
    }
}

static int *mmul(int size, int *A, int *B) {
    int *result = malloc(size * size * sizeof(int));
    assert(result != NULL);
    long long t = readTSC();
//    int B = 2;

//
//    for (int i = 0; i < size; i += B)
//        for (int j = 0; j < size; j += B)
//            for (int k = 0; k < size; k += B)
///* B x B mini matrix multiplications */
//                for (int i1 = i; i1 < i + B; i1++)
//                    for (int j1 = j; j1 < j + B; j1++)
//                        for (int k1 = k; k1 < k + B; k1++)
//                            result[i1 * size + j1] += a[i1 * size + k1] * b[k1 * size + j1];


    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int sum = 0;
            for (int k = 0; k < size; k++) {
                sum += A[i * size + k] * B[k * size + j];
            }
            result[i * size + j] = sum;
        }
    }

    long long resultcycles = readTSC();
    long long resultcycles2 = resultcycles - t;

    printf(" %d: %d \n", size, resultcycles2);
    return result;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("USAGE: mmul <matrix_size>\n");
        return 1;
    }
    int size = atoi(argv[1]);
    printf("SPCA Simple Matrix Multiplicator. Matrix size: %d\n", size);

    for (int i = 0; i < size; ++i) {


        int *A = randmatrix(i);
        //printmatrix(size, A, 'A');

        int *B = randmatrix(i);
        //printmatrix(size, B, 'B');


        int *C = mmul(i, A, B);


        // the result needs to be used, otherwise C will not be computed
        // redirect the output with ./mmul 5000 > output
        printmatrix(i, C, 'C');

        free(A);
        free(B);
        free(C);
    }

    return 0;
}
