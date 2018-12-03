#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
#include <emmintrin.h>

#define SM (CLS / sizeof (int))


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

static void simple(int *A, int *B, int *result, int size) {
    //Simple
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int sum = 0;
            for (int k = 0; k < size; k++) {
                sum += A[i * size + k] * B[k * size + j];
            }
            result[i * size + j] = sum;
        }
    }
}

static void blocktiling(int *A, int *B, int *result, int size) {
    int block = 2;
    //Block tiling
    for (int i = 0; i < size; i += block) {
        for (int j = 0; j < size; j += block) {
            for (int k = 0; k < size; k += block) {
                /* block x block mini matrix multiplications */
                for (int i1 = i; i1 < i + block; i1++) {
                    for (int j1 = j; j1 < j + block; j1++) {

                        for (int k1 = k; k1 < k + block; k1++) {
                            result[i1 * size + j1] += A[i1 * size + k1] * B[k1 * size + j1];
                        }
                    }
                }
            }
        }
    }
}

static void blocktilingunrolled(int *A, int *B, int *result, int size) {
    int block = 2;
//    Unrolled block tiling
    for (int i = 0; i < size; i += block) {
        for (int j = 0; j < size; j += block) {
            for (int k = 0; k < size; k += block) {
                /* block x block mini matrix multiplications */
                result[i * size + j] = A[i * size + k] * B[k * size + j] + A[i * size + k + 1] * B[(k + 1) * size + j]
                                       + result[i * size + j];
                result[(i + 1) * size + j] =
                        A[(i + 1) * size + k] * B[k * size + j] + A[(i + 1) * size + k + 1] * B[(k + 1) * size + j]
                        + result[(i + 1) * size + j];
                result[i * size + (j + 1)] =
                        A[i * size + k] * B[k * size + (j + 1)] + A[i * size + k + 1] * B[(k + 1) * size + (j + 1)]
                        + result[i * size + (j + 1)];
                result[(i + 1) * size + (j + 1)] = A[(i + 1) * size + k] * B[k * size + (j + 1)]
                                                   + A[(i + 1) * size + k + 1] * B[(k + 1) * size + (j + 1)] +
                                                   result[(i + 1) * size + (j + 1)];
            }
        }
    }
}

static void transposesimple(int *A, int *B, int *result, int size) {
    int *tmp = malloc(size * size * sizeof(int));


    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            tmp[i * size + j] = B[j * size + i];
        }
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            int sum = 0;
            for (int k = 0; k < size; ++k) {
                sum += A[i * size + k] * tmp[j * size + k];

            }
            result[i * size + j] += sum;
        }
    }

}

static void blockvectorize(int *A, int *B, int *result, int size) {
//    Not working

    int i, i2, j, j2, k, k2;
    double *restrict rres;
    double *restrict rmul1;
    double *restrict rmul2;
    for (i = 0; i < size; i += SM)
        for (j = 0; j < size; j += SM)
            for (k = 0; k < size; k += SM)
                for (i2 = 0, rres = &result[i * sizeof +j], rmul1 = &A[i * size + k]; i2 < SM;
                     ++i2, rres += size, rmul1 += size) {
                    _mm_prefetch (&rmul1[8], _MM_HINT_NTA);
                    for (k2 = 0, rmul2 = &B[k * size + j]; k2 < SM; ++k2, rmul2 += size) {
                        __m128d m1d = _mm_load_sd(&rmul1[k2]);
                        m1d = _mm_unpacklo_pd(m1d, m1d);
                        for (j2 = 0; j2 < SM; j2 += 2) {
                            __m128d m2 = _mm_load_pd(&rmul2[j2]);
                            __m128d r2 = _mm_load_pd(&rres[j2]);
                            _mm_store_pd(&rres[j2],
                                         _mm_add_pd(_mm_mul_pd(m2, m1d), r2));
                        }
                    }
                }
}


static int *mmul(int size, int *A, int *B) {
    int *result = malloc(size * size * sizeof(int));
    assert(result != NULL);

    long long t = readTSC();

    simple(A, B, result, size);


    double resultcycles = readTSC();
    double resultcycles2 = resultcycles - t;

    if (resultcycles2 < 0) {
        printf("0\n");
    } else {
        printf("%g\n", resultcycles2);
//        printf("%d\n", resultcycles2);
    }


    return result;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("USAGE: mmul <matrix_size>\n");
        return 1;
    }
    int size = atoi(argv[1]);
    printf("SPCA Simple Matrix Multiplicator. Matrix size: %d\n", size);

//    for (int i = 0; i < size; i += 1) {


    int *A = randmatrix(size);
    //printmatrix(size, A, 'A');

    int *B = randmatrix(size);
    //printmatrix(size, B, 'B');


    int *C = mmul(size, A, B);


    // the result needs to be used, otherwise C will not be computed
    // redirect the output with ./mmul 5000 > output
    printmatrix(size, C, 'C');

    free(A);
    free(B);
    free(C);
//    }

    return 0;
}
