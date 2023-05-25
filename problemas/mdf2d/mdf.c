#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define N 16

void solveMatrixSystem(double A[N][N], double b[N], double x[N]) {
    int i, j, k;
    float factor, sum;

    // Forward elimination
    for (k = 0; k < N - 1; k++) {
        for (i = k + 1; i < N; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k; j < N; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Backward substitution
    for (i = N - 1; i >= 0; i--) {
        sum = 0.0;
        for (j = i + 1; j < N; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
    for (int i=0; i < N; i++){
        printf("%f \n", x[i]);
    }
    printf("---");
}

int main(){
    printf("MDF!\n");

    int conect[16][4] = {
        {2,   0,   0,   5},
        {3,   1,   0,   6},
        {4,   2,   0,   7},
        {0,   3,   0,   8},
        {6,   0,   1,   9},
        {7,   5,   2,   10},
        {8,   6,   3,   11},
        {0,   7,   4,   12},
        {10,  0,   5,   13},
        {11,  9,   6,   14},
        {12,  10,  7,   15},
        {0,   11,  8,   16},
        {14,  0,   9,   0},
        {15,  13,  10,  0},
        {16,  14,  11,  0},
        {0,   15,  12,  0}
    };

    int cc[16][2] = {
        {1,   100},
        {1,   75},
        {1,   75},
        {1,   0},
        {1,   100},
        {0,   0},
        {0,   0},
        {1,   0},
        {1,   100},
        {0,   0},
        {0,   0},
        {1,   0},
        {1,   100},
        {1,   25},
        {1,   25},
        {1,   0}
    };

    int bloco[5] = {4, -1, -1, -1, -1};

    int n = 16;
    double A[16][16];
    memset(A, 0, sizeof(A));
    double b[16];
    memset(b, 0, sizeof(b));

    //Assembly
    for (int i=0; i<n; i++){
        A[i][i] = bloco[0];
        for (int j=0; j<4; j++){
            int col = conect[i][j];
            if (col != 0){
                if (cc[col,1] == 0){
                    A[i][col] = bloco[j+1];
                }
                else {
                    b[i] -= bloco[j+1]*cc[col][1];
                }
            }
        }
    }

    // bc
    for (int i=0; i<n; i++){
        if (cc[i][0] == 1){
            memset(A[i], 0, sizeof(double));
            A[i][i] = 1;
            b[i] = cc[i][1];
        }
    }
    // for (int i=0; i < n; i++){
    //     printf("%f \n", b[i]);
    // }

    double x[16];
    solveMatrixSystem(A, b, x);
    
}