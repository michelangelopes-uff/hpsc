#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ITERATION_LIMIT 1000

// TODO: ajustar tamanho das matrizes
int mdf(int connect[5][5], double cc[5][5], int nx, int ny, int nz){
    int nn = nx * ny * nz;
    double A[nn][nn];
    double b[nn];

    memset(A, 0, sizeof(A));
    memset(b, 0, sizeof(b));

    int bloco[] = {1, 1, 1, 1, 1, 1};

    printf("criando matriz A\n");

    int e_t; //top
    int e_d; //down
    int e_r; //right
    int e_l; //left
    int e_f; //forward
    int e_b; //backward

    for(int e = 0; e < nn; e++) {
        e_r = e + ny;
        e_t = e - 1;
        e_l = e - ny;
        e_d = e + 1;
        e_f = e - nn;
        e_b = e + nn;

        if(cc[e][0] == 0) {
            A[e][e] = -6;
            
            A[e][e_r] = 1;
            A[e][e_t] = 1;
            A[e][e_l] = 1;
            A[e][e_d] = 1;
            A[e][e_f] = 1;
            A[e][e_b] = 1;
        }
        else {
            A[e][e] = 1;
            b[e] = cc[e][1];
        }
    }

    printf("calculando gaussseidel");
    double x[nn];
    gaussseidel(A, b, x, nn); // TODO: usar ve
    printf("fim do calculo");

    
    return 0;
}

void gaussseidel(double** A, double* b, double* x, int nn) {
    memset(x, 0, sizeof(x));

    double x_new[sizeof(x)];
    double s1, s2;

    for(int it_count = 0; it_count < ITERATION_LIMIT; it_count++) {
        memset(x_new, 0, sizeof(x_new));

        for(int i = 0; i < nn; i++) {
            s1 = dot_product(A[i], x_new, nn, 0);    
            s2 = dot_product(A[i], x_new, nn, 1); 

            x_new[i] = (b[i] - (s1 + s2)) / A[i][i];
        }

        if(is_approx(x, x_new, nn)){
            break;
        }

        x = x_new;
    }
}

int is_approx(double *x, double* x_new, int nn) {
    double tolerance = 0.001;

    double sum = 0;
    for(int i = 0; i < nn; i++) {
        sum += (x[i] - x_new[i]);
    }

    if(sum > tolerance) {
        return 0;
    }

    return 1;
}

double dot_product(double* m, double* n, int nn, int aux){
    double res = 0;

    for(int i = aux; i < nn + (aux-1); i++) {
        res += m[i] * n[i];
    }
    return res;
}

void main() {
    int nx, ny, nz;
    int nel = nx * ny * nz;


    // mdf()
}

