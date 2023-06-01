#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "read_image.c"

/*
Para compilar: gcc mdf3d.c -o mdf3d -lm
Para rodar: ./mdf3d exemplos/teste_3x3
*/

// Variáveis Globais
// float* cor = NULL;
// int* connect = NULL;
// float* cc = NULL;

#define ITERATION_LIMIT 1000

double dot_product(double *m, double *n, int nn, int aux)
{   
    // printf("\nCHAMADA\n");
    // for (int ii = 0; ii < nn; ii++){
    //     printf("m[%i] = %f \n", ii, m[ii]);
    // }
    double res = 0;
    for (int i = aux; i < nn + (aux - 1); i++)
    {
        res += m[i] * n[i];
    }
    
    return res;
}

int is_approx(double *x, double *x_new, int nn)
{
    double tolerance = 0.001;

    double sum = 0;
    for (int i = 0; i < nn; i++)
    {
        sum += (x[i] - x_new[i]);
    }

    if (sum > tolerance)
    {
        return 0;
    }

    return 1;
}

int criaConnect(int nn, int connect[nn][6], float cc[nn][2], int cor[m_nmat], double dT, int nz)
{
    printf("\ncor 01 = %i e cor 02 = %i \n", cor[0], cor[1]);
    int layerAtual = nz;
    int dim = (m_nx-1) * (m_ny-1) * (m_nz-1);
    int xy = (m_nx-1) * (m_ny-1);
    int x = m_nx-1;
    bool isContorno = false;
    int indConnect = 0;
    int r;
    int slices = m_nz - 1;
	int rows = m_nx - 1;
	int cols = m_ny - 1;
    int index = 0;
    int index01 = 0;
    int material_map[dim];

    //turning material map in a vector
    for(int k = 0; k < slices; k++) {
		for(int i = 0; i < rows; i++){
			for(int j = 0; j < cols; j++){
				index = i + j * rows + k * rows * cols;
				material_map[index01] = elem_material_map[index];
                index01 += 1;
			}
			printf("\n");
		}

		if(k + 1 == slices - 1) {
			printf("\n\n");
		}
	}

    for (int e = 1; e <= dim; e++)
    {
        isContorno = false;

        // Direita
        // r = resto da divisão de (e-1) por x;
        r = (e - 1) % x;
        if (r == 0)
        {
            isContorno = true;
        }
        else
        {
            connect[indConnect][0] = e - 1;
        }

        // Baixo
        r = (e - 1) % xy;
        if (r >= xy - x)
        {
            isContorno = true;
        }
        else{
            connect[indConnect][1]= e + x;
        }

        // Cima
        r = (e - 1) % xy;
        if (r < x)
        {
            isContorno = true;
        }
        else
        {
            connect[indConnect][2] = e - x;
        }

        // Esquerda
        r = (e - 1) % x;
        if (r == (x - 1))
        {
            isContorno = true;
        }
        else
        {
            connect[indConnect][3] = e + 1;
        }

        int s = e - xy;
        if (s < 1)
        {
            isContorno = true;
        }
        else
        {
            connect[indConnect][4] = s;
        }

        int w = e + xy;
        if (w > dim)
        {
            isContorno = true;
        }
        else
        {
            connect[indConnect][5] = w;
        }

        if (isContorno)
        {
            r = (e - 1) % xy;
            if (r == 0 && e != 1)
            {
                layerAtual -= 1;
            }
            cc[e-1][0] = 1;
            cc[e-1][1] = cor[material_map[indConnect]] + floor(dT) * layerAtual;
        }
        indConnect += 1;
    }
}

void circshift(int *arr, int size, int shift)
{
    int temp[size];

    // Calculate the effective shift amount
    shift = shift % size;

    // Handling negative shifts
    if (shift < 0)
    {
        shift = size + shift;
    }

    // Copy the array into a temporary array
    for (int i = 0; i < size; i++)
    {
        temp[i] = arr[i];
    }
    // Perform the circular shift
    for (int i = 0; i < size; i++)
    {
        arr[i] = temp[(i + shift) % size];
    }
}

void criaCoords(int nn, int connect[nn][6], float cc[nn][2], float deltaT)
{
    int dim = (m_nx-1) * (m_ny-1) * (m_nz-1);
    int cor[m_nmat];
    printf("\nm_nmat = %i\n", m_nmat);
    memset(cor, 0, sizeof(cor));
    float nz;
    int count = 0;
    for (int m = 0; m < 256; m++)
    {
        if (props_keys[m] > 0)
        {
            cor[count] = m - 1;
            count += 1;
        }
    }

    int size = sizeof(cor) / sizeof(cor[0]);
    int shift = 1;

    if (cor[m_nmat-1] == 0)
    {   
        cor[0] +=1;
        circshift(cor, size, shift);
    }
    else{
        cor[1] += 1;
    }
    // To do: nz;
    nz = (m_nz) / 2;

    float dT = deltaT / nz;

    criaConnect(nn, connect, cc, cor, dT, nz);
}

void gaussseidel(int nn, double A[nn][nn], double b[nn], double x[nn]) {
    double x_new[nn];
    memset(x, 0, sizeof(x_new));
    double s1, s2;

    for(int it_count = 0; it_count < ITERATION_LIMIT; it_count++) {
        memset(x_new, 0, sizeof(x_new));

        for(int i = 0; i < nn; i++) {
            s1 = dot_product(A[i], x_new, nn, 0);

            s2 = dot_product(A[i], x, nn, 1); 
            //printf("\nb = %f  s1 = %f  s2 = %f \n", b[i], s1, s2);

            x_new[i] = (b[i] - (s1 + s2)) / A[i][i];
            printf("\n x_new[%i] = %f\n", i, x_new[i]);
        }
        if(is_approx(x, x_new, nn)){
            for (int ii = 0; ii < nn; ii++){
                x[ii] = x_new[ii];
                printf("\n x[%i] = %f ", ii, x[ii]);
            }
            break;
            }
    }
}

void outputRes(int nn, int nx, int ny, double x[nn]){
    // write results    
    FILE *fptr;
    fptr = fopen("outputImage.txt", "w");
    for (int ii=0; ii < nn; ii++){
        if (fptr == NULL) {
            printf("Error opening file!\n");
        }
        if(ii == nn - 1){
            fprintf(fptr, "%f.", x[ii]);
        }
        else{
            fprintf(fptr, "%f,", x[ii]);
        }
    }
    
    fclose(fptr);
}

// TODO: ajustar tamanho das matrizes
int mdf(int nn, int connect[nn][6], float cc[nn][2], double deltaT){
    double A[nn][nn];
    double b[nn];
    for (int ii = 0; ii < nn; ii++){
            printf(" cc[%i,0] = %f cc[%i,1] = %f \n", ii, cc[ii][0], ii, cc[ii][1]);
    }

    memset(A, 0, sizeof(A));
    memset(b, 0, sizeof(b));

    int bloco[] = {1, 1, 1, 1, 1, 1};

    printf("criando matriz A\n");

    for(int e = 0; e < nn; e++) {
        if(cc[e][0] == 0) {
            A[e][e] = -6;
            for(int i = 0; i < 6; i++) {
                int c_ei = connect[e][i];

                A[e][c_ei] = bloco[i];
            }
        }
        else {
            A[e][e] = 1;
            b[e] = cc[e][1];
        }
    }

    // int e_t; //top
    // int e_d; //down
    // int e_r; //right
    // int e_l; //left
    // int e_f; //forward
    // int e_b; //backward

    // for(int e = 0; e < nn; e++) {
    //     e_r = e + ny;
    //     e_t = e - 1;
    //     e_l = e - ny;
    //     e_d = e + 1;
    //     e_f = e - nn;
    //     e_b = e + nn;

    //     if(cc[e][0] == 0) {
    //         A[e][e] = -6;
            
    //         A[e][e_r] = 1;
    //         A[e][e_t] = 1;
    //         A[e][e_l] = 1;
    //         A[e][e_d] = 1;
    //         A[e][e_f] = 1;
    //         A[e][e_b] = 1;
    //     }
    //     else {
    //         A[e][e] = 1;
    //         b[e] = cc[e][1];
    //     }
    // }

    printf("\ncalculando gaussseidel");
    double x[nn];
    gaussseidel(nn, A, b, x); // TODO: usar ve
    printf("\nfim do calculo");
    for (int ii = 0; ii < nn; ii++){
        printf(" x[%i] = %f\n", ii, x[ii]);
    }

    float layers = (m_nz-1)/2;

    outputRes(nn, m_nx-1, m_ny-1, x);
    
    return 0;
}

void main(int argc, char** argv)
{
    printf("Construindo Modelo \n");

    char nfFilename[25];
    char rawFilename[25];

    strcat(strcpy(nfFilename, argv[1]), ".nf");
    strcat(strcpy(rawFilename, argv[1]), ".raw");

    readData(nfFilename);
    readMaterialMapRAW(rawFilename);
    // m_nx, m_ny e m_nz já podem ser usadas!
    printMaterialMap();

    int nx, ny, nz, nMat, cor;
    double deltaT = 10;
    int dim = (m_nx-1) * (m_ny-1) * (m_nz-1);
    int layers = (nz - 1) / 2;
    printf("\n x : %i \n", m_nx);

    int connect[dim][6];
    float cc[dim][2];
    memset(cc, 0, sizeof(cc));
    memset(connect, 0, sizeof(connect));

    printf("\nConstruindo Connect \n");
    criaCoords(dim, connect, cc, deltaT);
    printf("\nImprimindo Connect: \n");
    for (int ii = 0; ii < dim; ii++){
        for(int jj = 0; jj < 6; jj++){
            printf(" %i ", connect[ii][jj]);
        }
        printf(" \n ");
    }
    printf("\nRealizando análise \n");
    mdf(dim, connect, cc, deltaT);
}