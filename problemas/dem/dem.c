#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/////////////////////////////////////////////////////////////////////////////////////////
// Private Fundctions
/////////////////////////////////////////////////////////////////////////////////////////

double leapfrog(int N, double h, int ne, int ndofs, double F0[ne][2], double x0[ne], double y0[ne], int conect[ne][5], int restrs0[ne][2], double kspr, double raio, double mass) {

    
    double u[ndofs], v[ndofs], a[ndofs], fi[ndofs], res[N];
    memset(u, 0, sizeof(u));
    memset(v, 0, sizeof(v));
    memset(a, 0, sizeof(a));
    memset(fi, 0, sizeof(fi));
    memset(res, 0, sizeof(res));

    double xj, yj, xk, yk, dX, dY, di, d2, dx, dy;
    int k;
    double F[ndofs], restrs[ndofs];

    int ii = 0;
    for (int tt = 0; tt < ne; tt++) {
        for (int kk = 0; kk < 2; kk++) {
            F[ii]  = F0[tt][kk];
            restrs[ii]  = restrs0[tt][kk];
            ii = ii + 1;
        }
    }

    // a .= (F .- fi)./mass  
    for (int i = 0; i < ndofs; i++) {
        a[i] = (F[i] - fi[i])/mass;
    }

    for (int i = 0; i < N; i++) {

        // v .+= a .* (0.5*h)
        // u .+= v .* h
        for (int j = 0; j < ndofs; j++) {
            v[j] += a[j]*0.5*h;
            u[j] += v[j]*h;
            fi[j] = 0.0;
        }

        // Loop over elements
        for (int j = 1; j <= ne; j++) {
            if (restrs[2*j-2] == 1) 
                u[2*j-2] = 0.0;
            if (restrs[2*j-1] == 1) 
                u[2*j-1] = 0.0;
            xj = x0[j-1] + u[2*j-2];
            yj = y0[j-1] + u[2*j-1];
            // Check interaction
            for (int index = 1; index <= conect[j-1][0]; index++) {
                k = conect[j-1][index];
                xk = x0[k-1] + u[2*k-2];
                yk = y0[k-1] + u[2*k-1];
                dX = xj-xk;
                dY = yj-yk;
                di = sqrt(dX*dX+dY*dY);
                // printf("di = %.16f \n",di);
                d2 = (di - 2*raio);
                dx = d2*dX/di;
                dy = d2*dY/di;
                fi[2*j-2] += kspr*dx;
                fi[2*j-1] += + kspr*dy;
            }
        }
        // a .= (F .- fi)./mass 
        // v .+= a .* (0.5*h) 
        for (int j = 0; j < ndofs; j++) {
            a[j] = (F[j] - fi[j])/mass;
            v[j] = v[j] + a[j]*0.5*h;
        }
        // Plot
        res[i] = u[32];
        printf("u[32] = %.16f \n", u[32]);

    }
    
    return res[N];
}
//------------------------------------------------------------------------------
// readData(char *filename)
// {
// 	char str[STR_BUFFER_SIZE];
// 	unsigned short int mat_id;
// 	FILE *file;
// 	file = fopen(filename, "r");
// 	if (file)
// 	{
// 		while (fscanf(file, "%s", str) != EOF)
// 		{
// 			if (!strcmp(str, "%type_of_analysis"))
// 			{
// 				fscanf(file, "%i", &m_analysis_flag);
// 				if (m_analysis_flag != HMG_THERMAL)
// 				{
// 					fclose(file);
// 					return HMG_FALSE;
// 				}
// 			}
// 			else if (!strcmp(str, "%voxel_size"))
// 			{
// 				fscanf(file, "%lf", &m_elem_size);
// 			}
// 			else if (!strcmp(str, "%solver_tolerance"))
// 			{
// 				fscanf(file, "%lf", &m_tol);
// 			}
// 			else if (!strcmp(str, "%number_of_iterations"))
// 			{
// 				fscanf(file, "%i", &m_nmaxit);
// 			}
// 			else if (!strcmp(str, "%image_dimensions"))
// 			{
// 				fscanf(file, "%i %i %i", &m_nx, &m_ny, &m_nz);
// 				if (m_nz)
// 				{
// 					m_dim_flag = HMG_3D;
// 				}
// 				else
// 				{
// 					m_dim_flag = HMG_2D;
// 				}
// 			}
// 			else if (!strcmp(str, "%refinement"))
// 			{
// 				fscanf(file, "%i", &m_mesh_refinement);
// 				m_nx *= m_mesh_refinement;
// 				m_ny *= m_mesh_refinement;
// 				m_nz *= m_mesh_refinement;
// 			}
// 			else if (!strcmp(str, "%image_divisions"))
// 			{
// 				fscanf(file, "%i %i %i", &m_nSDx, &m_nSDy, &m_nSDz);
// 				if ( (m_nSDx%2 !=0) || (m_nSDy%2 !=0) || (m_nSDz%2 !=0) ){
// 					printf("\n\nAs it is, the program can not accept odd number of subdomains in any of the directions (Data race in the S-by-S computations).\n");
// 					return HMG_FALSE;
// 				}
// 				if (m_dim_flag == HMG_2D)
// 				{
// 					if (((m_nx % m_nSDx) != 0) || ((m_ny % m_nSDy) != 0))
// 					{
// 						fclose(file);
// 						printf("\nNumber of elements is not divisible by the number of subdomains!\n");
// 						return HMG_FALSE;
// 					}
// 				}
// 				if (m_dim_flag == HMG_3D)
// 				{
// 					if (((m_nSDx) == 0) || ((m_nSDy) == 0) || ((m_nSDz) == 0))
// 					{
// 						printf("\n\n Atention! Check number of subdomains!\n");
// 						sleep(2);
// 						return HMG_FALSE;
// 					}
// 					if (((m_nx % m_nSDx) != 0) || ((m_ny % m_nSDy) != 0) || ((m_nz % m_nSDz) != 0))
// 					{
// 						fclose(file);
// 						printf("\nNumber of elements is not divisible by the number of subdomains!\n");
// 						return HMG_FALSE;
// 					}
// 				}
// 			}
// 			else if (!strcmp(str, "%number_of_materials"))
// 			{
// 				fscanf(file, "%i", &m_nmat);
// 				if (m_nmat > MAX_NUM_MAT)
// 				{
// 					fclose(file);
// 					return HMG_FALSE;
// 				}
// 			}
// 			else if (!strcmp(str, "%properties_of_materials"))
// 			{
// 				for (unsigned int i = 0; i < m_nmat; i++)
// 				{
// 					if (m_analysis_flag == HMG_ELASTIC)
// 						fscanf(file, "%hi %lf %lf", &mat_id, &props[2 * i], &props[2 * i + 1]);
// 					else
// 						fscanf(file, "%hi %lf", &mat_id, &props[i]);
// 					if (mat_id >= MAX_NUM_MAT)
// 					{
// 						fclose(file);
// 						return HMG_FALSE;
// 					}
// 					props_keys[mat_id] = i;
// 				}
// 			}
// 		}
// 		fclose(file);
// 	}
// 	else
// 		return HMG_FALSE;
// 	return HMG_TRUE;
// }
// //------------------------------------------------------------------------------
// logical readMaterialMap(char *filename)
// {
// 	// Check file format before reading
// 	unsigned long int str_len = strlen(filename);
// 	if (str_len < 3)
// 		return HMG_FALSE;
// 	if (!strcmp(&filename[str_len - 3], ".nf"))
// 		return readMaterialMapNF(filename);
// 	if (!strcmp(&filename[str_len - 4], ".raw"))
// 		return readMaterialMapRAW(filename);
// 	return HMG_FALSE;
// }
// //------------------------------------------------------------------------------
// logical readMaterialMapNF(char *filename)
// {
// 	map_value buffer;
// 	unsigned int i;
// 	char str[STR_BUFFER_SIZE];
// 	FILE *file;
// 	file = fopen(filename, "r");
// 	if (file)
// 	{
// 		while (fscanf(file, "%s", str) != EOF)
// 		{
// 			if (!strcmp(str, "%ELEMENTS"))
// 				for (i = 0; i < m_nElems; i++)
// 				{
// 					fscanf(file, "%hhu", &buffer);
// 					elem_material_map[i] = props_keys[buffer];
// 				}
// 		}
// 		fclose(file);
// 	}
// 	else
// 		return HMG_FALSE;
// 	return HMG_TRUE;
// }
// //------------------------------------------------------------------------------
// logical readMaterialMapRAW(char *filename)
// {
// 	map_value buffer;
// 	unsigned int ii, jj, kk, i, j, k, SD, d_SD, c_SD, r_SD, i_SD, j_SD, k_SD, elem_SD;
// 	unsigned int slices, mesh_refinement_z;

// 	if (m_dim_flag == HMG_3D)
// 	{
// 		slices = m_nz;
// 		mesh_refinement_z = m_mesh_refinement;
// 	}
// 	else
// 	{
// 		slices = 1;
// 		m_nzSD = 1;
// 		mesh_refinement_z = 1;
// 	}

// 	char str[STR_BUFFER_SIZE];
// 	FILE *file;
// 	file = fopen(filename, "rb");
// 	if (file)
// 	{
// 		// Loops to transpose data. Raw file is line by line, our indexing is column by column
// 		for (k = 0; k < slices / mesh_refinement_z; k++)
// 		{
// 			for (i = 0; i < m_ny / m_mesh_refinement; i++)
// 			{
// 				for (j = 0; j < m_nx / m_mesh_refinement; j++)
// 				{
// 					if (fread(&buffer, sizeof(map_value), 1, file) != EOF)
// 					{
// 						for (kk = 0; kk < mesh_refinement_z; kk++)
// 						{
// 							for (jj = 0; jj < m_mesh_refinement; jj++)
// 							{
// 								for (ii = 0; ii < m_mesh_refinement; ii++)
// 								{
// 									r_SD = ((i * m_mesh_refinement) + ii) / m_nySD;
// 									c_SD = ((j * m_mesh_refinement) + jj) / m_nxSD;
// 									d_SD = ((k * m_mesh_refinement) + kk) / m_nzSD;
// 									i_SD = ((i * m_mesh_refinement) + ii) % m_nySD;
// 									j_SD = ((j * m_mesh_refinement) + jj) % m_nxSD;
// 									k_SD = ((k * m_mesh_refinement) + kk) % m_nzSD;
// 									SD = r_SD + c_SD * m_nSDy + d_SD * m_nSDx * m_nSDy;
// 									elem_SD = i_SD + j_SD * m_nySD + k_SD * m_nxSD * m_nySD;
// 									elem_material_map[SD * (m_nxSD * m_nySD * m_nzSD) + elem_SD] = props_keys[buffer];
// 								}
// 							}
// 						}
// 					}
// 					else
// 					{
// 						// If reached EOF before expected, return false
// 						fclose(file);
// 						return HMG_FALSE;
// 					}
// 				}
// 			}
// 		}
// 		fclose(file);
// 		return ;
// 	}
// 	return ;
// }
//------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////////////////////////
// Main Fundction
/////////////////////////////////////////////////////////////////////////////////////////
int main() {

    // Read Input
    // double *coords;
    // coords = (double *) malloc(sizeof(double)*N);
    // function_read_image(*coords) // Verifica onde tá o pixel, guarda o centroide do pixel como posicao da bolinha

    //Inicialization
    // Elements
    int ne = 18;
    //double x0[ne], y0[ne];
    double x0[18] = {1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 5.0, 5.0, 5.0, 7.0, 7.0, 7.0, 9.0, 9.0, 9.0, 11.0, 11.0, 11.0};
    double y0[18] = {1.0, 3.0, 5.0, 1.0, 3.0, 5.0, 1.0, 3.0, 5.0, 1.0, 3.0, 5.0, 1.0, 3.0, 5.0, 1.0, 3.0, 5.0};

    int ndofs = 2*ne;
    double raio = 1.0;
    double mass = 7850.0;
    double kspr = 210000000000.0;

    // Solver inicializations
    int N = 600;
    double h = 0.00004;

    double F[ne][2], res[N];

    // double *res = (double *) malloc(sizeof(double)*N);


    int restrs[ne][2];

    int conect[18][5] = { 
        {2,    2,    4,    0,    0},
        {3,    1,    3,    5,    0},
        {2,    2,    6,    0,    0},
        {3,    5,    1,    7,    0},
        {4,    4,    6,    2,    8},
        {3,    5,    3,    9,    0},
        {3,    8,    4,   10,    0},
        {4,    7,    9,    5,   11},
        {3,    8,    6,   12,    0},
        {3,   11,    7,   13,    0},
        {4,   10,   12,    8,   14},
        {3,   11,    9,   15,    0},
        {3,   14,   10,   16,    0},
        {4,   13,   15,   11,   17},
        {3,   14,   12,   18,    0},
        {2,   17,   13,    0,    0},
        {3,   16,   18,   14,    0},
        {2,   17,   15,    0,    0}
    };

    F[0][0] = 0.0; F[0][1] = 0.0;
    F[1][0] = 0.0; F[1][1] = 0.0;
    F[2][0] = 0.0; F[2][1] = 0.0;
    F[3][0] = 0.0; F[3][1] = 0.0;
    F[4][0] = 0.0; F[4][1] = 0.0;
    F[5][0] = 0.0; F[5][1] = 0.0;
    F[6][0] = 0.0; F[6][1] = 0.0;
    F[7][0] = 0.0; F[7][1] = 0.0;
    F[8][0] = 0.0; F[8][1] = 0.0;
    F[9][0] = 0.0; F[9][1] = 0.0;
    F[10][0] = 0.0; F[10][1] = 0.0;
    F[11][0] = 0.0; F[11][1] = 0.0;
    F[12][0] = 0.0; F[12][1] = 0.0;
    F[13][0] = 0.0; F[13][1] = 0.0;
    F[14][0] = 0.0; F[14][1] = 0.0;
    F[15][0] = -1000.0; F[15][1] = 0.0;
    F[16][0] = -1000.0; F[16][1] = 0.0;
    F[17][0] = -1000.0; F[17][1] = 0.0;

    restrs[0][0] = 1; restrs[0][1] = 1;
    restrs[1][0] = 1; restrs[1][1] = 1;
    restrs[2][0] = 1; restrs[2][1] = 1;
    restrs[3][0] = 0; restrs[3][1] = 0;
    restrs[4][0] = 0; restrs[4][1] = 0;
    restrs[5][0] = 0; restrs[5][1] = 0;
    restrs[6][0] = 0; restrs[6][1] = 0;
    restrs[7][0] = 0; restrs[7][1] = 0;
    restrs[8][0] = 0; restrs[8][1] = 0;
    restrs[9][0] = 0; restrs[9][1] = 0;
    restrs[10][0] = 0; restrs[10][1] = 0;
    restrs[11][0] = 0; restrs[11][1] = 0;
    restrs[12][0] = 0; restrs[12][1] = 0;
    restrs[13][0] = 0; restrs[13][1] = 0;
    restrs[14][0] = 0; restrs[14][1] = 0;
    restrs[15][0] = 0; restrs[15][1] = 0;
    restrs[16][0] = 0; restrs[16][1] = 0;
    restrs[17][0] = 0; restrs[17][1] = 0;

    // Call LeapFrog
    res[N] = leapfrog(N, h, ne, ndofs, F, x0, y0, conect, restrs, kspr, raio, mass);

    // write results    
    FILE *fptr;
    fptr = fopen("outputC.txt", "w");
    for (int zz=0; zz<N; zz++){
        if (fptr == NULL) {
            printf("Error opening file!\n");
            return 1;
        }
        fprintf(fptr, "%f,", res[zz]);
    }
    
    fclose(fptr);

    // printf("Hello, world!\n");
    return 0;
}