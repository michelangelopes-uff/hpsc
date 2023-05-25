#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HMG_FALSE 0
#define HMG_TRUE 1
#define HMG_2D 2
#define HMG_3D 3
#define HMG_THERMAL 0
#define HMG_ELASTIC 3 // VALOR ALEATÓRIO ESCOLHIDO
#define STR_BUFFER_SIZE 1000 // VALOR ALEATÓRIO ESCOLHIDO
#define MAX_NUM_MAT 1000 // VALOR ALEATÓRIO ESCOLHIDO

double m_elem_size = -1; // tamanho de um elemento
double m_tol = -1; // valor da tolerância da solução
int m_analysis_flag = -1; // flag do tipo de análise a ser realizada: HMG_THERMAL ou HMG_ELASTIC, sendo HMG = homogeneity
int m_dim_flag = -1; // flag da dimensão da malha: 3d ou 2d
int m_nelem = -1; // número de elementos
int m_nmat = -1; // número de materiais
int m_nmaxit = -1; // número máximo de iterações
int m_nx = -1, m_ny = -1, m_nz = -1; // valores em cada dimensão: x, y e z (3d); ou x e y (2d)
int m_mesh_refinement = -1; // valor do refinamento da malha
short int mat_id = -1; // identificador do material

double* props = NULL; // array relacionando o índice ordenado do material com o valor referente à propriedade do material
int* props_keys = NULL; // array relacionando o identificador do material com o índice ordenado do material
int* elem_material_map = NULL; // array relacionando o elemento com o identificador do seu material

typedef unsigned char map_value;

//------------------------------------------------------------------------------
int readData(char * filename){
	char str[STR_BUFFER_SIZE];
	unsigned short int mat_id;
	FILE * file;
	file = fopen( filename , "r");
	if (file) {
	    while (fscanf(file, "%s", str)!=EOF){
			if (!strcmp(str,"\%type_of_analysis")){
				fscanf(file, "%i", &m_analysis_flag);
				if (m_analysis_flag != HMG_THERMAL && m_analysis_flag != HMG_ELASTIC){
					fclose(file);
					return HMG_FALSE;
				}
			//} else if (!strcmp(str,"%type_of_solver")){

			//} else if (!strcmp(str,"%type_of_rhs")){

			} else if (!strcmp(str,"\%voxel_size")){
				fscanf(file, "%lf", &m_elem_size);
			} else if (!strcmp(str,"\%solver_tolerance")){
				fscanf(file, "%lf", &m_tol);
			} else if (!strcmp(str,"\%number_of_iterations")){
				fscanf(file, "%i", &m_nmaxit);
			} else if (!strcmp(str,"\%image_dimensions")){
	        	fscanf(file, "%i %i %i", &m_nx, &m_ny, &m_nz);
				if (m_nz){
					m_nelem = m_nx * m_ny * m_nz;
					m_dim_flag = HMG_3D;
					m_nz++;
				} else {
					m_nelem = m_nx * m_ny;
					m_dim_flag = HMG_2D;
				}
				m_nx++;m_ny++;

			} else if (!strcmp(str,"\%refinement")){
				fscanf(file, "%i", &m_mesh_refinement);
			} else if (!strcmp(str,"\%number_of_materials")){
				fscanf(file, "%i", &m_nmat);
				if (m_nmat > MAX_NUM_MAT){
					fclose(file);
					return HMG_FALSE;
				}
			} else if (!strcmp(str,"\%properties_of_materials")){
				props_keys = (int*) calloc(MAX_NUM_MAT, sizeof(int));
				props = (double*) calloc(m_nmat, sizeof(double));

				for (unsigned int i=0; i<m_nmat; i++){
					if (m_analysis_flag == HMG_ELASTIC)
						fscanf(file, "%hi %lf %lf", &mat_id, &props[2*i], &props[2*i+1]);
					else
						fscanf(file, "%hi %lf", &mat_id, &props[i]);
					if (mat_id >= MAX_NUM_MAT){
						fclose(file);
						return HMG_FALSE;
					}
					props_keys[mat_id] = i;
				}
			//} else if (!strcmp(str,"%volume_fraction")){

			//} else if (!strcmp(str,"%data_type")){

			}
		}
	    fclose(file);
	} else
		return HMG_FALSE;
	return HMG_TRUE;
}
//------------------------------------------------------------------------------
int readMaterialMapNF(char * filename){
	map_value buffer;
	unsigned int i;
	char str[STR_BUFFER_SIZE];
	FILE * file;
	file = fopen(filename,"r");
	if (file) {
        while (fscanf(file, "%s", str)!=EOF){
            if (!strcmp(str,"\%ELEMENTS"))
				if(elem_material_map == NULL) {
					elem_material_map = (int*) calloc(m_nelem, sizeof(int));
				}

                for (i=0;i<m_nelem;i++){
                    fscanf(file, "%hhu", &buffer);
					elem_material_map[i] = props_keys[buffer];
				}
        }
        fclose(file);
	} else
		return HMG_FALSE;
	return HMG_TRUE;
}
//------------------------------------------------------------------------------
int readMaterialMapRAW(char * filename){
	map_value buffer;
	unsigned int i, j, k;
	unsigned int rows = m_ny-1;
	unsigned int cols = m_nx-1;
	unsigned int slices;
	if (m_dim_flag == HMG_3D)
		slices = m_nz-1;
	else
		slices = 1;
	char str[STR_BUFFER_SIZE];
	FILE * file;
	file = fopen(filename,"rb");
	if (file) {
		// Loops to transpose data. Raw file is line by line, our indexing is
		// column by column
		if(elem_material_map == NULL) {
			elem_material_map = (int*) calloc(m_nelem, sizeof(int));
		}
		for (k = 0; k<slices; k++){
			for (i = 0; i<rows; i++){
				for (j = 0; j<cols; j++){
					if (fread(&buffer,sizeof(map_value),1,file)!=EOF){
						elem_material_map[i+j*rows+k*rows*cols] = props_keys[buffer];
					} else {
						// If reached EOF before expected, return false
						fclose(file);
						return HMG_FALSE;
					}
				}
			}
		}
        fclose(file);
		return HMG_TRUE;
	}
	return HMG_FALSE;
}
//------------------------------------------------------------------------------
int readMaterialMap(char * filename){
	// Check file format before reading
	unsigned long int str_len = strlen(filename);
	if (str_len < 3)
		return HMG_FALSE;
	if (!strcmp(&filename[str_len-3],".nf"))
		return readMaterialMapNF(filename);
	if (!strcmp(&filename[str_len-4],".raw"))
		return readMaterialMapRAW(filename);

	printf("\n\n HERE \n\n");
	return HMG_FALSE;
}
//------------------------------------------------------------------------------
void printMaterialMap() {
	int slices = m_nz - 1;
	int rows = m_nx - 1;
	int cols = m_ny - 1;
	int index = -1;

	for(int k = 0; k < slices; k++) {
		for(int i = 0; i < rows; i++){
			for(int j = 0; j < cols; j++){
				index = i + j * rows + k * rows * cols;
				printf("%d\t", elem_material_map[index]);
			}
			printf("\n");
		}

		if(k + 1 == slices - 1) {
			printf("\n\n");
		}
	}
}
//------------------------------------------------------------------------------
int main() {
	
	int a = readData("teste.nf");
	// int b = readMaterialMapNF("teste.nf");
	int c = readMaterialMapRAW("teste.raw");

	printMaterialMap();

	return 0;
}