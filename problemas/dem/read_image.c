
#include<stdio.h>
#include<string.h>

#define HMG_FALSE 0
#define HMG_TRUE 1
#define HMG_THERMAL 0
#define HMG_ELASTIC 3
#define STR_BUFFER_SIZE 1000
#define HMG_3D 3
#define HMG_2D 2
#define MAX_NUM_MAT 1000

int m_analysis_flag;
double m_elem_size;
double m_tol;
int m_nmaxit;
int m_nx, m_ny, m_nz;
int m_dim_flag;
int m_mesh_refinement;
int m_nmat;
short int mat_id;
double *props;
int* props_keys;
int m_nelem = 10;
int* elem_material_map;

typedef unsigned int map_value;

//------------------------------------------------------------------------------
int readData(char * filename){
	char str[STR_BUFFER_SIZE];
	unsigned short int mat_id;
	FILE * file;
	file = fopen( filename , "r");
	if (file) {
	    while (fscanf(file, "%s", str)!=EOF){
			if (!strcmp(str,"%type_of_analysis")){
				fscanf(file, "%i", &m_analysis_flag);
				if (m_analysis_flag != HMG_THERMAL && m_analysis_flag != HMG_ELASTIC){
					fclose(file);
					return HMG_FALSE;
				}
			//} else if (!strcmp(str,"%type_of_solver")){

			//} else if (!strcmp(str,"%type_of_rhs")){

			} else if (!strcmp(str,"%voxel_size")){
				fscanf(file, "%lf", &m_elem_size);
			} else if (!strcmp(str,"%solver_tolerance")){
				fscanf(file, "%lf", &m_tol);
			} else if (!strcmp(str,"%number_of_iterations")){
				fscanf(file, "%i", &m_nmaxit);
			} else if (!strcmp(str,"%image_dimensions")){
	        	fscanf(file, "%i %i %i", &m_nx, &m_ny, &m_nz);
				m_nx++;m_ny++;
				if (m_nz){
					m_dim_flag = HMG_3D;
					m_nz++;
				} else {
					m_dim_flag = HMG_2D;
				}
			} else if (!strcmp(str,"%refinement")){
				fscanf(file, "%i", &m_mesh_refinement);
			} else if (!strcmp(str,"%number_of_materials")){
				fscanf(file, "%i", &m_nmat);
				if (m_nmat > MAX_NUM_MAT){
					fclose(file);
					return HMG_FALSE;
				}
			} else if (!strcmp(str,"%properties_of_materials")){
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
            if (!strcmp(str,"%ELEMENTS"))
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

int main() {
	
	// int a = readData("teste.nf");
	// int b = readMaterialMapNF("teste.nf");
	int c = readMaterialMapRAW("teste.raw");

	// printf("%d %d %d", a, b, c);
	printf("%d", c);

	return 0;
}