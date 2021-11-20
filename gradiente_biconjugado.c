#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iohb.h>
#include <mpi.h>

#define kN_COLUMNS 5
#define kQTD_ARGS 2
#define mat(i_row, i_col) pd_mat_a[i_row * kN_COLUMNS + i_col]
#define kMAIN_PROC 0

// ---------------------------------------------------------------------------------------------------------------------------------------

typedef struct {
  //int i_size_ptr;
  //int i_size_indexes;
  //int i_size_values;
  int *pi_pointers;
  int *pi_indexes;
  double *pd_values;
} matrix_hb_t;

// ---------------------------------------------------------------------------------------------------------------------------------------

void read_matrix(char *pc_file, int *pi_m, int *pi_n, int *pi_non_zeros, int **ppi_colptr, int **ppi_rows, double **ppd_values);
matrix_hb_t* prepare_matrix(matrix_hb_t *ps_matrix, int i_non_zeros, int i_M);
void init_vector(double *pd_vector, int i_n);
double *aloc_vector(int i_size);
void mult_mat_vector(int i_n, double *pd_val, int *pi_col, int *pi_ptr, double *pd_vet, double *pd_res);
void mult_mat_row_vector(int i_n, double *pd_val, int *pd_row, int *pd_ptr, double *pd_vet, double *pd_res);
double mult_vector_row_vector(int i_n, double *pd_vector1, double *pd_vector2);
void sub_vector_vector(int i_n, double *pd_vector1, double *pd_vector2, double *pd_result);
void copy_vector(int i_n, double *pd_vector1, double *pd_result);
void print_vector_aux(int i_n, double *pd_vector);
void print_vector(int i_n, double *pd_vector);
int main(int argc, char **argv);

// ---------------------------------------------------------------------------------------------------------------------------------------

/* Parameters
   pc_file: file name
   pi_m: number of rows
   pi_n: number of columns
   pi_non_zeros: number of nonzeros in the matrix
   ppi_colptr: columns pointer
   ppi_rows: rows pointer
   ppd_values: values pointer
*/
void read_matrix(char *pc_file, int *pi_m, int *pi_n, int *pi_non_zeros, int **ppi_colptr, int **ppi_rows, double **ppd_values) {

	int i_result = 0, i_nrhs = 0;
	char *pc_type = NULL;

	i_result = readHB_info(pc_file, pi_m, pi_n, pi_non_zeros, &pc_type, &i_nrhs);

  if(i_result == 0) {
	  printf("Erro ao ler as informacoes da matriz\n");
		exit(-1);
	}

	printf("Linhas: %d \t Colunas: %d \t Nao Zeros: %d \n\n", *pi_m, *pi_n, *pi_non_zeros);

	*ppd_values = (double*) malloc(*pi_non_zeros * sizeof(double));
	*ppi_rows = (int*) malloc(*pi_non_zeros * sizeof(int));
	*ppi_colptr = (int*) malloc((*pi_n+1) * sizeof(int));

	i_result = readHB_mat_double(pc_file, *ppi_colptr, *ppi_rows, *ppd_values);

	if(i_result == 0) {
    printf("Erro ao ler os valores da matriz!\n");
    exit(-1);
  }
}

matrix_hb_t* prepare_matrix(matrix_hb_t *ps_matrix, int i_non_zeros, int i_M) {

  matrix_hb_t *ps_matrix_b = malloc(sizeof(matrix_hb_t));
  *ps_matrix_b = *ps_matrix;

  ps_matrix_b->pi_pointers = calloc(sizeof(int), i_M + 1);
  ps_matrix_b->pi_indexes = calloc(sizeof(int), i_non_zeros);
  ps_matrix_b->pd_values = calloc(sizeof(double), i_non_zeros);

  int i_n = i_non_zeros; 
  int i_nnz = i_non_zeros;
  int i_n_row = i_M + 1;
  int i_jj = 0, i_row = 0;
  int i_col = 0, i_dest = 0, i_temp = 0, i_cumsum = 0, i_last = 0;

  for(i_n=0 ; i_n<i_nnz ; i_n++) {
    ps_matrix_b->pi_pointers[ps_matrix->pi_indexes[i_n]-1]++;
  }
  
  for(i_col=0, i_cumsum=1 ; i_col<i_n_row ; i_col++) {
    i_temp = ps_matrix_b->pi_pointers[i_col];
    ps_matrix_b->pi_pointers[i_col] = i_cumsum;
    i_cumsum += i_temp;
  }

  for(i_row=0 ; i_row<i_n_row-1 ; i_row++) {
    for(i_jj=(ps_matrix->pi_pointers[i_row] - 1) ; i_jj<(ps_matrix->pi_pointers[i_row + 1] - 1) ; i_jj++) {
      i_col = ps_matrix->pi_indexes[i_jj] - 1;
      i_dest = ps_matrix_b->pi_pointers[i_col] - 1;

      ps_matrix_b->pi_indexes[i_dest] = i_row + 1;
      ps_matrix_b->pd_values[i_dest] = ps_matrix->pd_values[i_jj];

      ps_matrix_b->pi_pointers[i_col]++;
    }
  }

  for(i_col=0, i_last=1 ; i_col<=i_n_row-1 ; i_col++) {
    i_temp = ps_matrix_b->pi_pointers[i_col];
    ps_matrix_b->pi_pointers[i_col] = i_last;
    i_last = i_temp;
  }

  return ps_matrix_b;
}

void init_vector(double *pd_vector, int i_n) {

    for(int i_index=0; i_index<i_n; i_index++)  {
        pd_vector[i_index] = 1;
    }
}

double *aloc_vector(int i_size) {
  return calloc(i_size, sizeof(double));
}

void mult_mat_vector(int i_n, double *pd_val, int *pi_col, int *pi_ptr, double *pd_vet, double *pd_res) {

  int i_index = 0, i_j = 0;

  for(i_index=0; i_index<i_n ; i_index++) {

    pd_res[i_index] = 0;
    for(i_j=pi_ptr[i_index] - 1 ; i_j < pi_ptr[i_index + 1] - 1 ; i_j++) {
        pd_res[i_index] += pd_val[i_j] * pd_vet[pi_col[i_j] - 1];
    }
  }  
}

void mult_mat_row_vector(int i_n, double *pd_val, int *pi_row, int *pi_ptr, double *pd_vet, double *pd_res)
{
  int i_index = 0, i_j = 0;
  
  for(i_index=0 ; i_index<i_n ; i_index++) {

    pd_res[i_index] = 0;
    
    for(i_j=(pi_ptr[i_index]-1) ; i_j<(pi_ptr[i_index+1] - 1) ; i_j++) {
        pd_res[i_index] += pd_val[i_j] * pd_vet[pi_row[i_j] - 1];
      }
    }  
}

double mult_vector_row_vector(int i_n, double *pd_vector1, double *pd_vector2) {

  double d_result = 0;
  int i_index = 0;

  for(i_index=0 ; i_index<i_n ; i_index++) {
    d_result += pd_vector1[i_index] * pd_vector2[i_index];
  }

  return d_result;
}

void sub_vector_vector(int i_n, double *pd_vector1, double *pd_vector2, double *pd_result) {

  int i_index = 0;

  for(i_index=0 ; i_index<i_n ; i_index++) {

    pd_result[i_index] = pd_vector1[i_index] - pd_vector2[i_index];
  }
}

void copy_vector(int i_n, double *pd_vector1, double *pd_result) {

  int i_index = 0;

  for(i_index=0 ; i_index<i_n ; i_index++) {
    pd_result[i_index] = pd_vector1[i_index];
  }
}

void print_vector_aux(int i_n, double *pd_vector) {

    printf("\nVetor\t");

    for(int i_index=0 ; i_index<i_n ; i_index++) {
        printf("%f\t", pd_vector[i_index]);
    }

    printf("\n\n");
}

void print_vector(int i_n, double *pd_vector) {

	for(int i_index=0; i_index<i_n; i_index++) {
		printf("%.4f\n", pd_vector[i_index]);
	}

	printf("\n");
}

int main(int argc, char **argv) {

  //All processes has ----------------------------------------------------------------------------
  int i_imax = 1000000;
  int i_iteration = 1;
  int i_M = 0, i_N = 0, i_non_zeros = 0;
  int i_id = 0, i_n_proc = 0;

  double d_rho0;
  double d_beta;
  double d_alpha;
  double d_error = 0.00001;
  double d_calculated_error = 0;

  matrix_hb_t *ps_mat_csc = malloc(sizeof(matrix_hb_t));
  matrix_hb_t *ps_mat_csr = NULL;

  MPI_Status v_mpi_status = {0};

  //Send variables -------------------------------------------------------------------------------
  int i_M_send = 0, i_N_send = 0, i_non_zeros_send = 0;

  //CSC matrix to send
  int *pi_pointers_send;
  int *pi_indexes_send;
  double *pd_values_send;

  //Send variables -------------------------------------------------------------------------------

  //All processes has ----------------------------------------------------------------------------

  MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &i_id);
	MPI_Comm_size(MPI_COMM_WORLD, &i_n_proc);

  // Main process
  if(i_id==kMAIN_PROC) {

    if(argc != kQTD_ARGS) {
      printf("%s <file>\n", argv[0]);
      exit(-1);
    }

    read_matrix(argv[1], &i_M, &i_N, &i_non_zeros, &ps_mat_csc->pi_pointers, &ps_mat_csc->pi_indexes, &ps_mat_csc->pd_values);
  }
  //Another process
  else {
    
  }

  MPI_Bcast(&i_M, 1, MPI_INT, kMAIN_PROC, MPI_COMM_WORLD);
  MPI_Bcast(&i_N, 1, MPI_INT, kMAIN_PROC, MPI_COMM_WORLD);
  MPI_Bcast(&i_non_zeros, 1, MPI_INT, kMAIN_PROC, MPI_COMM_WORLD);

  printf("ID: %d Non zeros = %d\n", i_id, i_non_zeros);

  /*ps_mat_csr = prepare_matrix(ps_mat_csc, i_non_zeros, i_M);

  double *pd_vector_x = aloc_vector(i_N);
  double *pd_vector_p = aloc_vector(i_N);
  double *pd_vector_p2 = aloc_vector(i_N);

  double *pd_vector_r = aloc_vector(i_N);
  double *pd_vector_aux = aloc_vector(i_N);
  double *pd_vector_v = aloc_vector(i_N);
  double *pd_vector_r2 = aloc_vector(i_N);

  double *pd_vector_b = aloc_vector(i_N);
  init_vector(pd_vector_b, i_N);

  //r = b - A*x;
  mult_mat_vector(i_N, ps_mat_csr->pd_values, ps_mat_csr->pi_indexes, ps_mat_csr->pi_pointers, pd_vector_x, pd_vector_aux);
  sub_vector_vector(i_N, pd_vector_b, pd_vector_aux, pd_vector_r);

  //r2 = r;
  copy_vector(i_N, pd_vector_r, pd_vector_r2);

  //rho = 1;
  double d_rho = 1;

  while(i_iteration < i_imax) {

    //rho0 = rho;
    //rho = r2' * r;
    //beta = rho / rho0
    d_rho0 = d_rho;
    d_rho = mult_vector_row_vector(i_N, pd_vector_r2, pd_vector_r);
    d_beta = d_rho / d_rho0;

    //p = r + beta*p;
    //p2 = r2 + beta*p2;
    for(int i_index=0; i_index<i_N; i_index++) {

      pd_vector_p[i_index] = pd_vector_r[i_index] + d_beta * pd_vector_p[i_index];
      pd_vector_p2[i_index] = pd_vector_r2[i_index] + d_beta * pd_vector_p2[i_index];
    }

    //v = A * p;    
    mult_mat_vector(i_N, ps_mat_csr->pd_values, ps_mat_csr->pi_indexes, ps_mat_csr->pi_pointers, pd_vector_p, pd_vector_v);

    //alpha = rho/(p2'*v)
    d_alpha = d_rho / mult_vector_row_vector(i_N, pd_vector_p2, pd_vector_v);

    //x = x + alpha*p;
    for(int i_index=0 ; i_index<i_N ; i_index++) {
      pd_vector_x[i_index] = pd_vector_x[i_index] + d_alpha * pd_vector_p[i_index];
    }

    d_calculated_error = mult_vector_row_vector(i_N, pd_vector_r, pd_vector_r);

    if(d_calculated_error < (d_error*d_error)) {
      break;
    }

    //r = r - alpha * v;
    //r2 = r2 - alpha *A' * p2;
    mult_mat_row_vector(i_N, ps_mat_csc->pd_values, ps_mat_csc->pi_indexes, ps_mat_csc->pi_pointers, pd_vector_p2, pd_vector_aux);
    for(int i_index=0 ; i_index<i_N ; i_index++) {

      pd_vector_r[i_index] = pd_vector_r[i_index] - d_alpha * pd_vector_v[i_index];
      pd_vector_r2[i_index] = pd_vector_r2[i_index] - d_alpha * pd_vector_aux[i_index];
    }

    //print_vector_aux(i_n, pd_vector_r);
    //print_vector_aux(i_n, pd_vector_r2);

    i_iteration += 1;
  }

  printf("Iteracoes = %d\n\n", i_iteration);
  printf("Resposta:\n");
  print_vector(i_N, pd_vector_x);

  free(ps_mat_csr->pd_values);
  free(ps_mat_csr->pi_pointers);
  free(ps_mat_csr->pi_indexes);

  free(ps_mat_csc->pd_values);
  free(ps_mat_csc->pi_pointers);
  free(ps_mat_csc->pi_indexes);

  free(pd_vector_b);
  free(pd_vector_x);
  free(pd_vector_p);
  free(pd_vector_p2);
  free(pd_vector_r);
  free(pd_vector_aux);
  free(pd_vector_v);
  free(pd_vector_r2);*/

  MPI_Finalize();
}

// ---------------------------------------------------------------------------------------------------------------------------------------