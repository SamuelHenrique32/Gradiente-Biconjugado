#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iohb.h>

#define kN_COLUMNS 5
#define kQTD_ARGS 2
#define mat(i_row, i_col) pd_mat_a[i_row * kN_COLUMNS + i_col]

// ---------------------------------------------------------------------------------------------------------------------------------------

typedef struct
{
  int i_size_ptr;
  int i_size_indexes;
  int i_size_values;
  int *pi_pointers;
  int *pi_indexes;
  double *pd_values;
} matrix_hb_t;

// ---------------------------------------------------------------------------------------------------------------------------------------

void read_matrix(char *pc_file, int *pi_m, int *pi_n, int *pi_non_zeros, int **ppi_colptr, int **ppi_rows, double **ppd_values);
void init_mat(double *pd_mat_a);
void init_vector(double *pd_vector, int i_n);
double *aloc_vector(int i_size);
void mult_mat_vector(int i_n, double *pd_mat_a, double *pd_vector, double *pd_result);
void mult_mat_row_vector(int i_n, double *pd_mat_a, double *pd_vector, double *pd_result);
double mult_vector_row_vector(int i_n, double *pd_vector1, double *pd_vector2);
void sub_vector_vector(int i_n, double *pd_vector1, double *pd_vector2, double *pd_result);
void copy_vector(int i_n, double *pd_vector1, double *pd_result);
void print_vector_aux(int i_n, double *pd_vector);
void print_vector(int i_n, double *pd_vector);
int main(int argc, char **argv);

// ---------------------------------------------------------------------------------------------------------------------------------------

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

matrix_hb_t* prepare_matrix(matrix_hb_t *ps_matrix) {

  matrix_hb_t *ps_matrix_b = malloc(sizeof(matrix_hb_t));
  *ps_matrix_b = *ps_matrix;

  ps_matrix_b->pi_pointers = calloc(sizeof(int), ps_matrix->i_size_ptr);
  ps_matrix_b->pi_indexes = calloc(sizeof(int), ps_matrix->i_size_indexes);
  ps_matrix_b->pd_values = calloc(sizeof(double), ps_matrix->i_size_indexes);

  int i_n = ps_matrix->i_size_values; 
  int i_nnz = ps_matrix->i_size_values;
  int i_n_row = ps_matrix->i_size_ptr;
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

void init_mat(double *pd_mat_a) {

  mat(0, 0) = 5;
  mat(0, 1) = 0;
  mat(0, 2) = 0;
  mat(0, 3) = 2;
  mat(0, 4) = 1;

  mat(1, 0) = 0;
  mat(1, 1) = 5;
  mat(1, 2) = 0;
  mat(1, 3) = 0;
  mat(1, 4) = 2;

  mat(2, 0) = 0;
  mat(2, 1) = 0;
  mat(2, 2) = 5;
  mat(2, 3) = 0;
  mat(2, 4) = 0;

  mat(3, 0) = 2;
  mat(3, 1) = 0;
  mat(3, 2) = 0;
  mat(3, 3) = 5;
  mat(3, 4) = 0;

  mat(4, 0) = 1;
  mat(4, 1) = 2;
  mat(4, 2) = 0;
  mat(4, 3) = 0;
  mat(4, 4) = 5;
}

void init_vector(double *pd_vector, int i_n) {

    for(int i_index=0; i_index<i_n; i_index++)  {
        pd_vector[i_index] = 1;
    }
}

double *aloc_vector(int i_size) {
  return calloc(i_size, sizeof(double));
}

void mult_mat_vector(int i_n, double *pd_mat_a, double *pd_vector, double *pd_result) {

    for(int i_index=0 ; i_index<i_n ; i_index++) {

        pd_result[i_index] = 0;

        for(int i_index2=0 ; i_index2<i_n ; i_index2++) {
            pd_result[i_index] += mat(i_index, i_index2) * pd_vector[i_index2];
        }
    }
}

void mult_mat_row_vector(int i_n, double *pd_mat_a, double *pd_vector, double *pd_result) {

  for(int i_index=0 ; i_index<i_n ; i_index++) {

    pd_result[i_index] = 0;

    for(int i_index2=0; i_index2<i_n; i_index2++) {
      pd_result[i_index] += mat(i_index2, i_index) * pd_vector[i_index2];
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

  int i_imax = 1000000;
  int i_iteration = 1;
  int i_n = kN_COLUMNS;
  int i_M = 0, i_N = 0, i_non_zeros = 0;

  double *pd_mat_a = aloc_vector(i_n * i_n); //Vector as a matrix
  init_mat(pd_mat_a);

  double *pd_vector_b = aloc_vector(i_n);
  init_vector(pd_vector_b, i_n);

  //x = zeros(n,1);
  //p = zeros(n,1);
  //p2 = zeros(n,1)
  double *pd_vector_x = aloc_vector(i_n);
  double *pd_vector_p = aloc_vector(i_n);
  double *pd_vector_p2 = aloc_vector(i_n);

  double *pd_vector_r = aloc_vector(i_n);
  double *pd_vector_aux = aloc_vector(i_n);
  double *pd_vector_v = aloc_vector(i_n);
  double *pd_vector_r2 = aloc_vector(i_n);
  double d_rho0;
  double d_beta;
  double d_alpha;
  double d_error = 0.00001;
  double d_calculated_error = 0;

  matrix_hb_t *ps_mat_csc = malloc(sizeof(matrix_hb_t));
  matrix_hb_t *ps_mat_csr = NULL;

  if(argc != kQTD_ARGS) {
    printf("%s <file>\n", argv[0]);
    exit(-1);
  }

  read_matrix(argv[1], &i_M, &i_N, &i_non_zeros, &ps_mat_csc->pi_pointers, &ps_mat_csc->pi_indexes, &ps_mat_csc->pd_values);

  ps_mat_csc->i_size_values = i_non_zeros;
  ps_mat_csc->i_size_ptr = i_M + 1;
  ps_mat_csc->i_size_indexes = i_non_zeros;

  ps_mat_csr = prepare_matrix(ps_mat_csc);

  //r = b - A*x;
  mult_mat_vector(i_n, pd_mat_a, pd_vector_x, pd_vector_aux);
  sub_vector_vector(i_n, pd_vector_b, pd_vector_aux, pd_vector_r);

  //r2 = r;
  copy_vector(i_n, pd_vector_r, pd_vector_r2);

  //rho = 1;
  double d_rho = 1;

  while(i_iteration < i_imax) {

    //rho0 = rho;
    //rho = r2' * r;
    //beta = rho / rho0
    d_rho0 = d_rho;
    d_rho = mult_vector_row_vector(i_n, pd_vector_r2, pd_vector_r);
    d_beta = d_rho / d_rho0;

    //p = r + beta*p;
    //p2 = r2 + beta*p2;
    for(int i_index=0; i_index<i_n; i_index++) {

      pd_vector_p[i_index] = pd_vector_r[i_index] + d_beta * pd_vector_p[i_index];
      pd_vector_p2[i_index] = pd_vector_r2[i_index] + d_beta * pd_vector_p2[i_index];
    }

    //v = A * p;
    mult_mat_vector(i_n, pd_mat_a, pd_vector_p, pd_vector_v);

    //alpha = rho/(p2'*v)
    d_alpha = d_rho / mult_vector_row_vector(i_n, pd_vector_p2, pd_vector_v);

    //x = x + alpha*p;
    for(int i_index=0 ; i_index<i_n ; i_index++) {
      pd_vector_x[i_index] = pd_vector_x[i_index] + d_alpha * pd_vector_p[i_index];
    }

    d_calculated_error = mult_vector_row_vector(i_n, pd_vector_r, pd_vector_r);

    if(d_calculated_error < (d_error*d_error)) {
      break;
    }

    //r = r - alpha * v;
    //r2 = r2 - alpha *A' * p2;
    mult_mat_row_vector(i_n, pd_mat_a, pd_vector_p2, pd_vector_aux);
    for(int i_index=0 ; i_index<i_n ; i_index++) {

      pd_vector_r[i_index] = pd_vector_r[i_index] - d_alpha * pd_vector_v[i_index];
      pd_vector_r2[i_index] = pd_vector_r2[i_index] - d_alpha * pd_vector_aux[i_index];
    }

    //print_vector_aux(i_n, pd_vector_r);
    //print_vector_aux(i_n, pd_vector_r2);

    i_iteration += 1;
  }

  printf("Iteracoes = %d\n\n", i_iteration);
  printf("Resposta:\n");
  print_vector(i_n, pd_vector_x);
}

// ---------------------------------------------------------------------------------------------------------------------------------------