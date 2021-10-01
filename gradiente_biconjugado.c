#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define kN_COLUMNS 5
#define mat(i_row, i_col) pd_mat_a[i_row * kN_COLUMNS + i_col]

// ---------------------------------------------------------------------------------------------------------------------------------------

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

void init_mat(double *pd_mat_a) {

  mat(0, 0) = 4;
  mat(0, 1) = 0;
  mat(0, 2) = 0;
  mat(0, 3) = 0;
  mat(0, 4) = 0;

  mat(1, 0) = 1;
  mat(1, 1) = 4;
  mat(1, 2) = 0;
  mat(1, 3) = 0;
  mat(1, 4) = 0;

  mat(2, 0) = 0;
  mat(2, 1) = 1;
  mat(2, 2) = 4;
  mat(2, 3) = 0;
  mat(2, 4) = 0;

  mat(3, 0) = 0;
  mat(3, 1) = 0;
  mat(3, 2) = 1;
  mat(3, 3) = 4;
  mat(3, 4) = 0;

  mat(4, 0) = 5;
  mat(4, 1) = 0;
  mat(4, 2) = 0;
  mat(4, 3) = 1;
  mat(4, 4) = 4;
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

    print_vector_aux(i_n, pd_vector_r);
    print_vector_aux(i_n, pd_vector_r2);

    i_iteration += 1;
  }

  printf("Iteracao = %d\n", i_iteration);
  printf("x:\n");
  print_vector(i_n, pd_vector_x);
}

// ---------------------------------------------------------------------------------------------------------------------------------------