#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define N_COL 5
#define A(row, col) A[row * N_COL + col]

double *aloca_vetor(int n)
{
  return calloc(n, sizeof(double));
}

void inicializa_matriz(double *A)
{

  A(0, 0) = 4;
  A(0, 1) = 0;
  A(0, 2) = 0;
  A(0, 3) = 0;
  A(0, 4) = 0;

  A(1, 0) = 1;
  A(1, 1) = 4;
  A(1, 2) = 0;
  A(1, 3) = 0;
  A(1, 4) = 0;

  A(2, 0) = 0;
  A(2, 1) = 1;
  A(2, 2) = 4;
  A(2, 3) = 0;
  A(2, 4) = 0;

  A(3, 0) = 0;
  A(3, 1) = 0;
  A(3, 2) = 1;
  A(3, 3) = 4;
  A(3, 4) = 0;

  A(4, 0) = 5;
  A(4, 1) = 0;
  A(4, 2) = 0;
  A(4, 3) = 1;
  A(4, 4) = 4;
}

void inicializa_vetor(double *V, int n)
{
  for (int i = 0; i < n; i++)
  {
    V[i] = 1;
  }
}

void multiplicacao_matriz_vetor(int n, double *A, double *V, double *R)
{
  for (int i = 0; i < n; i++)
  {
    R[i] = 0;
    for (int j = 0; j < n; j++)
    {
      R[i] += A(i, j) * V[j];
    }
  }
}

void print_vetor(int n, double *V){
	int i;
	for (i = 0; i < n; i++){
		printf("%.4f\n", V[i]);
	}
	printf("\n");
}

void multiplicacao_matriz_linha_vetor(int n, double *A, double *V, double *R)
{
  for (int i = 0; i < n; i++)
  {
    R[i] = 0;
    for (int j = 0; j < n; j++)
    {
      R[i] += A(j, i) * V[j];
    }
  }
}

void subtracao_vetor_vetor(int n, double *V1, double *V2, double *R)
{
  int i = 0;
  for (i = 0; i < n; i++)
  {
    R[i] = V1[i] - V2[i];
  }
}

void copia_vetor(int n, double *V, double *R)
{
  int i = 0;
  for (i = 0; i < n; i++)
  {
    R[i] = V[i];
  }
}

double multiplicacao_vetor_linha_vetor(int n, double *Vl, double *V)
{
  double f = 0;
  int i = 0;
  for (i = 0; i < n; i++)
  {
    f += Vl[i] * V[i];
  }
  return f;
}

void printV(int n, double * v){
  printf("\nVetor\t");
  for(int i=0;i<n;i++){
    printf("%f\t", v[i]);
  }
  printf("\n\n");
}

int main()
{
  // imax = 1000;
  // erro = 0.0001;
  // i = 1;
  // n = length(A);
  int imax = 1000000;
  double erro = 0.00001;
  int i = 1;
  int n = N_COL;

  double *A = aloca_vetor(n * n); // Uso de vetor como matriz, usar v[ linha * N_COL + col ]
  inicializa_matriz(A);

  double *b = aloca_vetor(n);
  inicializa_vetor(b, n);

  // x = zeros(n,1);
  // p = zeros(n,1);
  // p2 = zeros(n,1)
  double *x = aloca_vetor(n);
  //printV(n, x);
  double *p = aloca_vetor(n);
  //printV(n, p);
  double *p2 = aloca_vetor(n);
  //printV(n, p2);

  double *r = aloca_vetor(n);
  double *rAux = aloca_vetor(n);
  double *vetAux = aloca_vetor(n);
  double *v = aloca_vetor(n);
  double *r2 = aloca_vetor(n);
  double rho0;
  double beta;
  double alpha;

  // r = b - A*x;
  multiplicacao_matriz_vetor(n, A, x, vetAux);
  subtracao_vetor_vetor(n, b, vetAux, r);
  //printV(n, r);

  // r2 = r;
  copia_vetor(n, r, r2);
  //printV(n, r);

  //rho = 1;
  double rho = 1;

  while (i < imax)
  {
    // rho0 = rho;
    // rho = r2' * r;
    // beta = rho / rho0
    rho0 = rho;
    //printf("rho0 = %f\n", rho0);
    rho = multiplicacao_vetor_linha_vetor(n, r2, r);
    //printf("rho = %f\n", rho);
    beta = rho / rho0;
    //printf("beta = %f\n", beta);

    //p = r + beta*p;
    //p2 = r2 + beta*p2;
    for (int i = 0; i < n; i++)
    {
      p[i] = r[i] + beta * p[i];
      p2[i] = r2[i] + beta * p2[i];
    }
    // printV(n, p);
    // printV(n, p2);

    //v = A * p;
    multiplicacao_matriz_vetor(n, A, p, v);
    //printV(n, v);

    // alpha = rho/(p2'*v)
    alpha = rho / multiplicacao_vetor_linha_vetor(n, p2, v);
    //printf("alpha = %f\n", alpha);

    //	x = x + alpha*p;
    for (int i = 0; i < n; i++)
    {
      x[i] = x[i] + alpha * p[i];
    }
    //printV(n, x);

    double error = multiplicacao_vetor_linha_vetor(n, r, r);
    // printf("Error=%f\n", error);
    if (error < erro * erro)
    {
      break;
    }

    // r = r - alpha * v;
    // r2 = r2 - alpha *A' * p2;
    multiplicacao_matriz_linha_vetor(n, A, p2, vetAux);
    for (int i = 0; i < n; i++)
    {
      r[i] = r[i] - alpha * v[i];
      r2[i] = r2[i] - alpha * vetAux[i];
    }
    printV(n, r);
    printV(n, r2);

    i = i + 1;
  }

  printf("i = %d\n", i);

  printf("x:\n");
  print_vetor(n, x);
}