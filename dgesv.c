#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <openblas/lapacke.h>
//#include <mkl_lapacke.h>

void printMatrix(double *mat, int size)
{
  for (int i=0; i<size; i++)
  {
    for (int j=0; j<size; j++)
    {
      printf("%f\t",mat[i*size+j]);
    }
    printf("\n");
  }
}

void printMatrixMatrix(double *m1, double *m2, int size)
{
  for (int i=0; i<size; i++)
  {
    for (int j=0; j<size; j++)
    {
      printf("%f\t", m1[i*size+j]);
    }
    for (int j=0; j<size; j++)
    {
      printf("%f\t",m2[i*size+j]);
    }
    printf("\n");
  }
}
void copyMatrix(double *original, double *copia, int size)
{
  for (int i=0; i<size; i++)
  {
    for (int j=0; j<size; j++)
    {
      copia[i*size+j]=original[i*size+j];
    }
  }     
}

double *generate_matrix(int size)
{
  int i;
  double *matrix = (double *) malloc(sizeof(double) * size * size);


  for (i = 0; i < size * size; i++) {
    matrix[i] = rand() % 100;
  }

  return matrix;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return abs(x - y) <= epsilon * abs(x);
  // see Knuth section 4.2.2 pages 217-218
}

int check_result(double *bref, double *b, int size)
{
  int i;

  for(i = 0; i < size*size; i++) {
    if (!is_nearly_equal(bref[i], b[i]))
      return 0;
  }

  return 1;
}

int my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{

  //Replace next line to use your own DGESV implementation
  //LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);

  //Este codigo es el mas eficiente de todos XD
  /*
  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n; j++)
    {
      if (i==j)	b[i*n+j]=1;
      else b[i*n+j]=0;
    }
    printf("\n");
  }
  */

  //A->L
  // printf ("De A a L\n");
  for (int i=0; i<n; i++) //fila que contiene al pivot
  {
    if (a[i*n+i]!=0) //Si de hecho la fila tiene un pivot
    {
      // Escalonamos la matriz ampliada (a|b)

      for (int j=i+1; j<n; j++) //fila a eliminar elementos debajo del pivot
      {
	double factor = a[j*n+i]/a[i*n+i]; //Constante a multiplicar  //Seria mas optimo haber creado la variable afuera?
	//printf("factor=%f\n",factor);
	for (int k=0; k<n; k++)
	 {
	   a[j*n+k]=a[j*n+k]-factor*a[i*n+k];
	   b[j*n+k]=b[j*n+k]-factor*b[i*n+k];
	 }
      }
      // printMatrixMatrix(a,b,n);
    }
    else
    {
      printf("Habria que buscar una permutacion que funcione!\n");
    }
  }
  //L->D
  //  printf("De L a D\n");
  for (int i=n-1; i>=0; i--) //filas pivot
  {
    if (a[i*n+i]!=0)
    {
      for (int j=0;j<i; j++) //fila a eliminar valores intrusos (?)
      {
	double factor = a[j*n+i]/a[i*n+i];
	//printf("Factor: %f\n",factor);
	for (int k=0;k<n; k++)
	{
	  a[j*n+k]=a[j*n+k]-a[i*n+k]*factor; 
          b[j*n+k]=b[j*n+k]-b[i*n+k]*factor;
	}
	//printMatrixMatrix(a,b,n);
      }
    }
    else
    {
      printf("Tenemos problemas. No se puede resolver.\n");
      //Hay que comunicar algo en la variable info
      return 1; //No vi qué valor devolver aún en este caso.
    }
  }

  //  printf("De D a I\n");
  for (int i=0; i<n; i++)
  {
    double factor =1/a[i*n+i];
    for (int j=0; j<n; j++)
    {
      
      a[i*n+j]*=factor;
      b[i*n+j]*=factor;
    }
  }
  // printMatrixMatrix(a,b,n);
  
}

void main(int argc, char *argv[])
{
  int size = atoi(argv[1]);

  double *a, *aref;
  double *b, *bref;

  srand(time(0)); //Pongo acá para que los números sean aleatorios realmente.

  
  aref = generate_matrix(size);
  bref = generate_matrix(size);
  a = generate_matrix(size); //Uso esto y el siguiente porque tiene malloc
  b = generate_matrix(size);
    
  //aref = generate_matrix(size);
  //bref = generate_matrix(size);
  copyMatrix(aref,a,size);
  copyMatrix(bref,b,size);

  
  // printf("a:\n"); printMatrix(a,size);
  //printf("aref:\n"); printMatrix(aref,size);
  //printf("b:\n"); printMatrix(b,size);
  //printf("bref:\n"); printMatrix(bref,size);
  

  
  // Using LAPACK dgesv OpenBLAS implementation to solve the system
  int n = size, nrhs = size, lda = size, ldb = size, info;
  int *ipiv = (int *) malloc(sizeof(int) * size);

  clock_t tStart = clock();
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
  printf("Time taken by OpenBLAS LAPACK: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  int *ipiv2 = (int *) malloc(sizeof(int) * size);
  
  
  tStart = clock();
  my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
  printf("Time taken by my implementation: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  //printf("a:\n"); printMatrix(a,size);
  //printf("aref:\n"); printMatrix(aref,size);
  //printf("b:\n"); printMatrix(b,size);
  //printf("bref:\n"); printMatrix(bref,size);



  if (check_result(bref, b, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");
}
