#include<stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define RK_ACCURACY 1e-08
#define GSL_ODEIV_STEP_TYPE gsl_odeiv_step_rkf45
//#define GSL_ODEIV_STEP_TYPE gsl_odeiv_step_rk8pd

#define INV_SQRT_TWO_PI 0.39894228040143270286 /* = 1.0/sqrt(2*pi) */

/* from lapack, to compute inverse matricies. */
//extern void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
//extern void dgetrs_(char *trans, int *n, int *nrhs, double *a, int *lda, 
//                    int *ipiv, double *b, int *ldb, int *info);
extern void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
extern void dpotri_(char *uplo, int *n, double *a, int *lda, int *info);

static int dim = 0; 
static int dim2 = 0; 
static int nfacet = 0; /* the number of facets*/
static int nfacet2 = 0; /* nfacet^2 */
static int rank = 0; /* rank = number of the non-empty faces */
static double *a = NULL; 
static double *b = NULL;
static double *a0 = NULL; 
static double *b0 = NULL;
static double *a1 = NULL; 
static double *b1 = NULL;
static double *dadt = NULL; 
static double *dbdt = NULL;
static double *alpha = NULL;
static double *inv_alpha = NULL;
static double *det_alpha = NULL;
static double *dgdb = NULL;

static void simplex(char **);
static void simplex_init(double *g);
static void simplex_rk(double *g);
static int simplex_function(double t, const double g[], double dg[], void *params);
static void simplex_update(double t);

static void simplex2(char **);
static void simplex2_rk(double *g);
static int simplex2_function(double t, const double g[], double dg[], void *params);
static void simplex2_update(double t);

static void simplicial_cone(char **);
static void simplicial_cone_init(double *g);
static void simplicial_cone_rk(double *g);
static int simplicial_cone_function(double t, const double g[], double dg[], void *params);
static void simplicial_cone_update(double t);

static void tukey(char **);
static void tukey_rk(double *g);
static int tukey_function(double t, const double g[], double dg[], void *params);
static void tukey2(char **);
static void tukey2_rk(double *g);
static int tukey2_function(double t, const double g[], double dg[], void *params);
static void tukey3(char **);
static void tukey3_rk(double *g);
static int tukey3_function(double t, const double g[], double dg[], void *params);

static void update_alpha(void);
static void update_inv_alpha(void);
static int inv_submat_alpha(int i, double *inv_submat);
static int inv_submat_alpha2(int i,int j, double *inv_submat, double *det_submat);
static int cardinality(int J, const int max);
static double get_del_aij(const double g[], int J, int i, int j);
static double get_del_bibj(const double g[], int J, int i, int j);
static double get_del_bj(const double g[], int J, int j);
static void get_prob(double *g);


static void print_vector(FILE *, int, double *, char *);
static void print_matrix(FILE *fp, int nrow, int ncol, double *matrix,const char *str);

static void numerical_experiment1(char **);
static void numerical_experiment2(char **);
static void numerical_experiment3(char **);
static void numerical_experiment4(char **);

static void 
simplex(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  printf("d=%d\n", d);

  int n = d+1;
  double a_loc[d*n], b_loc[n];
  int i,j;
  for ( i = 0; i < n; i++){
    for ( j = 0; j < d; j++)
      sscanf(*++p, "%lf", a_loc+j+i*d);
    sscanf(*++p, "%lf", b_loc+i);
  }
  printf("a,b:\n");
  for ( i = 0; i < n; i++){
    for ( j = 0; j < d; j++)
      printf("%10.6f ", a_loc[j+i*d]);
    printf("%10.6f\n", b_loc[i]);
  }

  /* Setting global variables */
  dim = d;
  nfacet = dim + 1;
  rank = (1 << nfacet) -1;

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[rank];
  //double g[(1<<nfacet)];
  //int size_a = dim2;
  int size_b = nfacet;
  int size_alpha=nfacet2;
  int size_inv_alpha=(1<<nfacet)*nfacet2;
  int size_det_alpha=(1<<nfacet);
  int size_dgdb=rank*dim;
  double work_space[size_b 
                   +size_alpha 
                   +size_inv_alpha 
                   +size_det_alpha
                   +size_dgdb];
  double *w = work_space;
  a = a_loc;           
  b = w;          w += size_b;
  a1 = a_loc;
  b1 = b_loc;
  dadt = NULL;
  dbdt = b_loc;
  alpha = w;      w += size_alpha;
  inv_alpha = w;  w += size_inv_alpha;
  det_alpha = w;  w += size_det_alpha;
  dgdb = w;  w += size_dgdb;

  printf("dim = %d\n", dim);
  print_matrix(stdout, dim, nfacet, a1, "The value of matrix a:");
  print_vector(stdout, nfacet, b1, "The value of vector b:\n");
  printf("rank = %d\n", rank);

  update_alpha();
  print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  update_inv_alpha();
  int J;
  printf("determinant: ");
  for(J = 0; J < rank; J++){
    //char s[10];
    //sprintf(s, "%d:", J);
    //print_matrix(stdout, nfacet, nfacet, inv_alpha+J*nfacet2, s);
    printf("%g ", det_alpha[J]);
  }
  printf("\n");

  simplex_init(g);
  print_vector(stdout, rank+1, g, "The initial value of g:\n");

  simplex_rk(g);
  print_vector(stdout, rank, g, "output of runge-kutta:\n");

  get_prob(g);
  print_vector(stdout, rank, g, "result:\n");

  printf("Probability = %g\n", g[0]);
  return;
}

static void 
simplex_init(double *g)
{
  int J;
  for (J = 0; J < rank; J++)
    g[J] = 0.0;

  int i;
  int m = 0;
    
  for (i = 0; i < nfacet; i++)
    m = m | (1<<i);
    
  for (i = 0; i < nfacet; i++){
    J = ~(1<<i)&m;
    g[J] = 1.0 / det_alpha[J];
  }
  return;
}

static void 
simplex_rk(double *g)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = GSL_ODEIV_STEP_TYPE;
  gsl_odeiv_step *s  = gsl_odeiv_step_alloc (T, rank);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (rank);
  gsl_odeiv_system sys = {simplex_function, NULL, rank, NULL};

  double t = 0.0;
  while (t < 1.0){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, 1.0, &h, g);
    if (status != GSL_SUCCESS)
      break;
#ifdef _DEBUG
    //fprintf(stderr, "%lf %lf\n",t, g[0]);
    //print_vector(stdout, rank+1, g, "\n");
#endif
  }
  /* free */
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return;
}

static int 
simplex_function(double t, const double g[], double dg[], void *params)
{
  simplex_update(t);
  int J, j;
  for (J = 0; J < rank; J++){
    dg[J] = 0.0;
    for (j = 0; j < nfacet; j++)
      dg[J] += dbdt[j] * get_del_bj(g, J, j);
  }
  return GSL_SUCCESS;
}

static void 
simplex_update(double t)
{
  int i;
  for ( i = 0; i < nfacet ; i++)
    b[i] = t * b1[i];
  return;
}


static void 
simplex2(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  printf("d=%d\n", d);

  int n = d+1;
  double a_loc[d*n], b_loc[n];
  int i,j;
  for ( i = 0; i < n; i++){
    for ( j = 0; j < d; j++)
      sscanf(*++p, "%lf", a_loc+j+i*d);
    sscanf(*++p, "%lf", b_loc+i);
  }
  printf("a,b:\n");
  for ( i = 0; i < n; i++){
    for ( j = 0; j < d; j++)
      printf("%10.6f ", a_loc[j+i*d]);
    printf("%10.6f\n", b_loc[i]);
  }

  /* Setting global variables */
  dim = d;
  nfacet = dim + 1;
  rank = (1 << nfacet) -1;

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[(1<<nfacet)];
  g[rank] = 0.0;
  int size_a = dim*nfacet;
  int size_b = nfacet;
  int size_alpha=nfacet2; //?
  int size_inv_alpha=(1<<nfacet)*nfacet2;
  int size_det_alpha=(1<<nfacet);
  int size_dgdb=rank*dim;
  double work_space[3*size_a
                   +2*size_b 
                   +size_alpha 
                   +size_inv_alpha 
                   +size_det_alpha
                   +size_dgdb];
  double *w = work_space;
  a = w;     w += size_a;
  b = w;     w += size_b;
  a0 = w;     w += size_a;
  b0 = w;     w += size_b;
  a1 = a_loc;
  b1 = b_loc;
  dadt = w;  w += size_a;
  dbdt = b_loc;
  alpha = w; w += size_alpha;
  inv_alpha = w; w += size_inv_alpha;
  det_alpha = w; w += size_det_alpha;
  dgdb = w;  w += size_dgdb;

  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++)
      a0[j+i*dim] = (i==j)? 1.0 : 0.0;
    b0[i] = 0.0;
  }
  for ( j = 0; j < dim; j++)
    a0[j+dim2] = -1.0;

  printf("dim = %d\n", dim);
  print_matrix(stdout, dim, nfacet, a0, "The value of matrix a0:");
  print_vector(stdout, nfacet, b0, "The value of vector b0:\n");
  print_matrix(stdout, dim, nfacet, a1, "The value of matrix a1:");
  print_vector(stdout, nfacet, b1, "The value of vector b1:\n");
  printf("rank = %d\n", rank);

  for ( i = 0; i < dim; i++)
    for ( j = 0; j < nfacet; j++)
      dadt[i+j*dim] = -a0[i+j*dim]+a1[i+j*dim];
  print_matrix(stdout, dim, nfacet, dadt, "dadt:");

  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++)
      a[j+i*dim] = a0[j+i*dim];
    b[i] = b0[i];
  }

  update_alpha();
  print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  update_inv_alpha();
  int J;
  printf("determinant: ");
  for(J = 0; J < rank; J++){
    //char s[10];
    //sprintf(s, "%d:", J);
    //print_matrix(stdout, nfacet, nfacet, inv_alpha+J*nfacet2, s);
    printf("%g ", det_alpha[J]);
  }
  printf("\n");

  simplex_init(g);
  print_vector(stdout, rank+1, g, "The initial value of g:\n");

  simplex2_rk(g);
  print_vector(stdout, rank, g, "output of runge-kutta:\n");

  get_prob(g);
  print_vector(stdout, rank, g, "result:\n");

  printf("Probability = %g\n", g[0]);
  return;
}

static void 
simplex2_rk(double *g)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = GSL_ODEIV_STEP_TYPE;
  gsl_odeiv_step *s  = gsl_odeiv_step_alloc (T, rank);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (rank);
  gsl_odeiv_system sys = {simplex2_function, NULL, rank, NULL};

  double t = 0.0;
  while (t < 1.0){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, 1.0, &h, g);
    if (status != GSL_SUCCESS)
      break;
    //fprintf(stdout, "r=%lf\n",t);
    //print_vector(stdout, rank+1, g, "g:");
  }
  /* free */
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return;
}

static int 
simplex2_function(double t, const double g[], double dg[], void *params)
{
  simplex2_update(t);
  int J, i, j;
  for (J = 0; J < rank; J++){//for (J = 0; J < rank; J++){
    dg[J] = 0.0;
    for (j = 0; j < nfacet; j++){
      dg[J] += dbdt[j] * get_del_bj(g, J, j);
    }
    for ( i = 0; i < dim; i++){
      for ( j = 0; j < nfacet; j++){
	if ( i != j){ /*???*/
	  dg[J] += dadt[i+j*dim] * get_del_aij(g, J, i,j) ;
	  //printf("J=%d, i=%d, j=%d, del_aij=%f\n", J,i,j,get_del_aij(g, J, i,j));
	}
      }
    }
  }
  //print_vector(stdout, rank, dg, "dg:");
  return GSL_SUCCESS;
}

static void 
simplex2_update(double t)
{
  int i,j;
  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++)
      a[j+i*dim] = (1.0 - t)*a0[j+i*dim] + t*a1[j+i*dim];
    b[i] = t*b1[i];//b[i] = (1-t)*b0[i] + t*b1[i];
  }
  update_alpha();
  update_inv_alpha();
}

static void
simplicial_cone(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  printf("d=%d\n", d);

  int n = d;
  double a_loc[d*n], b_loc[n];
  int i,j;
  for ( i = 0; i < n; i++){
    for ( j = 0; j < d; j++)
      sscanf(*++p, "%lf", a_loc+j+i*d);
    sscanf(*++p, "%lf", b_loc+i);
  }
  printf("a,b:\n");
  for ( i = 0; i < n; i++){
    for ( j = 0; j < d; j++)
      printf("%10.6f ", a_loc[j+i*d]);
    printf("%10.6f\n", b_loc[i]);
  }

  /* Setting global variables */
  dim = d;
  nfacet = d;
  rank = (1 << nfacet);

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[(1<<nfacet)];
  g[rank] = 0.0;
  int size_a = dim*nfacet;
  int size_b = nfacet;
  int size_alpha=nfacet2; 
  int size_inv_alpha=(1<<nfacet)*nfacet2;
  int size_det_alpha=(1<<nfacet);
  int size_dgdb=rank*dim;
  double work_space[3*size_a
                   +2*size_b 
                   +size_alpha 
                   +size_inv_alpha 
                   +size_det_alpha
                   +size_dgdb];
  double *w = work_space;
  a = w;     w += size_a;
  b = w;     w += size_b;
  a0 = w;     w += size_a;
  b0 = w;     w += size_b;
  a1 = a_loc;
  b1 = b_loc;
  dadt = w;  w += size_a;
  dbdt = b_loc;
  alpha = w; w += size_alpha;
  inv_alpha = w; w += size_inv_alpha;
  det_alpha = w; w += size_det_alpha;
  dgdb = w;  w += size_dgdb;

  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++)
      a0[j+i*dim] = (i==j)? 1.0 : 0.0;
    b0[i] = 0.0;
  }

  printf("dim = %d\n", dim);
  print_matrix(stdout, dim, nfacet, a0, "The value of matrix a0:");
  print_vector(stdout, nfacet, b0, "The value of vector b0:\n");
  print_matrix(stdout, dim, nfacet, a1, "The value of matrix a1:");
  print_vector(stdout, nfacet, b1, "The value of vector b1:\n");
  printf("rank = %d\n", rank);

  for ( i = 0; i < dim; i++)
    for ( j = 0; j < nfacet; j++)
      dadt[i+j*dim] = -a0[i+j*dim]+a1[i+j*dim];
  print_matrix(stdout, dim, nfacet, dadt, "dadt:");

  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++)
      a[j+i*dim] = a0[j+i*dim];
    b[i] = b0[i];
  }

  update_alpha();
  print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  update_inv_alpha();
  int J;
  printf("determinant: ");
  for(J = 0; J < rank; J++){
    //char s[10];
    //sprintf(s, "%d:", J);
    //print_matrix(stdout, nfacet, nfacet, inv_alpha+J*nfacet2, s);
    printf("%g ", det_alpha[J]);
  }
  printf("\n");

  simplicial_cone_init(g);
  print_vector(stdout, rank+1, g, "The initial value of g:\n");

  simplicial_cone_rk(g);
  print_vector(stdout, rank, g, "output of runge-kutta:\n");

  get_prob(g);
  print_vector(stdout, rank, g, "result:\n");

  printf("Probability = %g\n", g[0]);
  return;
}

static void 
simplicial_cone_init(double *g)
{
  const double c = sqrt(0.5*M_PI);
  int i, J;

  for (J = 0; J < rank; J++){
    g[J] = 1.0;
    for (i = 0; i < dim; i++)
      if ( J & (1<<i))
	g[J] *= 1.0 / a0[i*(dim+1)];
      else
	g[J] *= c;
  }
  return;
}

static void 
simplicial_cone_rk(double *g)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = GSL_ODEIV_STEP_TYPE;
  gsl_odeiv_step *s  = gsl_odeiv_step_alloc (T, rank);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (rank);
  gsl_odeiv_system sys = {simplicial_cone_function, NULL, rank, NULL};

  double t = 0.0;
  while (t < 1.0){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, 1.0, &h, g);
    if (status != GSL_SUCCESS)
      break;
    //fprintf(stdout, "r=%lf\n",t);
    //print_vector(stdout, rank+1, g, "g:");
  }
  /* free */
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return;
}

static int 
simplicial_cone_function(double t, const double g[], double dg[], void *params)
{
  simplicial_cone_update(t);
  int J, i, j;
  for (J = 0; J < rank; J++){
    dg[J] = 0.0;
    for (j = 0; j < nfacet; j++){
      dg[J] += dbdt[j] * get_del_bj(g, J, j);
    }
    for ( i = 0; i < dim; i++){
      for ( j = 0; j < nfacet; j++){
	if ( i != j){ /*???*/
	  dg[J] += dadt[i+j*dim] * get_del_aij(g, J, i,j) ;
	  //printf("J=%d, i=%d, j=%d, del_aij=%f\n", J,i,j,get_del_aij(g, J, i,j));
	}
      }
    }
  }
  //print_vector(stdout, rank, dg, "dg:");
  return GSL_SUCCESS;
}

static void 
simplicial_cone_update(double t)
{
  int i,j;
  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++)
      a[j+i*dim] = (1.0 - t)*a0[j+i*dim] + t*a1[j+i*dim];
    b[i] = t*b1[i];//b[i] = (1-t)*b0[i] + t*b1[i];
  }
  update_alpha();
  update_inv_alpha();
}

static void
tukey(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  //printf("d=%d\n", d);
  double c;
  sscanf(*++p, "%lf", &c);
  //printf("c=%f\n", c);

  /* Setting global variables */
  dim = d-1;
  nfacet = d;
  rank = (1<<nfacet)-1;

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[nfacet2];
  int size_a = dim*nfacet;
  int size_alpha=nfacet2;
  int size_b = 3*nfacet2; /* b <- C_ijk */
  double work_space[size_a
                   +size_alpha
                   +size_b];
  double *w = work_space;
  a = w;          w += size_a;
  alpha = w;      w += size_alpha;
  b = w; 

  /* a[j+i*dim] <- (transpose matrix of a in our paper) */
  int i,j,k;
  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++){
      a[j+i*dim] = 0.0;
    }
  }
  for ( i = 0; i < d-1; i++)
    a[i+i*dim] = sqrt(i+2);
  for ( i = 0; i < d-2; i++){
    a[i+(i+1)*dim] = -sqrt(i+1);
    a[i+(d-1)*dim] = -1.0/sqrt((i+1)*(i+2));
  }
  a[(d-2)+(d-1)*dim] = -sqrt((double) d/(d-1));
  //printf("dim = %d\n", dim);
  //print_matrix(stdout, dim, nfacet, a, "The value of matrix a:");

  /* alpha <- a^\top a*/
  update_alpha();
  //print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  /* fac_d <- d! */
  double fac_d = 1.0;
  for ( i = 0; i < d; i++){
    fac_d *= (i+1.0);
  }
  //printf("d!=%g\n", 1.0 * fac_d);

  /* b <- C_ijk, g <- initial value */
  for ( i = 0; i < size_b; i++) b[i] = 0.0;
  for (i = 0; i < nfacet2; i++) g[i] = 0.0;
  b[0] = c;
  for ( i = 0; i < d-1; i++)
    for ( j = i+1; j < d; j++){
      double *b_ij = b + 3*(j+i*nfacet);
      //fprintf(stderr, "done %d\n",nfacet2);
      double inv[nfacet2], det;
      //fprintf(stderr, "done\n");
      inv_submat_alpha2(i,j, inv, &det );      

      /* g <- initial value */
      if ( j == i+1){
        g[i+j*nfacet] = fac_d/det;
      }

      /* b <- C_ijk */
      b_ij[0] = -inv[(d-1)+(d-1)*nfacet];

      if ( i+1 != j){
        b_ij[1] = 0.0;
        b_ij[2] = 0.0;
        for ( k = 0; k < i; k++){
	  b_ij[1] += -inv[(d-1)+k*nfacet]*alpha[k+i*nfacet];
	  b_ij[2] += -inv[(d-1)+k*nfacet]*alpha[k+(j-1)*nfacet];
        }
        for ( k = j; k < d; k++){
	  b_ij[1] += -inv[(d-1)+k*nfacet]*alpha[k+i*nfacet];
	  b_ij[2] += -inv[(d-1)+k*nfacet]*alpha[k+(j-1)*nfacet];
        }
        b_ij[1] *= INV_SQRT_TWO_PI;
        b_ij[2] *= INV_SQRT_TWO_PI;
      } else {
        b_ij[1] = 0.0;
        b_ij[2] = 0.0;
      }
      //printf("i,j %d %d,", i+1, j+1);
      //printf(" %g, %g, %g\n", b_ij[0],b_ij[1],b_ij[2]);
      //printf(" %d, %d, %d\n", 0+3*(j+i*nfacet),1+3*(j+i*nfacet),2+3*(j+i*nfacet));
    }
  //print_matrix(stdout, nfacet, nfacet, g, "The initial value of g:");
  tukey_rk(g);
  //print_matrix(stdout, nfacet, nfacet, g, "output of runge-kutta:");

  //printf("Probability = %20.18g\n", 1-g[0]);
  printf("%20.18g\t", g[0]);
  return;
}

static void 
tukey_rk(double *g)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = GSL_ODEIV_STEP_TYPE;
  gsl_odeiv_step *s  = gsl_odeiv_step_alloc (T, nfacet2);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (nfacet2);
  gsl_odeiv_system sys = {tukey_function, NULL, nfacet2, NULL};

  double t = 0.0, t1 = b[0];
    while (t < t1){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, g);
    if (status != GSL_SUCCESS)
      break;
    //fprintf(stderr, "%lf %lf\n",t, g[0]);
    //print_vector(stdout, rank+1, g, "\n");
  }
  /* free */
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return;
}

static int 
tukey_function(double t, const double g[], double dg[], void *params)
{
  int d = nfacet;
  dg[0] = INV_SQRT_TWO_PI * g[(d-1)*nfacet];

  int i,j;
  double *b_ij;
  for ( i = 0; i < d-1; i++)
    for ( j = i+1; j < d; j++){
      b_ij  = b + 3*(j+i*nfacet);
      dg[i+j*nfacet]  = t* b_ij[0] * g[i+j*nfacet];
      if ( i+1 != j){
        dg[i+j*nfacet] +=    b_ij[1] * g[(i+1)+j*nfacet];
        dg[i+j*nfacet] +=    b_ij[2] * g[i+(j-1)*nfacet];
      }
    }
  return GSL_SUCCESS;
}

static void
tukey2(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  printf("d=%d\n", d);
  double c;
  sscanf(*++p, "%lf", &c);
  printf("c=%f\n", c);

  /* Setting global variables */
  dim = d-1;
  nfacet = d;
  rank = (1<<nfacet)-1;

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[nfacet2];
  int size_a = dim*nfacet;
  int size_alpha=nfacet2;
  int size_b = 3*nfacet2; /* b <- C_ijk */
  double work_space[size_a
                   +size_alpha
                   +size_b];
  double *w = work_space;
  a = w;          w += size_a;
  alpha = w;      w += size_alpha;
  b = w; 

  /* a[j+i*dim] <- (transpose matrix of a in our paper) */
  int i,j,k;
  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++){
      a[j+i*dim] = 0.0;
    }
  }
  for ( i = 0; i < d-1; i++)
    a[i+i*dim] = sqrt(i+2);
  for ( i = 0; i < d-2; i++){
    a[i+(i+1)*dim] = -sqrt(i+1);
    a[i+(d-1)*dim] = -1.0/sqrt((i+1)*(i+2));
  }
  a[(d-2)+(d-1)*dim] = -sqrt((double) d/(d-1));
  printf("dim = %d\n", dim);
  //print_matrix(stdout, dim, nfacet, a, "The value of matrix a:");

  /* alpha <- a^\top a*/
  update_alpha();
  //print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  for ( i = 0; i < nfacet*dim; i++) a[i] = 0.0;

  double inv[nfacet2], det;
  /* b <- C_ijk, g <- initial value */
  for ( i = 0; i < size_b; i++) b[i] = 0.0;
  for (i = 0; i < nfacet2; i++) g[i] = 0.0;

  for( i = d-2; i > -1; i--)  //for ( i = 0; i < d-1; i++)
    for ( j = i+1; j < d; j++){
      /* g <- initial value */
      if ( j == i+1){
        g[i+j*nfacet] = 1.0;
      }

      double *b_ij = b + 3*(j+i*nfacet);
      inv_submat_alpha2(i,j, inv, &det );
      a[i+j*dim] = det;
      //printf("det %d %d %lf %d\n", i+1,j+1,det,i+j*dim);

      /* b <- C_ijk */
      b_ij[0] = -inv[(d-1)+(d-1)*nfacet];

      if ( i+1 != j){
        b_ij[1] = 0.0;
        b_ij[2] = 0.0;
        for ( k = 0; k < i; k++){
	  b_ij[1] += -inv[(d-1)+k*nfacet]*alpha[k+i*nfacet];
	  b_ij[2] += -inv[(d-1)+k*nfacet]*alpha[k+(j-1)*nfacet];
        }
        for ( k = j; k < d; k++){
	  b_ij[1] += -inv[(d-1)+k*nfacet]*alpha[k+i*nfacet];
	  b_ij[2] += -inv[(d-1)+k*nfacet]*alpha[k+(j-1)*nfacet];
        }
        b_ij[1] *= (j-i)*INV_SQRT_TWO_PI*a[i+j*dim]/a[(i+1)+j*dim];
        b_ij[2] *= (j-i)*INV_SQRT_TWO_PI*a[i+j*dim]/a[i+(j-1)*dim];
      } else {
        b_ij[1] = 0.0;
        b_ij[2] = 0.0;
      }
      printf("i,j %d %d,", i+1, j+1);
      printf(" %g, %g, %g\n", b_ij[0],b_ij[1],b_ij[2]);
      printf(" %d, %d, %d\n", 0+3*(j+i*dim),1+3*(j+i*dim),2+3*(j+i*dim));
    }
  //for ( i = 0; i < nfacet*dim; i++) a[i] = i;
  //print_matrix(stdout, dim, nfacet, a, "determinants a:");
  b[0] = d*INV_SQRT_TWO_PI/sqrt(alpha[nfacet2-1]);
  //print_vector(stdout, 3*nfacet2, b, "b");

  //print_matrix(stdout, nfacet, nfacet, g, "The initial value of g:");
  a[0] = c;
  tukey2_rk(g);
  //print_matrix(stdout, nfacet, nfacet, g, "output of runge-kutta:");

  printf("err = %20.18g\n", 1-g[0]);
  //printf("Probability = %20.18g\n", g[0]);
  return;
}

static void 
tukey2_rk(double *g)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = GSL_ODEIV_STEP_TYPE;
  gsl_odeiv_step *s  = gsl_odeiv_step_alloc (T, nfacet2);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (nfacet2);
  gsl_odeiv_system sys = {tukey2_function, NULL, nfacet2, NULL};

  double t = 0.0, t1 = a[0];
    while (t < t1){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, g);
    if (status != GSL_SUCCESS)
      break;
    //fprintf(stderr, "%lf %lf\n",t, g[0]);
    //print_vector(stderr, nfacet2, g, "");
  }
  /* free */
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return;
}

static int 
tukey2_function(double t, const double g[], double dg[], void *params)
{
  int d = nfacet;
  dg[0] = b[0] * g[(d-1)*nfacet];

  int i,j;
  double *b_ij;
  for ( i = 0; i < d-1; i++)
    for ( j = i+1; j < d; j++){
      b_ij  = b + 3*(j+i*nfacet);
      dg[i+j*nfacet]  = t* b_ij[0] * g[i+j*nfacet];
      if ( i+1 != j){
        dg[i+j*nfacet] +=    b_ij[1] * g[(i+1)+j*nfacet];
        dg[i+j*nfacet] +=    b_ij[2] * g[i+(j-1)*nfacet];
      }
    }
  return GSL_SUCCESS;
}

static void
tukey3(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  printf("d=%d\n", d);
  double c;
  sscanf(*++p, "%lf", &c);
  printf("c=%f\n", c);

  /* Setting global variables */
  dim = d-1;
  nfacet = d;
  rank = (1<<nfacet)-1;

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[nfacet2];
  int size_a = nfacet2; //dim*nfacet;
  int size_alpha=nfacet2;
  int size_b = 5*nfacet2; /* b <- C_ijk */
  double work_space[size_a
                   +size_alpha
                   +size_b];
  double *w = work_space;
  a = w;          w += size_a;
  alpha = w;      w += size_alpha;
  b = w; 

  /* a[j+i*dim] <- (transpose matrix of a in our paper) */
  int i,j,k;
  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++){
      a[j+i*dim] = 0.0;
    }
  }
  for ( i = 0; i < d-1; i++)
    a[i+i*dim] = sqrt(i+2);
  for ( i = 0; i < d-2; i++){
    a[i+(i+1)*dim] = -sqrt(i+1);
    a[i+(d-1)*dim] = -1.0/sqrt((i+1)*(i+2));
  }
  a[(d-2)+(d-1)*dim] = -sqrt((double) d/(d-1));
  printf("dim = %d\n", dim);
  //print_matrix(stdout, dim, nfacet, a, "The value of matrix a:");

  /* alpha <- a^\top a*/
  update_alpha();
  //print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  //for ( i = 0; i < nfacet*dim; i++) a[i] = 0.0;

  /* b <- C_ijk, g <- initial value */
  for ( i = 0; i < size_b; i++) b[i] = 0.0;
  for (i = 0; i < nfacet2; i++) g[i] = 0.0;
  g[0] = 0.0;
  double inv[nfacet2], det, lambda[dim*nfacet];
  for( i = d-2; i > -1; i--)  //for ( i = 0; i < d-1; i++)
    for ( j = i+1; j < d; j++){
      /* g <- initial value */
      if ( j == i+1){
        g[i+j*nfacet] = 1.0;
      }

      double *b_ij = b + 5*(j+i*nfacet);
      inv_submat_alpha2(i,j, inv, &det );
      a[i+j*dim] = det;
      lambda[i+j*dim] = inv[(d-1)+(d-1)*nfacet];
      //printf("i=%d j=%d det=%lf %d ", i+1,j+1,det,i+j*dim);
      //printf("lambda=%lf %d\n", lambda[j+i*dim],j+i*dim);

      /* b <- C_ijk */
      b_ij[0] = 0.0;
      //b_ij[0] = -inv[(d-1)+(d-1)*nfacet] + lambda[i+j*dim];
      //printf("b_ij=%lf\n",b_ij[0]);

      if ( i+1 != j){
        b_ij[1] = 0.0;
        b_ij[2] = 0.0;
        for ( k = 0; k < i; k++){
	  b_ij[1] += -inv[(d-1)+k*nfacet]*alpha[k+i*nfacet];
	  b_ij[2] += -inv[(d-1)+k*nfacet]*alpha[k+(j-1)*nfacet];
        }
        for ( k = j; k < d; k++){
	  b_ij[1] += -inv[(d-1)+k*nfacet]*alpha[k+i*nfacet];
	  b_ij[2] += -inv[(d-1)+k*nfacet]*alpha[k+(j-1)*nfacet];
        }
        b_ij[1] *= (j-i-1)*INV_SQRT_TWO_PI*a[i+j*dim]/a[(i+1)+j*dim];
        b_ij[2] *= (j-i-1)*INV_SQRT_TWO_PI*a[i+j*dim]/a[i+(j-1)*dim];
        b_ij[3] = -0.5 * (lambda[i+1+j*dim] - lambda[i+j*dim]);
        b_ij[4] = -0.5 * (lambda[i+(j-1)*dim] - lambda[i+j*dim]);
      } else {
        b_ij[1] = 0.0;
        b_ij[2] = 0.0;
        b_ij[3] = 0.0;
        b_ij[4] = 0.0;
      }
      //printf("i,j:%2d %2d:", i+1, j+1);
      //printf("%g\t%5.3g\t%5.3g\t%5.3g\t%5.3g\n", b_ij[0],b_ij[1],b_ij[2],b_ij[3],b_ij[4]);
      //printf(" %d, %d, %d\n", 0+3*(j+i*dim),1+3*(j+i*dim),2+3*(j+i*dim));
    }
  b[0] = d*(d-1)*INV_SQRT_TWO_PI/sqrt(2.0); //sqrt(alpha[nfacet2-1]); /*2?*/
  //printf("alpha[nfacet2-1]=%lf\n",alpha[nfacet2-1]);
  b[1] = -0.25;
  //print_matrix(stdout, nfacet, nfacet, a, "The value of matrix a:");
  //print_matrix(stdout, nfacet, nfacet, g, "The initial value of g:");
  a[0] = c;
  tukey3_rk(g);
  print_matrix(stdout, nfacet, nfacet, g, "output of runge-kutta:");
  //print_vector(stdout, nfacet2, g, "output of runge-kutta:");

  printf("err = %20.18g\n", 1.0-g[0]);
  printf("Probability = %20.18g\n", g[0]);
  return;
}

static void 
tukey3_rk(double *g)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = GSL_ODEIV_STEP_TYPE;
  gsl_odeiv_step *s  = gsl_odeiv_step_alloc (T, nfacet2);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (nfacet2);
  gsl_odeiv_system sys = {tukey3_function, NULL, nfacet2, NULL};

  double t = 0.0, t1 = a[0];
    while (t < t1){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, g);
    if (status != GSL_SUCCESS)
      break;
    //fprintf(stderr, "%lf %lf\n",t, g[0]);
    fprintf(stderr, "%lf ",t);
    print_vector(stderr, nfacet2, g, "");
  }
  /* free */
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return;
}

static int 
tukey3_function(double t, const double g[], double dg[], void *params)
{
  int d = nfacet;
  double t2 = t*t;
  dg[0] = b[0] * exp(t2*b[1]) * g[(d-1)*nfacet];

  int i,j;
  double *b_ij;
  for ( i = 0; i < d-1; i++)
    for ( j = i+1; j < d; j++){
      b_ij  = b + 5*(j+i*nfacet);
      //dg[i+j*nfacet]  = t* b_ij[0] * g[i+j*nfacet];
      dg[i+j*nfacet]  = 0.0;
      if ( i+1 != j){
        dg[i+j*nfacet] += b_ij[1] * exp(b_ij[3]*t2) * g[(i+1)+j*nfacet];
        dg[i+j*nfacet] += b_ij[2] * exp(b_ij[4]*t2) * g[i+(j-1)*nfacet];
      }
    }
  return GSL_SUCCESS;
}

static void 
update_alpha(void)
{
  int i,j;
  for(i = 0; i< nfacet; i++)
    for(j = 0; j< nfacet; j++){
      int k;
      double sum = 0.0;
      for(k = 0; k < dim; k++)
        sum += a[k+i*dim] * a[k+j*dim];
      alpha[i + j*nfacet] = sum;
    }
  return;
}

static void 
update_inv_alpha(void)
{
  int J;
  for(J = 0; J < rank; J++){
    inv_submat_alpha(J, inv_alpha+J*nfacet2);
  }
  return;
}

static int
inv_submat_alpha(int I, double *inv)
{
  if (I==0){
    det_alpha[I] = 1.0;
    /* inv <- (unit matrix of size nfacet) */
    int i;
    for ( i = 0; i < nfacet*nfacet; i++)
      inv[i] = 0.0;
    for ( i = 0; i < nfacet; i++)
      inv[i*(nfacet+1)] = 1.0;
    return 0;
  }

  int size_of_submat = cardinality(I, nfacet);
  double mat[size_of_submat * size_of_submat];
  int i,j,k,l;
  k = 0;
  for(i=0; i<nfacet; i++)
    if(I&1<<i){
      l = 0;
      for(j=0; j<nfacet; j++)
	if(I&1<<j){
	  mat[l+k*size_of_submat] = alpha[i+j*nfacet];
	  l++;
	}
      k++;
    }

  /* mat <- cholesky decomposition of mat */
  int info, n, m;
  n = m = size_of_submat;
  dpotrf_("U", &n, mat, &m, &info);
  if (info != 0) {
    printf("info:I=%d\n", I);
    return info;
  }

  det_alpha[I] = 1.0;
  for ( k = 0; k < size_of_submat; k++)
    det_alpha[I] *= mat[k*(size_of_submat+1)];

  n = m = size_of_submat;
  dpotri_("U", &n, mat, &m, &info);
  if (info != 0) {
    printf("info:I=%d\n", I);
    return info;
  }

  for ( i = 0; i<nfacet2; i++) inv[i] = 0.0;

  k = 0;
  for(i=0; i<nfacet; i++)
    if(I&1<<i){
      l = k;
      for(j=i; j<nfacet; j++)
	if(I&1<<j){
	  inv[i+j*nfacet] = inv[j+i*nfacet] = mat[k+l*size_of_submat];
	  l++;
	}else if(i==j){
	  inv[i+j*nfacet] = 1.0;
	}else{
	  inv[i+j*nfacet] = 0.0;
	}
      k++;
    }else{
      for(j=0; j<nfacet; j++)
	if(i==j)
	  inv[i+j*nfacet] = 1.0;
	else
	  inv[i+j*nfacet] = 0.0;
    }
  return info;
}

static int
inv_submat_alpha2(int ii, int jj, double *inv, double *det)
{
  //int size_of_submat = cardinality(I, nfacet);
  int size_of_submat= ii - jj + nfacet;
  //printf("\t ii,jj,d %d %d %d %d\n",ii+1,jj+1,d,size_of_submat);
  
  double mat[size_of_submat * size_of_submat];
  int i,j,k,l;
  k = 0;
  for(i=0; i<nfacet; i++)
    if( i<ii || i>jj-1){
      l = 0;
      for(j=0; j<nfacet; j++)
	if( j<ii || j>jj-1){
	  mat[l+k*size_of_submat] = alpha[i+j*nfacet];
	  l++;
	}
      k++;
    }

  /* mat <- cholesky decomposition of mat */
  int info, n, m;
  n = m = size_of_submat;
  dpotrf_("U", &n, mat, &m, &info);
  if (info != 0) {
    printf("info:(i,j)=(%d,%d)\n", ii,jj);
    return info;
  }

  *det = 1.0;
  for ( k = 0; k < size_of_submat; k++)
    *det *= mat[k*(size_of_submat+1)];

  n = m = size_of_submat;
  dpotri_("U", &n, mat, &m, &info);
  if (info != 0) {
    printf("info:(i,j)=(%d,%d)\n", ii,jj);
    return info;
  }

  for ( i = 0; i<nfacet2; i++) inv[i] = 0.0;

  k = 0;
  for(i=0; i<nfacet; i++)
    if( i<ii || i>jj-1){
      l = k;
      for(j=i; j<nfacet; j++)
	if( j<ii || j>jj-1){
	  inv[i+j*nfacet] = inv[j+i*nfacet] = mat[k+l*size_of_submat];
	  l++;
	}else if(i==j){
	  inv[i+j*nfacet] = 1.0;
	}else{
	  inv[i+j*nfacet] = 0.0;
	}
      k++;
    }else{
      for(j=0; j<nfacet; j++)
	if(i==j)
	  inv[i+j*nfacet] = 1.0;
	else
	  inv[i+j*nfacet] = 0.0;
    }
  return info;
}

static int
cardinality(int J, int max)
{
  int s = 0, j;
  for ( j = 0; j < max; j++)
    if ( J & (1<<j) )
      s++;
  return s;
}

static double
get_del_aij(const double g[], int J, int i, int j)
{
  double s = 0.0;
  int k;
  int kdim=0;
  for ( k = 0; k < nfacet; k++){
    s += a[i+kdim] * get_del_bibj(g, J, k, j);
    kdim += dim;
  }
  return s;
}

static double 
get_del_bibj(const double g[], int J, int i, int j)
{
  double retv = 0.0;
  int shift_j = 1 << j;
  if ( J & shift_j ){
    double s1 = 0.0;
    int k;
    for ( k = 0; k < nfacet; k++)
      if ( J & ( 1<< k)){
	double s2 = 0.0;
	if ( k == i)
	  s2 += g[J];
	s2 += b[k] * get_del_bj(g, J, i);
	//s2 += b[k] * dgdb[i+J*dim];
	int l;
	for ( l = 0; l < nfacet; l++){
	  int shift_l = 1<<l;
	  if (!( J&shift_l ) ){
	    s2 += alpha[k+l*nfacet] * get_del_bj(g, J|shift_l, i);
	    //s2 += alpha[k+l*nfacet] * dgdb[i+(J|shift_l)*dim];
	    //printf("\t\t J=%d, i=%d, bj=%f\n", J|shift_l,i,get_del_bj(g, J|shift_l, i));
	  }
	}
	s1 -= inv_alpha[j+k*nfacet+J*nfacet2] * s2;
	//	printf("\t J=%d, i=%d, j=%d, k=%d, inv_alpha=%f, s2=%f\n", J,i,j,k, inv_alpha[j+k*nfacet+J*nfacet2], s2);
      }
    retv = s1;
  } else{
    retv = get_del_bj(g, J|shift_j, i);
    //retv = dgdb[i+(J|shift_j)*dim];
  }
  //printf("J=%d, i=%d, j=%d, bibj=%f\n", J,i,j,retv);
  return retv;
}

static double 
get_del_bj(const double g[], int J, int j)
{
  double retv;
  if ( J == rank) {
    retv = 0.0;
  }else if ( J & (1<<j)){
    int k, l;
    double s1 , s2;
    s1 = 0.0;
    for (k = 0; k < nfacet; k++)
      if ( J & (1<<k)){
        s2 = b[k] * g[J];
        for (l = 0; l < nfacet; l++)
          if ( !(J & (1<<l)) )
            s2 += alpha[k+l*nfacet] * g[J|(1<<l)];
        s1 -= inv_alpha[j+k*nfacet+J*nfacet2] * s2;
	//printf("\t J=%d, j=%d, k=%d, inv_alpha=%f\n", J,j,k, inv_alpha[j+k*nfacet+J*nfacet2]);
      }
    retv = s1;
  } else
    retv = g[J|(1<<j)];
  //dgdb[j+J*dim] = retv;
  return retv;
}

static void
get_prob(double *g)
{
  double normal_const[dim+1];
  const double sqrt_two_pi = sqrt(2*M_PI);
  int i;
  normal_const[0] = 1.0;
  for (i = 0; i < dim; i++)
    normal_const[i+1] = sqrt_two_pi * normal_const[i];

  int J;
  for ( J = 0; J < rank; J++){
    int d = dim;
    for ( i = 0; i < nfacet; i++)
      if ( J & (1<<i))
	d--;
    g[J] *= det_alpha[J]/normal_const[d];
  }
  return;
}

static void 
print_vector(FILE *fp, int length, double *v, char *str)
{
  int i;
  fprintf(fp, "%s", str);
  for(i=0; i<length; i++)
    fprintf(fp, "%15.10f ", v[i]);
  fprintf(fp, "\n");
}

static void 
print_matrix(FILE *fp, int nrow, int ncol, double *matrix, const char *str)
{
  int i, j;

  fprintf(fp,"%s\n", str);
  for(i=0; i<nrow; i++){
    for (j=0; j<ncol; j++)
      fprintf(fp,"%15.10lf ",matrix[i+j*nrow]);
    fprintf(fp,"\n");
  }
}    

static void 
numerical_experiment1(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  //printf("d=%d\n", d);

  int n = d+1;
  double a_loc[d*n], b_loc[n];
  int i,j;
  for ( i = 0; i < n*d; i++)
    a_loc[i] = 0.0;
  for ( j = 0; j < d; j++){
    a_loc[j+j*d] = 1.0;
    a_loc[j+d*d] = -1.0;
  }
  double h_sqrtd = 0.5 * sqrt(d);
  for ( i = 0; i < n; i++){
    b_loc[i] = h_sqrtd;
  }
  //for ( i = 0; i < n; i++){
  //  for ( j = 0; j < d; j++)
  //    sscanf(*++p, "%lf", a_loc+j+i*d);
  //  sscanf(*++p, "%lf", b_loc+i);
  //}
  //printf("a,b:\n");
  //for ( i = 0; i < n; i++){
  //  for ( j = 0; j < d; j++)
  //    printf("%10.6f ", a_loc[j+i*d]);
  //  printf("%10.6f\n", b_loc[i]);
  //}

  /* Setting global variables */
  dim = d;
  nfacet = dim + 1;
  rank = (1 << nfacet) -1;

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[rank];
  //double g[(1<<nfacet)];
  //int size_a = dim2;
  int size_b = nfacet;
  int size_alpha=nfacet2;
  int size_inv_alpha=(1<<nfacet)*nfacet2;
  int size_det_alpha=(1<<nfacet);
  int size_dgdb=rank*dim;
  double work_space[size_b 
                   +size_alpha 
                   +size_inv_alpha 
                   +size_det_alpha
                   +size_dgdb];
  double *w = work_space;
  a = a_loc;           
  b = w;          w += size_b;
  a1 = a_loc;
  b1 = b_loc;
  dadt = NULL;
  dbdt = b_loc;
  alpha = w;      w += size_alpha;
  inv_alpha = w;  w += size_inv_alpha;
  det_alpha = w;  w += size_det_alpha;
  dgdb = w;  w += size_dgdb;

  //printf("dim = %d\n", dim);
  //print_matrix(stdout, dim, nfacet, a1, "The value of matrix a:");
  //print_vector(stdout, nfacet, b1, "The value of vector b:\n");
  //printf("rank = %d\n", rank);

  update_alpha();
  //print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  update_inv_alpha();
  //int J;
  //printf("determinant: ");
  //for(J = 0; J < rank; J++){
    //char s[10];
    //sprintf(s, "%d:", J);
    //print_matrix(stdout, nfacet, nfacet, inv_alpha+J*nfacet2, s);
    //printf("%g ", det_alpha[J]);
  //}
  //printf("\n");

  simplex_init(g);
  //print_vector(stdout, rank+1, g, "The initial value of g:\n");

  simplex_rk(g);
  //print_vector(stdout, rank, g, "output of runge-kutta:\n");

  get_prob(g);
  //print_vector(stdout, rank, g, "result:\n");

  printf("probability= %20.18f\n", g[0]);
  return;
}

static void 
numerical_experiment2(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  //printf("d=%d\n", d);

  int n = d;
  double a_loc[d*n], b_loc[n];
  int i,j;
  for ( i = 0; i < n*d; i++)
    a_loc[i] = 0.0;
  for ( i = 0; i < n; i++){
    a_loc[i+i*d] = 1.0;
    for ( j = i+1; j < d; j++){
      a_loc[j+i*d] = 0.01 * (i+j+2);
    }
  }
  double h_sqrtd = 0.5 * sqrt(d);
  for ( i = 0; i < n; i++){
    b_loc[i] = h_sqrtd;
  }
  //printf("a,b:\n");
  //for ( i = 0; i < n; i++){
  //  for ( j = 0; j < d; j++)
      //printf("%10.6f ", a_loc[j+i*d]);
    //printf("%10.6f\n", b_loc[i]);
  // }

  /* Setting global variables */
  dim = d;
  nfacet = d;
  rank = (1 << nfacet);

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[(1<<nfacet)];
  g[rank] = 0.0;
  int size_a = dim*nfacet;
  int size_b = nfacet;
  int size_alpha=nfacet2; 
  int size_inv_alpha=(1<<nfacet)*nfacet2;
  int size_det_alpha=(1<<nfacet);
  int size_dgdb=rank*dim;
  double work_space[3*size_a
                   +2*size_b 
                   +size_alpha 
                   +size_inv_alpha 
                   +size_det_alpha
                   +size_dgdb];
  double *w = work_space;
  a = w;     w += size_a;
  b = w;     w += size_b;
  a0 = w;     w += size_a;
  b0 = w;     w += size_b;
  a1 = a_loc;
  b1 = b_loc;
  dadt = w;  w += size_a;
  dbdt = b_loc;
  alpha = w; w += size_alpha;
  inv_alpha = w; w += size_inv_alpha;
  det_alpha = w; w += size_det_alpha;
  dgdb = w;  w += size_dgdb;

  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++)
      a0[j+i*dim] = (i==j)? 1.0 : 0.0;
    b0[i] = 0.0;
  }

  //printf("dim = %d\n", dim);
  //print_matrix(stdout, dim, nfacet, a0, "The value of matrix a0:");
  //print_vector(stdout, nfacet, b0, "The value of vector b0:\n");
  //print_matrix(stdout, dim, nfacet, a1, "The value of matrix a1:");
  //print_vector(stdout, nfacet, b1, "The value of vector b1:\n");
  //printf("rank = %d\n", rank);

  for ( i = 0; i < dim; i++)
    for ( j = 0; j < nfacet; j++)
      dadt[i+j*dim] = -a0[i+j*dim]+a1[i+j*dim];
  //print_matrix(stdout, dim, nfacet, dadt, "dadt:");

  for ( i = 0; i < nfacet; i++){
    for ( j = 0; j < dim; j++)
      a[j+i*dim] = a0[j+i*dim];
    b[i] = b0[i];
  }

  update_alpha();
  //print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  update_inv_alpha();
  //int J;
  //printf("determinant: ");
  //for(J = 0; J < rank; J++){
    //char s[10];
    //sprintf(s, "%d:", J);
    //print_matrix(stdout, nfacet, nfacet, inv_alpha+J*nfacet2, s);
    //printf("%g ", det_alpha[J]);
  //}
  //printf("\n");

  simplicial_cone_init(g);
  //print_vector(stdout, rank+1, g, "The initial value of g:\n");

  simplicial_cone_rk(g);
  //print_vector(stdout, rank, g, "output of runge-kutta:\n");

  get_prob(g);
  //print_vector(stdout, rank, g, "result:\n");

  printf("probability= %20.18f\n", g[0]);
  return;
}

static void 
numerical_experiment3(char **p)
{
  /* Input */
  int d;
  sscanf(*++p,"%d", &d);
  //printf("d=%d\n", d);

  int n = d+1;
  double a_loc[d*n], b_loc[n];
  int i,j;
  for ( i = 0; i < n*d; i++)
    a_loc[i] = 0.0;
  for ( j = 0; j < d; j++){
    a_loc[j+j*d] = 1.0;
    a_loc[j+d*d] = -1.0;
  }
  double h_sqrtd = 0.5 * sqrt(d);
  for ( i = 0; i < n-1; i++){
    b_loc[i] = (-1.0) * h_sqrtd;
  }
  b_loc[n-1] = (d+0.5) * sqrt(d);
  //for ( i = 0; i < n; i++){
  //  for ( j = 0; j < d; j++)
  //    sscanf(*++p, "%lf", a_loc+j+i*d);
  //  sscanf(*++p, "%lf", b_loc+i);
  //}
  //printf("a,b:\n");
  //for ( i = 0; i < n; i++){
  //  for ( j = 0; j < d; j++)
  //    printf("%10.6f ", a_loc[j+i*d]);
  //  printf("%10.6f\n", b_loc[i]);
  //}

  /* Setting global variables */
  dim = d;
  nfacet = dim + 1;
  rank = (1 << nfacet) -1;

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[rank];
  //double g[(1<<nfacet)];
  //int size_a = dim2;
  int size_b = nfacet;
  int size_alpha=nfacet2;
  int size_inv_alpha=(1<<nfacet)*nfacet2;
  int size_det_alpha=(1<<nfacet);
  int size_dgdb=rank*dim;
  double work_space[size_b 
                   +size_alpha 
                   +size_inv_alpha 
                   +size_det_alpha
                   +size_dgdb];
  double *w = work_space;
  a = a_loc;           
  b = w;          w += size_b;
  a1 = a_loc;
  b1 = b_loc;
  dadt = NULL;
  dbdt = b_loc;
  alpha = w;      w += size_alpha;
  inv_alpha = w;  w += size_inv_alpha;
  det_alpha = w;  w += size_det_alpha;
  dgdb = w;  w += size_dgdb;

  //printf("dim = %d\n", dim);
  //print_matrix(stdout, dim, nfacet, a1, "The value of matrix a:");
  //print_vector(stdout, nfacet, b1, "The value of vector b:\n");
  //printf("rank = %d\n", rank);

  update_alpha();
  //print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  update_inv_alpha();
  //int J;
  //printf("determinant: ");
  //for(J = 0; J < rank; J++){
    //char s[10];
    //sprintf(s, "%d:", J);
    //print_matrix(stdout, nfacet, nfacet, inv_alpha+J*nfacet2, s);
    //printf("%g ", det_alpha[J]);
  //}
  //printf("\n");

  simplex_init(g);
  //print_vector(stdout, rank+1, g, "The initial value of g:\n");

  simplex_rk(g);
  //print_vector(stdout, rank, g, "output of runge-kutta:\n");

  get_prob(g);
  //print_vector(stdout, rank, g, "result:\n");

  printf("probability= %20.18e\n", g[0]);
  return;
}

static void 
numerical_experiment4(char **p)
{
  /* Input */
  int d;
  double c0;
  sscanf(*++p,"%d", &d);
  sscanf(*++p,"%lf", &c0);
  //printf("d=%d\n", d);
  d = d - 1;

  int n = d+1;
  double a_loc[d*n], b_loc[n];
  int i,j;
  for ( i = 0; i < n*d; i++)
    a_loc[i] = 0.0;
  a_loc[0] = sqrt(2.0);
  for ( i = 1; i < d; i++){
    a_loc[i-1+i*d] = -sqrt(i);
    a_loc[i+i*d] = sqrt(i+2);
  }
  for ( j = 0; j < d-1; j++){
    a_loc[j+d*d] = -1.0/sqrt((j+1)*(j+2));
  }
  a_loc[d-1+d*d] = -(d+1.0)/sqrt(d*(d+1));
  //double h_sqrtd = 0.5 * sqrt(d);
  for ( i = 0; i < n-1; i++){
    b_loc[i] = 0.0;
  }
  b_loc[n-1] = c0;
  //for ( i = 0; i < n; i++){
  //  for ( j = 0; j < d; j++)
  //    sscanf(*++p, "%lf", a_loc+j+i*d);
  //  sscanf(*++p, "%lf", b_loc+i);
  //}

  //printf("a,b:\n");
  //for ( i = 0; i < n; i++){
  //  for ( j = 0; j < d; j++)
  //    printf("%10.6f ", a_loc[j+i*d]);
  //  printf("%10.6f\n", b_loc[i]);
  //}

  /* Setting global variables */
  dim = d;
  nfacet = dim + 1;
  rank = (1 << nfacet) -1;

  dim2 = dim * dim;
  nfacet2 = nfacet*nfacet;

  double g[rank];
  //double g[(1<<nfacet)];
  //int size_a = dim2;
  int size_b = nfacet;
  int size_alpha=nfacet2;
  int size_inv_alpha=(1<<nfacet)*nfacet2;
  int size_det_alpha=(1<<nfacet);
  int size_dgdb=rank*dim;
  double work_space[size_b 
                   +size_alpha 
                   +size_inv_alpha 
                   +size_det_alpha
                   +size_dgdb];
  double *w = work_space;
  a = a_loc;           
  b = w;          w += size_b;
  a1 = a_loc;
  b1 = b_loc;
  dadt = NULL;
  dbdt = b_loc;
  alpha = w;      w += size_alpha;
  inv_alpha = w;  w += size_inv_alpha;
  det_alpha = w;  w += size_det_alpha;
  dgdb = w;  w += size_dgdb;

  //printf("dim = %d\n", dim);
  //print_matrix(stdout, dim, nfacet, a1, "The value of matrix a:");
  //print_vector(stdout, nfacet, b1, "The value of vector b:\n");
  //printf("rank = %d\n", rank);

  update_alpha();
  //print_matrix(stdout, nfacet, nfacet, alpha, "The value of matrix alpha:");

  update_inv_alpha();
  //int J;
  //printf("determinant: ");
  //for(J = 0; J < rank; J++){
    //char s[10];
    //sprintf(s, "%d:", J);
    //print_matrix(stdout, nfacet, nfacet, inv_alpha+J*nfacet2, s);
    //printf("%g ", det_alpha[J]);
  //}
  //printf("\n");

  simplex_init(g);
  //print_vector(stdout, rank+1, g, "The initial value of g:\n");

  double factorial_d=1.0;
  for ( i = 0; i < d+1; i++)
    factorial_d *= (i+1);
  for ( i = 0; i < rank; i++)
    g[i] *= factorial_d;
  //print_vector(stdout, rank+1, g, "The initial value of g:\n");

  simplex_rk(g);
  //print_vector(stdout, rank, g, "output of runge-kutta:\n");

  get_prob(g);
  //print_vector(stdout, rank, g, "result:\n");

  printf("probability= %20.18g\n", g[0]);
  return;
}

int 
main(int argc, char *argv[])
{
  char **p = argv;

  int c;
  sscanf(*++p, "%d", &c);

  switch(c){
  case 1:
    simplex(p);
    break;
  case 2:
    simplicial_cone(p);
    break;
  case 3:
    simplex2(p);
    break;
  case 4:
    tukey(p);
    break;
  case 5:
    tukey2(p);
    break;
  case 6:
    tukey3(p);
    break;
  case 11:
    numerical_experiment1(p);
    break;
  case 12:
    numerical_experiment2(p);
    break;
  case 13:
    numerical_experiment3(p);
    break;
  case 14:
    numerical_experiment4(p);
    break;
  default:
    fprintf(stderr, "Error\n");
    exit(EXIT_FAILURE);
    break;
  }
  return 0;
}
