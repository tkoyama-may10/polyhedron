#include<stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define RK_ACCURACY 1e-08
#define GSL_ODEIV_STEP_TYPE gsl_odeiv_step_rkf45

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

static void update_alpha(void);
static void update_inv_alpha(void);
static int inv_submat_alpha(int i, double *inv_submat);
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

  int info, n, m;
  n = m = size_of_submat;
  dpotrf_("U", &n, mat, &m, &info);
  if (info != 0) {
    printf("info:I=%d\n", I);
    return info;
  }

  //  /*
  det_alpha[I] = 1.0;
  for ( k = 0; k < size_of_submat; k++)
    det_alpha[I] *= mat[k*(size_of_submat+1)];
  // */

  n = m = size_of_submat;
  dpotri_("U", &n, mat, &m, &info);
  if (info != 0) {
    printf("info:I=%d\n", I);
    return info;
  }

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
  case 11:
    numerical_experiment1(p);
    break;
  case 12:
    numerical_experiment2(p);
    break;
  case 13:
    numerical_experiment3(p);
    break;
  default:
    fprintf(stderr, "Error\n");
    exit(EXIT_FAILURE);
    break;
  }
  return 0;
}
