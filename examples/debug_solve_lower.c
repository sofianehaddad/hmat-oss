/*
  HMat-OSS (HMatrix library, open source software)

  Copyright (C) 2014-2015 Airbus Group SAS

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

  http://github.com/jeromerobert/hmat-oss
*/

// Cylinder
#include <stdio.h>
#include <math.h>
#include "hmat/hmat.h"
#include "examples.h"
#include "common/chrono.h"

/** This is a simple example showing how to use the HMatrix library.

    In this example, we assemble and do a decomposition of a Matrix such that:
    \f[A_{ij} = \frac{e^{i\kappa |x_i - x_j|}}{4 \pi |x_i - x_j|}\f]
    with the points \f$(x_i)\f$ on a cylinder.
    In the real case we use 1 / r instead.
 */

/**
  Define interaction between 2 degrees of freedoms  (real case)
  Use the Matern 3/2 kernel with a range parameter rho = 1.359047 (sqrt(3) / theta), which is the value of the range parameter such that the correlation between two points at distance 1 is equal to 0.01.
*/
double interaction_real(double* points, int i, int j)
{
  if (i == j) return 1.0;
  double rho = 1.359047;
  double scaled_dx = rho * fabs(points[i] - points[j]);
  return (1.0 + scaled_dx) * exp(-scaled_dx);
}

inline static void PointsToFile(double *points, int size, int dim, const char *filename)
{
  int i, j;
  FILE *fp = fopen(filename, "w");
  for (i = 0; i < size; i++)
  {
    for (j = 0; j < dim; j++)
      fprintf(fp, "%e ", points[dim * i + j]);
    fprintf(fp, "\n");
    // fprintf(fp, "%e %e %e\n", points[3 * i], points[3 * i + 1], points[3 * i + 2]);
  }
  fclose(fp);
}

/**
  This structure contains data related to our problem.
*/
typedef struct {
  int type;
  int n;
  double* points;
} problem_data_t;

/**
  This structure is a glue between hmat library and application.
*/
typedef struct {
  int row_start;
  int col_start;
  int* row_hmat2client;
  int* col_hmat2client;
  problem_data_t* user_context;
} block_data_t;

/**
  Function to free our block_data_t structure.  As we only store pointers, there is no
*/
void
free_hmat(void *data)
{
  free((block_data_t*) data);
}

/**
  prepare_hmat is called by hmat library to prepare assembly of
  a cluster block.  We allocate a block_data_t, which is then
  passed to the compute_hmat function.  We have to store all
  datas needed by compute_hmat into this block_data_t structure.
*/
void
prepare_hmat(int row_start,
             int row_count,
             int col_start,
             int col_count,
             int *row_hmat2client,
             int *row_client2hmat,
             int *col_hmat2client,
             int *col_client2hmat,
             void *user_context,
             hmat_block_info_t * block_info)
{
  /* Silence C/C++ compiler warnings */
  (void) row_count; (void) col_count; (void) row_client2hmat; (void) col_client2hmat;

  block_info->user_data = calloc(1, sizeof(block_data_t));
  block_info->release_user_data = free_hmat;
  block_data_t* bdata = (block_data_t*) block_info->user_data;

  bdata->row_start = row_start;
  bdata->col_start = col_start;
  bdata->row_hmat2client = row_hmat2client;
  bdata->col_hmat2client = col_hmat2client;
  bdata->user_context = (problem_data_t*) user_context;
}

/**
  Compute all values of a block and store them into an array,
  which had already been allocated by hmat.  There is no
  padding, all values computed in this block are contiguous,
  and stored in a column-major order.
  This block is not necessarily identical to the one previously
  processed by prepare_hmat, it may be a sub-block.  In fact,
  it is either the full block, or a full column, or a full row.
*/
void
compute_hmat(void *data,
             int rowBlockBegin,
             int rowBlockCount,
             int colBlockBegin,
             int colBlockCount,
             void *values)
{
  int i, j;
  double *dValues = (double *) values;
  //float *fValues = (float *)values;
  block_data_t *bdata = (block_data_t *)data;

  int type = bdata->user_context->type;
  int pos = 0;
  for (j = 0; j < colBlockCount; ++j) {
      int col = bdata->col_hmat2client[j + colBlockBegin + bdata->col_start];
      switch (type) {
	  case HMAT_DOUBLE_PRECISION:
	    for (i = 0; i < rowBlockCount; ++i, ++pos)
	      dValues[pos] = interaction_real(bdata->user_context->points, bdata->row_hmat2client[i + rowBlockBegin + bdata->row_start], col);
    	break;
    case HMAT_SIMPLE_PRECISION:
      for (i = 0; i < rowBlockCount; ++i, ++pos)
        dValues[pos] = interaction_real(bdata->user_context->points, bdata->row_hmat2client[i + rowBlockBegin + bdata->row_start], col);
      break;
    default:
      break;
    }
  }
}

int main(int argc, char **argv) {
  double* points;
  hmat_settings_t settings;
  hmat_interface_t hmat;
  hmat_value_t scalar_type;
  hmat_info_t mat_info;
  int i, N, M;
  double tolerance;
  char arithmetic;
  hmat_clustering_algorithm_t* clustering, * clustering_algo;
  hmat_cluster_tree_t* cluster_tree;
  hmat_matrix_t* hmatrix;
  problem_data_t problem_data;
  int isTransposed;


  if (argc != 6)
  {
    fprintf(stdout, "All arguments weren't passed to executable!");
    //  Size of the Matrix in consideration:
    N = 6000;
    // Size of Matrices at leaf level:
    M = 200;
    // Tolerance of problem
    tolerance = pow(10, -12);
    scalar_type = HMAT_DOUBLE_PRECISION;
    isTransposed = 1;
    fprintf(stdout, "Using Default Arguments, solving linear traingular system with transposition = %d\n", isTransposed);
  }

  else
  {
    // Size of the Matrix in consideration:
    N = atoi(argv[1]);
    // Size of Matrices at leaf level:
    M = atoi(argv[2]);
    // Tolerance of problem
    tolerance = pow(10, -atoi(argv[3]));
    arithmetic = argv[4][0];
    switch (arithmetic)
    {
    case 'S':
      scalar_type = HMAT_SIMPLE_PRECISION;
      break;
    case 'D':
      scalar_type = HMAT_DOUBLE_PRECISION;
      break;
    default:
      fprintf(stderr, "Unknown arithmetic code %c, exiting...\n", arithmetic);
      return 1;
    }
    isTransposed = atoi(argv[5]);
    fprintf(stdout, "Transposition = %d\n", isTransposed);
  }

  hmat_init_default_interface(&hmat, scalar_type);
  fprintf(stdout, "Transposition = %d\n", isTransposed);

  hmat_get_parameters(&settings);
  settings.maxLeafSize = M;
  hmat_set_parameters(&settings);
  if (0 != hmat.init())
  {
    fprintf(stderr, "Unable to initialize HMat library\n");
    return 1;
  }

  points = malloc(N * sizeof(double));
  printf("Generating the point cloud...\n");
  for (i = 0; i < N; ++i)
      points[i] = -1.0 + 2.0 * i / (N - 1);
  printf("done.\n");
  PointsToFile(points, N, 1, "points.txt");

  problem_data.n = N;
  problem_data.points = points;
  problem_data.type = scalar_type;

  clustering_algo = hmat_create_clustering_geometric();
  clustering = hmat_create_clustering_max_dof(clustering_algo, M);
  cluster_tree = hmat_create_cluster_tree(points, 1, N, clustering);
  hmat_delete_clustering(clustering);
  hmat_delete_clustering(clustering_algo);
  printf("ClusterTree node count = %d\n", hmat_tree_nodes_count(cluster_tree));
  hmat_admissibility_t *admissibilityCondition = hmat_create_admissibility_standard(3.0);
  hmatrix = hmat.create_empty_hmatrix_admissibility(
      cluster_tree, cluster_tree, 1, admissibilityCondition);
  hmat.set_low_rank_epsilon(hmatrix, tolerance);
  hmat_delete_admissibility(admissibilityCondition);
  hmat.get_info(hmatrix, &mat_info);
  printf("HMatrix node count = %d\n", mat_info.nr_block_clusters);
  Time start = now();
  hmat_assemble_context_t ctx_assemble;
  hmat_assemble_context_init(&ctx_assemble);
  ctx_assemble.compression = hmat_create_compression_aca_random(tolerance);
  ctx_assemble.user_context = &problem_data;
  ctx_assemble.prepare = prepare_hmat;
  ctx_assemble.block_compute = compute_hmat;
  ctx_assemble.lower_symmetric = 1;
  int rc = hmat.assemble_generic(hmatrix, &ctx_assemble);
  if (0 != rc) {
      fprintf(stderr, "Error during assembly, aborting\n");
      exit(rc);
  }
  Time end = now();
  fprintf(stdout, "elapsed time = %f\n", time_diff(start, end));

  hmat_delete_compression(ctx_assemble.compression);

  start = now();
  hmat_factorization_context_t ctx_facto;
  hmat_factorization_context_init(&ctx_facto);
  ctx_facto.factorization = hmat_factorization_llt;
  rc = hmat.factorize_generic(hmatrix, &ctx_facto);
  end = now();
  fprintf(stdout, "elapsed time = %f\n", time_diff(start, end));
   if (0 != rc)
  {
    fprintf(stderr, "Error during factorization, aborting\n");
    exit(rc);
  }

  // solve lineat system
  double* rhs = malloc(N * sizeof(double));
  double* sol = malloc(N * sizeof(double));
  // isTransposed = 0 --> no problem
  // isTransposed = 1 --> bad result
  for (i = 0; i < N; ++i)
  {
    rhs[i] = 1.0;
    sol[i] = 1.0;
  }
  hmat.solve_lower_triangular(hmatrix, isTransposed, sol, 1);

  // check the solution
  // isTransposed, alpha, 
  double* y = malloc(N * sizeof(double));
  double alpha = 1.0;
  double beta = 0.0;
  hmat.gemv(isTransposed, &alpha, hmatrix, sol, &beta, &y[0], 1);
  double error = 0.0;
  for (i = 0; i < N; ++i)
    error += fabs(y[i] - rhs[i]);
  fprintf(stdout, "error = %e\n", error);
  hmat.destroy(hmatrix);
  hmat_delete_cluster_tree(cluster_tree);
  hmat.finalize();
  free(points);
  free(rhs);
  free(sol);
  free(y);
  return 0;
}
