/*
 * This source code is the property of Airbus Group S.A.S. No part of it shall
 * be reproduced or transmitted without the express prior written authorization
 * of Airbus Group S.A.S., and its contents shall not be disclosed.
 * Copyright Airbus Group S.A.S.
 */

#include <stdio.h>
#include <math.h>

#include <hmat/hmat.h>
#include <string.h> /* strcmp */
#include "common/chrono.h"

#ifdef _MSC_VER
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

int global =0;
typedef struct {
  int n;
  int dim;
  double* points;
  double l;
} problem_data_t;

/**
  Define interaction between 2 degrees of freedoms  (real case)
 */
void interaction_real(void* data, int i, int j, void* result)
{
  problem_data_t* pdata = (problem_data_t*) data;
  double* points = pdata->points;
  int dim = pdata->dim;
  int k;
  double r = 0.0;
  if (i == j)
  {
    *((double *)result) = 1.0;
  }
  else
  {
    for (k = 0; k < dim; ++k)
    {
      r += (points[dim * i + k] - points[dim * j + k]) * (points[dim * i + k] - points[dim * j + k]);
    }
    *((double *)result) = exp(-r / pdata->l);
  }
  global = global + 1;
}

/**
  Define interaction between 2 degrees of freedoms (real case) using Matern 5/2
 */
void interaction_real_matern_5_2(void* data, int i, int j, void* result)
{
  problem_data_t* pdata = (problem_data_t*) data;
  double* points = pdata->points;
  int dim = pdata->dim;
  int k;
  double r = 0.0;
  double prod = 1.0;
  if (i == j)
  {
    *((double *)result) = 1.0;
  }
  else
  {
    for (k = 0; k < dim; ++k)
    {
      // r over rho
      r = fabs(points[dim * i + k] - points[dim * j + k]) / pdata->l;
      prod *= ((1.0 + sqrt(5) * r + 5/3 * (r * r)) * exp(-sqrt(5) * r));
    }
    *((double *)result) = prod;
  }
  global = global + 1;
}

/**
  Define interaction between 2 degrees of freedoms (real case) using Matern 3/2
 */
void interaction_real_matern_3_2(void* data, int i, int j, void* result)
{
  problem_data_t* pdata = (problem_data_t*) data;
  double* points = pdata->points;
  int dim = pdata->dim;
  int k;
  double r = 0.0;
  double prod = 1.0;
  double sqrt3 = sqrt(3);
  if (i == j)
  {
    *((double *)result) = 1.0;
  }
  else
  {
    for (k = 0; k < dim; ++k)
    {
      // r over rho
      r = fabs(points[dim * i + k] - points[dim * j + k]) / pdata->l;
      prod *= ((1.0 + sqrt3 * r) * exp(-sqrt3 * r));
    }
    *((double *)result) = prod;
  }
  global = global + 1;
}

void generate_data(void* data)
{
  problem_data_t* pdata = (problem_data_t*) data;
  int dim = pdata->dim;
  int N = pdata->n;
  double *points = malloc(N * dim * sizeof(double));
  int i;
  if (dim == 1)
  {
    for (i = 0; i < N; ++i) points[i] = -1.0 + 2.0 * i / (N-1);
  }

  if (dim == 2)
  {
    int sqrt_N = sqrt(N);
    int k = 0;
    for (int i = 0; i < sqrt_N; ++i)
    {
      for (int j = 0; j < sqrt_N; j++)
      {
        points[k] = -1.0 + 2.0 * i /(sqrt_N-1);
        points[k + 1] = -1.0 + 2.0 * j /(sqrt_N-1);
        k = k + 2;
      }
    }
  }
  pdata->points = points;
}

int main(int argc, char **argv) {

  int N, M, dim;
  int kLowerSymmetric = 1; /* =0 if not Symmetric */
  double tolerance, l;
  char arithmetic;

  hmat_value_t type;
  hmat_interface_t hmat;
  hmat_settings_t settings;
  hmat_info_t mat_info;
  hmat_clustering_algorithm_t *clustering;
  hmat_cluster_tree_t *cluster_tree;
  hmat_matrix_t *hmatrix;
  problem_data_t problem_data;

  if (argc != 6)
  {
    // fprintf(stdout, "All arguments weren't passed to executable!");
    // fprintf(stdout, "Using Default Arguments" );
    //  Size of the Matrix in consideration:
    N = 900;
    // Size of Matrices at leaf level:
    M = 250;
    // Dimensionality of the problem:
    dim = 2;
    // Tolerance of problem
    tolerance = pow(10, -12);
    type = HMAT_DOUBLE_PRECISION;
    }

    else
    {
      // Size of the Matrix in consideration:
      N          = atoi(argv[1]);
      // Size of Matrices at leaf level:
      M          = atoi(argv[2]);
      // Dimensionality of the problem:
      dim        = atoi(argv[3]);
      // Tolerance of problem
      tolerance  = pow(10, -atoi(argv[4]));
      arithmetic = argv[5][0];
      switch (arithmetic) 
      {
        case 'S':
          type = HMAT_SIMPLE_PRECISION;
          break;
        case 'D':
          type = HMAT_DOUBLE_PRECISION;
          break;
        default:
          fprintf(stderr, "Unknown arithmetic code %c, exiting...\n", arithmetic);
          return 1;
       }
    }
    hmat_init_default_interface(&hmat, type);

    hmat_get_parameters(&settings);
    settings.maxLeafSize = M;
    hmat_set_parameters(&settings);
    if (0 != hmat.init())
    {
      fprintf(stderr, "Unable to initialize HMat library\n");
      return 1;
  }


  l = 1.0;
  printf("correlationLength = %le\n", l);
  problem_data.n = N;
  problem_data.dim = dim;
  generate_data(&problem_data);
  //problem_data.points = points;
  problem_data.l = l;

  clustering = hmat_create_clustering_geometric();
  // Set MaxDof
  hmat_clustering_algorithm_t *algodof = hmat_create_clustering_max_dof(clustering, M);
  cluster_tree = hmat_create_cluster_tree(problem_data.points, dim, 1 * N, clustering);
  printf("ClusterTree node count = %d\n", hmat_tree_nodes_count(cluster_tree));
  hmat_admissibility_t * admissibilityCondition = hmat_create_admissibility_hodlr();
  hmatrix = hmat.create_empty_hmatrix_admissibility(cluster_tree, cluster_tree,
                                                    kLowerSymmetric, admissibilityCondition);
  hmat.set_low_rank_epsilon(hmatrix, tolerance);
  hmat_delete_admissibility(admissibilityCondition);
  hmat_delete_clustering(algodof);

  hmat.get_info(hmatrix, &mat_info);
  printf("HMatrix node count = %d\n", mat_info.nr_block_clusters);
  hmat.dump_info(hmatrix, "rbf");

  fprintf(stdout,"Assembly...");
  Time start = now();
  hmat_assemble_context_t ctx_assemble;
  hmat_assemble_context_init(&ctx_assemble);
  ctx_assemble.compression = hmat_create_compression_aca_full(tolerance);
  ctx_assemble.user_context = &problem_data;
  ctx_assemble.simple_compute = interaction_real_matern_5_2;
  ctx_assemble.lower_symmetric = kLowerSymmetric;
  hmat.assemble_generic(hmatrix, &ctx_assemble);
  hmat_delete_compression(ctx_assemble.compression);
  Time end = now();
  fprintf(stdout, "done.\n");
  fprintf(stdout, "elapsed time = %f\n", time_diff(start, end));

  fprintf(stdout,"Factorisation...");
  start = now();
  hmat_factorization_context_t ctx_facto;
  hmat_factorization_context_init(&ctx_facto);
  ctx_facto.factorization = hmat_factorization_hodlrsym;
  hmat.factorize_generic(hmatrix, &ctx_facto);
  end = now();
  fprintf(stdout, "done.\n");
  fprintf(stdout, "elapsed time = %f\n", time_diff(start, end));
  fprintf(stdout, "Number of total calls to exp function = %d\n", global);

  free(problem_data.points);
  hmat_delete_cluster_tree(cluster_tree);
  hmat.finalize();
  return 0;

}


