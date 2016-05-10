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

#ifndef _C_WRAPPING_HPP
#define _C_WRAPPING_HPP

#include <string>
#include <cstring>

#include "common/context.hpp"
#include "common/my_assert.h"
#include "full_matrix.hpp"
#include "h_matrix.hpp"
#include "uncompressed_values.hpp"

namespace
{
template<typename T, template <typename> class E>
hmat_matrix_t * create_empty_hmatrix(hmat_cluster_tree_t* rows_tree, hmat_cluster_tree_t* cols_tree, int lower_sym)
{
    hmat::SymmetryFlag sym = (lower_sym ? hmat::kLowerSymmetric : hmat::kNotSymmetric);
    return (hmat_matrix_t*) new hmat::HMatInterface<T, E>(
            (hmat::ClusterTree*)rows_tree,
            (hmat::ClusterTree*)cols_tree, sym);
}

template<typename T, template <typename> class E>
hmat_matrix_t * create_empty_hmatrix_admissibility(
  hmat_cluster_tree_t* rows_tree,
  hmat_cluster_tree_t* cols_tree, int lower_sym,
  hmat_admissibility_t* condition)
{
    hmat::SymmetryFlag sym = lower_sym ? hmat::kLowerSymmetric : hmat::kNotSymmetric;
    return (hmat_matrix_t*) new hmat::HMatInterface<T, E>(
            static_cast<hmat::ClusterTree*>(static_cast<void*>(rows_tree)),
            static_cast<hmat::ClusterTree*>(static_cast<void*>(cols_tree)),
            sym, (hmat::AdmissibilityCondition*)condition);
}


template<typename T>
class SimpleCAssemblyFunction : public hmat::SimpleAssemblyFunction<T> {
private:
  hmat_interaction_func_t functor;
  void* functor_extra_args;

  /** Constructor.

      \param _mat The FullMatrix<T> the values are taken from.
   */
public:
  SimpleCAssemblyFunction(void* user_context, hmat_interaction_func_t &f)
    : hmat::SimpleAssemblyFunction<T>(), functor(f), functor_extra_args(user_context) {}

  typename hmat::Types<T>::dp interaction(int i, int j) const {
    typename hmat::Types<T>::dp result;
    (*functor)(functor_extra_args, i, j, &result);
    return result;
  }
};

template<typename T, template <typename> class E>
void assemble_generic(hmat_matrix_t* matrix, hmat_assemble_context_t * ctx) {
    DECLARE_CONTEXT;
    hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*)matrix;
    bool assembleOnly = ctx->factorization == hmat_factorization_none;
    hmat::SymmetryFlag sf = ctx->lower_symmetric ? hmat::kLowerSymmetric : hmat::kNotSymmetric;
    if(ctx->assembly != NULL) {
        HMAT_ASSERT(ctx->block_compute == NULL && ctx->simple_compute == NULL);
        hmat::Assembly<T> * cppAssembly = (hmat::Assembly<T> *)ctx->assembly;
        hmat->assemble(*cppAssembly, sf, assembleOnly, ctx->progress);
        if(!assembleOnly)
            hmat->factorize(ctx->factorization, ctx->progress);
    } else if(ctx->block_compute != NULL) {
        HMAT_ASSERT(ctx->simple_compute == NULL && ctx->assembly == NULL);
        hmat::BlockAssemblyFunction<T> f(hmat->rows(), hmat->cols(), ctx->user_context,
                                         ctx->prepare,  ctx->block_compute);
        hmat->assemble(f, sf, assembleOnly, ctx->progress);
        if(!assembleOnly)
            hmat->factorize(ctx->factorization, ctx->progress);
    } else {
        HMAT_ASSERT(ctx->block_compute == NULL && ctx->assembly == NULL);
        SimpleCAssemblyFunction<T> f(ctx->user_context, ctx->simple_compute);
        hmat->assemble(f, sf, assembleOnly, ctx->progress);
        if(!assembleOnly)
            hmat->factorize(ctx->factorization, ctx->progress);
    }
}

template<typename T, template <typename> class E>
int assemble(hmat_matrix_t * holder,
                              void* user_context,
                              hmat_prepare_func_t prepare,
                              hmat_compute_func_t compute,
                              int lower_symmetric) {
    hmat_assemble_context_t ctx;
    hmat_assemble_context_init(&ctx);
    ctx.user_context = user_context;
    ctx.prepare = prepare;
    ctx.block_compute = compute;
    ctx.lower_symmetric = lower_symmetric;
    assemble_generic<T, E>(holder, &ctx);
    return 0;
}

template<typename T, template <typename> class E>
int assemble_factor(hmat_matrix_t * holder,
                              void* user_context,
                              hmat_prepare_func_t prepare,
                              hmat_compute_func_t compute,
                              int lower_symmetric, hmat_factorization_t f_type) {
    hmat_assemble_context_t ctx;
    hmat_assemble_context_init(&ctx);
    ctx.user_context = user_context;
    ctx.prepare = prepare;
    ctx.block_compute = compute;
    ctx.lower_symmetric = lower_symmetric;
    ctx.factorization = f_type;
    assemble_generic<T, E>(holder, &ctx);
    return 0;
}

template<typename T, template <typename> class E>
int assemble_simple_interaction(hmat_matrix_t * holder,
                          void* user_context,
                          hmat_interaction_func_t compute,
                          int lower_symmetric) {
    hmat_assemble_context_t ctx;
    hmat_assemble_context_init(&ctx);
    ctx.user_context = user_context;
    ctx.simple_compute = compute;
    ctx.lower_symmetric = lower_symmetric;
    assemble_generic<T, E>(holder, &ctx);
    return 0;
}

template<typename T, template <typename> class E>
hmat_matrix_t* copy(hmat_matrix_t* holder) {
  DECLARE_CONTEXT;
  return (hmat_matrix_t*) ((hmat::HMatInterface<T, E>*) holder)->copy();
}

template<typename T, template <typename> class E>
int destroy(hmat_matrix_t* holder) {
  delete (hmat::HMatInterface<T, E>*)(holder);
  return 0;
}

template<typename T, template <typename> class E>
void factorize_generic(hmat_matrix_t* holder, hmat_factorization_context_t * ctx) {
    DECLARE_CONTEXT;
    hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*) holder;
    hmat->factorize(ctx->factorization, ctx->progress);
}

template<typename T, template <typename> class E>
int factor(hmat_matrix_t* holder, hmat_factorization_t t) {
    hmat_factorization_context_t ctx;
    hmat_factorization_context_init(&ctx);
    ctx.factorization = t;
    factorize_generic<T, E>(holder, &ctx);
    return 0;
}

template<typename T, template <typename> class E>
int finalize() {
  hmat::HMatInterface<T, E>::finalize();
  return 0;
}

template<typename T, template <typename> class E>
int full_gemm(char transA, char transB, int mc, int nc, void* c,
                             void* alpha, void* a, hmat_matrix_t * holder, void* beta) {
  DECLARE_CONTEXT;

  const hmat::HMatInterface<T, E>* b = (hmat::HMatInterface<T, E>*)holder;
  hmat::FullMatrix<T> matC((T*)c, mc, nc);
  hmat::FullMatrix<T>* matA = NULL;
  if (transA == 'N') {
    matA = new hmat::FullMatrix<T>((T*)a, mc, transB == 'N' ? b->rows()->size()
                             : b->cols()->size());
  } else {
    matA = new hmat::FullMatrix<T>((T*)a, transB == 'N' ? b->rows()->size()
                             : b->cols()->size(), mc);
  }
  hmat::HMatInterface<T, E>::gemm(matC, transA, transB, *((T*)alpha), *matA, *b, *((T*)beta));
  delete matA;
  return 0;
}

template<typename T, template <typename> class E>
int gemm(char trans_a, char trans_b, void *alpha, hmat_matrix_t * holder,
                   hmat_matrix_t * holder_b, void *beta, hmat_matrix_t * holder_c) {
  DECLARE_CONTEXT;
  hmat::HMatInterface<T, E>* hmat_a = (hmat::HMatInterface<T, E>*)holder;
  hmat::HMatInterface<T, E>* hmat_b = (hmat::HMatInterface<T, E>*)holder_b;
  hmat::HMatInterface<T, E>* hmat_c = (hmat::HMatInterface<T, E>*)holder_c;
  hmat_c->gemm(trans_a, trans_b, *((T*)alpha), hmat_a, hmat_b, *((T*)beta));
  return 0;
}

template<typename T, template <typename> class E>
int gemv(char trans_a, void* alpha, hmat_matrix_t * holder, void* vec_b,
                   void* beta, void* vec_c, int nrhs) {
  DECLARE_CONTEXT;
  hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*)holder;
  const hmat::ClusterData* bData = (trans_a == 'N' ? hmat->cols(): hmat->rows());
  const hmat::ClusterData* cData = (trans_a == 'N' ? hmat->rows(): hmat->cols());
  hmat::FullMatrix<T> mb((T*) vec_b, bData->size(), nrhs);
  hmat::FullMatrix<T> mc((T*) vec_c, cData->size(), nrhs);
  hmat->gemv(trans_a, *((T*)alpha), mb, *((T*)beta), mc);
  return 0;
}

template<typename T, template <typename> class E>
int add_identity(hmat_matrix_t* holder, void *alpha) {
  DECLARE_CONTEXT;
  hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*)holder;
  hmat->addIdentity(*((T*)alpha));
  return 0;
}

template<typename T, template <typename> class E>
int init() {
    return hmat::HMatInterface<T, E>::init();
}

template<typename T, template <typename> class E>
double norm(hmat_matrix_t* holder) {
  DECLARE_CONTEXT;
  return ((hmat::HMatInterface<T, E>*)holder)->norm();
}

template<typename T, template <typename> class E>
int scale(void *alpha, hmat_matrix_t* holder) {
  DECLARE_CONTEXT;
  ((hmat::HMatInterface<T, E>*)holder)->scale(*((T*)alpha));
  return 0;
}

template<typename T, template <typename> class E>
int solve_mat(hmat_matrix_t* hmat, hmat_matrix_t* hmatB) {
  ((hmat::HMatInterface<T, E>*)hmat)->solve(*(hmat::HMatInterface<T, E>*)hmatB);
    return 0;
}

template<typename T, template <typename> class E>
int solve_systems(hmat_matrix_t* holder, void* b, int nrhs) {
  DECLARE_CONTEXT;
  hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*)holder;
  hmat::FullMatrix<T> mb((T*) b, hmat->cols()->size(), nrhs);
  hmat->solve(mb);
  return 0;
}

template<typename T, template <typename> class E>
int transpose(hmat_matrix_t* hmat) {
  DECLARE_CONTEXT;
  ((hmat::HMatInterface<T, E>*)hmat)->transpose();
  return 0;
}

template<typename T, template <typename> class E>
int hmat_get_info(hmat_matrix_t* holder, hmat_info_t* info) {
  DECLARE_CONTEXT;
  hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*) holder;
  hmat->info(*info);
  return 0;
}

template<typename T, template <typename> class E>
int hmat_dump_info(hmat_matrix_t* holder, char* prefix) {
  DECLARE_CONTEXT;
  hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*) holder;
  std::string fileps(prefix);
  fileps += ".ps";
  std::string filejson(prefix);
  filejson += ".json";
  hmat->createPostcriptFile( fileps );
  hmat->dumpTreeToFile( filejson );
  return 0;
}


template<typename T, template <typename> class E>
int set_cluster_trees(hmat_matrix_t* holder, hmat_cluster_tree_t * rows, hmat_cluster_tree_t * cols) {
  hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*) holder;
  hmat->engine().hmat->setClusterTrees((hmat::ClusterTree*)rows, (hmat::ClusterTree*)cols);
  return 0;
}

template<typename T, template <typename> class E>
int extract_diagonal(hmat_matrix_t* holder, void* diag, int size)
{
  DECLARE_CONTEXT;
  (void)size; //for API compatibility
  hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*) holder;
  hmat->engine().hmat->extractDiagonal(static_cast<T*>(diag));
  hmat::FullMatrix<T> permutedDiagonal(static_cast<T*>(diag), hmat->cols()->size(), 1);
  hmat::restoreVectorOrder(&permutedDiagonal, hmat->cols()->indices());
  return 0;
}

template<typename T, template <typename> class E>
int solve_lower_triangular(hmat_matrix_t* holder, int transpose, void* b, int nrhs)
{
  DECLARE_CONTEXT;
  hmat::HMatInterface<T, E>* hmat = (hmat::HMatInterface<T, E>*)holder;
  hmat::FullMatrix<T> mb((T*) b, hmat->cols()->size(), nrhs);
  hmat->solveLower(mb, transpose);
  return 0;
}

template <typename T, template <typename> class E>
int get_block(struct hmat_get_values_context_t *ctx) {
    hmat::HMatInterface<T, E> *hmat = (hmat::HMatInterface<T, E> *)ctx->matrix;
    hmat::IndexSet rows(ctx->row_offset, ctx->row_size);
    hmat::IndexSet cols(ctx->col_offset, ctx->col_size);
    typename E<T>::UncompressedBlock view;
    view.uncompress(hmat->engine().data(), rows, cols, (T*)ctx->values);
    if (ctx->renumber_rows)
        view.renumberRows();
    ctx->col_indices = view.colsNumbering();
    ctx->row_indices= view.rowsNumbering();
    return 0;
}

template <typename T, template <typename> class E>
int get_values(struct hmat_get_values_context_t *ctx) {
    hmat::HMatInterface<T, E> *hmat = (hmat::HMatInterface<T, E> *)ctx->matrix;
    typename E<T>::UncompressedValues view;
    view.uncompress(hmat->engine().data(),
                    ctx->row_indices, ctx->row_size,
                    ctx->col_indices, ctx->col_size,
                    (T*)ctx->values);
    return 0;
}

template <typename T, template <typename> class E>
int walk(hmat_matrix_t* holder, hmat_procedure_t* proc) {
    hmat::HMatInterface<T, E> *hmat = (hmat::HMatInterface<T, E> *) holder;
    hmat::TreeProcedure *functor = (hmat::TreeProcedure *) proc;
    hmat->walk(functor);
    return 0;
}

}  // end anonymous namespace

namespace hmat {

template<typename T, template <typename> class E>
static void createCInterface(hmat_interface_t * i)
{
    i->assemble = assemble<T, E>;
    i->assemble_factorize = assemble_factor<T, E>;
    i->assemble_simple_interaction = assemble_simple_interaction<T, E>;
    i->copy = copy<T, E>;
    i->create_empty_hmatrix = create_empty_hmatrix<T, E>;
    i->create_empty_hmatrix_admissibility = create_empty_hmatrix_admissibility<T, E>;
    i->destroy = destroy<T, E>;
    i->factorize = factor<T, E>;
    i->finalize = finalize<T, E>;
    i->full_gemm = full_gemm<T, E>;
    i->gemm = gemm<T, E>;
    i->gemv = gemv<T, E>;
    i->add_identity = add_identity<T, E>;
    i->init = init<T, E>;
    i->norm = norm<T, E>;
    i->scale = scale<T, E>;
    i->solve_mat = solve_mat<T, E>;
    i->solve_systems = solve_systems<T, E>;
    i->transpose = transpose<T, E>;
    i->internal = NULL;
    i->get_info  = hmat_get_info<T, E>;
    i->dump_info = hmat_dump_info<T, E>;
    i->set_cluster_trees = set_cluster_trees<T, E>;
    i->extract_diagonal = extract_diagonal<T, E>;
    i->solve_lower_triangular = solve_lower_triangular<T, E>;
    i->assemble_generic = assemble_generic<T, E>;
    i->factorize_generic = factorize_generic<T, E>;
    i->get_values = get_values<T, E>;
    i->get_block = get_block<T, E>;
    i->walk = walk<T, E>;
}

}  // end namespace hmat

#endif  // _C_WRAPPING_HPP
