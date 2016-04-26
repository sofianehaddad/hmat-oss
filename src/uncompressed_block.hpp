/*
  HMat-OSS (HMatrix library, open source software)

  Copyright (C) 2014-2016 Airbus Group SAS

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
#pragma once
#include "h_matrix.hpp"
#include "cluster_tree.hpp"
#include "rk_matrix.hpp"
#include "common/my_assert.h"

namespace hmat {

/**
 * Base class for uncompressing a block
 * T is the scalar type, M is the matrix type (HMatrix),
 * I this class sub type when subclassing (CRTP)
 * Subclasses must provide: isLeaf, finish, matrix, getValues (see UncompressedBlock)
 */
template <typename T, template <typename> class M, typename I> class UncompressedBlockBase {

    void getValuesImpl() {
        if (rowIndexSet_.size() == 0 || colIndexSet_.size() == 0)
            return;
        if (me()->isLeaf()) {
            me()->getValues();
        } else {
            for (int i = 0; i < matrix_->nbChild(); i++) {
                I view;
                view.init(me());
                view.uncompress(matrix_->getChild(i), *me());
            }
        }
    }

    // https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
    I* me() {
        return static_cast<I*>(this);
    }

    void uncompress(const M<T> * matrix, const UncompressedBlockBase & o) {
        matrix_ = matrix;
        rowIndexSet_.intersection(o.rowIndexSet_, *me()->matrix().rows());
        colIndexSet_.intersection(o.colIndexSet_, *me()->matrix().cols());
        lDim_ = o.ld();
        int drow = rowIndexSet_.offset() - o.rowIndexSet_.offset();
        int dcol = colIndexSet_.offset() - o.colIndexSet_.offset();
        values_ = o.values_ + drow + ((size_t)lDim_) * dcol;
        getValuesImpl();
    }

  protected:
    const M<T> * matrix_;
    IndexSet rowIndexSet_, colIndexSet_;
    // TODO replace by a FullMatrix ?
    T *values_;
    int lDim_;

  public:

    /** The rows of the block to decompress */
    const IndexSet & rows() const {
        return rowIndexSet_;
    }

    /** The columns of the block to decompress */
    const IndexSet & cols() const {
        return colIndexSet_;
    }

    /** The leading dimension of the target buffer */
    int ld() const {
        return lDim_;
    }

    /** Actually uncompress */
    void uncompress(const M<T> * matrix, const IndexSet & rows, const IndexSet & cols,
                    T *values, int ld = -1) {
        matrix_ = matrix;
        rowIndexSet_ = rows;
        colIndexSet_ = cols;
        values_ = values;
        lDim_ = ld;
        if(lDim_ == -1)
            lDim_ = rows.size();
        me()->init(NULL);
        getValuesImpl();
        me()->finish();
    }

    /** The column numbering of the uncompressed block */
    int * colsNumbering() {
        return me()->matrix().cols()->indices() + colIndexSet_.offset();
    }

    /** The row numbering of the uncompressed block */
    int * rowsNumbering() {
        return me()->matrix().rows()->indices() + rowIndexSet_.offset();
    }

    /** Renumber rows. Only available for full heigh blocks. */
    void renumberRows() {
        HMAT_ASSERT_MSG(matrix_->father == NULL && rowIndexSet_ == *me()->matrix().rows(),
                        "Cannot renumber");
        FullMatrix<T> fm(values_, rowIndexSet_.size(), colIndexSet_.size(), ld());
        restoreVectorOrder(&fm, rowsNumbering());
    }
};

/**
 * Specialisation of UncompressedBlockBase for HMatrix
 */
template <typename T> class UncompressedBlock:
        public UncompressedBlockBase<T, HMatrix, UncompressedBlock<T> > {

    void getNullValues() {
        T *toFill = this->values_;
        for (int col = 0; col < this->colIndexSet_.size(); col++) {
            for (int row = 0; row < this->rowIndexSet_.size(); row++) {
                toFill[row] = Constants<T>::zero;
            }
            toFill += this->ld();
        }
    }

    void getFullValues() {
        int nr = this->rows().size();
        int nc = this->cols().size();
        FullMatrix<T> target(this->values_, nr, nc, this->ld());
        int localRowOffset = this->rowIndexSet_.offset() - matrix().rows()->offset();
        int localColOffset = this->colIndexSet_.offset() - matrix().cols()->offset();
        assert(localRowOffset >= 0);
        assert(localColOffset >= 0);
        T *sa = matrix().full()->m + localRowOffset +
                ((size_t)matrix().full()->lda) * localColOffset;
        FullMatrix<T> source(sa, nr, nc, matrix().full()->lda);
        target.copyMatrixAtOffset(&source, 0, 0);
    }

    void getRkValues() {
        // TODO may be we could ask the caller to provide a clean target array
        // to avoid this call.
        getNullValues();
        int nr = this->rowIndexSet_.size();
        int nc = this->colIndexSet_.size();
        int k = matrix().rank();
        FullMatrix<T> result(this->values_, nr, nc, this->ld());
        T *sa = matrix().rk()->a->m + this->rowIndexSet_.offset() - matrix().rows()->offset();
        FullMatrix<T> a(sa, nr, k, matrix().rows()->size());
        T *sb = matrix().rk()->b->m + this->colIndexSet_.offset() - matrix().cols()->offset();
        FullMatrix<T> b(sb, nc, k, matrix().cols()->size());
        result.gemm('N', 'T', Constants<T>::pone, &a, &b, Constants<T>::zero);
    }

  public:
    void init(UncompressedBlock *) {}

    bool isLeaf() {
        return matrix().isLeaf();
    }

    void finish(){}

    const HMatrix<T> & matrix() const {
        return *this->matrix_;
    }

    void getValues() {
         if (matrix().isNull()) {
            getNullValues();
        } else if (matrix().isRkMatrix()) {
            getRkValues();
        } else if (matrix().isFullMatrix()) {
            getFullValues();
        } else {
            assert(false);
        }
    }
};
}