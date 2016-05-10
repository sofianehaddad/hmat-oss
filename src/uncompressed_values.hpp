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

#include "h_matrix.hpp"
#include "cluster_tree.hpp"
#include <limits>
#include <algorithm>

namespace hmat {


template <typename T, template <typename> class M, typename I> class UncompressedValuesBase {
  protected:
    const M<T> * matrix_;
    T *values_;
    int valuesLd_;
    typedef std::vector<std::pair<int, int> >::iterator IndiceIt;
    /**
     * The first element of the pair a number of element to get, using hmat numbering.
     * The second element of the pair is the position of the element in the original query.
     */
    IndiceIt rowStart_, rowEnd_, colStart_, colEnd_;

    void compatibleQuery(const IndexSet & clusterData, IndiceIt & begin, IndiceIt & end) {
        int lb = clusterData.offset();
        int ub = lb + clusterData.size() - 1;
        std::pair<int, int> lbP(lb, 0);
        // use max int to ensure that upper_bound->first will be greater than ubP
        std::pair<int, int> ubP(ub, std::numeric_limits<int>::max());
        IndiceIt newBegin = std::lower_bound(begin, end, lbP);
        if(newBegin == end) {
            // empty intersection
            begin = newBegin;
            return;
        }
        assert(newBegin->first >= lb);
        IndiceIt newEnd = std::upper_bound(begin, end, ubP);
        assert((newEnd-1)->first <= ub);
        begin = newBegin;
        end = newEnd;
    }

    /**
     * @brief createQuery Convert the C API query to a vector<pair<>> where each pair is
     * <hmat id, original query id>
     * @param query, querySize the C API query
     * @param indices the result
     */
    void createQuery(const ClusterData & clusterData, int * query, int querySize, std::vector<std::pair<int, int> > & indices) {
        indices.resize(querySize);
        for(int i = 0; i < querySize; i++) {
            indices[i].first = clusterData.indices_rev()[query[i] - 1];
            indices[i].second = i;
        }
        std::sort(indices.begin(), indices.end());
    }

    void getValues() {
        if (rowStart_ == rowEnd_ || colStart_ == colEnd_)
            return;
        if (!matrix_->isLeaf()) {
            for (int i = 0; i < matrix_->nrChild(); i++) {
                I view;
                M<T> * child = matrix_->getChild(i);
                view.matrix_ = child;
                view.values_ = values_;
                view.valuesLd_ = valuesLd_;
                view.compatibleQuery(*child->rows(), rowStart_, rowEnd_);
                view.compatibleQuery(*child->cols(), colStart_, colEnd_);
                view.getValues();
            }
        } else
            static_cast<I*>(this)->getLeafValues();
    }
  public:
    /** Actually uncompress */
    void uncompress(const M<T> * matrix, int * rows, int rowSize, int * cols, int colSize, T *values, int ld = -1)
    {
        assert(matrix->father == NULL);
        matrix_ = matrix;
        values_ = values;
        valuesLd_ = ld == -1 ? rowSize : ld;
        std::vector<std::pair<int, int> > rowsIndices, colsIndices;
        createQuery(*matrix->rows(), rows, rowSize, rowsIndices);
        rowStart_ = rowsIndices.begin();
        rowEnd_ = rowsIndices.end();
        createQuery(*matrix->cols(), cols, colSize, colsIndices);
        colStart_ = colsIndices.begin();
        colEnd_ = colsIndices.end();
        getValues();
    }
};

template <typename T> class UncompressedValues: public UncompressedValuesBase<T, HMatrix, UncompressedValues<T> > {
    typedef std::vector<std::pair<int, int> >::iterator IndiceIt;
    friend class UncompressedValuesBase<T, HMatrix, UncompressedValues<T> >;
    void getValue(IndiceIt r, IndiceIt c, T v) {
        this->values_[r->second + ((size_t)this->valuesLd_) * c->second] = v;
    }

    void getNullValues() {
        for(IndiceIt r = this->rowStart_; r != this->rowEnd_; ++r) {
            for(IndiceIt c = this->colStart_; c != this->colEnd_; ++c) {
                getValue(r, c, Constants<T>::zero);
            }
        }
    }

    void getFullValues() {
        const HMatrix<T> & m = *this->matrix_;
        // Check for not supported cases
        assert(m.full()->pivots == NULL);
        assert(m.full()->diagonal == NULL);
        int ro = m.rows()->offset();
        int co = m.cols()->offset();
        for(IndiceIt r = this->rowStart_; r != this->rowEnd_; ++r) {
            for(IndiceIt c = this->colStart_; c != this->colEnd_; ++c) {
                getValue(r, c, m.full()->get(r->first - ro, c->first - co));
            }
        }
    }

    void getRkValues();

    void getLeafValues() {
        if (this->matrix_->isNull()) {
            getNullValues();
        } else if (this->matrix_->isRkMatrix()) {
            getRkValues();
        } else if (this->matrix_->isFullMatrix()) {
            getFullValues();
        } else {
            assert(false);
        }
    }
};
}
