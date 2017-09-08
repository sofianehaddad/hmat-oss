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

#ifndef CUSTOM_TRUNCATE_HPP_
#define CUSTOM_TRUNCATE_HPP_

#include <vector>
#include <algorithm>
#include <utility>

namespace hmat {
template<typename T> class RkMatrix;
template<typename T> class ScalarArray;
template<typename T> class Vector;

template<typename T> void mGSTruncate(RkMatrix<T> *Rk, double prec);
}// end namespace hmat

#endif /* CUSTOM_TRUNCATE_HPP_ */