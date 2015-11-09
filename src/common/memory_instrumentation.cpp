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

/*! \file
  \ingroup HMatrix
  \brief Memory Allocation tracking.
*/
#include "memory_instrumentation.hpp"
#include "common/my_assert.h"
#include <algorithm>

#ifdef __GLIBC__
#include <malloc.h>

#if __GLIBC_MINOR__ >= 10
#define HAVE_MALLOC_INFO
#endif
#if defined(HMAT_FORCE_MALLOC_INFO) && !defined(HAVE_MALLOC_INFO)
#define HAVE_MALLOC_INFO
extern "C" {
int __attribute__((weak)) malloc_info(int options, FILE *fp);
}
#endif

#if defined(HMAT_HEAPTOP)
extern "C" {
size_t __attribute__((weak)) heaptop_allocated();
}
#endif

// Do not care about thread safety. This is an acceptable approximation.
static struct mallinfo global_mallinfo;
static size_t system_current, aspace_total, aspace_mprotect;
static int mallinfo_counter;
static int mallinfo_sampling;
static int write_counter;
static int write_sampling;
static int malloc_info_counter;
static int malloc_info_sampling;
#endif

namespace hmat {

static size_t get_res_mem(void *)
{
    size_t resident = 0;
#ifdef __linux__
    FILE * statm_file = fopen("/proc/self/statm", "r");
    fscanf(statm_file, "%*s %lu", &resident);
    fclose(statm_file);
#endif
    return resident * 4096;
}

MemoryInstrumenter::MemoryInstrumenter(): enabled_(false) {
    addType("Time", false);
#if __GNUC__
    addType("FullMatrix", false);
#else
    addType("FullMatrix", true);
#endif
#ifdef __GLIBC__
    // Same as executable maps + arena so not needed when MALLOC_ARENA_MAX=1
    addType("RSS", false, get_res_mem, NULL);

    // Only work export MALLOC_ARENA_MAX=1
    addType("Non-mmapped space (arena)", false);

    //Not needed
    //addType("ordblks", false);
    //addType("smblks", false);
    //addType("hblks", false);

    //Always almost zero with hmat so disabled
    addType("Space in mmapped regions (hblkhd)", false);

    //Not needed
    //addType("usmblks", false);
    //addType("fsmblks", false);
    addType("Total allocated space (uordblks)", false);

    // Same as arena - uordblks
    // addType("Total free space (fordblks)", false);
    addType("Top-most, releasable (keepcost)", false);
#endif
    char * ws = getenv("HMAT_MEMINSTR_WS");
    write_sampling = ws ? atoi(ws) : 1;
    char * mi = getenv("HMAT_MEMINSTR_MI");
    mallinfo_sampling = mi ? atoi(mi) : 100;
    mi = getenv("HMAT_MEMINSTR_MI2");
    malloc_info_sampling = mi ? atoi(mi) : 100;
}

void MemoryInstrumenter::setFile(const std::string & filename) {
#ifdef HMAT_MEM_INSTR
    filename_ = filename;
    output_ = fopen(filename.c_str(), "w+");
    HMAT_ASSERT_MSG(output_ != NULL, "Cannot open %s", filename.c_str());
    start_ = now();
    fullMatrixMem_ = 0;
    mallinfo_counter = 0;

    FILE * labelsf = fopen((filename_+".labels").c_str(), "w");
    for(int i = 0; i < labels_.size(); i++) {
        fputs(labels_[i].c_str(), labelsf);
        fputc('\n', labelsf);
    }
    fclose(labelsf);

    enabled_ = true;
#endif
}

char MemoryInstrumenter::addType(const std::string & label, bool cumul,
                                 HookFunction hook, void * param) {
    HMAT_ASSERT_MSG(output_ == NULL, "Cannot call addType after setFile");
    HMAT_ASSERT_MSG(write_sampling == 1 || !cumul,
                    "Cannot use write sub sampling with cumulative records.");
    cumulatives_.push_back(cumul);
    labels_.push_back(label);
    hooks_.push_back(hook);
    hookParams_.push_back(param);
    return labels_.size() - 1;
}

static int previous_xml_tag(int tag_number, int pos, char * buffer) {
    int nb_gt = 0;
    // 70s way of parsing XML
    for(; pos >= 0; pos--) {
        if(buffer[pos] == '>')
            nb_gt ++;
        if(nb_gt == tag_number)
            break;
    }
    return pos;
}

static size_t parse_xml_tag(int i, char * buffer) {
    i -= 3;
    for(; i >= 0; i--) {
        if(buffer[i] == '"')
            break;
    }
    i++;
    return strtoll(buffer+i, NULL, 10);
}

static void parse_malloc_info() {
#ifdef HAVE_MALLOC_INFO
    int pos = 0;
    char * buffer;
    size_t buffer_size;
    FILE * fb = open_memstream(&buffer, &buffer_size);
    int r = malloc_info(0, fb);
    fclose(fb);
    assert(r == 0);
    pos = previous_xml_tag(2, buffer_size - 1, buffer);
    aspace_mprotect = parse_xml_tag(pos, buffer);

    pos = previous_xml_tag(1, pos, buffer);
    aspace_total = parse_xml_tag(pos, buffer);

    pos = previous_xml_tag(2, pos, buffer);
    system_current = parse_xml_tag(pos, buffer);

    free(buffer);
#endif
}

void MemoryInstrumenter::allocImpl(mem_t size, char type) {
    if(enabled_) {
        std::vector<mem_t> buffer(labels_.size());
        assert(output_ != NULL);
        assert(type < buffer.size() - 1);
        std::fill(buffer.begin(), buffer.end(), 0);
        buffer[0] = nanoTime();
#ifdef __GNUC__
        if(type == FULL_MATRIX) {
            buffer[FULL_MATRIX] = __sync_add_and_fetch(&fullMatrixMem_, size);
        } else
#endif
        if(type > 0)
            buffer[type] = size;

#ifdef __GLIBC__
        mallinfo_counter ++;
        if(mallinfo_counter >= mallinfo_sampling && mallinfo_sampling > 0) {
            global_mallinfo = mallinfo();
            mallinfo_counter = 0;
        }

        malloc_info_counter ++;
        if(malloc_info_counter >= malloc_info_sampling && malloc_info_sampling > 0) {
            malloc_info_counter = 0;
            parse_malloc_info();
        }

        int k = 3;
#ifdef HAVE_MALLOC_INFO
        buffer[k++] = system_current;
#else
        buffer[k++] = global_mallinfo.arena;
#endif

        //buffer[k++] = global_mallinfo.ordblks;
        //buffer[k++] = global_mallinfo.smblks;
        //buffer[k++] = global_mallinfo.hblks;
        buffer[k++] = global_mallinfo.hblkhd;
        //buffer[k++] = global_mallinfo.usmblks;
        //buffer[k++] = global_mallinfo.fsmblks;


#ifdef HMAT_HEAPTOP
        buffer[k++] = heaptop_allocated();
#elif defined(HAVE_MALLOC_INFO)
        // mallinfo does not support value greater than 2GiB so we overwrite
        // the result with a value taken from malloc_info
        buffer[k++] = aspace_total;
#else
        buffer[k++] = global_mallinfo.uordblks;
#endif

        //buffer[k++] = global_mallinfo.fordblks;
        buffer[k++] = global_mallinfo.keepcost;
#endif
        for(int i = 0; i < hooks_.size(); i++) {
            if(hooks_[i]) {
                assert(type != i);
                buffer[i] = hooks_[i](hookParams_[i]);
            }
        }
        assert(buffer[0] > 0);
        write_counter ++;
        if(write_counter >= write_sampling) {
            fwrite(buffer.data(), sizeof(size_t), buffer.size(), output_);
            fflush(output_);
            write_counter = 0;
        }
    }
}

void MemoryInstrumenter::freeImpl(mem_t size, char type) {
    allocImpl(-size, type);
}

void MemoryInstrumenter::finish() {
    if(output_ == NULL)
        return;
    int n = cumulatives_.size();
    std::vector<mem_t> cumulator(n);
    std::vector<mem_t> buffer_(n);
    rewind(output_);
    fpos_t pos;
    if(std::find(cumulatives_.begin(), cumulatives_.end(), true) != cumulatives_.end()) {
        while(true) {
            fgetpos(output_, &pos);
            int r = fread(buffer_.data(), sizeof(size_t) * buffer_.size(), 1, output_);
            if(!r)
                break;
            assert(r == 1);
            assert(buffer_[0] > 0);
            for(int i = 0; i < buffer_.size(); i++) {
                if(cumulatives_[i])
                    cumulator[i] += buffer_[i];
                else
                    cumulator[i] = buffer_[i];
            }
            assert(cumulator[0] > 0);
            fsetpos(output_, &pos);
            r = fwrite(cumulator.data(), sizeof(size_t) * cumulator.size(), 1, output_);
            assert(r == 1);
        }
    }
    fclose(output_);
    output_ = NULL;
    enabled_ = false;
}

MemoryInstrumenter::~MemoryInstrumenter() {
    finish();
}

void MemoryInstrumenter::enable() {
    enabled_ = true;
}

void MemoryInstrumenter::disable() {
    enabled_ = false;
}

size_t MemoryInstrumenter::nanoTime() {
    return time_diff_in_nanos(start_, now());
}
}
