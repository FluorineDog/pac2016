#include "LE_SymSprsMatAux.h"
#include <atomic>
#include <cfloat>
#include <condition_variable>
#include <iostream>
#include <memory>
#include <mutex>
#include <pthread.h>
#include <thread>
// #include <hbwmalloc.h>
using std::cout;
using std::cerr;
using std::endl;

#include <cassert>
#include <map>
#include <set>
#include <vector>

#define dog_calloc(size, type)                                                 \
  (type *)aligned_alloc(64, ((size) * sizeof(type) + 63) & ~63);
#define dog_free(ptr) free(ptr);

constexpr int BLOCK = 1;
constexpr int THREAD_NUM = 0; // magic !! don't modify !!
typedef double __attribute((aligned(64))) aligned_double;
using su_t = SprsUMatRealStru;
void AdditionLU_SymbolicSymG(su_t *pFU) {
  int *rs_u = pFU->uMax.rs_u;
  int *j_u = pFU->uMax.j_u;
  int iDim = pFU->uMax.iDim;
  std::vector<std::set<int>> refTable(iDim + 1);
  for (int i = 1; i < iDim + 1; ++i) {
    int kbeg = rs_u[i];
    int kend = rs_u[i + 1];
    for (int k = kbeg; k < kend; ++k) {
      int j = j_u[k];
      refTable[j].insert(i);
    }
  }
  std::vector<std::set<int>> seqPart(iDim + 1);
  std::vector<std::vector<int>> paraPart(iDim + 1);
  for (int basei = 1; basei < iDim + 1; ++basei) {
    int iend = std::min(basei + BLOCK, iDim + 1);
    for (int i = basei; i < iend; ++i) {
      const auto &ori = refTable[i];
      for (auto x : ori) {
        if (x < basei) {
          paraPart[i].push_back(x);
        } else {
          seqPart[i].insert(x);
          assert(x != i);
          seqPart[i].insert(seqPart[x].begin(), seqPart[x].end());
        }
      }
    }
  }
  int lines_sum = 0;
  for (int i = 1; i < iDim + 1; ++i) {
    lines_sum += (paraPart[i].size() + 7) / 8;
    lines_sum += (seqPart[i].size() + 7) / 8;
  }
  pFU->dogUMat_upper.iDim = iDim;
  pFU->dogUMat_upper.alloc_size = lines_sum * 8;
  dog_free(pFU->dogUMat_upper.columns);
  dog_free(pFU->dogUMat_upper.SeqRanges);
  dog_free(pFU->dogUMat_upper.paraRanges);
  dog_free(pFU->values_upper);

  pFU->dogUMat_upper.SeqRanges = dog_calloc(iDim + 1, pair_ii_t);
  pFU->dogUMat_upper.paraRanges = dog_calloc(iDim + 1, pair_ii_t);
  pFU->dogUMat_upper.columns = dog_calloc(lines_sum * 8, int);
  pFU->values_upper = dog_calloc(lines_sum * 8, double);

  int col_index = 0;
  for (int i = 1; i < iDim + 1; ++i) {
    col_index = (col_index + 7) & ~7;
    pFU->dogUMat_upper.paraRanges[i].beg = col_index;
    for (auto x : paraPart[i]) {
      pFU->dogUMat_upper.columns[col_index++] = x;
    }
    pFU->dogUMat_upper.paraRanges[i].end = col_index;

    col_index = (col_index + 7) & ~7;
    pFU->dogUMat_upper.SeqRanges[i].beg = col_index;
    for (auto x : seqPart[i]) {
      pFU->dogUMat_upper.columns[col_index++] = x;
    }
    pFU->dogUMat_upper.SeqRanges[i].end = col_index;
  }
}

void AdditionLU_NumericSymG(SprsUMatRealStru *pFU) {
  // init
  int *rs_u = pFU->uMax.rs_u;
  int *j_u = pFU->uMax.j_u;
  double *u_u = pFU->u_u;
  auto paraRanges = pFU->dogUMat_upper.paraRanges;
  auto seqRanges = pFU->dogUMat_upper.SeqRanges;
  auto columns = pFU->dogUMat_upper.columns;
  auto values_upper = pFU->values_upper;
  memset(values_upper, 0, pFU->dogUMat_upper.alloc_size * sizeof(double));
  auto iDim = pFU->dogUMat_upper.iDim;
  int debug_i = -1, debug_j = -1;

  std::vector<std::map<int, double>> refTable(iDim + 1);
  for (int i = 1; i < iDim + 1; ++i) {
    int kbeg = rs_u[i];
    int kend = rs_u[i + 1];
    for (int k = kbeg; k < kend; ++k) {
      int j = j_u[k];
      refTable[j][i] = u_u[k];
    }
  }

  for (int basei = iDim; basei > 0; basei -= BLOCK) {
    int iend = std::max(basei - BLOCK, 0);
    // (i,j) -> col_index;
    std::map<int, std::map<int, int>> mapping;
    for (int i = basei; i > iend; i--) {

      for (int k = paraRanges[i].beg; k < paraRanges[i].end; ++k) {
        int j = columns[k];
        mapping[i][j] = k;
      }

      for (int k = seqRanges[i].beg; k < seqRanges[i].end; ++k) {
        int j = columns[k];
        mapping[i][j] = k;
      }

      // int kbeg = rs_u[i];
      // int kend = rs_u[i + 1];
      // for (int k = kbeg; k < kend; k++) {
      for (auto pp:refTable[i]) {
        int j = pp.first;
        int k = mapping[i][j]; 
        int col_id = mapping[i][j];
        __float128 coef = -pp.second;
        values_upper[col_id] += coef;
        if (j < basei) {
          continue;
        }
        for (auto pair : mapping[j]) {
          int cur_j = pair.first;
          if (cur_j < basei) {
            break;
          }
          int cur_col_id = pair.second;
          // god bless me !!!!
          values_upper[mapping[i][cur_j]] += coef * values_upper[cur_col_id];
        }
      }
    }
  }
}
