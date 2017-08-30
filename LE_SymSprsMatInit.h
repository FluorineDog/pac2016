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

constexpr int BLOCK = 528;
constexpr int THREAD_NUM = 0; // magic !! don't modify !!
typedef double __attribute((aligned(64))) aligned_double;
using su_t = SprsUMatRealStru;
//
// inline void* __attribute((malloc)) dog_calloc_sub_(const int alignment,
// size_t size){
// size = (size+63)&~63;
// void* ptr;
// int ret = hbw_posix_memalign(&ptr, alignment, size);
// return ptr;
// }

#define dog_calloc(size, type)                                                 \
  (type *)aligned_alloc(64, ((size) * sizeof(type) + 63) & ~63);

#define dog_free(ptr) free(ptr);

// 描    述:          // 稀疏实数向量内存初始化。数目、指针变量置零
void initMem_VecReal(VecRealStru *V) {
  // 内存初始化。数目、指针变量置零
  V->iNy = 0;
  V->pdVal = NULL;
}

// 描    述:          // 实数向量分配维数
// 输入参数:          // VecRealStru *V：实数向量
// 其    他:          // 调用该函数前向量元素数目V->iNy已经确定
//////////////////////////////////////////////////////////////////////
void allocate_VecReal(VecRealStru *V) {
  int m = 0;
  m = V->iNy + 1;
  V->pdVal = dog_calloc(m, double);
  //(double *)calloc(m, sizeof(double));
}

// 描    述:          // 指针变量内存释放
// 被deallocateNet()调用
void deallocate_VecReal(VecRealStru *V) {
  // 指针变量内存释放
  dog_free(V->pdVal);
}
void AdditionLU_SymbolicSymG(SprsUMatRealStru *pFU) {
  int *rs_u = pFU->uMax.rs_u;
  int *j_u = pFU->uMax.j_u;
  int iDim = pFU->uMax.iDim;
  // rs_u ---- store beg and end for each line
  // j_u  ---- colIndex flat storage
  // rowIndex -> [colIndex]
  std::vector<std::set<int>> refTable(iDim + 1);
  for (int i = 1; i < iDim + 1; ++i) {
    int kbeg = rs_u[i];
    int kend = rs_u[i + 1];
    refTable[i].insert(j_u + kbeg, j_u + kend);
  }

  std::vector<std::set<int>> seqPart(iDim + 1);
  std::vector<std::vector<int>> paraPart(iDim + 1);

  for (int basei = iDim; basei > 0; basei -= BLOCK) {
    int iend = std::max(basei - BLOCK, 0);
    for (int i = basei; i > iend; i--) {
      const auto &ori = refTable[i];
      for (auto x : ori) {
        if (x > basei) {
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
  pFU->dogUMat.iDim = iDim;
  pFU->dogUMat.alloc_size = lines_sum * 8;
  dog_free(pFU->dogUMat.columns);
  dog_free(pFU->dogUMat.SeqRanges);
  dog_free(pFU->dogUMat.paraRanges);
  dog_free(pFU->values);
  dog_free(pFU->rd_u);

  pFU->dogUMat.SeqRanges = dog_calloc(iDim + 1, pair_ii_t);
  pFU->dogUMat.paraRanges = dog_calloc(iDim + 1, pair_ii_t);
  pFU->dogUMat.columns = dog_calloc(lines_sum * 8, int);
  pFU->values = dog_calloc(lines_sum * 8, double);
  pFU->rd_u = dog_calloc(iDim + 1, double);

  int col_index = 0;
  for (int i = iDim; i > 0; --i) {
    col_index = (col_index + 7) & ~7;
    pFU->dogUMat.paraRanges[i].beg = col_index;
    for (auto x : paraPart[i]) {
      pFU->dogUMat.columns[col_index++] = x;
    }
    pFU->dogUMat.paraRanges[i].end = col_index;

    col_index = (col_index + 7) & ~7;
    pFU->dogUMat.SeqRanges[i].beg = col_index;
    for (auto x : seqPart[i]) {
      pFU->dogUMat.columns[col_index++] = x;
    }
    pFU->dogUMat.SeqRanges[i].end = col_index;
  }
}

static void dog_init_task(su_t *U);
using std::pair;
using std::make_pair;
void AdditionLU_NumericSymG(SprsUMatRealStru *pFU) {
  // init
  int *rs_u = pFU->uMax.rs_u;
  int *j_u = pFU->uMax.j_u;
  double *u_u = pFU->u_u;
  auto paraRanges = pFU->dogUMat.paraRanges;
  auto seqRanges = pFU->dogUMat.SeqRanges;
  auto columns = pFU->dogUMat.columns;
  auto values = pFU->values;
  memset(values, 0, pFU->dogUMat.alloc_size * sizeof(double));
  auto iDim = pFU->dogUMat.iDim;
  auto rd_u = pFU->rd_u;
  rd_u[0] = 0;
  for (int i = 1; i < iDim + 1; ++i) {
    rd_u[i] = 1 / pFU->d_u[i];
  }
  for (int i = iDim + 1; i < ((iDim + 1 + 7) & ~7); ++i) {
    rd_u[i] = 0;
  }

  int debug_i = -1, debug_j = -1;
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

      int kbeg = rs_u[i];
      int kend = rs_u[i + 1];
      for (int k = kbeg; k < kend; k++) {
        int j = j_u[k];
        int col_id = mapping[i][j];
        __float128 coef = -u_u[k];
        values[col_id] += coef;
        if (j > basei) {
          continue;
        }
        for (auto pair : mapping[j]) {
          int cur_j = pair.first;
          if (cur_j > basei) {
            break;
          }
          int cur_col_id = pair.second;
          // god bless me !!!!
          values[mapping[i][cur_j]] += coef * values[cur_col_id];
        }
      }
    }
  }
  dog_init_task(pFU);
}

template <int n_> class Barrier {
public:
  Barrier() : nwait_(0), step_(0) {}

  bool Wait() {
    unsigned int step = step_;

    if (nwait_.fetch_add(1) == n_ - 1) {
      // OK, last thread to come.
      // step_.fetch_add(1);
      ++step_;
      nwait_.store(0); // XXX: maybe can use relaxed ordering here ??
                       // yes, just use different barriers intersectually
      return true;
    } else {
      // Run in circles and scream like a little girl.
      while (step_ == step)
        ;
      return false;
    }
  }

protected:
  /* Number of synchronized threads. */
  /* Number of threads currently spinning.  */
  std::atomic<unsigned int> nwait_;

  /* Number of barrier syncronizations completed so far,
   * it's OK to wrap.  */
  // std::atomic<unsigned int> step_;
  volatile unsigned int step_;
};

// class dogBarrier {
// public:
//   static constexpr int barrier_thread = 16;
//   dogBarrier() {
//     assert(!THREAD_NUM || THREAD_NUM == barrier_thread);
//     for (auto &line : lock) {
//       for (auto &lock : line) {
//         lock = 1;
//       }
//     }
//   };
//   void Wait(int id) {
//     // cerr << "Damn";
//     lock[0][id ^ 0x1] = 0;
//     lock[1][id ^ 0x2] = 0;
//     lock[2][id ^ 0x4] = 0;
//     lock[3][id ^ 0x8] = 0;
//     // lock[4][id ^ 0x10] = 0;
//     // lock[5][id ^ 0x20] = 0;
//     while (lock[0][id]);
//     while (lock[1][id]);
//     while (lock[2][id]);
//     while (lock[3][id]);
//     // while(lock[4][id]);
//     // while(lock[5][id]);
//     lock[0][id] = 1;
//     lock[1][id] = 1;
//     lock[2][id] = 1;
//     lock[3][id] = 1;
//     // lock[4][id] = 1;
//     // lock[5][id] = 1;
//   }

// protected:
//   // volatile int lock[3][8];
//   volatile int lock[4][THREAD_NUM]; // barrier_thread
//   // volatile int lock[2][THREAD_NUM]; // barrier_thread
//   // volatile int lock[6][64];
// };
