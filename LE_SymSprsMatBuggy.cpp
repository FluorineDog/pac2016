
#include "LE_SymSprsMatAux.h"
#include <atomic>
#include <cfloat>
#include <condition_variable>
#include <iostream>
#include <memory>
#include <mutex>
#include <pthread.h>
#include <thread>
using std::cout;
using std::cerr;
using std::endl;

#include <cassert>
#include <map>
#include <set>
#include <vector>

typedef double __attribute((aligned(64))) aligned_double;
constexpr int BLOCK = 1100;
constexpr int THREAD_NUM = 16; // magic !! don't modify !!
using su_t = SprsUMatRealStru;
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

  pFU->dogUMat.SeqRanges = dog_calloc(iDim + 1, pair_ii_t);
  pFU->dogUMat.paraRanges = dog_calloc(iDim + 1, pair_ii_t);
  pFU->dogUMat.columns = dog_calloc(lines_sum * 8, int);
  pFU->values = dog_calloc(lines_sum * 8, double);

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

#pragma float_control(fast, on, push)
using std::unique_lock;
using std::lock_guard;

class dogBarrier {
public:
  static constexpr int barrier_thread = 16;
  dogBarrier() {
    assert(!THREAD_NUM || THREAD_NUM == barrier_thread);
    for (auto &line : lock) {
      for (auto &lock : line) {
        lock = 1;
      }
    }
  };
  void Wait(int id) {
    // cerr << "Damn";
    lock[0][id ^ 0x1] = 0;
    lock[1][id ^ 0x2] = 0;
    lock[2][id ^ 0x4] = 0;
    lock[3][id ^ 0x8] = 0;
    // lock[4][id ^ 0x10] = 0;
    // lock[5][id ^ 0x20] = 0;
    while (lock[0][id])
      ;
    while (lock[1][id])
      ;
    while (lock[2][id])
      ;
    while (lock[3][id])
      ;
    // while(lock[4][id]);
    // while(lock[5][id]);
    lock[0][id] = 1;
    lock[1][id] = 1;
    lock[2][id] = 1;
    lock[3][id] = 1;
    // lock[4][id] = 1;
    // lock[5][id] = 1;
  }

protected:
  // volatile int lock[3][8];
  volatile int lock[4][THREAD_NUM]; // barrier_thread
  // volatile int lock[2][THREAD_NUM]; // barrier_thread
  // volatile int lock[6][64];
};

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
// static Barrier<THREAD_NUM + 1> barrierhead[2];
// static Barrier<THREAD_NUM> barrier[50];
// static Barrier<THREAD_NUM> barrier[2];
static dogBarrier barrier[3];

static std::vector<std::thread> threads;
// static std::atomic<int> glo_id[2];
static std::atomic<bool> die;
static aligned_double *glo_b;
static aligned_double *glo_x;

static aligned_double *tempx;

// static std::atomic<int> count[2];
static __attribute__((hot)) void workload(const su_t *__restrict pFU,
                                          int private_id) {
  // double *d_u = U->d_u;
  // double *u_u = U->u_u;
  // int *rs_u = U->uMax.rs_u;
  // int *j_u = U->uMax.j_u;
  int iDim = pFU->uMax.iDim;
  const auto paraRanges = pFU->dogUMat.paraRanges;
  const auto seqRanges = pFU->dogUMat.SeqRanges;
  const auto columns = pFU->dogUMat.columns;
  const auto values = pFU->values;
  while (true) {
    cerr << "damn";
    barrier[2].Wait(private_id);
    cerr << "die";
    if (die) {
      return;
    }
    int ia = iDim - private_id;
    int ib = iDim - private_id;
    for (int basei = iDim; basei > 0; basei -= BLOCK) {
      int iend = std::max(basei - BLOCK, 0);
      // this loop is ready for parallel
      // for (int i = basei; i > iend; i--)
      while (ia > iend) {
        int i = ia;
        double xc = glo_x[i];
        int kbeg = paraRanges[i].beg;
        int kend = paraRanges[i].end;
#pragma vector aligned
#pragma ivdep
        for (int k = kbeg; k < kend; ++k) {
          int j = columns[k];
          xc += values[k] * glo_x[j];
        }
        tempx[basei - i] = xc;
        // ia = glo_id[0].fetch_sub(1);
        ia -= THREAD_NUM;
      }
      // // barrier
      barrier[0].Wait(private_id);
      // this loop is ready for parallel
      // for (int i = basei; i > iend; i--) {
      while (ib > iend) {
        int i = ib;
        double xc = tempx[basei - i];
        int kbeg = seqRanges[i].beg;
        int kend = seqRanges[i].end;
#pragma vector aligned
#pragma ivdep
        for (int k = kbeg; k < kend; ++k) {
          int j = columns[k];
          xc += values[k] * tempx[basei - j];
        }
        glo_x[i] = xc;
        // ib = glo_id[1].fetch_sub(1);
        ib -= THREAD_NUM;
      }
      barrier[1].Wait(private_id);
    }
    // barrierhead[1].Wait();
  }
}

static void dog_init_task(su_t *U) {
  die.store(false);
  // 1 is for master mode
  for (int i = 1; i < THREAD_NUM; ++i) {
    auto th = std::thread(workload, U, i);
    threads.push_back(std::move(th));
  }
  tempx = dog_calloc(U->dogUMat.iDim + 1, double);
}

static void finalize_task(su_t *U) {
  die = true;
  if (threads.size()) {
    barrier[2].Wait(0);
  }
  for (auto &th : threads) {
    cerr << "recyclcing" << endl;
    cerr << "dying" << die;
    th.join();
    cerr << "recyclced" << endl;
  }
  cerr << "recy done" << endl;
}

// 描    述:          //对称矩阵前推回代方法求解方程组
// 输入参数:          // U阵结构及U阵值，右端项b
// 输出参数:          // 右端项x，维数为pU的维数（解向量）
// 其    他:          // 不会影响U阵中矩阵的任何信息,created by xdc 2014/6/16

void __attribute__((hot))
LE_FBackwardSym(SprsUMatRealStru *pFU, aligned_double *b, aligned_double *x) {
  // int iDim;
  // int *rs_u, *j_u;
  // double *d_u, *u_u;
  // double *d_u = pFU->d_u;
  // double *u_u = pFU->u_u;
  // int *rs_u = pFU->uMax.rs_u;
  // int *j_u = pFU->uMax.j_u;
  int iDim = pFU->uMax.iDim;
  auto paraRanges = pFU->dogUMat.paraRanges;
  auto seqRanges = pFU->dogUMat.SeqRanges;
  auto columns = pFU->dogUMat.columns;
  auto values = pFU->values;

  // for (int i = 1; i <= iDim; i++) {
  // b[i] -= i + 1;
  // x[i] = b[i];
  // }
  memcpy(x, b, (iDim + 1) * 8);

  // for (int i = 1; i <= iDim; i++) {
  //   double xc = x[i];
  //   int ks = rs_u[i];
  //   int ke = rs_u[i + 1];

  //   for (int k = ks; k < ke; k++) {
  //     int j = j_u[k];
  //     x[j] -= u_u[k] * xc;
  //     assert(u_u[k] == u_u[k]);
  //   }
  // }

  // for (int i = 1; i <= iDim; i++)
  //   x[i] /= d_u[i];

  // for (int i = iDim - 1; i >= 1; i--) {
  //   int ks = rs_u[i];
  //   int ke = rs_u[i + 1] - 1;
  //   double xc = x[i];

  //   for (int k = ke; k >= ks; k--) {
  //     int j = j_u[k];
  //     xc -= u_u[k] * x[j];
  //   }
  //   x[i] = xc;
  // }

  glo_x = x;
  glo_b = b;
  // count[0] = iDim;
  // count[1] = iDim;
  // barrierhead[0].Wait();
  // barrierhead[1].Wait();
  // cerr << "fuck";
  constexpr int private_id = 0;
  cerr << "fuck";
  // while (true) {
  barrier[2].Wait(private_id);
  int ia = iDim - private_id;
  int ib = iDim - private_id;
  for (int basei = iDim; basei > 0; basei -= BLOCK) {
    int iend = std::max(basei - BLOCK, 0);
    // this loop is ready for parallel
    // for (int i = basei; i > iend; i--)
    while (ia > iend) {
      int i = ia;
      double xc = glo_x[i];
      int kbeg = paraRanges[i].beg;
      int kend = paraRanges[i].end;
#pragma vector aligned
#pragma ivdep
      for (int k = kbeg; k < kend; ++k) {
        int j = columns[k];
        xc += values[k] * glo_x[j];
      }
      tempx[basei - i] = xc;
      // ia = glo_id[0].fetch_sub(1);
      ia -= THREAD_NUM;
    }
    // // barrier
    barrier[0].Wait(private_id);
    // this loop is ready for parallel
    // for (int i = basei; i > iend; i--) {
    while (ib > iend) {
      int i = ib;
      double xc = tempx[basei - i];
      int kbeg = seqRanges[i].beg;
      int kend = seqRanges[i].end;
#pragma vector aligned
#pragma ivdep
      for (int k = kbeg; k < kend; ++k) {
        int j = columns[k];
        xc += values[k] * tempx[basei - j];
      }
      glo_x[i] = xc;
      // ib = glo_id[1].fetch_sub(1);
      ib -= THREAD_NUM;
    }
    barrier[1].Wait(private_id);
  }
  // barrierhead[1].Wait();
  // }

  //   double tempx[BLOCK + 1];
  //   for (int basei = iDim; basei > 0; basei -= BLOCK) {
  //     int iend = std::max(basei - BLOCK, 0);
  //     // this loop is ready for parallel
  //     for (int i = basei; i > iend; i--) {
  //       double xc = x[i];
  //       int kbeg = paraRanges[i].beg;
  //       int kend = paraRanges[i].end;
  // // #pragma simd vectorlength(8), reduction(+ : xc)
  // #pragma vector aligned
  // #pragma ivdep
  //       for (int k = kbeg; k < kend; ++k) {
  //         int j = columns[k];
  //         xc += values[k] * x[j];
  //       }
  //       tempx[basei - i] = xc;
  //     }
  //     // barrier
  //     // this loop is ready for parallel
  //     for (int i = basei; i > iend; i--) {
  //       double xc = tempx[basei - i];
  //       int kbeg = seqRanges[i].beg;
  //       int kend = seqRanges[i].end;
  // // #pragma simd vectorlength(8), reduction(+:xc)
  // #pragma vector aligned
  // #pragma ivdep
  //       for (int k = kbeg; k < kend; ++k) {
  //         int j = columns[k];
  //         xc += values[k] * tempx[basei - j];
  //       }
  //       x[i] = xc;
  //     }
  //   }
}

// 描    述:
//内存初始化。数目、指针变量置零
void initMem_UMatReal(SprsUMatRealStru *U) {
  U->d_u = NULL;
  U->u_u = NULL;
  U->nzs = NULL;
  U->work = NULL;
  U->uMax.cs_u = NULL;
  U->uMax.rs_u = NULL;
  U->uMax.r_u = NULL;
  U->uMax.j_u = NULL;
  U->uMax.iDim = 0;
  U->uMax.iNzs = 0;
  // dog implement
  U->dogUMat.columns = nullptr;
  U->dogUMat.SeqRanges = nullptr;
  U->dogUMat.paraRanges = nullptr;
  U->values = nullptr;
}

// 函 数 名:          //deallocate_UMatReal
// 描    述:          //指针变量内存释放
void deallocate_UMatReal(SprsUMatRealStru *U) {
  free(U->d_u);
  free(U->u_u);
  free(U->nzs);
  free(U->work);
  free(U->uMax.cs_u);
  free(U->uMax.rs_u);
  free(U->uMax.r_u);
  free(U->uMax.j_u);
  dog_free(U->dogUMat.columns);
  dog_free(U->dogUMat.SeqRanges);
  dog_free(U->dogUMat.paraRanges);
  dog_free(U->values);
  initMem_UMatReal(U);
  finalize_task(U);
}

#pragma float_control(pop)