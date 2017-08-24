//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, 中国电力科学研究院
// 文 件 名: LE_SymSprsMat.c
// 版本: V1.0       时间: 2014.12.31
// 描    述: 线性方程组模块
//           包含LE_SymSprsMatFunc.h声明的所有子程序
// 其    他：用 #include <filename.h>
// 格式来引用标准库的头文件（编译器将从标准库目录开始搜索）
//           用 #include "filename.h"
//           格式来引用非标准库的头文件（编译器将从用户的工作目录开始搜索）
// 修改记录:       // 历史修改记录
// 1. 时间:
//    作者:
//    修改内容:
// 2. ...
//////////////////////////////////////////////////////////////////////
#include <atomic>
#include <condition_variable>
#include <iostream>
#include <math.h>
#include <memory>
#include <mutex>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using std::cout;
using std::cerr;
using std::endl;
#include "LE_SymSprsMatAuxFunc.h"
#define dog_calloc(size, type) (type *)aligned_alloc(64, (size) * sizeof(type))
// 描    述:          // 稀疏实数对称矩阵分配维数
// 输入参数:          // SprsMatRealStru *pA：稀疏实数对阵矩阵或上三角矩阵
// 其    他:          // 调用该函数前稀疏实数对称矩阵的维数A->Mat.iDim已经确定
// 因可能还有注入元引入，稀疏实数对阵矩阵A的元素实际数目不确定，故分配维数时一次将维数扩够
void allocate_MatReal(SprsMatRealStru *A) {
  int dimension = 0;
  int n = 0;

  dimension = A->Mat.iDim + 1;
  A->Mat.iNymax = dimension * (dimension + 1) / 2; // 矩阵元素最大数目
  n = A->Mat.iNymax + 1;
  // Guess: use rank 1 as the first element instead of rank 0
  // Fuck mathematicians
  A->Mat.piJno = (int *)calloc(n, sizeof(int));
  A->Mat.piIstart = (int *)calloc(dimension + 1, sizeof(int));
  A->Mat.piIdiag = (int *)calloc(dimension, sizeof(int));
  A->Mat.piLinkp = (int *)calloc(n, sizeof(int));
  A->Mat.piLinkn = (int *)calloc(n, sizeof(int));
  A->pdVal = (double *)calloc(n, sizeof(double));
}

// 描    述:          // 实数向量分配维数
// 输入参数:          // VecRealStru *V：实数向量
// 其    他:          // 调用该函数前向量元素数目V->iNy已经确定
void allocate_VecReal(VecRealStru *V) {
  int m = 0;
  m = V->iNy + 1;

  V->pdVal = (double *)calloc(m, sizeof(double));
}

// 描    述:          // 指针变量内存释放
// 被deallocateNet()调用
void deallocate_MatReal(SprsMatRealStru *A) {
  free(A->Mat.piJno);
  free(A->Mat.piIstart);
  free(A->Mat.piIdiag);
  free(A->Mat.piLinkn);
  free(A->Mat.piLinkp);
}

// 描    述:          // 指针变量内存释放
// 被deallocateNet()调用
// 输入参数:          // VecRealStru *V：
void deallocate_VecReal(VecRealStru *V) {
  // 指针变量内存释放
  free(V->pdVal);
}
// // CPFSO
// class Barrier {
// private:
//   std::mutex _mutex;
//   std::condition_variable _cv;
//   std::size_t _count;

// public:
//   explicit Barrier(std::size_t count) : _count{count} {}
//   void Wait() {
//     std::unique_lock<std::mutex> lock{_mutex};
//     if (--_count == 0) {
//       _cv.notify_all();
//     } else {
//       _cv.wait(lock, [this] { return _count == 0; });
//     }
//   }
// };

#include <cassert>
#include <map>
#include <set>
#include <vector>
constexpr int BLOCK = 1;
using su_t = SprsUMatRealStru;
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
  // cerr <<"fuck"<< refTable[iDim-2].size() << "ed" << endl;

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
  free(pFU->dogUMat.columns);
  free(pFU->dogUMat.SeqRanges);
  free(pFU->dogUMat.paraRanges);

  pFU->dogUMat.SeqRanges = dog_calloc(iDim + 1, pair_ii_t);
  pFU->dogUMat.paraRanges = dog_calloc(iDim + 1, pair_ii_t);
  pFU->dogUMat.columns = dog_calloc(lines_sum * 8, int);
  pFU->values = dog_calloc(lines_sum * 8, double);

  int col_index = 0;
  for (int i = iDim; i > 0; --i) {
    col_index = (col_index + 7) &~7;
    pFU->dogUMat.paraRanges[i].beg = col_index;
    for (auto x : paraPart[i]) {
      pFU->dogUMat.columns[col_index++] = x;
    }
    pFU->dogUMat.paraRanges[i].end = col_index;

    col_index = (col_index + 7) &~7;
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
  auto paraRange = pFU->dogUMat.paraRanges;
  auto seqRange = pFU->dogUMat.SeqRanges;
  auto columns = pFU->dogUMat.columns;
  auto values = pFU->values;
  // cerr << values;
  memset(values, 0, pFU->dogUMat.alloc_size);
  auto iDim = pFU->dogUMat.iDim;
  for (int basei = iDim; basei > 0; basei -= BLOCK) {
    int iend = std::max(basei - BLOCK, 0);
    // (i,j) -> col_index;
    std::map<int, std::map<int, int>> mapping;
    for (int i = basei; i > iend; i--) {
      int kbeg = rs_u[i];
      int kend = rs_u[i + 1];
      for (int k = paraRange[i].beg; k < paraRange[i].end; ++k) {
        int j = j_u[k];
        mapping[i][j] = columns[k];
      }
      for (int k = seqRange[i].beg; k < seqRange[i].end; ++k) {
        int j = j_u[k];
        mapping[i][j] = columns[k];
      }
      for (int k = kbeg; k < kend; k++) {
        int j = j_u[k];
        int col_id = mapping[i][j];
        double coef = -u_u[k];
        values[col_id] += coef;
        if (j > basei) {
          continue;
        }
        for (auto pair : mapping[j]) {
          int cur_j = pair.first;
          int cur_col_id = pair.second;
          // god bless me !!!!
          values[mapping[i][cur_j]] += coef * values[cur_col_id];
        }
      }
    }
  }
  dog_init_task(pFU);
}

using std::mutex;
using std::atomic;
using std::unique_lock;
using std::lock_guard;

static std::vector<std::thread> threads;
static std::condition_variable global_cv;
static std::mutex barrier_mutex;
static std::atomic<int> glo_id;
static std::atomic<int> glo_id2;
static double *glo_b;
static double *glo_x;
static std::atomic<bool> die(false);
static std::atomic<int> readiness(0);
static std::atomic<int> done(0);

constexpr int THREAD_NUM = 2;
static void workload(su_t *U, int private_id) {
  // double *d_u = U->d_u;
  // double *u_u = U->u_u;
  // int *rs_u = U->uMax.rs_u;
  // int *j_u = U->uMax.j_u;
  int iDim = U->uMax.iDim;
  while (true) {
    {
      unique_lock<mutex> ulk(barrier_mutex);
      global_cv.wait(ulk, [=]() -> bool { return readiness || die; });
      --readiness;
    }
    if (die) {
      return;
    }
    int row = glo_id.fetch_sub(1);
    while (row > 0) {
      glo_b[row] += row;
      row = glo_id.fetch_sub(1);
    }
    ++done;
    while (done == THREAD_NUM) {

    }
    row = glo_id2.fetch_sub(1);
    while (row > 0) {
      glo_b[row]++;
      row = glo_id2.fetch_sub(1);
    }
    ++done;
    while (done == 2 * THREAD_NUM) {
    }
  }
}

static void dog_init_task(su_t *U) {
  for (int i = 0; i < THREAD_NUM; ++i) {
    auto th = std::thread(workload, U, i);
    threads.push_back(std::move(th));
  }
}

static void finalize_task(su_t *U) {
  {
    lock_guard<mutex> lck(barrier_mutex);
    die.store(true);
    // readiness = THREAD_NUM;
  }
  cerr << "recy" << endl;
  global_cv.notify_all();
  for (auto &th : threads) {
    cerr << "recyclcing" << endl;
    th.join();
    cerr << "recyclced" << endl;
  }
}

// 输入参数:          // U阵结构及U阵值，右端项b
// 输出参数:          // 右端项x，维数为pU的维数（解向量）j
void LE_FBackwardSym(SprsUMatRealStru *pFU, double *__restrict__ b,
                     double *__restrict__ x) {

  // factorized A = (U^T)*D*L
  // solving         A x = b
  // i.e.      (U^T)*D*U x = b
  double *d_u = pFU->d_u;
  double *u_u = pFU->u_u;
  int *rs_u = pFU->uMax.rs_u;
  int *j_u = pFU->uMax.j_u;
  int iDim = pFU->uMax.iDim;
  auto paraRange = pFU->dogUMat.paraRanges;
  auto seqRange = pFU->dogUMat.SeqRanges;
  auto columns = pFU->dogUMat.columns;
  auto values = pFU->values;
  {
    lock_guard<mutex> lck(barrier_mutex);
    glo_b = b;
    glo_x = b;
    glo_id.store(iDim);
    glo_id2.store(iDim);
    done.store(0);
    readiness.store(THREAD_NUM);
  }
  global_cv.notify_all();

  while (done == 2 * THREAD_NUM)
    ;
  // let x0 = b
  for (int i = 1; i <= iDim; i++) {
    x[i] = b[i] - 1 - i;
  }

  // Solve U^T x1 = x0
  for (int i = 1; i <= iDim; i++) {
    double xc = x[i];
    int kbeg = rs_u[i];
    int kend = rs_u[i + 1];

    for (int k = kbeg; k < kend; k++) {
      int j = j_u[k];
      x[j] -= u_u[k] * xc;
    }
  }

  // Solve D x2 = x1
  for (int i = 1; i <= iDim; i++) {
    x[i] *= d_u[i]; // auto vectorized
  }

  // Solve U x3 = x2
  // x = x3
  // for (int i = iDim - 1; i >= 1; i--) {
  //   int kbeg = rs_u[i];
  //   int kend = rs_u[i + 1];
  //   double xc = x[i];
  //   // for (int k = kend - 1; k >= kbeg; k--) {
  //   for(int k = kbeg; k < kend; ++k) {
  //     int j = j_u[k];
  //     xc -= u_u[k] * x[j];
  //   }
  //   x[i] = xc;
  // }

  double tempx[BLOCK+1];
  for (int basei = iDim; basei > 0; basei -= BLOCK) {
    int iend = std::max(basei - BLOCK, 0);
    // para first
    for (int i = basei; i > iend; i--) {
      double xc = x[i];
      int kbeg = paraRange[i].beg;
      int kend = paraRange[i].end;
      int _kbeg = rs_u[i];
      int _kend = rs_u[i+1];
      // cerr << kend -kbeg <<  "-";
      // cerr << _kend -_kbeg << endl;
      // assert(kend-kbeg == _kend - _kbeg); 
      for (int k = kbeg; k < kend; ++k) {
        int _k = _kbeg++;
        int _j = j_u[_k];

        int j = columns[k];
        cerr << "k=" << k << endl;
        cerr << "row=" << i << endl;
        cerr << _j << "*"<< j << endl;
        xc += values[k] * x[j];
        assert(_j == j);
        cerr << values[k] << "--" << d_u[k] << endl;
        assert(values[k] == d_u[k]);
      }
      tempx[basei-i] = xc;
    } 
    // seq next
  //   for (int i = basei; i > iend; i--) {
  //     double xc = tempx[basei-i];
  //     int kbeg = seqRange[i].beg;
  //     int kend = seqRange[i].end;
  //     for (int k = kbeg; k < kend; ++k) {
  //       int j = columns[k];
  //       xc += values[k] * tempx[basei-j];
  //     }
  //     x[i] = xc;
  //   }
  }
}

// 描    述:          //内存初始化。数目、指针变量置零
// 输入参数:          // U
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

// 描    述:          //指针变量内存释放
// 输入参数:          // U
void deallocate_UMatReal(SprsUMatRealStru *U) {
  free(U->d_u);
  free(U->u_u);
  free(U->nzs);
  free(U->work);
  free(U->uMax.cs_u);
  free(U->uMax.rs_u);
  free(U->uMax.r_u);
  free(U->uMax.j_u);
  initMem_UMatReal(U);
  finalize_task(U);
}
