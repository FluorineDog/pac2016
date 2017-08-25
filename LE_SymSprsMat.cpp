
#include "LE_SymSprsMatAux.h"
#include <atomic>
#include <thread>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <iostream> 
using std::cout;
using std::cerr;
using std::endl;
#define dog_calloc(size, type) (type *)aligned_alloc(64, (size) * sizeof(type))
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
  V->pdVal = (double *)calloc(m, sizeof(double));
}

// 描    述:          // 指针变量内存释放
// 被deallocateNet()调用
void deallocate_VecReal(VecRealStru *V) {
  // 指针变量内存释放
  free(V->pdVal);
}
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
  auto paraRange = pFU->dogUMat.paraRanges;
  auto seqRange = pFU->dogUMat.SeqRanges;
  auto columns = pFU->dogUMat.columns;
  auto values = pFU->values;
  // cerr << values;
  memset(values, 0, pFU->dogUMat.alloc_size * sizeof(double));
  // cerr << "fuck" << values[0] << endl;
  auto iDim = pFU->dogUMat.iDim;
  int debug_i = -1, debug_j = -1;
  for (int basei = iDim; basei > 0; basei -= BLOCK) {
    int iend = std::max(basei - BLOCK, 0);
    // (i,j) -> col_index;
    std::map<int, std::map<int, int>> mapping;
    for (int i = basei; i > iend; i--) {
      for (int k = paraRange[i].beg; k < paraRange[i].end; ++k) {
        int j = columns[k];
        mapping[i][j] = k;
      }
      for (int k = seqRange[i].beg; k < seqRange[i].end; ++k) {
        int j = columns[k];
        mapping[i][j] = k;
      }
      int kbeg = rs_u[i];
      int kend = rs_u[i + 1];
      for (int k = kbeg; k < kend; k++) {
        int j = j_u[k];
        int col_id = mapping[i][j];
        // cerr << i << " "<< j;
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
static std::atomic<int> die(0);
static std::atomic<int> readiness(0);
static std::atomic<int> done(0);

constexpr int THREAD_NUM = 1;
static void workload(su_t *U, int private_id) {
  // double *d_u = U->d_u;
  // double *u_u = U->u_u;
  // int *rs_u = U->uMax.rs_u;
  // int *j_u = U->uMax.j_u;
  int iDim = U->uMax.iDim;
  while (true) {
    {
      unique_lock<mutex> ulk(barrier_mutex);
      global_cv.wait(ulk, [&]() -> bool { return readiness || die; });
      --readiness;
    }
    if (die) {
      return;
    }
    int row = glo_id.fetch_sub(1);
    while (row > 0) {
      // glo_b[row] += row;
      row = glo_id.fetch_sub(1);
    }
    ++done;
    while (done != THREAD_NUM) {
    }
    row = glo_id2.fetch_sub(1);
    while (row > 0) {
      // glo_b[row]++;
      row = glo_id2.fetch_sub(1);
    }
    ++done;
    while (done != 2 * THREAD_NUM) {
    }
    
  }
}

static void dog_init_task(su_t *U) {
  die.store(false);
  for (int i = 0; i < THREAD_NUM; ++i) {
    auto th = std::thread(workload, U, i);
    threads.push_back(std::move(th));
  }
}

static void finalize_task(su_t *U) {
  {
    lock_guard<mutex> lck(barrier_mutex);
    die.store(true);
    cerr << "dying" << die;
    // readiness = THREAD_NUM;
  }
  cerr << "recy" << endl;
  global_cv.notify_all();
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
void LE_FBackwardSym(SprsUMatRealStru *pFU, double b[], double x[]) {
  int i, j, k;
  int ks, ke;
  int iDim;
  int *rs_u, *j_u;
  double *d_u, *u_u;
  double xc;

  d_u = pFU->d_u;
  u_u = pFU->u_u;
  rs_u = pFU->uMax.rs_u;
  j_u = pFU->uMax.j_u;
  iDim = pFU->uMax.iDim;
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
  while (done != 2 * THREAD_NUM){}
  for (i = 1; i <= iDim; i++){
    // b[i] -= i + 1;
    x[i] = b[i];    
  }

  for (i = 1; i <= iDim; i++) {
    xc = x[i];
    ks = rs_u[i];
    ke = rs_u[i + 1];

    for (k = ks; k < ke; k++) {
      j = j_u[k];
      x[j] -= u_u[k] * xc;
      assert(u_u[k] == u_u[k]);
    }
  }

  for (i = 1; i <= iDim; i++)
    x[i] /= d_u[i];

  for (i = iDim - 1; i >= 1; i--) {
    ks = rs_u[i];
    ke = rs_u[i + 1] - 1;
    xc = x[i];

    for (k = ke; k >= ks; k--) {
      j = j_u[k];
      xc -= u_u[k] * x[j];
    }
    x[i] = xc;
  }
}

// 描    述:          //内存初始化。数目、指针变量置零
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
  initMem_UMatReal(U);
  finalize_task(U);
}
