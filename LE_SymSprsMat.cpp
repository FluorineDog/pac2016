constexpr int BLOCK = 1100;
constexpr int THREAD_NUM = 16; // magic !! don't modify !!
#include "LE_SymSprsMatInit.h"
#pragma float_control(fast, on, push)

static Barrier<THREAD_NUM + 1> barrierhead[2];
static Barrier<THREAD_NUM> barrier[2];

static std::vector<std::thread> threads;
static std::atomic<bool> die;
static aligned_double * __restrict__ glo_b;
static aligned_double * __restrict__ glo_x;

static aligned_double * __restrict__ tempx;

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
    barrierhead[0].Wait();
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
#pragma simd vectorlength(16)
        for (int k = kbeg; k < kend; ++k) {
          int j = columns[k];
          xc += values[k] * glo_x[j];
        }
        tempx[basei - i] = xc;
        // ia = glo_id[0].fetch_sub(1);
        ia -= THREAD_NUM;
      }
      // // barrier
      barrier[0].Wait();
      // this loop is ready for parallel
      // for (int i = basei; i > iend; i--) {
      while (ib > iend) {
        int i = ib;
        double xc = tempx[basei - i];
        int kbeg = seqRanges[i].beg;
        int kend = seqRanges[i].end;
#pragma vector aligned
#pragma ivdep
#pragma simd vectorlength(16)
        for (int k = kbeg; k < kend; ++k) {
          int j = columns[k];
          xc += values[k] * tempx[basei - j];
        }
        glo_x[i] = xc;
        // ib = glo_id[1].fetch_sub(1);
        ib -= THREAD_NUM;
      }
      barrier[1].Wait();
    }
    barrierhead[1].Wait();
  }
}

static void dog_init_task(su_t *U) {
  die.store(false);
  for (int i = 0; i < THREAD_NUM; ++i) {
    auto th = std::thread(workload, U, i);
    threads.push_back(std::move(th));
  }
  tempx = dog_calloc(U->dogUMat.iDim + 1, double);
}

static void finalize_task(su_t *U) {
  die = true;
  if (threads.size()) {
    barrierhead[0].Wait();
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
  double *d_u = pFU->d_u;
  double *u_u = pFU->u_u;
  int *rs_u = pFU->uMax.rs_u;
  int *j_u = pFU->uMax.j_u;
  int iDim = pFU->uMax.iDim;
  auto paraRanges = pFU->dogUMat.paraRanges;
  auto seqRanges = pFU->dogUMat.SeqRanges;
  auto columns = pFU->dogUMat.columns;
  auto values = pFU->values;

  // memcpy(x, b, (iDim + 1) * 8);
  for(int i = 0; i < iDim+1; ++i){
    x[i] = b[i];
  }

  for (int i = 1; i <= iDim; i++) {
    double xc = x[i];
    int ks = rs_u[i];
    int ke = rs_u[i + 1];

    for (int k = ks; k < ke; k++) {
      int j = j_u[k];
      x[j] -= u_u[k] * xc;
      assert(u_u[k] == u_u[k]);
    }
  }

  for (int i = 1; i <= iDim; i++)
    x[i] /= d_u[i];

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
  barrierhead[0].Wait();
  barrierhead[1].Wait();

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