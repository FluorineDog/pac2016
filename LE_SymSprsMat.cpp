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
#include <iostream>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tbb/tbb.h>
using std::cout;
using std::cerr;
using std::endl;
#include "LE_SymSprsMatAuxFunc.h"

// 描    述:          // 稀疏实数对称矩阵分配维数
// 输入参数:          // SprsMatRealStru *pA：稀疏实数对阵矩阵或上三角矩阵
// 其    他:          // 调用该函数前稀疏实数对称矩阵的维数A->Mat.iDim已经确定
// 因可能还有注入元引入，稀疏实数对阵矩阵A的元素实际数目不确定，故分配维数时一次将维数扩够
void allocate_MatReal(SprsMatRealStru *A) {
  int dimension = 0;
  int n = 0;
  MKL_INT opt = MKL_DSS_DEFAULTS;

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
  MKL_INT opt = MKL_DSS_DEFAULTS;

  V->pdVal = (double *)calloc(m, sizeof(double));
}

//////////////////////////////////////////////////////////////////////
// 描    述:          // 指针变量内存释放
// 被deallocateNet()调用
//////////////////////////////////////////////////////////////////////
void deallocate_MatReal(SprsMatRealStru *A) {
  free(A->Mat.piJno);
  free(A->Mat.piIstart);
  free(A->Mat.piIdiag);
  free(A->Mat.piLinkn);
  free(A->Mat.piLinkp);
}

// 函 数 名:          // deallocate_VecReal
// 描    述:          // 指针变量内存释放
// 被deallocateNet()调用
// 输入参数:          // VecRealStru *V：
void deallocate_VecReal(VecRealStru *V) {
  // 指针变量内存释放
  free(V->pdVal);
}

#include <map>


#include <cassert>
#include <set>
#include <vector>
constexpr int BLOCK = 50;
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
  int sumWork = 0;
  int seqWork = 0;
  int paraWork = 0;
  using Tp = std::tuple<int, int, int>;
  std::map<Tp, int> fuck;
  for (int i = iDim; i > 0; --i) {
    // cerr << "para: " << paraPart[i].size();
    // cerr << "   seq:" << seqPart[i].size();
    // cerr << endl;
    ++fuck[std::make_tuple(seqPart[i].size() + paraPart[i].size(), seqPart[i].size(), paraPart[i].size())];
    seqWork += seqPart[i].size();
    paraWork += paraPart[i].size();
  }
  int sumCount = 0;
  for(auto pp:fuck){
    cerr << "\t\tsumCount:" << (sumCount+=pp.second);
    cerr << "\t\tsum:" << std::get<0>(pp.first);
    cerr << "\t\tseq:" << std::get<1>(pp.first);
    cerr << "\t\tpara:" << std::get<2>(pp.first);
    cerr << "\t\tcount:" << pp.second;
    cerr << endl;
  }
  sumWork = seqWork + paraWork;

  cerr << endl;
  cerr << "paraBlock = " << BLOCK << endl;
  cerr << "seqWork = " << seqWork << endl;
  cerr << "paraWork = " << paraWork << endl;
  cerr << "sumWork = " << sumWork << endl;
  cerr << endl;

  cout << "paraBlock = " << BLOCK << endl;
  cout << "seqWork = " << seqWork << endl;
  cout << "paraWork = " << paraWork << endl;
  cout << "sumWork = " << sumWork << endl;
  cout << endl;
}

void AdditionLU_NumericSymG(SprsUMatRealStru *pFU) {}

// 输入参数:          // U阵结构及U阵值，右端项b
// 输出参数:          // 右端项x，维数为pU的维数（解向量）
void LE_FBackwardSym(SprsUMatRealStru *pFU, double *__restrict__ b,
                     double *__restrict__ x) {
  int iDim;
  int *rs_u, *j_u;
  double *d_u, *u_u;

  // factorized A = L*D*L^T
  // solving         A x = b
  // i.e.      L*D*L^T x = b
  d_u = pFU->d_u;
  u_u = pFU->u_u;
  rs_u = pFU->uMax.rs_u;
  j_u = pFU->uMax.j_u;
  iDim = pFU->uMax.iDim;
  // let x0 = b
  for (int i = 1; i <= iDim; i++)
    x[i] = b[i];

  // Solve L^T x1 = x0
  // for (int i = 1; i <= iDim; i++) {
  //   double xc = x[i];
  //   int kbeg = rs_u[i];
  //   int kend = rs_u[i + 1];

  //   for (int k = kbeg; k < kend; k++) {
  //     int j = j_u[k];
  //     x[j] -= u_u[k] * xc;
  //   }
  // }

  // Solve D x2 = x1
  for (int i = 1; i <= iDim; i++) {
    x[i] *= d_u[i]; // auto vectorized
  }

  // Solve L x3 = x2
  // x = x3
  for (int i = iDim - 1; i >= 1; i--) {
    int kbeg = rs_u[i];
    int kend = rs_u[i + 1] - 1;
    double xc = x[i];

    for (int k = kend; k >= kbeg; k--) {
      int j = j_u[k];
      xc -= u_u[k] * x[j];
    }
    x[i] = xc;
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
}
