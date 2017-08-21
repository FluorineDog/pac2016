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

void AdditionLU_SymbolicSymG(SprsUMatRealStru *pFU){
     
}
void AdditionLU_NumericSymG(SprsUMatRealStru *pFU){

}
// 输入参数:          // U阵结构及U阵值，右端项b
// 输出参数:          // 右端项x，维数为pU的维数（解向量）
void LE_FBackwardSym(SprsUMatRealStru *pFU, double b[], double x[]) {
  // int kbeg, kend;
  int iDim;
  int *rs_u, *j_u;
  double *d_u, *u_u;
  // double xc;

  d_u = pFU->d_u;
  u_u = pFU->u_u;
  rs_u = pFU->uMax.rs_u;
  j_u = pFU->uMax.j_u;
  iDim = pFU->uMax.iDim;

  for (int i = 1; i <= iDim; i++)
    x[i] = b[i];

  for (int i = 1; i <= iDim; i++) {
    double xc = x[i];
    int kbeg = rs_u[i];
    int kend = rs_u[i + 1];

    for (int k = kbeg; k < kend; k++) {
      int j = j_u[k];
      x[j] -= u_u[k] * xc;     
    }
  }

#pragma omp parallel for
#pragma vector always
#pragma ivdep
  for (int i = 1; i <= iDim; i++){
    x[i] *= d_u[i]; // auto vectorized
  }

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
