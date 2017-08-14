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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using std::cout;
using std::cerr;
using std::endl;

#include "LE_SymSprsMatDef.h"
#include "LE_SymSprsMatFunc.h"

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // initMem_MatReal
// 描    述:          // 稀疏全实数矩阵内存初始化。数目、指针变量置零
// 输入参数:          // 输入参数说明，包括每个参数的作
// 用、取值说明及参数间关系。
// 如参数非常复杂应举例说明。
// 输出参数:          // 对输出参数的说明。
// 返 回 值:          // 函数返回值的说明
// 其    他:          // 其它说明
//////////////////////////////////////////////////////////////////////

void initMem_MatReal(SprsMatRealStru *A) {
  // 内存初始化。数目、指针变量置零
  A->Mat.iDim = 0;
  A->Mat.iNy = 0;
  A->Mat.piIdiag = NULL;
  A->Mat.piIstart = NULL;
  A->Mat.piJno = NULL;
  A->Mat.piLinkn = NULL;
  A->Mat.piLinkp = NULL;
}

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // initMem_VecReal
// 描    述:          // 稀疏实数向量内存初始化。数目、指针变量置零
// 输入参数:          // 输入参数说明，包括每个参数的作
// 用、取值说明及参数间关系。
// 如参数非常复杂应举例说明。
// 输出参数:          // 对输出参数的说明。
// 返 回 值:          // 函数返回值的说明
// 其    他:          // 其它说明
//////////////////////////////////////////////////////////////////////
void initMem_VecReal(VecRealStru *V) {
  // 内存初始化。数目、指针变量置零
  V->iNy = 0;
  V->pdVal = NULL;
}

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // allocate_MatReal
// 描    述:          // 稀疏实数对称矩阵分配维数
// 输入参数:          // SprsMatRealStru *pA：稀疏实数对阵矩阵或上三角矩阵
// 输出参数:          // 无
// 返 回 值:          // 无
// 其    他:          // 调用该函数前稀疏实数对称矩阵的维数A->Mat.iDim已经确定
// 因可能还有注入元引入，稀疏实数对阵矩阵A的元素实际数目不确定，故分配维数时一次将维数扩够
//////////////////////////////////////////////////////////////////////
void allocate_MatReal(SprsMatRealStru *A) {
  int dimension = 0;
  int n = 0;
  cerr << "iDim: " << A->Mat.iDim << endl;
  cerr << "iNymax: " << A->Mat.iNymax << endl;
  cerr << "iNy: " << A->Mat.iNy << endl;
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

//////////////////////////////////////////////////////////////////////
// 描    述:          // 实数向量分配维数
// 输入参数:          // VecRealStru *V：实数向量
// 其    他:          // 调用该函数前向量元素数目V->iNy已经确定
//////////////////////////////////////////////////////////////////////
void allocate_VecReal(VecRealStru *V) {
  int m = 0;
  m = V->iNy + 1;
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

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // deallocate_VecReal
// 描    述:          // 指针变量内存释放
// 被deallocateNet()调用
// 输入参数:          // VecRealStru *V：
//////////////////////////////////////////////////////////////////////
void deallocate_VecReal(VecRealStru *V) {
  // 指针变量内存释放
  free(V->pdVal);
}

// 输入参数:          // G阵结构
// 输出参数:          // U阵结构,并申请U阵内存，包括工作相量
void LU_SymbolicSymG(SprsMatRealStru *pG, SprsUMatRealStru *pFU) {}

// 输入参数:          // G阵结构及G阵值
// 输出参数:          // U阵中的值
void LU_NumbericSymG(SprsMatRealStru *pG, SprsUMatRealStru *pFU) {}

// 输入参数:          // U阵结构及U阵值，右端项b
// 输出参数:          // 右端项x，维数为pU的维数（解向量）
void LE_FBackwardSym(SprsUMatRealStru *pFU, double b[], double x[]) {}

// 描    述:          //内存初始化。数目、指针变量置零
// 输入参数:          // U
void initMem_UMatReal(SprsUMatRealStru *U) {
  /* -------------------------------------------------------------------- */
  /* .. Setup Pardiso control parameters. */
  /* -------------------------------------------------------------------- */
  for (int i = 0; i < 64; i++) {
    U->iparm[i] = 0;
  }
  U->iparm[0] = 1;   /* No solver default */
  U->iparm[1] = 2;   /* Fill-in reordering from METIS */
  U->iparm[3] = 0;   /* No iterative-direct algorithm */
  U->iparm[4] = 0;   /* No user fill-in reducing permutation */
  U->iparm[5] = 0;   /* Write solution into x */
  U->iparm[6] = 0;   /* Not in use */
  U->iparm[7] = 2;   /* Max numbers of iterative refinement steps */
  U->iparm[8] = 0;   /* Not in use */
  U->iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
  U->iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
  U->iparm[11] = 0;  /* Not in use */
  U->iparm[12] = 0;  /* Maximum weighted matching algorithm is switched-off
                        (default for symmetric). Try iparm[12] = 1 in case of
                        inappropriate accuracy */
  U->iparm[13] = 0;  /* Output: Number of perturbed pivots */
  U->iparm[14] = 0;  /* Not in use */
  U->iparm[15] = 0;  /* Not in use */
  U->iparm[16] = 0;  /* Not in use */
  U->iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  U->iparm[18] = -1; /* Output: Mflops for LU factorization */
  U->iparm[19] = 0;  /* Output: Numbers of CG Iterations */
  U->iparm[34] = 1;  // zero based indexing
  U->mtype = -2;      /* Real symmetric matrix */
  U->nrhs = 1;     /* Number of right hand sides. */
  U->maxfct = 1;     /* Maximum number of numerical factorizations. */
  U->mnum = 1;       /* Which factorization to use. */
  U->msglvl = 1;     /* Print statistical information in file */
  U->error = 0;      /* Initialize error flag */
  /* -------------------------------------------------------------------- */
  /*  Initialize the internal solver memory pointer. This is only       */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  for (int i = 0; i < 64; i++) {
    U->pt[i] = 0;
  }
}

//////////////////////////////////////////////////////////////////////
// 描    述:          //指针变量内存释放
// 输入参数:          // U
//////////////////////////////////////////////////////////////////////
void deallocate_UMatReal(SprsUMatRealStru *U) {
  int phase = -1; // termination and release 
  PARDISO(U->pt, &U->maxfct, &U->mnum, &U->mtype, &phase, &U->n, &U->ddum,
          &U->idum, &U->idum, &U->idum, &U->nrhs, U->iparm, &U->msglvl,
          &U->ddum, &U->ddum, &U->error);
  );
}
