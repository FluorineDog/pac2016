//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, 中国电科院
// 文 件 名: LE_SymSprsMatDef.h
// 描    述: 线性方程组之稀疏求解法
//           包含：结构体及全局变量声明
// 其    他：用 #include <filename.h>
// 格式来引用标准库的头文件（编译器将从标准库目录开始搜索）
//           用 #include "filename.h"
//           格式来引用非标准库的头文件（编译器将从用户的工作目录开始搜索）
// 修改记录:     // 修改历史记录列表，每条修改记录应包括修改日期、修改
// 者及修改内容简述
// 1. 时间:
//    作者:
//    修改内容:
// 2. ...
//////////////////////////////////////////////////////////////////////

#ifndef _LE_SYMSPRSMATDEF_H__
#define _LE_SYMSPRSMATDEF_H__

#include <stdio.h>

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>
// 稀疏二维全矩阵结构
typedef struct {

  int iDim; // 矩阵维数。（MAXN）
  int iNy;  // 矩阵元素实际数目。（NY1、NU）
  int iNymax; // 矩阵元素最大数目,iNy = iDim*(iDim+1)/2，分配维数时使用
  int *piJno; // 每个矩阵元素列号，维数iNymax+1。（JNOY1、JNOU）
  int *piIstart; // 每行矩阵元素在iJno中的起始位置，维数iDim+2。（IYD1、IYDU）
  int *piIdiag;  // 对角线元素的位置。没有对角线元素时，
                 //指向Aij，j>i且最接近,维数iDim+1。（IYD1、IYDU）
  int *piLinkp;  // 位置纵向链,
                 // 矩阵中同一列下一个元素的位置，矩阵最后一行元素的linkp指向第一个元素,维数iNymax+1（LP1、LPU）
  int *piLinkn;  // 行号纵向链，下个元素的行号，
                 //矩阵最后一行元素的linkn为第一元素行号,维数iNymax+1。（LR1、LRU）

} SprsMatStru;

// 实数向量
typedef struct {
  int iNy;       // 向量元素数目
  double *pdVal; // 向量元素值，维数iNy。（I0）
} VecRealStru;

// 稀疏全实数矩阵
typedef struct {
  SprsMatStru Mat; // 矩阵结构
  double *pdVal;   // 矩阵元素值，维数iNy
} SprsMatRealStru;

// U'DU分解后的U阵结构描述
typedef struct {
  int iDim;  //矩阵维数
  int iNzs;  //上三角（不包括对角元)的非零元素个数
  int *rs_u; //上三角行向开始位置,iDim+2
  int *cs_u; //上三角列向开始位置,iDim+2
  int *r_u;  //上三角列向元素行号,iNzs+1
  int *j_u;  //上三角行向元素列号,iNzs+1
} SprsUMatStru;

struct pair_ii_t{
  int beg;
  int end;
};
struct UMat_t {
  int iDim;
  int alloc_size;            // sum of cache lines of double
  pair_ii_t *SeqRanges;  // seq part [beg, end)
  pair_ii_t *paraRanges; // para part [beg, end)
  int *columns;
};
// LU分解后的U阵值描述
typedef struct {
  SprsUMatStru uMax; //矩阵结构
  double *d_u;       //分解后对角元值,iDim+1维
  double *u_u;       //上三角元素值(行向),已归一化,iNzs+1维

  //以下两个数组为LU数值分解时所需要的工作数组
  int *nzs;     //工作数组 iDim+1
  double *work; //工作数组 iDim+1

  // dog implement
  UMat_t dogUMat;
  double *values;
} SprsUMatRealStru;

#endif
