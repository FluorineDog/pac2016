//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, 中国电科院
// 文 件 名: LE_SymSprsMatFunc.h
// 作    者:       版本:        时间: // 作者、版本及完成日期
// 描    述: 线性方程组
//           包含：外部函数声明
// 其    他：用 #include <filename.h> 格式来引用标准库的头文件（编译器将从标准库目录开始搜索）
//           用 #include "filename.h" 格式来引用非标准库的头文件（编译器将从用户的工作目录开始搜索）
// 修改记录:     // 修改历史记录列表，每条修改记录应包括修改日期、修改
                 // 者及修改内容简述  
// 1. 时间:
//    作者:
//    修改内容:
// 2. ...
//////////////////////////////////////////////////////////////////////

#ifndef _LE_SYMSPRSMATFUNC_H__
#define _LE_SYMSPRSMATFUNC_H__

#include <stdio.h>

#include "LE_SymSprsMatDef.h"

//-----------外部函数声明--------------------------------------------
// 稀疏全实数矩阵相关函数
//// 内存初始化。数目、指针变量置零
extern void initMem_MatReal(SprsMatRealStru *A);
//// 指针变量内存申请
extern void allocate_MatReal(SprsMatRealStru *A);
//// 指针变量内存释放
extern void deallocate_MatReal(SprsMatRealStru *A);

// 稀疏实数向量相关函数
//// 内存初始化。数目、指针变量置零
extern void initMem_VecReal(VecRealStru *V);
//// 指针变量内存申请
extern void allocate_VecReal(VecRealStru *V);
//// 指针变量内存释放
extern void deallocate_VecReal(VecRealStru *V);

// 稀疏矩阵加链
// extern void SparseMatrix_adlink(SprsMatRealStru *pA);
// 内存初始化。数目、指针变量置零
extern void initMem_UMatReal(SprsUMatRealStru *U);
// 指针变量内存释放
extern void deallocate_UMatReal(SprsUMatRealStru *U);

// 第二种LU分解和网络求解：xdc提供
//cmt：初始形成G阵时，需调用符号因子分解LU_SymbolicSymG，为U阵分配空间
//cmt：每一时步，如果G阵结构未变，只有值变化，可以只调用值分解LU_NumbericSymG
//稀疏实数对称矩阵符号因子分解，确定U阵结构和元素个数
extern void LU_SymbolicSymG(SprsMatRealStru *pG,SprsUMatRealStru *pFU);
//稀疏实数对称矩阵数值因子分解，确定U阵元素值
extern void LU_NumbericSymG(SprsMatRealStru *pG,SprsUMatRealStru *pFU);
//对称矩阵前推回代方法求解方程组
extern void LE_FBackwardSym(SprsUMatRealStru *pFU,double b[],double x[]);
/////////end of xdc added 2014-6-1cod6

#endif

