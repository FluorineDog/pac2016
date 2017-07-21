//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, XXXX公司
// 文 件 名: Pub_AllocFunc.h
// 作    者:       版本:        时间: // 作者、版本及完成日期
// 描    述: 二维及以上数组分配空间接口函数声明
// 其    他：用 #include <filename.h> 格式来引用标准库的头文件（编译器将从标准库目录开始搜索）
//           用 #include "filename.h" 格式来引用非标准库的头文件（编译器将从用户的工作目录开始搜索）
// 修改记录:     // 修改历史记录列表，每条修改记录应包括修改日期、修改
                 // 者及修改内容简述  
// 1. 时间:
//    作者:
//    修改内容:
// 2. ...
//////////////////////////////////////////////////////////////////////
#pragma once

#ifndef _PUB_ALLOCFUNC_H__
#define _PUB_ALLOCFUNC_H__

#include "Pub_Def.h"

//-----------函数声明--------------------------------------------
// 分配int型二维数组a[m][n]内存
extern void allocateI2(int m,int n,int ***a);
// 分配double型二维数组a[m][n]内存
extern void allocateD2(int m,int n,double ***a);
// 为YElmStru型二维数组a[m][n]分配内存
extern void allocYElmStru2(int m,int n,YElmStru ***a);

#endif