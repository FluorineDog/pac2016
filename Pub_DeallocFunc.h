//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, XXXX公司
// 文 件 名: dealloctefunc.h
// 作    者:        版本:       时间:
// 描    述: 二维及以上数组空间释放相关接口函数声明 
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

#ifndef _PUB_DEALLOCFUNC_H__
#define _PUB_DEALLOCFUNC_H__

#include "Pub_Def.h"
//-----------函数声明--------------------------------------------
// 释放int型二维数组a内存
extern void deallocateI2(int **ppiArr);
// 释放double型二维数组a内存
extern void deallocateD2(double **ppdArr);
// 为YElmStru型二维数组a[m][n]分配内存
extern void deallocateYElm2(YElmStru **ppiArr);

#endif