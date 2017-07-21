//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, XXXX公司
// 文 件 名: Pub_DbgDefFunc.h
// 作    者:       版本:        时间: // 作者、版本及完成日期
// 描    述: 程序debug相关常量定义：包括打印提示信息类型、输出位置
//           程序debug相关接口函数声明
// 其    他：用 #include <filename.h> 格式来引用标准库的头文件（编译器将从标准库目录开始搜索）
//           用 #include "filename.h" 格式来引用非标准库的头文件（编译器将从用户的工作目录开始搜索）
// 修改记录:     // 修改历史记录列表，每条修改记录应包括修改日期、修改
                 // 者及修改内容简述  
// 1. 时间:2016.02.19
//    作者:彭红英
//    修改内容:修改printSimMsg()函数声明
// 2. ...
//////////////////////////////////////////////////////////////////////
#ifndef _PUB_DBGDEFFUNC_H__
#define _PUB_DBGDEFFUNC_H__

#include <stdio.h>
#include <stdlib.h>

#include "Pub_Def.h"

// 常量定义
// 打印提示信息类型：错误、警告、提示信息，用于计算中的提示
# define WS_FATAL     1		// 致命错误
# define WS_ERROR     2		// 错误
# define WS_WARNING   3		// 警告
# define WS_INFO      4		// 提示

// 输出的位置
# define DOUT_SCREEN  1		// 屏幕打印
# define DOUT_FILE    2		// 文件
# define DOUT_NET     3		// 网络



//-----------外部函数声明--------------------------------------------
// Debug文件的打开和关闭
extern void DebugFileOpenClose(const char *pcProjOut,int iIsOpen);
// 根据类型iType打印输入信息str至位置iOut
//extern void printString(FixedString str,int iType,int iOut);
// 打印double型一维数组arr[iNB]
extern void printDArray(int iNB,const double *pdArr,const char *str,int iOut);
// 打印int型一维数组arr[iNB]
extern void printIArray(int iNB,const int *piArr,const char *str,int iOut);
// 打印double型二维数组arr[iNRow][iNCol]
extern void printD2Array(int iNRow,int iNCol,double **ppdArr,const char *str,int iOut);
// 打印仿真信息
//extern void printSimMsg(int iECTyp,int iECNo,int iENo,int iOut,FixedString str);
//extern void printSimMsg(char *str);
extern void printSimMsg(const char *str,...);
// 处理错误信息
//extern void dealSimMsg(int iECTyp,int iECNo,int iENo,int iType,int iOut,FixedString str);


#endif