//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, XXXX公司
// 文 件 名: Pub_ErrDefFunc.h
// 作    者:       版本:        时间: // 作者、版本及完成日期
// 描    述:     // 用于详细说明此程序文件完成的主要功能，与其他模块
                 // 或函数的接口，输出值、取值范围、含义及参数间的控
                 // 制、顺序、独立或依赖等关系
// 其    他：用 #include <filename.h> 格式来引用标准库的头文件（编译器将从标准库目录开始搜索）
//           用 #include "filename.h" 格式来引用非标准库的头文件（编译器将从用户的工作目录开始搜索）
// 修改记录:     // 修改历史记录列表，每条修改记录应包括修改日期、修改
                 // 者及修改内容简述  
// 1. 时间:2016.02.19
//    作者:彭红英
//    修改内容:错误信息结构ErrorStru增加变量iELevel
// 2. ...
//////////////////////////////////////////////////////////////////////
#ifndef _PUB_ERRDEFFUNC_H__
#define _PUB_ERRDEFFUNC_H__


// 可处理最大错误数
#define MAXERRNUM 100
//----------系统仿真错误号定义---------------------------------------
// 系统仿真错误号定义
// 必备文件不存在或打开出错，致命错误，程序退出
# define ERROR_FILEOPEN 1
// 必备文件存在但为空，致命错误，程序退出
# define ERROR_FILEEMPTY 2
// 元件参数填写有误,非致命的，可由程序修正
# define ERROR_PARA 3
// 元件参数有误，致命错误，程序退出
# define ERROR_PARAFATAL 4
// 元件类型不存在
# define ERROR_CTYP 5
// 输出变量号不存在
# define ERROR_VARNO 6
// 元件不存在输出变量
# define ERROR_NOVAROUT 7
// 内存分配失败
# define ERROR_ALLOC 8
// 内存重新分配失败
# define ERROR_REALLOC 9
// 创建文件失败
# define ERROR_FILECREAT 10
// 数目超过定义的最大维数，致命错误，程序退出
# define ERROR_OUTOFMAX 11
// 建模错误
# define ERROR_MODELING 12
// 子程序运行出错
# define ERROR_FUNC 13
// 除数为0
# define ERROR_DIV 14
// 矩阵求逆失败
# define ERROR_INVMAT 15
// 元件两端接地
# define ERROR_GNDBOTH 16
// 矩阵奇异，LU分解失败
# define ERROR_LUMAT 17
// 存在除基准、有名、标幺以外的输入输出变量类型
# define ERROR_IOCALTYP 18
// 自定义计算失败
# define ERROR_UDMCAL 19
// 插值算法自然通断插值时间计算失败
# define ERROR_INTERPTIMECAL 20
// 计算不收敛
# define ERROR_NOTCONVERGENT 21

//----------存储错误信息数组---------------------------------------
// 错误信息结构
typedef struct
{
	int iENo;	// 错误号
	int iECTyp;	// 错误元件类型
	int iECNo;	// 错误元件号
	int iELevel;// 错误等级：
}ErrorStru;


//-----------外部函数声明--------------------------------------------
// 仿真错误信息数组值清零
extern void clearSimError();
// 将输入的仿真错误赋值给错误信息数组中第一个错误号为0的位置
extern void setSimError(ErrorStru *Err);
// 统一处理仿真错误
extern int dealSimError();
// 内存分配成功与否检查，并打印相应信息
extern void checkPoint(void *p,char *str);

#endif