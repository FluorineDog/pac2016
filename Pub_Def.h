//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, XXXX公司
// 文 件 名: Pub_Def.h
// 作    者:       版本:        时间: // 作者、版本及完成日期
// 描    述: 全局宏定义
//           全局常量定义
//           常用自定义类型
//           通过结构体定义
// 其    他：用 #include <filename.h> 格式来引用标准库的头文件（编译器将从标准库目录开始搜索）
//           用 #include "filename.h" 格式来引用非标准库的头文件（编译器将从用户的工作目录开始搜索）
// 修改记录:     // 修改历史记录列表，每条修改记录应包括修改日期、修改
                 // 者及修改内容简述  
// 1. 时间:2014.8.4
//    作者:彭红英
//    修改内容:MAXSTRLEN 256 修改为 MAXSTRLEN 512
//    原因：*.TB2文件一行记录字符串长度越限
// 2. 时间:2014.12.29
//    作者:彭红英
//    修改内容:MAXSTRLEN 512 修改为 MAXSTRLEN 2048
//    原因：*.GEN文件一行记录字符串长度越限
// 3. 时间:2015.4.26
//    作者:彭红英
//    修改内容:EPSINON 1.0e-10修改为 EPSINON 1.0e-14
//    原因：当元件参数较小时，如作为小开关的线路，阻抗矩阵的行列式为较小值，如1e-12，但不应该认为矩阵奇异
// 4. 时间:2015.4.29
//    作者:彭红英
//    修改内容:MAXNODE 2000修改为 MAXNODE 4000
//             MAXNY 32000修改为  MAXNY 128000
// 5. 时间:2015.8.14
//    作者:彭红英
//    修改内容:EPSINON 1.0e-14修改为 EPSINON 1.0e-16
//    原因：在step_GOV1()中判断dAVZ<=dEK时，dAVZ为1.0e-15，dEK=0.0时，判断错误
//////////////////////////////////////////////////////////////////////
#ifndef _PUB_DEF_H__
#define _PUB_DEF_H__

// 全局宏定义
// 求最小值
#ifndef min
	#define min(a, b)  (((a) < (b)) ? (a) : (b))
#endif
// 求最大值
#ifndef max
	#define max(a, b)  (((a) > (b)) ? (a) : (b))
#endif
// 求最近整数（四舍五入）
#define NINT(a) ((a >= 0)? (int)(a+0.5) : (int)(a-0.5)) // 同Fortran语言IDNINT

// 全局常量定义
// 允许的误差（即精度）
#define EPSINON 1.0e-16 // 1.0e-10-->1.0e-14-->1.0e-16
// 数值上限
#define VALUPLIMIT 1.0e+10
// 圆周率
#define PI 3.1415926535897932
// 中间变量最大维数
#define MAXMIDVAR 5000
// 节点最大数目
#define MAXNODE 4000
// 上三角导纳阵非零元数目, 100个节点系统, 0 元占96%,500个节点,为99.2%
#define MAXNY 128000
// 节点最大出线数
#define MAXOUTLINE 1000
// 子网最大数目
#define MAXSUB 10
// 字符串最大长度
#define MAXSTRLEN 2048 
// 文件允许最大列数
#define MAXCOL 100
// LU分解后，U矩阵最大数目
#define MAXNU 50000
// 真
#define TRUE 1
// 假
#define FALSE 0


// 常用自定义类型
// 固定大小的字符数组
typedef char FixedString[MAXSTRLEN];

// 通用结构定义
// 相量结构
typedef struct
{
	double dRms;      // 有效值
	double dTheta;	  // 相角
}PhasorStru;


// 发电机潮流结构
//typedef struct
//{
//	double dP0;      // 有功功率初值
//	double dQ0;	     // 无功功率初值
//}GenFLStru;

// 元件变量结构
typedef struct
{
	int iCTyp;	// 元件类型
	int iCNo;	// 元件序号
	int iVNo;	// 变量序号
}CVarIndexStru;

// 可转换为三个单相元件的三相元件结构
typedef struct
{
	int	iNum;		// 元件数目
	int	*piIndexNo;// 三相元件转换为三个单相元件时，其A相在单相元件中的编号
}ThrPhaseConvtStru;

// 触发脉冲信息结构
typedef struct  
{
	int iPulse;	    // 触发脉冲(01信号)
	double dTTag;	// 脉冲分度
}TrigPulseStru;

// 导纳阵元素信息结构
typedef struct  
{
	int iNo;	    // 行号或列号
	double dVal;	// 值
}YElmStru;

#endif