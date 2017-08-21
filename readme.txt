基于Intel众核平台的对称稀疏线性方程组求解并行优化
1、程序介绍：
程序提供两个源文件，LUSolve.cpp,LE_SymSprsMats.cpp以及一个Makefile文件：
其中：
a)LUSolve.cpp是主程序代码，提供数据读写，调用计算函数，输出结果。本文件代码不是优化内容部分。
b)LE_SymSprsMats.cpp是稀疏矩阵求解计算函数，是基于众核平台优化内容。
直接make可以生成可执行文件test1，直接运行可执行文件得出程序计时。
编译：
$make
icc -c LUSolve.cpp
icc -c LE_SymSprsmats.cpp
icc LUSolve.o  LE_SymSprsmats.o -o test1
运行：
./test1
The program elapsed :   程序执行时间单位秒

2、比赛规则:
 a) 优化内容不包括矩阵和右端项的读写磁盘操作；
 b) 优化内容限定于某一个右端项B分量的方程组求解部分，不得展开循环将所有循环合并统一并行优化;
 c) 不得修改LU_SymbolicSymG，LU_NumbericSymG和LE_FBackwardSym三个函数接口参数定义和格式；
 d) 优化平台限定在Intel 众核平台上。