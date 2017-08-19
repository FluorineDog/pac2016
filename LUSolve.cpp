#include "LUSolve.h"
#include <iostream> 
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, const char *argv[]) {
  unsigned long long timeStart; // 仿真开始时间
  unsigned long long timeEnd;   // 仿真结束时间
  double delapseTime = 0;       // 仿真过程花费时间
  // AX=B
  SprsMatRealStru A;
  SprsUMatRealStru U;
  VecRealStru *X;
  VecRealStru *B;
  VecRealStru *Reference;
  int nsize, i, j;

  initMem_MatReal(&A);
  initMem_UMatReal(&U);
  // add by wsh

  // add end by wsh 20170621
  //添加构建A阵的读写函数
  printf("Begin ReadMatrixA...\n");
  ReadMatrixA(&A, "A.txt");
  printf("ReadMatrixA finish!\n");
  ////////////////////////

  //测试写AA////////////////
  // WriteMatrixA(&A,"AA.txt");
  //////////////////////////
  nsize = 10000;
  //  B = new VecRealStru[nsize];
  printf("Begin ReadVectorB...\n");
  //添加构建B向量的函数
  ReadVectorB(&B, nsize, "B.txt");
  printf("ReadVectorB finish!\n");
  ////////////////////////
  // printf("%22.15e\n", B[1].pdVal[1]);
  // for (i = 0; i < nsize;i++)
  // WriteVectorB(&B[i], nsize, "BB.txt");

  //   WriteVectorX(B, nsize, "X.txt");
  //初始化X数组
  X = new VecRealStru[nsize];
  for (i = 0; i < nsize; i++) {
    initMem_VecReal(&X[i]);
    X[i].iNy = B[i].iNy;
    allocate_VecReal(&X[i]);
  }

  ////////////////////////////////

  printf("The Program is Running...\n");
  // 获取仿真开始时间
  //	timeStart = clock();//modify by wsh 2017.7.31
  timeStart = rdtsc();
  {
    LU_SymbolicSymG(&A, &U);
    LU_NumbericSymG(&A, &U);
    cerr << "h";
    for (i = 0; i < nsize; i++) {
      for (j = 0; j < 1; j++) {
        cerr << "fuck";
        LE_FBackwardSym(&U, B[i].pdVal, X[i].pdVal);
      }
    }

    // 获取仿真结束时间
    //	timeEnd = clock();  //modify by wsh 2017.7.31
  }
  timeEnd = rdtsc();
  delapseTime = (double)(timeEnd - timeStart) / (F * Time);

  printf("The program elapsed %13.8f s\n", delapseTime);

  //添加打印结果相量的代码
  printf("Begin Print Result...(disabled)\n");
  // WriteVectorX(X, nsize, "X1.txt");
  printf("Print Result finish!\n");
  ////////////////////////
  //初始化X数组
  Reference = new VecRealStru[nsize];
  for (i = 0; i < nsize; i++) {
    initMem_VecReal(&Reference[i]);
    Reference[i].iNy = B[i].iNy;
    allocate_VecReal(&Reference[i]);
  }

  CompareVectorX(Reference, nsize, X);
  deallocate_MatReal(&A);
  deallocate_UMatReal(&U);
  for (i = 0; i < nsize; i++) {
    deallocate_VecReal(&B[i]);
    deallocate_VecReal(&X[i]);
  }
  delete[] X;
  delete[] B;
  delete[] Reference;
  // getchar();
  return 0;
}
