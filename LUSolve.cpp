#include "optLUSolve.h"
#include <algorithm>
#include <iostream>
#include <memory>
using std::cout;
using std::cerr;
using std::endl;
void dog_checker(const SprsMatRealStru &A, const double x[],
                 const double std_b[]) {
  auto a = A.pdVal;
  auto ia = A.Mat.piIstart;
	auto ja = A.Mat.piJno;
	auto n = A.Mat.iDim + 1;
	std::unique_ptr<double[]> b(new double[n]);
  for (int i = 0; i < n; ++i) {
    double sum = 0;
    for (int jLoc = ia[i]; jLoc < ia[i + 1]; ++jLoc) {
			int j = ja[jLoc];
      b[i] += x[j]*a[jLoc];
      if(j != i){
        b[j] += x[i]*a[jLoc];
      }
		}	
  }
  double maxdiff = 0;
  for(int i = 0; i < n; ++i){
    b[i] -= std_b[i];
    maxdiff = std::max(b[i], maxdiff);
  }
  cerr << "maxdiff: " << maxdiff << endl;
  auto count = std::count_if(b.get(), b.get() + n, [=](double x){return x > maxdiff/10;});
  cerr << " and " << count << " numbers is in the same magnitude" << endl;
}
							

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

  //添加构建A阵的读写函数
  printf("Begin ReadMatrixA...\n");
  ReadMatrixA(&A, "A.txt");
  printf("ReadMatrixA finish!\n");
  ////////////////////////

  nsize = 10000;
  printf("Begin ReadVectorB...\n");
  //添加构建B向量的函数
  ReadVectorB(&B, nsize, "B.txt");
  printf("ReadVectorB finish!\n");
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
  timeStart = rdtsc();
  {
    LU_SymbolicSymG(&A, &U);
    LU_NumbericSymG(&A, &U);
    // cerr << "h";
    for (i = 0; i < nsize; i++) {
      for (j = 0; j < 20; j++) {
        // cerr << "fuck";
        LE_FBackwardSym(&U, B[i].pdVal, X[i].pdVal);
      }
    }
    // 获取仿真结束时间
  }
  timeEnd = rdtsc();
  delapseTime = (double)(timeEnd - timeStart) / (F * Time);

  printf("The program elapsed %13.8f s\n", delapseTime);

  //添加打印结果相量的代码
  printf("Begin Print Result...(disabled)\n");
  WriteVectorX(X, 1, "X1.txt");
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
  // cerr << "intel trash" << endl;
  // dog_checker(A, X[0].pdVal, B[0].pdVal);
  // cerr << "--------------------------" << endl;
  // cerr << "std trash" << endl;
  // dog_checker(A, Reference[0].pdVal, B[0].pdVal);

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
