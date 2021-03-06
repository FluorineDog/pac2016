#ifndef OPTLUSOLVE_H_
#define OPTLUSOLVE_H_

#ifdef LUSOLVE_H_
#error "never include LUSolve.h"
#endif
// LUSolve.cpp : Defines the entry point for the console application.
//
#include "LE_SymSprsMatDef.h"
#include "LE_SymSprsMatFunc.h"
#include <cassert>
#include <cmath>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// add by wsh 2017.7.31
#include <linux/types.h>
#define F 1.4E3
#define Time 1E6
inline unsigned long long rdtsc(void) {
  unsigned long hi = 0, lo = 0;

  __asm__ __volatile__("lfence;rdtsc" : "=a"(lo), "=d"(hi));

  return (((unsigned long long)lo)) | (((unsigned long long)hi) << 32);
}
// add end by wsh 2017.7.31
//从文件中读取，形成A矩阵,0-成功，其他-失败
int ReadMatrixA(SprsMatRealStru *A, const char *filename) {

  FILE *fp = NULL;
  // char LineBuffer[1024];
  // const char *line;
  // int i, first, last;
  if (!A) {
    fprintf(stderr, "ReadMatrixA():NULL arguments!\n");
    return -1;
  }
  if (strcmp(filename, "A.txt")) {
    fprintf(stderr, "ReadMatrixA():OPen file failed!\n");
    return -1;
  }
  fp = fopen("A.dat", "rb");
  if (!fp) {
    fprintf(stderr, "ReadMatrixA():OPen file failed!\n");
    return -2;
  }

  // fgets(LineBuffer, 1024, fp);

  // if (fscanf(fp, "%d %d %d\n", &(A->Mat.iDim), &(A->Mat.iNy),
  //  &(A->Mat.iNymax)) != 3) {
  // fprintf(stderr, "ReadMatrixA():read error...!\n");
  // return -3;
  // }
  fread(&A->Mat.iDim, sizeof(int), 1, fp);
  fread(&A->Mat.iNy, sizeof(int), 1, fp);
  fread(&A->Mat.iNymax, sizeof(int), 1, fp);

  if (A->Mat.iDim < 0 || A->Mat.iNy < 0 || A->Mat.iNymax < 0) {
    fprintf(stderr, "ReadMatrixA():read error!\n");
    return -4;
  }
  allocate_MatReal(A);
  // for (int i = 0; i < A->Mat.iNymax + 1; i++) {
  // if (fscanf(fp, "%d\n", &(A->Mat.piJno[i])) != 1) {
  if (fread(A->Mat.piJno, sizeof(int), A->Mat.iNymax + 1, fp) !=
      A->Mat.iNymax + 1) {
    fprintf(stderr, "ReadMatrixA():read error!\n");
    return -5;
  }

  // for (i = 0; i < A->Mat.iDim + 2; i++) {
  // if (fscanf(fp, "%d\n", &(A->Mat.piIstart[i])) != 1) {
  if (fread(A->Mat.piIstart, sizeof(int), A->Mat.iDim + 2, fp) !=
      A->Mat.iDim + 2) {
    fprintf(stderr, "ReadMatrixA():read error!\n");
    return -5;
  }
  // for (i = 0; i < A->Mat.iDim + 1; i++) {
  // if (fscanf(fp, "%d\n", &(A->Mat.piIdiag[i])) != 1) {
  if (fread(A->Mat.piIdiag, sizeof(int), A->Mat.iDim + 1, fp) !=
      A->Mat.iDim + 1) {
    fprintf(stderr, "ReadMatrixA():read error!\n");
    return -5;
  }
  // for (i = 0; i < A->Mat.iNymax + 1; i++) {
  // if (fscanf(fp, "%d\n", &(A->Mat.piLinkp[i])) != 1) {
  if (fread(A->Mat.piLinkp, sizeof(int), A->Mat.iNymax + 1, fp) !=
      A->Mat.iNymax + 1) {
    fprintf(stderr, "ReadMatrixA():read error!\n");
    return -5;
  }

  // for (i = 0; i < A->Mat.iNymax + 1; i++) {
  // if (fscanf(fp, "%d\n", &(A->Mat.piLinkn[i])) != 1) {
  if (fread(A->Mat.piLinkn, sizeof(int), A->Mat.iNymax + 1, fp) !=
      A->Mat.iNymax + 1) {
    fprintf(stderr, "ReadMatrixA():read error!\n");
    return -5;
  }
  // for (i = 0; i < A->Mat.iNymax + 1; i++) {
  // if (fscanf(fp, "%lf\n", &(A->pdVal[i])) != 1) {
  if (fread(A->pdVal, sizeof(double), A->Mat.iNymax + 1, fp) !=
      A->Mat.iNymax + 1) {
    fprintf(stderr, "ReadMatrixA():read error!\n");
    return -5;
  }
  // int WriteMatrixA(SprsMatRealStru * A, const char *filename);
  // WriteMatrixA(A, "A.dat");
  return 0;
}
//将A矩阵写入文件,0-成功，其他-失败

int WriteMatrixA(SprsMatRealStru *A, const char *filename) {
  FILE *fp = NULL;
  int i;
  // char LineBuffer[1024];

  if (!A) {
    fprintf(stderr, "WriteMatrixA:NULL arguments!\n");
    return -1;
  }

  fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "WriteMatrixA:open failed\n");
    return -2;
  }
  // strcpy(
  // LineBuffer,
  // "SparsMatRealStruA.idm,SparsMatRealStruA.iNy,SparsMatRealStruA.iNymax");
  // fprintf(fp, "%s\n", LineBuffer);
  if (A->Mat.iDim < 0 || A->Mat.iNy < 0 || A->Mat.iNymax < 0) {
    fprintf(stderr, "WriteMatrixA: input error\n");
    return -2;
  }

  // fprintf(fp, "%d %d %d\n", A->Mat.iDim, A->Mat.iNy, A->Mat.iNymax);
  fwrite(&A->Mat.iDim, sizeof(int), 1, fp);
  fwrite(&A->Mat.iNy, sizeof(int), 1, fp);
  fwrite(&A->Mat.iNymax, sizeof(int), 1, fp);

  // memset(LineBuffer, 0, 1024);
  // sprintf(LineBuffer, "SparsMatRealStruA->piJno,size=%d", A->Mat.iNymax + 1);
  // fprintf(fp, "%s\n", LineBuffer);
  if (A->Mat.piJno != NULL) {
    // for (i = 0; i < A->Mat.iNymax + 1; i++) {
    // fprintf(fp, "%d\n", A->Mat.piJno[i]);
    // }
    fwrite(A->Mat.piJno, sizeof(int), A->Mat.iNymax + 1, fp);
  } else {
    fprintf(stderr, "WriteMatrixA:input value NULL!\n");
    return -3;
  }
  // memset(LineBuffer, 0, 1024);
  // sprintf(LineBuffer, "SparsMatRealStruA->piIstart,size=%d", A->Mat.iDim +
  // 2);
  // fprintf(fp, "%s\n", LineBuffer);
  if (A->Mat.piIstart != NULL) {
    // for (i = 0; i < A->Mat.iDim + 2; i++) {
    // fprintf(fp, "%d\n", A->Mat.piIstart[i]);
    // }
    fwrite(A->Mat.piIstart, sizeof(int), A->Mat.iDim + 2, fp);
  } else
    return -4;
  // memset(LineBuffer, 0, 1024);
  // sprintf(LineBuffer, "SparsMatRealStruA->piIdiag,size=%d", A->Mat.iDim + 1);
  // fprintf(fp, "%s\n", LineBuffer);
  if (A->Mat.piIdiag != NULL) {
    // for (i = 0; i < A->Mat.iDim + 1; i++) {
    // fprintf(fp, "%d\n", A->Mat.piIdiag[i]);
    // }
    fwrite(A->Mat.piIdiag, sizeof(int), A->Mat.iDim + 1, fp);
  } else
    return -5;
  // memset(LineBuffer, 0, 1024);
  // sprintf(LineBuffer, "SparsMatRealStruA->piLinkp,size=%d", A->Mat.iNymax +
  // 1);
  // fprintf(fp, "%s\n", LineBuffer);
  if (A->Mat.piLinkp != NULL) {
    fwrite(A->Mat.piLinkp, sizeof(int), A->Mat.iNymax + 1, fp);
  } else
    return -6;
  // memset(LineBuffer, 0, 1024);
  // sprintf(LineBuffer, "SparsMatRealStruA->piLinkn,size=%d", A->Mat.iNymax +
  // 1);
  // fprintf(fp, "%s\n", LineBuffer);
  if (A->Mat.piLinkn != NULL) {
    fwrite(A->Mat.piLinkn, sizeof(int), A->Mat.iNymax + 1, fp);
  } else
    return -7;
  // memset(LineBuffer, 0, 1024);
  // sprintf(LineBuffer, "SparsMatRealStruA->pdval,size=%d", A->Mat.iNymax + 1);
  // fprintf(fp, "%s\n", LineBuffer);
  if (A->pdVal != NULL) {
    // for (i = 0; i < A->Mat.iNymax + 1; i++) {
    // fprintf(fp, "%22.15e\n", A->pdVal[i]);
    // }
    fwrite(A->pdVal, sizeof(double), A->Mat.iNymax + 1, fp);
  } else
    return -8;

  fclose(fp);
  return 0;
}
// 从文件中读取，形成B向量数组,0-成功，其他-失败
// B是指针数组，传入NULL指针，函数执行完成后，指针非NULL
// nsize是指针数组的维数，传入0或者任意值，传出指针数组维数
// int ReadVectorB(VecRealStru **B, int &nsize, const char *filename) {
//   FILE *fp = NULL;
//   int i, j, m, tnsize;
//   double val;
//   int tNy;
//   char LineBuffer[1024];

//   fp = fopen(filename, "rb");
//   if (!fp) {
//     fprintf(stderr, "ReadVectorB:open error!\n");
//     return -1;
//   }
//   *B = new VecRealStru[nsize];
//   for (i = 0; i < nsize; i++) {
//     memset(LineBuffer, 0, 1024);
//     fgets(LineBuffer, 1024, fp);
//     if (fscanf(fp, "%d\n", &tnsize) != 1) {
//       fprintf(stderr, "ReadVectorB:read error!\n");
//       return -2;
//     }
//     if (nsize != tnsize) {
//       fprintf(stderr, "ReadVectorB: nsize error!\n");
//       return -3;
//     }
//     memset(LineBuffer, 0, 1024);
//     fgets(LineBuffer, 1024, fp);
//     if (fscanf(fp, "%d\n", &tNy) != 1) {
//       fprintf(stderr, "ReadVectorB:read error!\n");
//       return -2;
//     }

//     ((*B) + i)->iNy = tNy;
//     ((*B) + i)->pdVal = new double[((*B) + i)->iNy + 1];

//     memset(LineBuffer, 0, 1024);
//     fgets(LineBuffer, 1024, fp);
//     for (j = 0; j < ((*B) + i)->iNy + 1; j++) {
//       if (fscanf(fp, "%lf\n", &val) != 1) {
//         fprintf(stderr, "ReadVectorB:read error!\n");
//         return -2;
//       }
//       ((*B) + i)->pdVal[j] = val;
//       val = 0;
//     }
//   }
//   fclose(fp);
//   int WriteVectorB(VecRealStru *B, int nsize, const char *filename);
//   WriteVectorB(*B, nsize, "B.dat");
//   return 0;
// }
int ReadVectorB(VecRealStru **BB, int &nsize, const char *filename) {
  FILE *fp = NULL;
  // int i, j, m;
  int tnsize;
  double val;
  int tNy;
  // char LineBuffer[1024];

  if (!strcmp(filename, "B.txt")) {
    fp = fopen("B.dat", "rb");
  } else if (!strcmp(filename, "X.txt")) {
    fp = fopen("X.dat", "rb");
  } else {
    fprintf(stderr, "ReadMatrixA():OPen file failed!\n");
    return -1;
  }
  if (!fp) {
    fprintf(stderr, "ReadVectorB:open error!\n");
    return -1;
  }
  if (!strcmp(filename, "B.txt")) {
    *BB = new VecRealStru[nsize];
  }
  auto B = *BB;
  if (1 != fread(&tnsize, sizeof(int), 1, fp)) {
    fprintf(stderr, "ReadVectorB:read error!\n");
    return -2;
  }
  // fprintf(stderr, "%d %d", nsize, tnsize);
  if (nsize != tnsize) {
    fprintf(stderr, "ReadVectorB: nsizexxxx error!\n");
    return -3;
  }
  for (int i = 0; i < nsize; i++) {
    // memset(LineBuffer, 0, 1024);
    // fgets(LineBuffer, 1024, fp);
    if (1 != fread(&tNy, sizeof(int), 1, fp)) {
      fprintf(stderr, "ReadVectorB:tNy error!\n");
      return -2;
    }

    (B + i)->iNy = tNy;
    (B + i)->pdVal = new double[tNy + 1];

    // for (int j = 0; j < (B + i)->iNy + 1; j++) {
    // if (fscanf(fp, "%lf\n", &val) != 1) {
    if (fread((B + i)->pdVal, sizeof(double), tNy + 1, fp) != tNy + 1) {
      fprintf(stderr, "ReadVectorB: pdVal error!\n");
      return -4;
    }
  }
  fclose(fp);
  // int WriteVectorB(VecRealStru * B, int nsize, const char *filename);
  // WriteVectorB(B, nsize, "B.dat");
  return 0;
}
//将B向量数组写入文件,0-成功，其他-失败
// B是指针数组，传入非NULL指针，函数执行过程中不得修改B的内容
// nsize是指针数组的维数，传入维数，函数执行过程中不得修改
int WriteVectorB(VecRealStru *B, int nsize, const char *filename) {
  FILE *fp = NULL;
  char LineBuffer[1024];
  fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "WriteVectorB:open failed!\n");
    return -1;
  }
  fprintf(stderr, "fuck%dfuck", nsize);
  fwrite(&nsize, sizeof(int), 1, fp);
  for (int i = 0; i < nsize; ++i) {
    fwrite(&B[i].iNy, sizeof(int), 1, fp);
    fwrite(B[i].pdVal, sizeof(double), B->iNy + 1, fp);
  }
  fclose(fp);
  return 0;
}
// X向量数组读取文件，0成功，其他失败
// X输入向量数组，nsize是指针数组的维数，filename是文件名字
int ReadVectorX(VecRealStru *X, int &nsize, const char *filename) {
  FILE *fp = NULL;
  int i, j, m, tnsize;
  double val;
  int tNy;
  char LineBuffer[1024];

  fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "ReadVectorX:open error!\n");
    return -1;
  }
  fgets(LineBuffer, 1024, fp);
  if (fscanf(fp, "%d\n", &tnsize) != 1) {
    fprintf(stderr, "ReadVectorX:read error!\n");
    return -2;
  }
  if (nsize != tnsize) {
    fprintf(stderr, "ReadVectorX:read size err!\n");
    return -3;
  }
  memset(LineBuffer, 0, 1024);
  fgets(LineBuffer, 1024, fp);
  for (i = 0; i < tnsize; i++)
    if (fscanf(fp, "%d\n", &(X[i].iNy)) != 1) {
      fprintf(stderr, "ReadVectorX:read error!\n");
      return -2;
    }
  //	memset(LineBuffer,0,1024);
  //	fgets(LineBuffer,1024,fp);
  for (i = 0; i < nsize; i++) {
    memset(LineBuffer, 0, 1024);
    fgets(LineBuffer, 1024, fp);
    for (j = 0; j < X[i].iNy + 1; j++) {
      fscanf(fp, "%lf\n", &(X[i].pdVal[j]));
    }
  }
  fclose(fp);
  WriteVectorB(X, nsize, "X.dat");
  return 0;
}
//读取X.txt并比较执行算例后的结果与从X.txt的结果的误差是否符合要求，不符合要求，记录log.txt中
int CompareVectorX(VecRealStru *X, int &nsize, VecRealStru *result) {
  int i, j, num;
  FILE *fp;
  fp = fopen("log.txt", "wb");
  num = 0;

  if (0 > ReadVectorB(&X, nsize, "X.txt")) {
    return -1;
  }

  int tmp = 2 * nsize;
  for (i = 0; i < nsize; i++) {
    for (j = 0; j < X[i].iNy + 1; j++) {
      assert(!std::isnan(result[i].pdVal[j]));
      if ((X[i].pdVal[j] - result[i].pdVal[j] > 1e-9) ||
          (result[i].pdVal[j] - X[i].pdVal[j] > 1e-9)) {
        if (i > tmp) {
          // IO is too slow:
          // only record the first vector
          continue;
        } else {
          tmp = i;
          fprintf(fp, "matrix dimensions=%dColumn number matrix =%d,reference "
                      "results=%22.15e,running result:%22.15e\n",
                  i, j, X[i].pdVal[j], result[i].pdVal[j]);
        }
        num++;
      }
    }
  }

  if (num > 0) {
    printf(
        " result exceed the set reference value  %d,please see the log.txt\n",
        num);
    return -1;
  }
  printf("the result up to standard\n");
  fclose(fp);
  return 0;
}

int WriteVectorX(VecRealStru *X, int nsize, const char *filename) {
  FILE *fp = NULL;
  int i, j;
  char LineBuffer[1024];
  fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "WriteVectorB:open failed!\n");
    return -1;
  }

  sprintf(LineBuffer, "Vector nsize");
  fprintf(fp, "%s\n", LineBuffer);
  fprintf(fp, "%d\n", nsize);

  memset(LineBuffer, 0, 1024);
  sprintf(LineBuffer, "VecReadStruN->iNy");
  fprintf(fp, "%s\n", LineBuffer);
  for (i = 0; i < nsize; i++)
    fprintf(fp, "%d\n", X[i].iNy);

  for (i = 0; i < nsize; i++) {
    memset(LineBuffer, 0, 1024);
    sprintf(LineBuffer, "VecReadStruN->pdval,size=%d", X[i].iNy + 1);
    fprintf(fp, "%s\n", LineBuffer);
    for (j = 0; j < X[i].iNy + 1; j++) {
      fprintf(fp, "%22.15e\n", X[i].pdVal[j]);
    }
  }
  fclose(fp);
  return 0;
}

#endif