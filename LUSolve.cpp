// LUSolve.cpp : Defines the entry point for the console application.
//
#include "LE_SymSprsMatDef.h"
#include "LE_SymSprsMatFunc.h"
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include <time.h>
//add by wsh 2017.7.31
#include<linux/types.h>
#define F 1.4E3
#define Time 1E6
inline unsigned long long rdtsc(void)
{
    unsigned long hi = 0, lo = 0;

    __asm__ __volatile__ ("lfence;rdtsc" : "=a"(lo), "=d"(hi));

    return (((unsigned long long)lo))|(((unsigned long long)hi)<<32);
}
//add end by wsh 2017.7.31
//从文件中读取，形成A矩阵,0-成功，其他-失败
int ReadMatrixA(SprsMatRealStru *A, char* filename)
{
    FILE * fp=NULL;
    char LineBuffer[1024];
    char *line;
    int i,first,last;
    if (!A)
    {
        fprintf(stderr,"ReadMatrixA():NULL arguments!\n");
        return -1;
    }

    fp = fopen(filename,"rb");
    if (!fp)
    {
        fprintf(stderr,"ReadMatrixA():OPen file failed!\n");
        return -2;
    }

    fgets(LineBuffer, 1024, fp);

    if (fscanf(fp, "%d %d %d\n", &(A->Mat.iDim), &(A->Mat.iNy), &(A->Mat.iNymax)) != 3)
    {
        fprintf(stderr,"ReadMatrixA():read error...!\n");
        return -3;
    }

    if (A->Mat.iDim < 0 || A->Mat.iNy < 0 || A->Mat.iNymax < 0)
    {
        fprintf(stderr, "ReadMatrixA():read error!\n");
        return -4;
    }
    allocate_MatReal(A);
    memset(LineBuffer,0,1024);
    fgets(LineBuffer, 1024, fp);
    for (i = 0; i < A->Mat.iNymax+1; i++)
    {
        if (fscanf(fp, "%d\n", &(A->Mat.piJno[i])) != 1)
        {
            fprintf(stderr, "ReadMatrixA():read error!\n");
            return -5;
        }
    }
    memset(LineBuffer, 0, 1024);
    fgets(LineBuffer, 1024, fp);
    for (i = 0; i < A->Mat.iDim+2; i++)
    {
        if (fscanf(fp, "%d\n", &(A->Mat.piIstart[i])) != 1)
        {
            fprintf(stderr, "ReadMatrixA():read error!\n");
            return -5;
        }
    }
    memset(LineBuffer, 0, 1024);
    fgets(LineBuffer, 1024, fp);
    for (i = 0; i < A->Mat.iDim + 1; i++)
    {
        if (fscanf(fp, "%d\n", &(A->Mat.piIdiag[i])) != 1)
        {
            fprintf(stderr, "ReadMatrixA():read error!\n");
            return -5;
        }
    }
    memset(LineBuffer, 0, 1024);
    fgets(LineBuffer, 1024, fp);
    for (i = 0; i < A->Mat.iNymax + 1; i++)
    {
        if (fscanf(fp, "%d\n", &(A->Mat.piLinkp[i])) != 1)
        {
            fprintf(stderr, "ReadMatrixA():read error!\n");
            return -5;
        }
    }
    memset(LineBuffer, 0, 1024);
    fgets(LineBuffer, 1024, fp);
    for (i = 0; i < A->Mat.iNymax + 1; i++)
    {
        if (fscanf(fp, "%d\n", &(A->Mat.piLinkn[i])) != 1)
        {
            fprintf(stderr, "ReadMatrixA():read error!\n");
            return -5;
        }
    }
    memset(LineBuffer, 0, 1024);
    fgets(LineBuffer, 1024, fp);
    for (i = 0; i < A->Mat.iNymax + 1; i++)
    {
        if (fscanf(fp, "%lf\n", &(A->pdVal[i])) != 1)
        {
            fprintf(stderr, "ReadMatrixA():read error!\n");
            return -5;
        }
    }
    return 0;
}
//将A矩阵写入文件,0-成功，其他-失败
int WriteMatrixA(SprsMatRealStru *A, char* filename)
{
    FILE *fp=NULL;
    int i;
    char LineBuffer[1024];

    if (!A)
    {
        fprintf(stderr,"WriteMatrixA:NULL arguments!\n");
        return -1;
    }

    fp = fopen(filename, "wb");
    if (!fp)
    {
        fprintf(stderr, "WriteMatrixA:open failed\n");
        return -2;
    }
    strcpy(LineBuffer, "SparsMatRealStruA.idm,SparsMatRealStruA.iNy,SparsMatRealStruA.iNymax");
    fprintf(fp,"%s\n",LineBuffer);
    if (A->Mat.iDim < 0 || A->Mat.iNy < 0 || A->Mat.iNymax < 0)
    {
        fprintf(stderr, "WriteMatrixA: input error\n");
        return -2;
    }

    fprintf(fp,"%d %d %d\n",A->Mat.iDim,A->Mat.iNy,A->Mat.iNymax);

    memset(LineBuffer,0,1024);
    sprintf(LineBuffer,"SparsMatRealStruA->piJno,size=%d",A->Mat.iNymax+1);
    fprintf(fp,"%s\n",LineBuffer);
    if (A->Mat.piJno != NULL)
    {
        for (i = 0; i < A->Mat.iNymax + 1; i++)
        {
            fprintf(fp, "%d\n", A->Mat.piJno[i]);
        }
    }
    else{
        fprintf(stderr,"WriteMatrixA:input value NULL!\n");
        return -3;
    }
    memset(LineBuffer, 0, 1024);
    sprintf(LineBuffer, "SparsMatRealStruA->piIstart,size=%d", A->Mat.iDim + 2);
    fprintf(fp, "%s\n", LineBuffer);
    if (A->Mat.piIstart != NULL)
    {
        for (i = 0; i < A->Mat.iDim + 2; i++)
        {
            fprintf(fp, "%d\n", A->Mat.piIstart[i]);
        }
    }else
        return -4;
    memset(LineBuffer, 0, 1024);
    sprintf(LineBuffer, "SparsMatRealStruA->piIdiag,size=%d", A->Mat.iDim + 1);
    fprintf(fp, "%s\n", LineBuffer);
    if (A->Mat.piIdiag != NULL)
    {
        for (i = 0; i < A->Mat.iDim + 1; i++)
        {
            fprintf(fp, "%d\n", A->Mat.piIdiag[i]);
        }
    }else
        return -5;
    memset(LineBuffer, 0, 1024);
    sprintf(LineBuffer, "SparsMatRealStruA->piLinkp,size=%d", A->Mat.iNymax + 1);
    fprintf(fp, "%s\n", LineBuffer);
    if (A->Mat.piLinkp != NULL)
    {
        for (i = 0; i < A->Mat.iNymax + 1; i++)
        {
            fprintf(fp, "%d\n", A->Mat.piLinkp[i]);
        }
    }else
        return -6;
    memset(LineBuffer, 0, 1024);
    sprintf(LineBuffer, "SparsMatRealStruA->piLinkn,size=%d", A->Mat.iNymax + 1);
    fprintf(fp, "%s\n", LineBuffer);
    if (A->Mat.piLinkn != NULL)
    {
        for (i = 0; i < A->Mat.iNymax + 1; i++)
        {
            fprintf(fp, "%d\n", A->Mat.piLinkn[i]);
        }
    }else
        return -7;
    memset(LineBuffer, 0, 1024);
    sprintf(LineBuffer, "SparsMatRealStruA->pdval,size=%d", A->Mat.iNymax + 1);
    fprintf(fp, "%s\n", LineBuffer);
    if (A->pdVal != NULL)
    {
        for (i = 0; i < A->Mat.iNymax + 1; i++)
        {
            fprintf(fp, "%22.15e\n", A->pdVal[i]);
        }
    }

    fclose(fp);
    return 0;
}
//从文件中读取，形成B向量数组,0-成功，其他-失败
//B是指针数组，传入NULL指针，函数执行完成后，指针非NULL
//nsize是指针数组的维数，传入0或者任意值，传出指针数组维数
int ReadVectorB(VecRealStru **B, int &nsize, char* filename)
{
    FILE * fp=NULL;
    int i,j,m,tnsize;
    double val;
    int tNy;
    char LineBuffer[1024];

    fp = fopen(filename,"rb");
    if (!fp)
    {
        fprintf(stderr,"ReadVectorB:open error!\n");
        return -1;
    }
    *B = new VecRealStru[nsize];
    for (i = 0; i < nsize; i++)
    {
        memset(LineBuffer,0,1024);
        fgets(LineBuffer,1024,fp);
        if (fscanf(fp, "%d\n", &tnsize) != 1)
        {
            fprintf(stderr,"ReadVectorB:read error!\n");
            return -2;
        }
        if (nsize != tnsize)
        {
            fprintf(stderr, "ReadVectorB: nsize error!\n");
            return -3;
        }
        memset(LineBuffer, 0, 1024);
        fgets(LineBuffer, 1024, fp);
        if (fscanf(fp, "%d\n", &tNy) != 1)
        {
            fprintf(stderr, "ReadVectorB:read error!\n");
            return -2;
        }

        ((*B)+i)->iNy = tNy;
        ((*B)+i)->pdVal = new double[((*B)+i)->iNy + 1];

        memset(LineBuffer, 0, 1024);
        fgets(LineBuffer, 1024, fp);
        for (j = 0; j < ((*B)+i)->iNy+1; j++)
        {
            if (fscanf(fp, "%lf\n", &val) != 1)
            {
                fprintf(stderr, "ReadVectorB:read error!\n");
                return -2;
            }
            ((*B)+i)->pdVal[j] = val;
            val = 0;

        }

    }
    fclose(fp);
    return 0;


}
//将B向量数组写入文件,0-成功，其他-失败
//B是指针数组，传入非NULL指针，函数执行过程中不得修改B的内容
//nsize是指针数组的维数，传入维数，函数执行过程中不得修改
int WriteVectorB(VecRealStru *B, int nsize, char* filename)
{
    FILE *fp = NULL;
    int i, j;
    char LineBuffer[1024];
    fp = fopen(filename, "ab+");
    if (!fp)
    {
        fprintf(stderr,"WriteVectorB:open failed!\n");
        return -1;
    }

    sprintf(LineBuffer,"Vector nsize");
    fprintf(fp, "%s\n",LineBuffer);
    fprintf(fp,"%d\n",nsize);

    memset(LineBuffer,0,1024);
    sprintf(LineBuffer,"VecReadStruN->iNy");
    fprintf(fp,"%s\n",LineBuffer);
    fprintf(fp, "%d\n", B->iNy);

    memset(LineBuffer, 0, 1024);
    sprintf(LineBuffer, "VecReadStruN->pdval,size=%d",B->iNy+1);
    fprintf(fp, "%s\n", LineBuffer);
    for (j = 0; j < B->iNy + 1; j++)
    {
        fprintf(fp, "%22.15e\n", B->pdVal[j]);
    }

    fclose(fp);
    return 0;
}
//X向量数组读取文件，0成功，其他失败
//X输入向量数组，nsize是指针数组的维数，filename是文件名字
int ReadVectorX(VecRealStru *X,int &nsize,char *filename)
{
    FILE * fp=NULL;
    int i,j,m,tnsize;
    double val;
    int tNy;
    char LineBuffer[1024];

    fp = fopen(filename,"rb");
    if (!fp)
    {
        fprintf(stderr,"ReadVectorX:open error!\n");
        return -1;
    }
    fgets(LineBuffer,1024,fp);
    if (fscanf(fp, "%d\n", &tnsize) != 1)
    {
        fprintf(stderr, "ReadVectorX:read error!\n");
        return -2;
    }
    if(nsize!=tnsize)
    {
        fprintf(stderr, "ReadVectorX:read size err!\n");
        return -3;
    }
    memset(LineBuffer,0,1024);
    fgets(LineBuffer,1024,fp);
    for (i = 0; i < tnsize;i++)
        if(fscanf(fp, "%d\n", &(X[i].iNy))!=1)
        {
            fprintf(stderr, "ReadVectorX:read error!\n");
            return -2;
        }
    //	memset(LineBuffer,0,1024);
    //	fgets(LineBuffer,1024,fp);
    for (i = 0; i < nsize; i++)
    {
        memset(LineBuffer, 0, 1024);
        fgets(LineBuffer,1024,fp);
        for (j = 0; j < X[i].iNy + 1; j++)
        {
            fscanf(fp, "%lf\n", &(X[i].pdVal[j]));
        }
    }
    fclose(fp);

    return 0;
}
//读取X.txt并比较执行算例后的结果与从X.txt的结果的误差是否符合要求，不符合要求，记录log.txt中
int CompareVectorX(VecRealStru *X,int &nsize,VecRealStru *result)
{
    int i,j,num;
    FILE *fp;
    fp=fopen("log.txt","wb");
    num=0;

    if(0>ReadVectorX(X,nsize,"X.txt"))
    {
        return -1;
    }

    for(i=0;i<nsize;i++)
    {
        for(j=0;j<X[i].iNy+1;j++)
        {
            if((X[i].pdVal[j]-result[i].pdVal[j]>1e-9)||(result[i].pdVal[j]-X[i].pdVal[j]>1e-9))
            {
                fprintf(fp,"matrix dimensions=%dColumn number matrix =%d,reference results=%22.15e,running result:%22.15e\n",i,j,X[i].pdVal[j],result[i].pdVal[j]);
                num++;
            }

        }

    }
    if(num>0)
    {
        printf(" result exceed the set reference value  %d,please see the log.txt\n",num);
        return -1;
    }
    printf("the result up to standard\n");
    fclose(fp);
    return 0;
}

int WriteVectorX(VecRealStru *X, int nsize, char * filename)
{
    FILE *fp = NULL;
    int i, j;
    char LineBuffer[1024];
    fp = fopen(filename, "wb");
    if (!fp)
    {
        fprintf(stderr, "WriteVectorB:open failed!\n");
        return -1;
    }

    sprintf(LineBuffer, "Vector nsize");
    fprintf(fp, "%s\n", LineBuffer);
    fprintf(fp, "%d\n", nsize);

    memset(LineBuffer, 0, 1024);
    sprintf(LineBuffer, "VecReadStruN->iNy");
    fprintf(fp, "%s\n", LineBuffer);
    for (i = 0; i < nsize;i++)
        fprintf(fp, "%d\n", X[i].iNy);


    for (i = 0; i < nsize; i++)
    {
        memset(LineBuffer, 0, 1024);
        sprintf(LineBuffer, "VecReadStruN->pdval,size=%d", X[i].iNy + 1);
        fprintf(fp, "%s\n", LineBuffer);
        for (j = 0; j < X[i].iNy + 1; j++)
        {
            fprintf(fp, "%22.15e\n", X[i].pdVal[j]);
        }
    }
    fclose(fp);
    return 0;
}
int main(int argc, char* argv[])
{
    unsigned long long timeStart;   // 仿真开始时间
    unsigned long long timeEnd;     // 仿真结束时间
    double delapseTime = 0; // 仿真过程花费时间
    //AX=B
    SprsMatRealStru A;
    SprsUMatRealStru U;
    VecRealStru *X;
    VecRealStru *B;
    VecRealStru *Reference;
    int nsize,i,j;

    initMem_MatReal(&A);
    initMem_UMatReal(&U);
    //add by wsh

    //add end by wsh 20170621
    //添加构建A阵的读写函数
    printf("Begin ReadMatrixA...\n");
    ReadMatrixA(&A,"A.txt");
    printf("ReadMatrixA finish!\n");
    ////////////////////////

    //测试写AA////////////////
    //WriteMatrixA(&A,"AA.txt");
    //////////////////////////
    nsize = 10000;
    //  B = new VecRealStru[nsize];
    printf("Begin ReadVectorB...\n");
    //添加构建B向量的函数
    ReadVectorB(&B,nsize,"B.txt");
    printf("ReadVectorB finish!\n");
    ////////////////////////
    //printf("%22.15e\n", B[1].pdVal[1]);
    //for (i = 0; i < nsize;i++)
    // WriteVectorB(&B[i], nsize, "BB.txt");

    //   WriteVectorX(B, nsize, "X.txt");
    //初始化X数组
    X=new VecRealStru[nsize];
    for(i=0; i<nsize;i++)
    {
        initMem_VecReal(&X[i]);
        X[i].iNy=B[i].iNy;
        allocate_VecReal(&X[i]);
    }


    ////////////////////////////////

    printf("The Program is Running...\n");
    // 获取仿真开始时间
    //	timeStart = clock();//modify by wsh 2017.7.31
    timeStart=rdtsc();

    LU_SymbolicSymG(&A,&U);
    LU_NumbericSymG(&A,&U);

    for(i=0; i<nsize; i++)
    {
        for(j=0;j < 10;j++)
        {
            LE_FBackwardSym(&U,B[i].pdVal,X[i].pdVal);
        }
    }

    // 获取仿真结束时间
    //	timeEnd = clock();  //modify by wsh 2017.7.31
    timeEnd=rdtsc();
    // 仿真过程花费时间
    //	delapseTime = (double)(timeEnd - timeStart)/CLOCKS_PER_SEC;//modify by wsh 2017.7.31
    delapseTime=(double)(timeEnd - timeStart)/(F*Time);

    printf("The program elapsed %13.8f s\n",delapseTime);

    //添加打印结果相量的代码
    printf("Begin Print Result...\n");
    WriteVectorX(X,nsize,"X1.txt");
    printf("Print Result finish!\n");
    ////////////////////////
    //初始化X数组
    Reference=new VecRealStru[nsize];
    for(i=0; i<nsize;i++)
    {
        initMem_VecReal(&Reference[i]);
        Reference[i].iNy=B[i].iNy;
        allocate_VecReal(&Reference[i]);
    }

    CompareVectorX(Reference,nsize,X);
    deallocate_MatReal(&A);
    deallocate_UMatReal(&U);
    for(i=0; i<nsize; i++)
    {
        deallocate_VecReal(&B[i]);
        deallocate_VecReal(&X[i]);
    }
    delete[] X;
    delete[] B;
    delete[] Reference;
    //getchar();
    return 0;
}

