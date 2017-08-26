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
//���ļ��ж�ȡ���γ�A����,0-�ɹ�������-ʧ��
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
//��A����д���ļ�,0-�ɹ�������-ʧ��
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
//���ļ��ж�ȡ���γ�B��������,0-�ɹ�������-ʧ��
//B��ָ�����飬����NULLָ�룬����ִ����ɺ�ָ���NULL
//nsize��ָ�������ά��������0��������ֵ������ָ������ά��
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
//��B��������д���ļ�,0-�ɹ�������-ʧ��
//B��ָ�����飬�����NULLָ�룬����ִ�й����в����޸�B������
//nsize��ָ�������ά��������ά��������ִ�й����в����޸�
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
//X���������ȡ�ļ���0�ɹ�������ʧ��
//X�����������飬nsize��ָ�������ά����filename���ļ�����
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
//��ȡX.txt���Ƚ�ִ��������Ľ�����X.txt�Ľ��������Ƿ����Ҫ�󣬲�����Ҫ�󣬼�¼log.txt��
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
    unsigned long long timeStart;   // ���濪ʼʱ��
    unsigned long long timeEnd;     // �������ʱ��
    double delapseTime = 0; // ������̻���ʱ��
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
    //��ӹ���A��Ķ�д����
    printf("Begin ReadMatrixA...\n");
    ReadMatrixA(&A,"A.txt");
    printf("ReadMatrixA finish!\n");
    ////////////////////////

    //����дAA////////////////
    //WriteMatrixA(&A,"AA.txt");
    //////////////////////////
    nsize = 10000;
    //  B = new VecRealStru[nsize];
    printf("Begin ReadVectorB...\n");
    //��ӹ���B�����ĺ���
    ReadVectorB(&B,nsize,"B.txt");
    printf("ReadVectorB finish!\n");
    ////////////////////////
    //printf("%22.15e\n", B[1].pdVal[1]);
    //for (i = 0; i < nsize;i++)
    // WriteVectorB(&B[i], nsize, "BB.txt");

    //   WriteVectorX(B, nsize, "X.txt");
    //��ʼ��X����
    X=new VecRealStru[nsize];
    for(i=0; i<nsize;i++)
    {
        initMem_VecReal(&X[i]);
        X[i].iNy=B[i].iNy;
        allocate_VecReal(&X[i]);
    }


    ////////////////////////////////

    printf("The Program is Running...\n");
    // ��ȡ���濪ʼʱ��
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

    // ��ȡ�������ʱ��
    //	timeEnd = clock();  //modify by wsh 2017.7.31
    timeEnd=rdtsc();
    // ������̻���ʱ��
    //	delapseTime = (double)(timeEnd - timeStart)/CLOCKS_PER_SEC;//modify by wsh 2017.7.31
    delapseTime=(double)(timeEnd - timeStart)/(F*Time);

    printf("The program elapsed %13.8f s\n",delapseTime);

    //��Ӵ�ӡ��������Ĵ���
    printf("Begin Print Result...\n");
    WriteVectorX(X,nsize,"X1.txt");
    printf("Print Result finish!\n");
    ////////////////////////
    //��ʼ��X����
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

