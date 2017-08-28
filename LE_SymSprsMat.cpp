//////////////////////////////////////////////////////////////////////
// ��Ȩ (C), 1988-1999, �й�������ѧ�о�Ժ
// �� �� ��: LE_SymSprsMat.c
// �汾: V1.0       ʱ��: 2014.12.31
// ��    ��: ���Է�����ģ��
//           ����LE_SymSprsMatFunc.h�����������ӳ���
// ��    ������ #include <filename.h> ��ʽ�����ñ�׼���ͷ�ļ������������ӱ�׼��Ŀ¼��ʼ������
//           �� #include "filename.h" ��ʽ�����÷Ǳ�׼���ͷ�ļ��������������û��Ĺ���Ŀ¼��ʼ������
// �޸ļ�¼:       // ��ʷ�޸ļ�¼
// 1. ʱ��:
//    ����:
//    �޸�����:
// 2. ...
//////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include<stdlib.h>
#include "omp.h"

#include "LE_SymSprsMatFunc.h"
#include "LE_SymSprsMatDef.h"

double **UT_trans, **U_trans;
int **UT_trans_in, **U_trans_in;
double *tempX;

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // initMem_MatReal
// ��    ��:          // ϡ��ȫʵ�������ڴ��ʼ������Ŀ��ָ���������
// �������:          // �������˵��������ÿ����������
// �á�ȡֵ˵�����������ϵ��
// ������ǳ�����Ӧ����˵����
// �������:          // �����������˵����
// �� �� ֵ:          // ��������ֵ��˵��
// ��    ��:          // ����˵��
//////////////////////////////////////////////////////////////////////

void initMem_MatReal(SprsMatRealStru *A)
{
    // �ڴ��ʼ������Ŀ��ָ���������
    A->Mat.iDim = 0;
    A->Mat.iNy = 0;
    A->Mat.piIdiag = NULL;
    A->Mat.piIstart = NULL;
    A->Mat.piJno = NULL;
    A->Mat.piLinkn = NULL;
    A->Mat.piLinkp = NULL;
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // initMem_VecReal
// ��    ��:          // ϡ��ʵ�������ڴ��ʼ������Ŀ��ָ���������
// �������:          // �������˵��������ÿ����������
// �á�ȡֵ˵�����������ϵ��
// ������ǳ�����Ӧ����˵����
// �������:          // �����������˵����
// �� �� ֵ:          // ��������ֵ��˵��
// ��    ��:          // ����˵��
//////////////////////////////////////////////////////////////////////
void initMem_VecReal(VecRealStru *V)
{
    // �ڴ��ʼ������Ŀ��ָ���������
    V->iNy = 0;
    V->pdVal = NULL;
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // allocate_MatReal
// ��    ��:          // ϡ��ʵ���Գƾ������ά��
// �������:          // SprsMatRealStru *pA��ϡ��ʵ���������������Ǿ���
// �������:          // ��
// �� �� ֵ:          // ��
// ��    ��:          // ���øú���ǰϡ��ʵ���Գƾ����ά��A->Mat.iDim�Ѿ�ȷ��
// ����ܻ���ע��Ԫ���룬ϡ��ʵ���������A��Ԫ��ʵ����Ŀ��ȷ�����ʷ���ά��ʱһ�ν�ά������
//////////////////////////////////////////////////////////////////////
void allocate_MatReal(SprsMatRealStru *A)
{
    int m = 0;
    int n = 0;

    m = A->Mat.iDim + 1;
    A->Mat.iNymax = m * (m + 1)/2;  // ����Ԫ�������Ŀ
    n = A->Mat.iNymax  + 1;

    A->Mat.piJno = (int *)calloc(n,sizeof(int));
    A->Mat.piIstart = (int *)calloc(m+1,sizeof(int));
    A->Mat.piIdiag = (int *)calloc(m,sizeof(int));
    A->Mat.piLinkp = (int *)calloc(n,sizeof(int));
    A->Mat.piLinkn = (int *)calloc(n,sizeof(int));
    A->pdVal = (double *)calloc(n,sizeof(double));
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // allocate_VecReal
// ��    ��:          // ʵ����������ά��
// �������:          // VecRealStru *V��ʵ������
// �������:          // ��
// �� �� ֵ:          // ��
// ��    ��:          // ���øú���ǰ����Ԫ����ĿV->iNy�Ѿ�ȷ��
//////////////////////////////////////////////////////////////////////
void allocate_VecReal(VecRealStru *V)
{
    int m = 0;
    m = V->iNy + 1;
    V->pdVal = (double *)calloc(m,sizeof(double));
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // deallocate_MatReal
// ��    ��:          // ָ������ڴ��ͷ�
// ��deallocateNet()����
// �������:          // ��
// �������:          // ��
// �� �� ֵ:          // ��
// ��    ��:          // ����˵��
//////////////////////////////////////////////////////////////////////
void deallocate_MatReal(SprsMatRealStru *A)
{
    free(A->Mat.piJno);
    free(A->Mat.piIstart);
    free(A->Mat.piIdiag);
    free(A->Mat.piLinkn );
    free(A->Mat.piLinkp);
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // deallocate_VecReal
// ��    ��:          // ָ������ڴ��ͷ�
// ��deallocateNet()����
// �������:          // VecRealStru *V��
// �������:          // ��
// �� �� ֵ:          // ��
// ��    ��:          // ����˵��
//////////////////////////////////////////////////////////////////////
void deallocate_VecReal(VecRealStru *V)
{
    // ָ������ڴ��ͷ�
    free(V->pdVal);
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // SparseMatrix_adlink
// ��    ��:          // ϡ������������
// �������:          // SprsMatRealStru *pA��ϡ��ʵ���������������Ǿ���
// �������:          // ��
// �� �� ֵ:          // ��
// ��    ��:          // ����˵��,writted by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void SparseMatrix_adlink(SprsMatRealStru *pA)
{
    int iDim;
    int i,j,k,m,r;
    int jmn,irow,ny;

    int *wb;
    int *istart,*jno;
    int *linkp,*linkn;
    double *pdval;

    iDim=pA->Mat.iDim;
    istart=pA->Mat.piIstart;
    jno=pA->Mat.piJno;
    pdval=pA->pdVal;
    ny=pA->Mat.iNy;

    linkp=pA->Mat.piLinkp;
    linkn=pA->Mat.piLinkn;

    wb=(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(wb,"SparseMatrix_adlink:wb");

    for(i=1; i<=ny; i++)
    {
        linkp[i]=0;
        linkn[i]=0;
    }

    for(i=1; i<=iDim; i++)
    {
        jmn=istart[i];
        wb[i]=jmn;
        linkp[jmn]=jmn;
        linkn[jmn]=i;
    }

    irow=1;
    for(i=1;i<=ny;i++)
    {
        if(i>=istart[irow+1])
            irow++;

        j=jno[i];

        if(irow!=j)
        {
            k = wb[j];
            m = linkp[k];
            r = linkn[k];
            linkp[i] = m;
            linkn[i] = r;
            linkp[k] = i;
            linkn[k] = irow;
            wb[j] = i;
        }
    }

    for(i=1; i<=iDim; i++)
    {
        jmn=istart[i];
        if(jmn==linkp[jmn])
        {
            linkn[jmn]=0;
            linkp[jmn]=0;
        }
    }

    free(wb);
    wb = NULL;
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // LU_EliminationTreeG
// ��    ��:          // ���ݵ��ɾ���G�Ľṹ�����ȥ��(G�Գ�)
// LU_SymbolicSymG����
// �������:          // �����о���ṹ��Ϣ��G����
// �������:          // ��ȥ����Ϣ��һά����
// �� �� ֵ:          // ��
// ��    ��:          // ����Ӱ��G���о�����κ���Ϣ,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void LU_EliminationTreeG(SprsMatRealStru *pG, int* pParent)
{
    int jmn;
    int i,j,t;
    int kp,kn,iDim;

    int *vParent;
    int *linkp,*linkn,*istart;

    iDim=pG->Mat.iDim;
    vParent=(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(vParent,"LU_EliminationTreeG:vParent");

    linkp=pG->Mat.piLinkp;
    linkn=pG->Mat.piLinkn;
    istart=pG->Mat.piIstart;

    for(i=1; i<=iDim; i++)
    {
        pParent[i] = 0;
        vParent[i] = 0;

        jmn = istart[i];
        kn = linkn[jmn];
        kp = linkp[jmn];
        if(kp==0) continue;

        while(kp!=jmn)
        {
            j = kn;
            while(vParent[j]!=0
                    && vParent[j] < i)
            {
                t = vParent[j];
                vParent[j] = i;
                j = t;
            }

            if(vParent[j]==0)
            {
                vParent[j] = i;
                pParent[j] = i;
            }

            kn=linkn[kp];
            kp=linkp[kp];
        }
    }

    free(vParent);
    vParent=NULL;
    return;
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // LU_SetUMatECountsG
// ��    ��:          // ȷ�����Ƿֽ����������Ԫ�ظ���
// LU_SymbolicSymG����
// �������:          // ��ȥ��pParent,��G��ṹ
// �������:          // U��ṹ��ά������ֵ
// �� �� ֵ:          // ��
// ��    ��:          // ����Ӱ��G���о�����κ���Ϣ,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void LU_SetUMatECountsG(SprsMatRealStru *pG,int* pParent,SprsUMatRealStru *pU)
{
    int i,j;
    int jmn,kp,kn,iDim;
    int *mark,*colcnt,*rowcnt;
    int *istart,*linkp,*linkn;
    int *r_u,*j_u,*rs_u,*cs_u;
    double *d_u,*u_u;

    iDim=pG->Mat.iDim;
    pU->uMax.iDim=iDim;

    mark =(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(mark,"LU_SetUMatECountsG:mark");
    colcnt=(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(colcnt,"LU_SetUMatECountsG:colcnt");
    rowcnt=(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(rowcnt,"LU_SetUMatECountsG:rowcnt");

    rs_u=(int *)calloc(iDim+2,sizeof(int));
    //checkPoint(rs_u,"LU_SetUMatECountsG:rs_u");
    pU->uMax.rs_u=rs_u;

    cs_u=(int *)calloc(iDim+2,sizeof(int));
    //checkPoint(cs_u,"LU_SetUMatECountsG:cs_u");
    pU->uMax.cs_u=cs_u;

    istart=pG->Mat.piIstart;
    linkp=pG->Mat.piLinkp;
    linkn=pG->Mat.piLinkn;

    for(i=1; i<=iDim; i++)
        colcnt[i]=0;

    for(i = 1; i <= iDim; i++)
    {
        mark[i] = i;
        rowcnt[i] = 0;

        jmn = istart[i];
        kn = linkn[jmn];
        kp = linkp[jmn];
        if(kp==0) continue;

        while(kp != jmn)
        {
            j = kn;
            while(mark[j] != i)
            {
                rowcnt[i] ++;
                colcnt[j] ++;
                mark[j] = i;
                j = pParent[j];
            }

            kn=linkn[kp];
            kp=linkp[kp];
        }
    }

    rs_u[1] = 1;
    cs_u[1] = 1;
    pU->uMax.iNzs=0;
    for(i=2; i<=iDim+1; i++)
    {
        pU->uMax.iNzs=pU->uMax.iNzs+colcnt[i-1];
        rs_u[i] = rs_u[i-1]+colcnt[i-1];
        cs_u[i] = cs_u[i-1]+rowcnt[i-1];
    }


    d_u=(double *)calloc(iDim+1,sizeof(double));
    //checkPoint(d_u,"LU_SetUMatECountsG:d_u");
    pU->d_u=d_u;

    u_u=(double *)calloc(pU->uMax.iNzs+1,sizeof(double));
    //checkPoint(u_u,"LU_SetUMatECountsG:u_u");
    pU->u_u=u_u;

    r_u=(int *)calloc(pU->uMax.iNzs+1,sizeof(int));
    //checkPoint(r_u,"LU_SetUMatECountsG:r_u");
    pU->uMax.r_u=r_u;

    j_u=(int *)calloc(pU->uMax.iNzs+1,sizeof(int));
    //checkPoint(j_u,"LU_SetUMatECountsG:j_u");
    pU->uMax.j_u=j_u;

    free(mark);
    free(colcnt);
    free(rowcnt);
    return;
}


//////////////////////////////////////////////////////////////////////
// �� �� ��:          // LU_SymbolicSymG
// ��    ��:          // ϡ��ʵ���Գƾ���������ӷֽ⣬ȷ��U��ṹ��Ԫ�ظ���
//
// �������:          // G��ṹ
// �������:          // U��ṹ,������U���ڴ棬������������
// �� �� ֵ:          // ��
// ��    ��:          // ����Ӱ��G���о�����κ���Ϣ,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void LU_SymbolicSymG(SprsMatRealStru *pG,SprsUMatRealStru *pFU)
{
    int p;
    int i,j,k;
    int jmn,kp,kn;
    int len,top,iDim;
    int *pParent;
    int *Lnz,*Flag,*Pattern;
    int *istart,*linkp,*linkn;
    int *rs_u,*j_u,*r_u;

    deallocate_UMatReal(pFU); //!!xdc!!20150513!!���ӣ�����ڴ�й¶����

    iDim=pG->Mat.iDim;
    pParent=(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(pParent,"LU_SymbolicSymG:pParent");

    Flag=(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(Flag,"LU_SymbolicSymG:Flag");

    Pattern=(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(Pattern,"LU_SymbolicSymG:Pattern");

    Lnz=(int *)calloc(iDim+1,sizeof(int));
    //checkPoint(Lnz,"LU_SymbolicSymG:Lnz");

    LU_EliminationTreeG(pG,pParent);
    LU_SetUMatECountsG(pG,pParent,pFU);

    istart=pG->Mat.piIstart;
    linkp=pG->Mat.piLinkp;
    linkn=pG->Mat.piLinkn;

    rs_u=pFU->uMax.rs_u;
    r_u=pFU->uMax.r_u;
    j_u=pFU->uMax.j_u;

    j = 1;
    for(k = 1; k <= iDim; k++)
    {
        jmn = istart[k];
        kn = linkn[jmn];
        kp = linkp[jmn];

        top = iDim+1;
        Flag[k] = k;
        Lnz[k] = 0;

        if(kp==0) continue;
        while(kp != jmn)
        {
            i = kn;
            for(len=1; Flag[i] != k; i=pParent[i])
            {
                Pattern[len++] = i;
                Flag[i] = k;
            }
            while(len>1)
            {
                Pattern[--top] = Pattern[--len];
            }

            kn=linkn[kp];
            kp=linkp[kp];
        }

        for (; top <= iDim; top++)
        {
            i = Pattern[top];
            p = rs_u[i]+Lnz[i];
            j_u[p] = k;
            Lnz[i]++;

            r_u[j] = i;
            j ++;
        }
    }

    pFU->nzs=Lnz; //ֱ�ӽ�Lnz���óɹ�������
    pFU->work=(double *)calloc(iDim+1,sizeof(double));
    //checkPoint(pFU->work,"LU_SymbolicSymG:work");

    free(Flag);
    free(Pattern);
    free(pParent);
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          // LU_NumbericSymG
// ��    ��:          //ϡ��ʵ���Գƾ�����ֵ���ӷֽ⣬ȷ��U��Ԫ��ֵ
//
// �������:          // G��ṹ��G��ֵ
// �������:          // U���е�ֵ
// �� �� ֵ:          // ��
// ��    ��:          // ����Ӱ��G���о�����κ���Ϣ,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void LU_NumbericSymG(SprsMatRealStru *pG,SprsUMatRealStru *pFU)
{
    int i,j,k;
    int p,m,n;
    int kn,kp;
    int jmn,jmx;
    int iDim;

    double dj,dk;
    double yc,yr;

    int *nzs;
    int *istart,*linkp,*linkn;
    double *g;

    int *rs_u,*cs_u,*r_u,*j_u;
    double *d_u,*u_u;
    double *work;

    nzs = pFU->nzs;
    work = pFU->work;
    iDim=pG->Mat.iDim;

    istart=pG->Mat.piIstart;
    linkp=pG->Mat.piLinkp;
    linkn=pG->Mat.piLinkn;
    g=pG->pdVal;

    rs_u=pFU->uMax.rs_u;
    cs_u=pFU->uMax.cs_u;
    r_u=pFU->uMax.r_u;
    j_u=pFU->uMax.j_u;
    d_u=pFU->d_u;
    u_u=pFU->u_u;

    for(i = 1; i <= iDim; i++)
    {
        nzs[i] = 0;
        work[i]=0.0;

        jmn = istart[i];
        kn = linkn[jmn];
        kp = linkp[jmn];

        if(kp==0)
        {
            dk = g[jmn];
            if(dk >= 0.0 && fabs(dk) < 1.0e-20)	dk = 1.0e-20;
            if(dk < 0.0 && fabs(dk) < 1.0e-20) dk = -1.0e-20;
            d_u[i] = dk;
            continue;
        }
        while(kp != jmn)
        {
            work[kn] = g[kp];
            kn = linkn[kp];
            kp = linkp[kp];
        }

        dk = g[jmn];

        jmn = cs_u[i];
        jmx = cs_u[i+1];
        for(k = jmn; k < jmx; k++)
        {
            j = r_u[k];
            dj = d_u[j];

            yr = work[j];
            work[j] = 0.0;

            n = rs_u[j] + nzs[j];
            for(p = rs_u[j]; p < n; p++)
            {
                m = j_u[p];
                work[m] -= u_u[p]*yr;
            }

            yc = yr/dj;
            u_u[p] = yc;
            dk -= yr*yc;
            nzs[j]++ ;
        }

        if(dk >= 0.0 && fabs(dk) < 1.0e-20)	dk = 1.0e-20;
        if(dk < 0.0 && fabs(dk) < 1.0e-20) dk = -1.0e-20;

        d_u[i] = dk;
    }


    /// ------------------------>> addition by Sixu Hu

    /// 1. remove zeros
    //printf("%d\n", pFU->uMax.iNzs);
    j = 0;                          // the new array counter
    p = rs_u[1];                    // init rs_u[n+1]
    rs_u[1] = 0;
    for(i = 1; i <= iDim; ++i) {
        m = p;
        n = rs_u[i + 1];                      // n <- p (last rs_u[n+1])

        for(k = m; k < n; ++k) {
            u_u[j] = u_u[k];
            j_u[j] = j_u[k];
            if(u_u[j] != 0.0)
                j++;
        }
        p = rs_u[i + 1];            // store rs_u[i+1] -> p
        rs_u[i + 1] = j;
    }

    /// 2. get UT^-1
    UT_trans = (double **)calloc(iDim + 1, sizeof(double *));
    UT_trans_in = (int **)calloc(iDim + 1, sizeof(int *));
    for(int i=0; i<iDim + 1; ++i){
        UT_trans[i] = (double *)calloc(iDim + 1, sizeof(double));
        memset(UT_trans[i], 0, sizeof(double) * (iDim + 1));

        UT_trans_in[i] = (int *)calloc(iDim + 1, sizeof(int));
    }

    for(i=1; i<=iDim; ++i){
        kn = rs_u[i];
        kp = rs_u[i+1];
        for(k = kn; k < kp; ++k){
            for(j = 1; j < i; ++j){
                if(UT_trans[i][j] != 0.0)
                    UT_trans[j_u[k]][j] -= u_u[k] * UT_trans[i][j];
            }
            UT_trans[j_u[k]][i] = -u_u[k];
        }
    }

    // remove zeros
    for(int i=1; i <= iDim; ++i){
        j = 1;
        for(k = 1; k <= i; ++k){
            if(UT_trans[i][k] != 0.0){
                UT_trans[i][j] = UT_trans[i][k];
                UT_trans_in[i][j] = k;
                j++;
            }
        }
        UT_trans_in[i][j] = 0;
    }

    //// <!-- debug -->
    //for(i = 1; i <= iDim; ++i){
    //    printf("%4d: ", i);
    //    for(j = 1; UT_trans_in[i][j] != 0; ++j)
    //        printf("%15le %d  ", UT_trans[i][j], UT_trans_in[i][j]);
    //    putchar('\n');
    //}
    //exit(0);

    /// 3. get U^-1
    U_trans = (double **)calloc(iDim + 1, sizeof(double *));
    U_trans_in = (int **)calloc(iDim + 1, sizeof(int *));
    for(int i=0; i<iDim + 1; ++i){
        U_trans[i] = (double *)calloc(iDim + 1, sizeof(double));
        memset(U_trans[i], 0, sizeof(double) * (iDim + 1));

        U_trans_in[i] = (int *)calloc(iDim + 1, sizeof(int));
    }

    for(i = iDim; i >= 1; --i){
        kn = rs_u[i];
        kp = rs_u[i+1];
        for(k = kp - 1; k >= kn; --k){
            for(j = j_u[k] + 1; j <= iDim; ++j){
                if(U_trans[j_u[k]][j] != 0.0)
                    U_trans[i][j] -= u_u[k] * U_trans[j_u[k]][j];
            }
            U_trans[i][j_u[k]] = -u_u[k];
        }
    }

    // remove zeros
    for(int i=1; i <= iDim; ++i){
        j = 1;
        for(k = i; k <= iDim; ++k){
            if(U_trans[i][k] != 0.0){
                U_trans[i][j] = U_trans[i][k];
                U_trans_in[i][j] = k;
                j++;
            }
        }
        U_trans_in[i][j] = 0;
    }

    // <!-- debug -->
    //for(i = 1; i <= iDim; ++i){
    //    printf("%4d: ", i);
    //    for(j = 1; U_trans_in[i][j] != 0; ++j)
    //        printf("%15le %d  ", U_trans[i][j], U_trans_in[i][j]);
    //    putchar('\n');
    //}
    //exit(0);

    tempX = (double *)calloc(iDim + 1, sizeof(double));
    //omp_set_num_threads(3);
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          //LE_FBackwardSym
// ��    ��:          //�Գƾ���ǰ�ƻش�������ⷽ����
//
// �������:          // U��ṹ��U��ֵ���Ҷ���b
// �������:          // �Ҷ���x��ά��ΪpU��ά������������
// �� �� ֵ:          // ��
// ��    ��:          // ����Ӱ��U���о�����κ���Ϣ,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void LE_FBackwardSym(SprsUMatRealStru *pFU,double b[],double x[])
{
    int i,j,k, r;
    int ks,ke;
    int iDim;
    int *rs_u,*j_u;
    double *d_u,*u_u;
    double xc;

    d_u = pFU->d_u;               // �Խ�Ԫ
    u_u = pFU->u_u;               // �����һ��
    rs_u = pFU->uMax.rs_u;        // ÿһ�Ф������һ��Ԫ����u_u�е�λ��
    j_u = pFU->uMax.j_u;          // ����Ԫ���к�
    iDim = pFU->uMax.iDim;        // ά��

    //for(int i=0; i<iDim; ++i){
    //    ks = rs_u[i];
    //    ke = rs_u[i+1];

    //    printf("%d", ks);
    //    for(k = ks; k < ke; k ++)
    //        printf("%15le %d", u_u[k], j_u[k]);
    //}
    // ��Ux=b
    #pragma omp parallel for private(i) shared(x, b)
    for(i = 1; i <= iDim; i++)  // x <- b
        x[i] = b[i];

    //static int temp = 0;
    //printf("%d\n", temp);
    //temp ++;

    // ->> modified
    #pragma omp parallel for private(i, j) shared(UT_trans, UT_trans_in, x, b) schedule(dynamic, 3)
    for(i = 1; i <= iDim; ++i){
        for(j = 1; UT_trans_in[i][j] != 0; ++j){
            x[i] += UT_trans[i][j] * b[UT_trans_in[i][j]];
        }
    }

    //// ->> original
    //for(i = 1; i <= iDim; i++)  //
    //{
    //    xc = x[i];
    //    ks = rs_u[i];
    //    ke = rs_u[i+1];

    //    for(k = ks; k < ke; k ++)
    //    {
    //        x[j_u[k]] -= u_u[k] * xc;
    //    }
    //}

    #pragma omp parallel for private(i) shared(x, d_u)
    for(i=1; i<=iDim; i++)
        x[i] /= d_u[i];


    // ->> modified
    #pragma omp parallel for private(i) shared(x, b)
    for(i = 1; i <= iDim; i++)  // b <- x
        tempX[i] = x[i];

    #pragma omp parallel for private(i, j) shared(U_trans, U_trans_in, x, b) schedule(dynamic, 3)
    for(i = iDim; i >= 1; --i){
        for(j = 1; U_trans_in[i][j] != 0; ++j){
            x[i] += U_trans[i][j] * tempX[U_trans_in[i][j]];
        }
    }

    // ->> original
    //for(i=iDim-1; i>=1; i--)
    //{
    //    ks = rs_u[i];
    //    ke = rs_u[i+1] - 1;
    //    xc = x[i];

    //    for(k=ke; k>=ks; k--)
    //    {
    //        xc -= u_u[k] * x[j_u[k]];
    //    }
    //    x[i] = xc;
    //}
}


//////////////////////////////////////////////////////////////////////
// �� �� ��:          //initMem_UMatReal
// ��    ��:          //�ڴ��ʼ������Ŀ��ָ���������
//
// �������:          // U
// �������:          //
// �� �� ֵ:          // ��
// ��    ��:          // writed by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void initMem_UMatReal(SprsUMatRealStru *U)
{
    U->d_u=NULL;
    U->u_u=NULL;
    U->nzs=NULL;
    U->work=NULL;
    U->uMax.cs_u=NULL;
    U->uMax.rs_u=NULL;
    U->uMax.r_u=NULL;
    U->uMax.j_u=NULL;
    U->uMax.iDim=0;
    U->uMax.iNzs=0;
}

//////////////////////////////////////////////////////////////////////
// �� �� ��:          //deallocate_UMatReal
// ��    ��:          //ָ������ڴ��ͷ�
//
// �������:          // U
// �������:          //
// �� �� ֵ:          // ��
// ��    ��:          // writed by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void deallocate_UMatReal(SprsUMatRealStru *U)
{
    free(U->d_u);
    free(U->u_u);
    free(U->nzs);
    free(U->work);
    free(U->uMax.cs_u);
    free(U->uMax.rs_u);
    free(U->uMax.r_u);
    free(U->uMax.j_u);
    initMem_UMatReal(U);
}

