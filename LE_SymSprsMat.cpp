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
#include "omp.h"

#include "LE_SymSprsMatFunc.h"
#include "LE_SymSprsMatDef.h"

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) < (y) ? (y) : (x))

double **UT_mu, **U_mu, **UT_su, **U_su;
int **UT_min, **U_min, **UT_sin, **U_sin;
int block_size;

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

    /// additions by husixu

    // 1. remove the zeros
    //printf("%d\n", pFU->uMax.iNzs);
    j = 1;                          // the new array counter
    p = rs_u[1];                    // init rs_u[n+1]
    rs_u[1] = 1;
    for(i = 1; i <= iDim; ++i) {
        m = p;
        n = rs_u[i + 1];                      // n <- p (last rs_u[n+1])

        for(k = m; k < n; ++k) {
            u_u[j] = u_u[k];                // TODO: can modified to if(u_u[k] == 0){...}else{...}
            j_u[j] = j_u[k];
            if(u_u[j] != 0.0)
                j++;
        }
        p = rs_u[i + 1];            // store rs_u[i+1] -> p
        rs_u[i + 1] = j;
    }

    UT_mu  = (double **)calloc(iDim + 1, sizeof(double *));
    UT_min  = (int **)calloc(iDim + 1, sizeof(int *));
    U_mu   = (double **)calloc(iDim + 1, sizeof(double *));
    U_min = (int **)calloc(iDim + 1, sizeof(int *));
    UT_su  = (double **)calloc(iDim + 1, sizeof(double *));
    UT_sin  = (int **)calloc(iDim + 1, sizeof(int *));
    U_su   = (double **)calloc(iDim + 1, sizeof(double *));
    U_sin = (int **)calloc(iDim + 1, sizeof(int *));
    for(int i = 1; i < iDim; ++i){

    }

    for(int i = 0; i < iDim + 1; ++i){
        UT_mu[i] = (double *)calloc(iDim + 1, sizeof(double));
        memset(UT_mu[i], 0, (iDim + 1) * sizeof(double));
        UT_min[i] = (int *)calloc(iDim + 1, sizeof(int));
        memset(UT_min[i], 0, (iDim + 1) * sizeof(int));
        U_mu[i] = (double *)calloc(iDim + 1, sizeof(double));
        memset(U_mu[i], 0, (iDim + 1) * sizeof(double));
        U_min[i] = (int *)calloc(iDim + 1, sizeof(int));
        memset(U_min[i], 0, (iDim + 1) * sizeof(int));

        UT_su[i] = (double *)calloc(iDim + 1, sizeof(double));
        memset(UT_su[i], 0, (iDim + 1) * sizeof(double));
        UT_sin[i] = (int *)calloc(iDim + 1, sizeof(int));
        memset(UT_sin[i], 0, (iDim + 1) * sizeof(int));
        U_su[i] = (double *)calloc(iDim + 1, sizeof(double));
        memset(U_su[i], 0, (iDim + 1) * sizeof(double));
        U_sin[i] = (int *)calloc(iDim + 1, sizeof(int));
        memset(U_sin[i], 0, (iDim + 1) * sizeof(int));
    }

    // 2. do elementary row operation in U^T
    block_size = 2;

    // get U^-1 directly, not percise enough
    //for(i = 1; i <= iDim; i++) {
    //    kn = rs_u[i];
    //    kp = rs_u[i+1];
    //    kp = rs_u[i+1];
    //    for(k = kn; k < kp; k ++){
    //        for(j = 1; j < i; ++j)
    //            UT_mu[j_u[k]][j] -= (u_u[k] * UT_mu[i][j]);
    //        UT_mu[j_u[k]][i] = -u_u[k];
    //    }
    //}

    // partitioned into blocks
    int blockStart, blockEnd;
    for(i = 1; i <= iDim; i++){
        kn = rs_u[i];
        kp = rs_u[i+1];
        blockStart = block_size * ((i-1) / block_size) + 1;
        blockEnd = block_size * ((i-1) / block_size + 1) + 1;

        for(k = kn; k < kp && j_u[k] < blockEnd; k++){
            for(j = blockStart; j < i; ++j)
                UT_mu[j_u[k]][j] -= (u_u[k] * UT_mu[i][j]);
            UT_mu[j_u[k]][i] = -u_u[k];

            for(j = 1; j < blockStart; ++j){
                if(UT_su[i][j] != 0)
                    UT_su[j_u[k]][j] -= (u_u[k] * UT_su[i][j]);
            }
        }
        for(; k < kp; ++k){
            UT_su[j_u[k]][i] = u_u[k];
        }
    }

    //// <!-- debug -->
    //for(int i=1; i <= iDim; ++i){
    //    printf("%4d: ", i);
    //    for(j=1; j <= iDim; ++j){
    //        if(UT_mu[i][j] != 0.0)
    //            printf("%15le %4d ", UT_mu[i][j], j);
    //    }
    //    putchar('\n');
    //}

    // remove zeros
    for(i = 1; i <= iDim; i++) {
        blockStart = block_size * ((i-1) / block_size) + 1;

        for(k = blockStart, j = 1; k <= i; ++j, ++k){
            for(; k <= i && UT_mu[i][k] == 0.0; ++k);
            if(k > i)
                break;
            UT_mu[i][j] = UT_mu[i][k];
            UT_min[i][j] = k;
        }
        UT_min[i][j] = 0;

        for(k=1, j=1; k < blockStart; ++j, ++k){
            for(; k < blockStart && UT_su[i][k] == 0.0; ++k);
            if(k >= blockStart)
                break;
            UT_su[i][j] = UT_su[i][k];
            UT_sin[i][j] = k;
        }
        UT_sin[i][j] = 0;
    }

    //<!-- debug -->
    //for(int i=1; i <= iDim; ++i){
    //    printf("%4d: ", i);
    //    for(j=1; UT_min[i][j] != 0 && j <= iDim; ++j)
    //        printf("%15le %4d   ", UT_mu[i][j], UT_min[i][j]);
    //    putchar('\n');
    //}
    //printf("-----\n");
    //for(int i=1; i <= iDim; ++i){
    //    printf("%4d: ", i);
    //    for(j=1; UT_sin[i][j] != 0 && j <= iDim; ++j)
    //        printf("%15le %4d   ", UT_su[i][j], UT_sin[i][j]);
    //    putchar('\n');
    //}

    // 2. do elementary column operation
    //for(i = 1; i <= iDim; i++){

    //}
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

    //// <!-- debug -->
    //static bool flag = true;
    //if(flag)
    //for(int i=0; i<=iDim; ++i){
    //    ks = rs_u[i];
    //    ke = rs_u[i+1];

    //    printf("%d, %d, %d: ", i, ks, ke);
    //    for(k = ks; k < ke; k ++)
    //        printf("%15le %d  ", u_u[k], j_u[k]);
    //    putchar('\n');
    //}
    //flag = false;

    #pragma omp parallel for private(i) shared(x, b)
    for(i = 1; i <= iDim; i++)  // x <- b
        x[i] = b[i];

    // ->> modified
    int blockEnd;
    #pragma omp parallel for private(i, j) shared(x, b, UT_mu, UT_min)
    for(i = 1; i <= iDim; ++i){
        for(j = 1; UT_min[i][j] != 0; ++j)
            x[i] += UT_mu[i][j] * b[UT_min[i][j]];
    }

    for(i = 1; i <= iDim; i += block_size){
        blockEnd = min(block_size * ((i-1) / block_size + 1) + 1, iDim + 1);

        // partition parallel
        //#pragma omp parallel for private(k, j) shared(blockEnd, i, x, UT_sin, UT_su)
        for(k = i; k < blockEnd; ++k){
            for(j = 1; UT_sin[k][j] != 0; ++j)
                x[k] -= UT_su[k][j] * x[UT_sin[k][j]];
        }
    }

    //// ->> original
    //for(i = 1; i <= iDim; i++)
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

    for(i=iDim-1; i>=1; i--)
    {
        ks = rs_u[i];
        ke = rs_u[i+1] - 1;
        xc = x[i];

        //#pragma omp parallel for reduction(-:xc) private(k) shared(u_u, j_u, x)
        for(k=ke; k>=ks; k--)
        {
            xc -= u_u[k] * x[j_u[k]];
        }
        x[i] = xc;
    }
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

