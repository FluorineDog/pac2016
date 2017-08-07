//////////////////////////////////////////////////////////////////////
// 版权 (C), 1988-1999, 中国电力科学研究院
// 文 件 名: LE_SymSprsMat.c
// 版本: V1.0       时间: 2014.12.31
// 描    述: 线性方程组模块  
//           包含LE_SymSprsMatFunc.h声明的所有子程序
// 其    他：用 #include <filename.h> 格式来引用标准库的头文件（编译器将从标准库目录开始搜索）
//           用 #include "filename.h" 格式来引用非标准库的头文件（编译器将从用户的工作目录开始搜索）
// 修改记录:       // 历史修改记录
// 1. 时间:
//    作者:
//    修改内容:
// 2. ...
//////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "LE_SymSprsMatFunc.h"
#include "LE_SymSprsMatDef.h"

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // initMem_MatReal
// 描    述:          // 稀疏全实数矩阵内存初始化。数目、指针变量置零
// 输入参数:          // 输入参数说明，包括每个参数的作
                      // 用、取值说明及参数间关系。
                      // 如参数非常复杂应举例说明。
// 输出参数:          // 对输出参数的说明。
// 返 回 值:          // 函数返回值的说明
// 其    他:          // 其它说明
//////////////////////////////////////////////////////////////////////

void initMem_MatReal(SprsMatRealStru *A)
{
	// 内存初始化。数目、指针变量置零
	A->Mat.iDim = 0;
	A->Mat.iNy = 0;
	A->Mat.piIdiag = NULL;
	A->Mat.piIstart = NULL;
	A->Mat.piJno = NULL;
	A->Mat.piLinkn = NULL;
	A->Mat.piLinkp = NULL;
}

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // initMem_VecReal
// 描    述:          // 稀疏实数向量内存初始化。数目、指针变量置零
// 输入参数:          // 输入参数说明，包括每个参数的作
                      // 用、取值说明及参数间关系。
                      // 如参数非常复杂应举例说明。
// 输出参数:          // 对输出参数的说明。
// 返 回 值:          // 函数返回值的说明
// 其    他:          // 其它说明
//////////////////////////////////////////////////////////////////////
void initMem_VecReal(VecRealStru *V)
{
	// 内存初始化。数目、指针变量置零
	V->iNy = 0;
	V->pdVal = NULL;
}

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // allocate_MatReal
// 描    述:          // 稀疏实数对称矩阵分配维数
// 输入参数:          // SprsMatRealStru *pA：稀疏实数对阵矩阵或上三角矩阵
// 输出参数:          // 无
// 返 回 值:          // 无
// 其    他:          // 调用该函数前稀疏实数对称矩阵的维数A->Mat.iDim已经确定
                      // 因可能还有注入元引入，稀疏实数对阵矩阵A的元素实际数目不确定，故分配维数时一次将维数扩够
//////////////////////////////////////////////////////////////////////
void allocate_MatReal(SprsMatRealStru *A)
{
	int dimension = 0;
	int n = 0;
  cerr << "iDim: " <<  A->Mat.iDim << endl;
  cerr << "iNymax: " <<  A->Mat.iNymax << endl;
  cerr << "iNy: " <<  A->Mat.iNy << endl;
	dimension = A->Mat.iDim + 1;
	A->Mat.iNymax = dimension * (dimension + 1)/2;  // 矩阵元素最大数目
	n = A->Mat.iNymax  + 1;
  // Guess: use rank 1 as the first element instead of rank 0
  // Fuck mathematicians
	A->Mat.piJno = (int *)calloc(n,sizeof(int));
	A->Mat.piIstart = (int *)calloc(dimension+1,sizeof(int));
	A->Mat.piIdiag = (int *)calloc(dimension,sizeof(int));
	A->Mat.piLinkp = (int *)calloc(n,sizeof(int));	
	A->Mat.piLinkn = (int *)calloc(n,sizeof(int));
	A->pdVal = (double *)calloc(n,sizeof(double));
 }

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // allocate_VecReal
// 描    述:          // 实数向量分配维数
// 输入参数:          // VecRealStru *V：实数向量
// 输出参数:          // 无
// 返 回 值:          // 无
// 其    他:          // 调用该函数前向量元素数目V->iNy已经确定
//////////////////////////////////////////////////////////////////////
void allocate_VecReal(VecRealStru *V)
{
	int m = 0;
	m = V->iNy + 1;
	V->pdVal = (double *)calloc(m,sizeof(double));
}

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // deallocate_MatReal
// 描    述:          // 指针变量内存释放
                      // 被deallocateNet()调用
// 输入参数:          // 无
// 输出参数:          // 无
// 返 回 值:          // 无
// 其    他:          // 其它说明
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
// 函 数 名:          // deallocate_VecReal
// 描    述:          // 指针变量内存释放
                      // 被deallocateNet()调用
// 输入参数:          // VecRealStru *V：
// 输出参数:          // 无
// 返 回 值:          // 无
// 其    他:          // 其它说明
//////////////////////////////////////////////////////////////////////
void deallocate_VecReal(VecRealStru *V)
{
	// 指针变量内存释放
	free(V->pdVal);
}

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // SparseMatrix_adlink
// 描    述:          // 稀疏矩阵加链程序
// 输入参数:          // SprsMatRealStru *pA：稀疏实数对阵矩阵或上三角矩阵
// 输出参数:          // 无
// 返 回 值:          // 无
// 其    他:          // 其它说明,writted by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
static void SparseMatrix_adlink(SprsMatRealStru *pA)
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
// 函 数 名:          // LU_EliminationTreeG
// 描    述:          // 根据导纳矩阵G的结构获得消去树(G对称)
                      // LU_SymbolicSymG调用
// 输入参数:          // 保存有矩阵结构信息的G矩阵
// 输出参数:          // 消去树信息，一维数组
// 返 回 值:          // 无
// 其    他:          // 不会影响G阵中矩阵的任何信息,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
static void LU_EliminationTreeG(SprsMatRealStru *pG, int* pParent)
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
// 函 数 名:          // LU_SetUMatECountsG
// 描    述:          // 确定三角分解后上三角阵元素个数
                      // LU_SymbolicSymG调用
// 输入参数:          // 消去树pParent,及G阵结构
// 输出参数:          // U阵结构中维数及其值
// 返 回 值:          // 无
// 其    他:          // 不会影响G阵中矩阵的任何信息,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
static void LU_SetUMatECountsG(SprsMatRealStru *pG,int* pParent,SprsUMatRealStru *pU)
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
// 函 数 名:          // LU_SymbolicSymG
// 描    述:          // 稀疏实数对称矩阵符号因子分解，确定U阵结构和元素个数
                      //
// 输入参数:          // G阵结构
// 输出参数:          // U阵结构,并申请U阵内存，包括工作相量
// 返 回 值:          // 无
// 其    他:          // 不会影响G阵中矩阵的任何信息,created by xdc 2014/6/16
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

   deallocate_UMatReal(pFU); //!!xdc!!20150513!!增加：解决内存泄露问题

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

   pFU->nzs=Lnz; //直接将Lnz设置成工作数组
   pFU->work=(double *)calloc(iDim+1,sizeof(double));
   //checkPoint(pFU->work,"LU_SymbolicSymG:work");

   free(Flag);
   free(Pattern);
   free(pParent);
}

//////////////////////////////////////////////////////////////////////
// 函 数 名:          // LU_NumbericSymG
// 描    述:          //稀疏实数对称矩阵数值因子分解，确定U阵元素值
                      //
// 输入参数:          // G阵结构及G阵值
// 输出参数:          // U阵中的值
// 返 回 值:          // 无
// 其    他:          // 不会影响G阵中矩阵的任何信息,created by xdc 2014/6/16
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
}

//////////////////////////////////////////////////////////////////////
// 函 数 名:          //LE_FBackwardSym
// 描    述:          //对称矩阵前推回代方法求解方程组
                      //
// 输入参数:          // U阵结构及U阵值，右端项b
// 输出参数:          // 右端项x，维数为pU的维数（解向量）
// 返 回 值:          // 无
// 其    他:          // 不会影响U阵中矩阵的任何信息,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void LE_FBackwardSym(SprsUMatRealStru *pFU,double b[],double x[])
{
   int i,j,k;
	int ks,ke;
   int iDim;
   int *rs_u,*j_u;
   double *d_u,*u_u;
   double xc;

   d_u=pFU->d_u;
   u_u=pFU->u_u;
   rs_u=pFU->uMax.rs_u;
   j_u=pFU->uMax.j_u;
   iDim=pFU->uMax.iDim;

   for(i = 1; i <= iDim; i++)
      x[i] = b[i];

	for(i = 1; i <= iDim; i++)
	{
		xc = x[i];
		ks = rs_u[i];
		ke = rs_u[i+1];

		for(k = ks; k < ke; k ++)
		{
			j = j_u[k];
			x[j] -= u_u[k]*xc;
		}
	}

   for(i=1; i<=iDim; i++)
      x[i] /= d_u[i];

   for(i=iDim-1; i>=1; i--)
   {
		ks = rs_u[i];
		ke = rs_u[i+1] - 1;
      xc = x[i];

      for(k=ke; k>=ks; k--)
      {
			j = j_u[k];
			xc -= u_u[k]*x[j];
		}
      x[i] = xc;
	}
}


//////////////////////////////////////////////////////////////////////
// 函 数 名:          //initMem_UMatReal
// 描    述:          //内存初始化。数目、指针变量置零
                      //
// 输入参数:          // U
// 输出参数:          // 
// 返 回 值:          // 无
// 其    他:          // writed by xdc 2014/6/16
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
// 函 数 名:          //deallocate_UMatReal
// 描    述:          //指针变量内存释放
                      //
// 输入参数:          // U
// 输出参数:          // 
// 返 回 值:          // 无
// 其    他:          // writed by xdc 2014/6/16
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

