#ifndef LE_SYM_SPRS_MAT_AUX_FUNC_H_
#define LE_SYM_SPRS_MAT_AUX_FUNC_H_
#include <algorithm>
#include <math.h>
#include "LE_SymSprsMatDef.h"
#include "LE_SymSprsMatFunc.h"

// 函 数 名:          // initMem_MatReal
// 描    述:          // 稀疏全实数矩阵内存初始化。数目、指针变量置零
// 如参数非常复杂应举例说明。
void initMem_MatReal(SprsMatRealStru *A) {
  // 内存初始化。数目、指针变量置零
  A->Mat.iDim = 0;
  A->Mat.iNy = 0;
  A->Mat.piIdiag = NULL;
  A->Mat.piIstart = NULL;
  A->Mat.piJno = NULL;
  A->Mat.piLinkn = NULL;
  A->Mat.piLinkp = NULL;
}

void aluss(SprsUMatRealStru *pFU) {
#ifdef TRASH_CODE
  constexpr int BLOCK = 36;
  int *rs_u = pFU->uMax.rs_u;
  int *j_u = pFU->uMax.j_u;
  int iDim = pFU->uMax.iDim;
  int barrier = 1;
  std::map<int, int> countTable;
  for (int i = 1; i <= iDim; ++i) {
    int kbeg = rs_u[i];
    int kend = rs_u[i + 1];
    if (i - barrier >= BLOCK) {
      barrier += BLOCK;
    }
    cerr << i << ": ";
    if (kbeg < kend) {
      int j = j_u[kbeg];
      // cerr << j - i << ": ";
      ++countTable[j - i];
    }
    for (int k = kbeg; k < kend; k++) {
      int j = j_u[k];
      cerr << j - i << " ";
    }
    cerr << endl;
  }
  int sum100 = 0;
  constexpr int barrier2 = 33;
  for (auto p : countTable) {
    if (p.first < barrier2) {
      cerr << p.first << ":: " << p.second << endl;
    } else {
      ++sum100;
    }
  }
  cerr << barrier2 << "+:: " << sum100 << endl;





  int sumWork = 0;
  int seqWork = 0;
  int paraWork = 0;
  using Tp = std::tuple<int, int, int>;
  std::map<Tp, int> fuck;
  for (int i = iDim; i > 0; --i) {
    // cerr << "para: " << paraPart[i].size();
    // cerr << "   seq:" << seqPart[i].size();
    // cerr << endl;
    ++fuck[std::make_tuple(seqPart[i].size() + paraPart[i].size(), seqPart[i].size(), paraPart[i].size())];
    seqWork += seqPart[i].size();
    paraWork += paraPart[i].size();
  }

  int sumCount = 0;
  for(auto pp:fuck){
    cerr << "\t\tsumCount:" << (sumCount+=pp.second);
    cerr << "\t\tsum:" << std::get<0>(pp.first);
    cerr << "\t\tseq:" << std::get<1>(pp.first);
    cerr << "\t\tpara:" << std::get<2>(pp.first);
    cerr << "\t\tcount:" << pp.second;
    cerr << endl;
  }
  sumWork = seqWork + paraWork;

  cerr << endl;
  cerr << "paraBlock = " << BLOCK << endl;
  cerr << "seqWork = " << seqWork << endl;
  cerr << "paraWork = " << paraWork << endl;
  cerr << "sumWork = " << sumWork << endl;
  cerr << endl;

  cout << "paraBlock = " << BLOCK << endl;
  cout << "seqWork = " << seqWork << endl;
  cout << "paraWork = " << paraWork << endl;
  cout << "sumWork = " << sumWork << endl;
  cout << endl;
#endif
}
// 函 数 名:          // initMem_VecReal
// 描    述:          // 稀疏实数向量内存初始化。数目、指针变量置零
// 输入参数:          // 输入参数说明，包括每个参数的作
// 用、取值说明及参数间关系。
// 如参数非常复杂应举例说明。
void initMem_VecReal(VecRealStru *V) {
  // 内存初始化。数目、指针变量置零
  V->iNy = 0;
  V->pdVal = NULL;
}

// 描    述:          // 稀疏矩阵加链程序
// 输入参数:          // SprsMatRealStru *pA：稀疏实数对阵矩阵或上三角矩阵
static void SparseMatrix_adlink(SprsMatRealStru *pA) {
  int iDim;
  int i, j, k, m, r;
  int jmn, irow, ny;

  int *wb;
  int *istart, *jno;
  int *linkp, *linkn;
  double *pdval;

  iDim = pA->Mat.iDim;
  istart = pA->Mat.piIstart;
  jno = pA->Mat.piJno;
  pdval = pA->pdVal;
  ny = pA->Mat.iNy;

  linkp = pA->Mat.piLinkp;
  linkn = pA->Mat.piLinkn;

  wb = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(wb,"SparseMatrix_adlink:wb");

  for (i = 1; i <= ny; i++) {
    linkp[i] = 0;
    linkn[i] = 0;
  }

  for (i = 1; i <= iDim; i++) {
    jmn = istart[i];
    wb[i] = jmn;
    linkp[jmn] = jmn;
    linkn[jmn] = i;
  }

  irow = 1;
  for (i = 1; i <= ny; i++) {
    if (i >= istart[irow + 1])
      irow++;

    j = jno[i];

    if (irow != j) {
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

  for (i = 1; i <= iDim; i++) {
    jmn = istart[i];
    if (jmn == linkp[jmn]) {
      linkn[jmn] = 0;
      linkp[jmn] = 0;
    }
  }

  free(wb);
  wb = NULL;
}

// 描    述:          // 根据导纳矩阵G的结构获得消去树(G对称)
// LU_SymbolicSymG调用
// 输入参数:          // 保存有矩阵结构信息的G矩阵
// 输出参数:          // 消去树信息，一维数组
// 其    他:          // 不会影响G阵中矩阵的任何信息,created by xdc 2014/6/16
static void LU_EliminationTreeG(SprsMatRealStru *pG, int *pParent) {
  int jmn;
  int i, j, t;
  int kp, kn, iDim;

  int *vParent;
  int *linkp, *linkn, *istart;

  iDim = pG->Mat.iDim;
  vParent = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(vParent,"LU_EliminationTreeG:vParent");

  linkp = pG->Mat.piLinkp;
  linkn = pG->Mat.piLinkn;
  istart = pG->Mat.piIstart;

  for (i = 1; i <= iDim; i++) {
    pParent[i] = 0;
    vParent[i] = 0;

    jmn = istart[i];
    kn = linkn[jmn];
    kp = linkp[jmn];
    if (kp == 0)
      continue;

    while (kp != jmn) {
      j = kn;
      while (vParent[j] != 0 && vParent[j] < i) {
        t = vParent[j];
        vParent[j] = i;
        j = t;
      }

      if (vParent[j] == 0) {
        vParent[j] = i;
        pParent[j] = i;
      }

      kn = linkn[kp];
      kp = linkp[kp];
    }
  }

  free(vParent);
  vParent = NULL;
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
static void LU_SetUMatECountsG(SprsMatRealStru *pG, int *pParent,
                               SprsUMatRealStru *pU) {
  int i, j;
  int jmn, kp, kn, iDim;
  int *mark, *colcnt, *rowcnt;
  int *istart, *linkp, *linkn;
  int *r_u, *j_u, *rs_u, *cs_u;
  double *d_u, *u_u;

  iDim = pG->Mat.iDim;
  pU->uMax.iDim = iDim;

  mark = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(mark,"LU_SetUMatECountsG:mark");
  colcnt = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(colcnt,"LU_SetUMatECountsG:colcnt");
  rowcnt = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(rowcnt,"LU_SetUMatECountsG:rowcnt");

  rs_u = (int *)calloc(iDim + 2, sizeof(int));
  // checkPoint(rs_u,"LU_SetUMatECountsG:rs_u");
  pU->uMax.rs_u = rs_u;

  cs_u = (int *)calloc(iDim + 2, sizeof(int));
  // checkPoint(cs_u,"LU_SetUMatECountsG:cs_u");
  pU->uMax.cs_u = cs_u;

  istart = pG->Mat.piIstart;
  linkp = pG->Mat.piLinkp;
  linkn = pG->Mat.piLinkn;

  for (i = 1; i <= iDim; i++)
    colcnt[i] = 0;

  for (i = 1; i <= iDim; i++) {
    mark[i] = i;
    rowcnt[i] = 0;

    jmn = istart[i];
    kn = linkn[jmn];
    kp = linkp[jmn];
    if (kp == 0)
      continue;

    while (kp != jmn) {
      j = kn;
      while (mark[j] != i) {
        rowcnt[i]++;
        colcnt[j]++;
        mark[j] = i;
        j = pParent[j];
      }

      kn = linkn[kp];
      kp = linkp[kp];
    }
  }

  rs_u[1] = 1;
  cs_u[1] = 1;
  pU->uMax.iNzs = 0;
  for (i = 2; i <= iDim + 1; i++) {
    pU->uMax.iNzs = pU->uMax.iNzs + colcnt[i - 1];
    rs_u[i] = rs_u[i - 1] + colcnt[i - 1];
    cs_u[i] = cs_u[i - 1] + rowcnt[i - 1];
  }

  d_u = (double *)calloc(iDim + 1, sizeof(double));
  // checkPoint(d_u,"LU_SetUMatECountsG:d_u");
  pU->d_u = d_u;

  u_u = (double *)calloc(pU->uMax.iNzs + 1, sizeof(double));
  // checkPoint(u_u,"LU_SetUMatECountsG:u_u");
  pU->u_u = u_u;

  r_u = (int *)calloc(pU->uMax.iNzs + 1, sizeof(int));
  // checkPoint(r_u,"LU_SetUMatECountsG:r_u");
  pU->uMax.r_u = r_u;

  j_u = (int *)calloc(pU->uMax.iNzs + 1, sizeof(int));
  // checkPoint(j_u,"LU_SetUMatECountsG:j_u");
  pU->uMax.j_u = j_u;

  free(mark);
  free(colcnt);
  free(rowcnt);
  return;
}

// 输入参数:          // G阵结构
// 输出参数:          // U阵结构,并申请U阵内存，包括工作相量
void LU_SymbolicSymG(SprsMatRealStru *pG, SprsUMatRealStru *pFU) {
  int p;
  int i, j, k;
  int jmn, kp, kn;
  int len, top, iDim;
  int *pParent;
  int *Lnz, *Flag, *Pattern;
  int *istart, *linkp, *linkn;
  int *rs_u, *j_u, *r_u;

  deallocate_UMatReal(pFU); //!!xdc!!20150513!!增加：解决内存泄露问题

  iDim = pG->Mat.iDim;
  pParent = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(pParent,"LU_SymbolicSymG:pParent");

  Flag = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(Flag,"LU_SymbolicSymG:Flag");

  Pattern = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(Pattern,"LU_SymbolicSymG:Pattern");

  Lnz = (int *)calloc(iDim + 1, sizeof(int));
  // checkPoint(Lnz,"LU_SymbolicSymG:Lnz");

  LU_EliminationTreeG(pG, pParent);
  LU_SetUMatECountsG(pG, pParent, pFU);

  istart = pG->Mat.piIstart;
  linkp = pG->Mat.piLinkp;
  linkn = pG->Mat.piLinkn;

  rs_u = pFU->uMax.rs_u;
  r_u = pFU->uMax.r_u;
  j_u = pFU->uMax.j_u;

  j = 1;
  for (k = 1; k <= iDim; k++) {
    jmn = istart[k];
    kn = linkn[jmn];
    kp = linkp[jmn];

    top = iDim + 1;
    Flag[k] = k;
    Lnz[k] = 0;

    if (kp == 0)
      continue;
    while (kp != jmn) {
      i = kn;
      for (len = 1; Flag[i] != k; i = pParent[i]) {
        Pattern[len++] = i;
        Flag[i] = k;
      }
      while (len > 1) {
        Pattern[--top] = Pattern[--len];
      }

      kn = linkn[kp];
      kp = linkp[kp];
    }

    for (; top <= iDim; top++) {
      i = Pattern[top];
      p = rs_u[i] + Lnz[i];
      j_u[p] = k;
      Lnz[i]++;

      r_u[j] = i;
      j++;
    }
  }

  pFU->nzs = Lnz; //直接将Lnz设置成工作数组
  pFU->work = (double *)calloc(iDim + 1, sizeof(double));

  free(Flag);
  free(Pattern);
	free(pParent);
	void AdditionLU_SymbolicSymG(SprsUMatRealStru *pFU);
	AdditionLU_SymbolicSymG(pFU);
}

// 输入参数:          // G阵结构及G阵值
// 输出参数:          // U阵中的值
void LU_NumbericSymG(SprsMatRealStru *pG, SprsUMatRealStru *pFU) {
  int i, j, k;
  int p, m, n;
  int kn, kp;
  int jmn, jmx;
  int iDim;

  double dj, dk;
  double yc, yr;

  int *nzs;
  int *istart, *linkp, *linkn;
  double *g;

  int *rs_u, *cs_u, *r_u, *j_u;
  double *d_u, *u_u;
  double *work;

  nzs = pFU->nzs;
  work = pFU->work;
  iDim = pG->Mat.iDim;

  istart = pG->Mat.piIstart;
  linkp = pG->Mat.piLinkp;
  linkn = pG->Mat.piLinkn;
  g = pG->pdVal;

  rs_u = pFU->uMax.rs_u;
  cs_u = pFU->uMax.cs_u;
  r_u = pFU->uMax.r_u;
  j_u = pFU->uMax.j_u;
  d_u = pFU->d_u;
  u_u = pFU->u_u;

  for (i = 1; i <= iDim; i++) {
    nzs[i] = 0;
    work[i] = 0.0;

    jmn = istart[i];
    kn = linkn[jmn];
    kp = linkp[jmn];

    if (kp == 0) {
      dk = g[jmn];
      if (false && dk >= 0.0 && fabs(dk) < 1.0e-20)
        dk = 1.0e-20;
      if (false && dk < 0.0 && fabs(dk) < 1.0e-20)
        dk = -1.0e-20;
      d_u[i] = dk;
      continue;
    }
    while (kp != jmn) {
      work[kn] = g[kp];
      kn = linkn[kp];
      kp = linkp[kp];
    }

    dk = g[jmn];

    jmn = cs_u[i];
    jmx = cs_u[i + 1];
    for (k = jmn; k < jmx; k++) {
      j = r_u[k];
      dj = d_u[j];

      yr = work[j];
      work[j] = 0.0;

      n = rs_u[j] + nzs[j];
      for (p = rs_u[j]; p < n; p++) {
        m = j_u[p];
        work[m] -= u_u[p] * yr;
      }

      yc = yr / dj;
      u_u[p] = yc;
      dk -= yr * yc;
      nzs[j]++;
    }

    if (false && dk >= 0.0 && fabs(dk) < 1.0e-20)
      dk = 1.0e-20;
    if (false && dk < 0.0 && fabs(dk) < 1.0e-20)
      dk = -1.0e-20;
    d_u[i] = dk;
  }
  for(int i = 1; i <= iDim; ++i){
    d_u[i] = 1 / d_u[i];
  }
	void AdditionLU_NumericSymG(SprsUMatRealStru *pFU);
	AdditionLU_NumericSymG(pFU);
}
#endif