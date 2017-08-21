//////////////////////////////////////////////////////////////////////
// ��Ȩ (C), 1988-1999, �й����Ժ
// �� �� ��: LE_SymSprsMatFunc.h
// ��    ��:       �汾:        ʱ��: // ���ߡ��汾���������
// ��    ��: ���Է�����
//           �������ⲿ��������
// ��    ������ #include <filename.h> ��ʽ�����ñ�׼���ͷ�ļ������������ӱ�׼��Ŀ¼��ʼ������
//           �� #include "filename.h" ��ʽ�����÷Ǳ�׼���ͷ�ļ��������������û��Ĺ���Ŀ¼��ʼ������
// �޸ļ�¼:     // �޸���ʷ��¼�б�ÿ���޸ļ�¼Ӧ�����޸����ڡ��޸�
// �߼��޸����ݼ���
// 1. ʱ��:
//    ����:
//    �޸�����:
// 2. ...
//////////////////////////////////////////////////////////////////////

#ifndef _LE_SYMSPRSMATFUNC_H__
#define _LE_SYMSPRSMATFUNC_H__

#include <stdio.h>

#include "LE_SymSprsMatDef.h"

//-----------�ⲿ��������--------------------------------------------
// ϡ��ȫʵ��������غ���
//// �ڴ��ʼ������Ŀ��ָ���������
extern void initMem_MatReal(SprsMatRealStru *A);
//// ָ������ڴ�����
extern void allocate_MatReal(SprsMatRealStru *A);
//// ָ������ڴ��ͷ�
extern void deallocate_MatReal(SprsMatRealStru *A);

// ϡ��ʵ��������غ���
//// �ڴ��ʼ������Ŀ��ָ���������
extern void initMem_VecReal(VecRealStru *V);
//// ָ������ڴ�����
extern void allocate_VecReal(VecRealStru *V);
//// ָ������ڴ��ͷ�
extern void deallocate_VecReal(VecRealStru *V);

// ϡ��������
extern void SparseMatrix_adlink(SprsMatRealStru *pA);
// �ڴ��ʼ������Ŀ��ָ���������
extern void initMem_UMatReal(SprsUMatRealStru *U);
// ָ������ڴ��ͷ�
extern void deallocate_UMatReal(SprsUMatRealStru *U);

// �ڶ���LU�ֽ��������⣺xdc�ṩ
//cmt����ʼ�γ�G��ʱ������÷������ӷֽ�LU_SymbolicSymG��ΪU�����ռ�
//cmt��ÿһʱ�������G��ṹδ�䣬ֻ��ֵ�仯������ֻ����ֵ�ֽ�LU_NumbericSymG
//ϡ��ʵ���Գƾ���������ӷֽ⣬ȷ��U��ṹ��Ԫ�ظ���
extern void LU_SymbolicSymG(SprsMatRealStru *pG,SprsUMatRealStru *pFU);
//ϡ��ʵ���Գƾ�����ֵ���ӷֽ⣬ȷ��U��Ԫ��ֵ
extern void LU_NumbericSymG(SprsMatRealStru *pG,SprsUMatRealStru *pFU);
//�Գƾ���ǰ�ƻش�������ⷽ����
extern void LE_FBackwardSym(SprsUMatRealStru *pFU,double b[],double x[]);
/////////end of xdc added 2014-6-16

#endif

