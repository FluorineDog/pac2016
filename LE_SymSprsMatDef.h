//////////////////////////////////////////////////////////////////////
// ��Ȩ (C), 1988-1999, �й����Ժ
// �� �� ��: LE_SymSprsMatDef.h
// ��    ��: ���Է�����֮ϡ����ⷨ
//           �������ṹ�弰ȫ�ֱ�������
// ��    ������ #include <filename.h> ��ʽ�����ñ�׼���ͷ�ļ������������ӱ�׼��Ŀ¼��ʼ������
//           �� #include "filename.h" ��ʽ�����÷Ǳ�׼���ͷ�ļ��������������û��Ĺ���Ŀ¼��ʼ������
// �޸ļ�¼:     // �޸���ʷ��¼�б�ÿ���޸ļ�¼Ӧ�����޸����ڡ��޸�
// �߼��޸����ݼ���
// 1. ʱ��:
//    ����:
//    �޸�����:
// 2. ...
//////////////////////////////////////////////////////////////////////

#ifndef _LE_SYMSPRSMATDEF_H__
#define _LE_SYMSPRSMATDEF_H__

#include <stdio.h>

// ϡ���άȫ����ṹ
typedef struct
{
    int iDim       ; // ����ά������MAXN��
    int iNy        ; // ����Ԫ��ʵ����Ŀ����NY1��NU��
    int iNymax     ; // ����Ԫ�������Ŀ,iNy = iDim*(iDim+1)/2������ά��ʱʹ�� ��husixu: ʵ������(iDim+1)*(iDim+2)/2��
    int *piJno     ; // ÿ������Ԫ���кţ�ά��iNymax+1����JNOY1��JNOU��
    int *piIstart  ; // ÿ�о���Ԫ����iJno�е���ʼλ�ã�ά��iDim+2����IYD1��IYDU��
    int *piIdiag   ; // �Խ���Ԫ�ص�λ�á�û�жԽ���Ԫ��ʱ��ָ��Aij��j>i����ӽ�,ά��iDim+1����IYD1��IYDU��
    int *piLinkp   ; // λ��������, ������ͬһ����һ��Ԫ�ص�λ�ã��������һ��Ԫ�ص�linkpָ���һ��Ԫ��,ά��iNymax+1��LP1��LPU��
    int *piLinkn   ; // �к����������¸�Ԫ�ص��кţ��������һ��Ԫ�ص�linknΪ��һԪ���к�,ά��iNymax+1����LR1��LRU��
}SprsMatStru;

// ʵ������
typedef struct
{
    int    iNy;        // ����Ԫ����Ŀ
    double *pdVal;     // ����Ԫ��ֵ��ά��iNy����I0��
}VecRealStru;

// ϡ��ȫʵ������
typedef struct
{
    SprsMatStru Mat;       // ����ṹ
    double      *pdVal;    // ����Ԫ��ֵ��ά��iNy
}SprsMatRealStru;

//U'DU�ֽ���U��ṹ����
typedef struct
{
    int iDim;           //����ά��
    int iNzs;           //�����ǣ��������Խ�Ԫ)�ķ���Ԫ�ظ���
    int *rs_u;          //����������ʼλ��,iDim+2
    int *cs_u;          //����������ʼλ��,iDim+2
    int *r_u;           //����������Ԫ���к�,iNzs+1
    int *j_u;           //����������Ԫ���к�,iNzs+1
} SprsUMatStru;

//LU�ֽ���U��ֵ����
typedef struct
{
    SprsUMatStru uMax;  //����ṹ
    double *d_u;        //�ֽ��Խ�Ԫֵ,iDim+1ά
    double *u_u;        //������Ԫ��ֵ(����),�ѹ�һ��,iNzs+1ά

    //������������ΪLU��ֵ�ֽ�ʱ����Ҫ�Ĺ�������
    int    *nzs;        //�������� iDim+1
    double *work;       //�������� iDim+1
} SprsUMatRealStru;

#endif

