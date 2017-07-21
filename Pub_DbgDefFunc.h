//////////////////////////////////////////////////////////////////////
// ��Ȩ (C), 1988-1999, XXXX��˾
// �� �� ��: Pub_DbgDefFunc.h
// ��    ��:       �汾:        ʱ��: // ���ߡ��汾���������
// ��    ��: ����debug��س������壺������ӡ��ʾ��Ϣ���͡����λ��
//           ����debug��ؽӿں�������
// ��    ������ #include <filename.h> ��ʽ�����ñ�׼���ͷ�ļ������������ӱ�׼��Ŀ¼��ʼ������
//           �� #include "filename.h" ��ʽ�����÷Ǳ�׼���ͷ�ļ��������������û��Ĺ���Ŀ¼��ʼ������
// �޸ļ�¼:     // �޸���ʷ��¼�б�ÿ���޸ļ�¼Ӧ�����޸����ڡ��޸�
                 // �߼��޸����ݼ���  
// 1. ʱ��:2016.02.19
//    ����:���Ӣ
//    �޸�����:�޸�printSimMsg()��������
// 2. ...
//////////////////////////////////////////////////////////////////////
#ifndef _PUB_DBGDEFFUNC_H__
#define _PUB_DBGDEFFUNC_H__

#include <stdio.h>
#include <stdlib.h>

#include "Pub_Def.h"

// ��������
// ��ӡ��ʾ��Ϣ���ͣ����󡢾��桢��ʾ��Ϣ�����ڼ����е���ʾ
# define WS_FATAL     1		// ��������
# define WS_ERROR     2		// ����
# define WS_WARNING   3		// ����
# define WS_INFO      4		// ��ʾ

// �����λ��
# define DOUT_SCREEN  1		// ��Ļ��ӡ
# define DOUT_FILE    2		// �ļ�
# define DOUT_NET     3		// ����



//-----------�ⲿ��������--------------------------------------------
// Debug�ļ��Ĵ򿪺͹ر�
extern void DebugFileOpenClose(const char *pcProjOut,int iIsOpen);
// ��������iType��ӡ������Ϣstr��λ��iOut
//extern void printString(FixedString str,int iType,int iOut);
// ��ӡdouble��һά����arr[iNB]
extern void printDArray(int iNB,const double *pdArr,const char *str,int iOut);
// ��ӡint��һά����arr[iNB]
extern void printIArray(int iNB,const int *piArr,const char *str,int iOut);
// ��ӡdouble�Ͷ�ά����arr[iNRow][iNCol]
extern void printD2Array(int iNRow,int iNCol,double **ppdArr,const char *str,int iOut);
// ��ӡ������Ϣ
//extern void printSimMsg(int iECTyp,int iECNo,int iENo,int iOut,FixedString str);
//extern void printSimMsg(char *str);
extern void printSimMsg(const char *str,...);
// ���������Ϣ
//extern void dealSimMsg(int iECTyp,int iECNo,int iENo,int iType,int iOut,FixedString str);


#endif