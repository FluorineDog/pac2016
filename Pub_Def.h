//////////////////////////////////////////////////////////////////////
// ��Ȩ (C), 1988-1999, XXXX��˾
// �� �� ��: Pub_Def.h
// ��    ��:       �汾:        ʱ��: // ���ߡ��汾���������
// ��    ��: ȫ�ֺ궨��
//           ȫ�ֳ�������
//           �����Զ�������
//           ͨ���ṹ�嶨��
// ��    ������ #include <filename.h> ��ʽ�����ñ�׼���ͷ�ļ������������ӱ�׼��Ŀ¼��ʼ������
//           �� #include "filename.h" ��ʽ�����÷Ǳ�׼���ͷ�ļ��������������û��Ĺ���Ŀ¼��ʼ������
// �޸ļ�¼:     // �޸���ʷ��¼�б�ÿ���޸ļ�¼Ӧ�����޸����ڡ��޸�
                 // �߼��޸����ݼ���  
// 1. ʱ��:2014.8.4
//    ����:���Ӣ
//    �޸�����:MAXSTRLEN 256 �޸�Ϊ MAXSTRLEN 512
//    ԭ��*.TB2�ļ�һ�м�¼�ַ�������Խ��
// 2. ʱ��:2014.12.29
//    ����:���Ӣ
//    �޸�����:MAXSTRLEN 512 �޸�Ϊ MAXSTRLEN 2048
//    ԭ��*.GEN�ļ�һ�м�¼�ַ�������Խ��
// 3. ʱ��:2015.4.26
//    ����:���Ӣ
//    �޸�����:EPSINON 1.0e-10�޸�Ϊ EPSINON 1.0e-14
//    ԭ�򣺵�Ԫ��������Сʱ������ΪС���ص���·���迹���������ʽΪ��Сֵ����1e-12������Ӧ����Ϊ��������
// 4. ʱ��:2015.4.29
//    ����:���Ӣ
//    �޸�����:MAXNODE 2000�޸�Ϊ MAXNODE 4000
//             MAXNY 32000�޸�Ϊ  MAXNY 128000
// 5. ʱ��:2015.8.14
//    ����:���Ӣ
//    �޸�����:EPSINON 1.0e-14�޸�Ϊ EPSINON 1.0e-16
//    ԭ����step_GOV1()���ж�dAVZ<=dEKʱ��dAVZΪ1.0e-15��dEK=0.0ʱ���жϴ���
//////////////////////////////////////////////////////////////////////
#ifndef _PUB_DEF_H__
#define _PUB_DEF_H__

// ȫ�ֺ궨��
// ����Сֵ
#ifndef min
	#define min(a, b)  (((a) < (b)) ? (a) : (b))
#endif
// �����ֵ
#ifndef max
	#define max(a, b)  (((a) > (b)) ? (a) : (b))
#endif
// ������������������룩
#define NINT(a) ((a >= 0)? (int)(a+0.5) : (int)(a-0.5)) // ͬFortran����IDNINT

// ȫ�ֳ�������
// ������������ȣ�
#define EPSINON 1.0e-16 // 1.0e-10-->1.0e-14-->1.0e-16
// ��ֵ����
#define VALUPLIMIT 1.0e+10
// Բ����
#define PI 3.1415926535897932
// �м�������ά��
#define MAXMIDVAR 5000
// �ڵ������Ŀ
#define MAXNODE 4000
// �����ǵ��������Ԫ��Ŀ, 100���ڵ�ϵͳ, 0 Ԫռ96%,500���ڵ�,Ϊ99.2%
#define MAXNY 128000
// �ڵ���������
#define MAXOUTLINE 1000
// ���������Ŀ
#define MAXSUB 10
// �ַ�����󳤶�
#define MAXSTRLEN 2048 
// �ļ������������
#define MAXCOL 100
// LU�ֽ��U���������Ŀ
#define MAXNU 50000
// ��
#define TRUE 1
// ��
#define FALSE 0


// �����Զ�������
// �̶���С���ַ�����
typedef char FixedString[MAXSTRLEN];

// ͨ�ýṹ����
// �����ṹ
typedef struct
{
	double dRms;      // ��Чֵ
	double dTheta;	  // ���
}PhasorStru;


// ����������ṹ
//typedef struct
//{
//	double dP0;      // �й����ʳ�ֵ
//	double dQ0;	     // �޹����ʳ�ֵ
//}GenFLStru;

// Ԫ�������ṹ
typedef struct
{
	int iCTyp;	// Ԫ������
	int iCNo;	// Ԫ�����
	int iVNo;	// �������
}CVarIndexStru;

// ��ת��Ϊ��������Ԫ��������Ԫ���ṹ
typedef struct
{
	int	iNum;		// Ԫ����Ŀ
	int	*piIndexNo;// ����Ԫ��ת��Ϊ��������Ԫ��ʱ����A���ڵ���Ԫ���еı��
}ThrPhaseConvtStru;

// ����������Ϣ�ṹ
typedef struct  
{
	int iPulse;	    // ��������(01�ź�)
	double dTTag;	// ����ֶ�
}TrigPulseStru;

// ������Ԫ����Ϣ�ṹ
typedef struct  
{
	int iNo;	    // �кŻ��к�
	double dVal;	// ֵ
}YElmStru;

#endif