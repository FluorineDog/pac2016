//////////////////////////////////////////////////////////////////////
// ��Ȩ (C), 1988-1999, XXXX��˾
// �� �� ��: Pub_ErrDefFunc.h
// ��    ��:       �汾:        ʱ��: // ���ߡ��汾���������
// ��    ��:     // ������ϸ˵���˳����ļ���ɵ���Ҫ���ܣ�������ģ��
                 // �����Ľӿڣ����ֵ��ȡֵ��Χ�����弰������Ŀ�
                 // �ơ�˳�򡢶����������ȹ�ϵ
// ��    ������ #include <filename.h> ��ʽ�����ñ�׼���ͷ�ļ������������ӱ�׼��Ŀ¼��ʼ������
//           �� #include "filename.h" ��ʽ�����÷Ǳ�׼���ͷ�ļ��������������û��Ĺ���Ŀ¼��ʼ������
// �޸ļ�¼:     // �޸���ʷ��¼�б�ÿ���޸ļ�¼Ӧ�����޸����ڡ��޸�
                 // �߼��޸����ݼ���  
// 1. ʱ��:2016.02.19
//    ����:���Ӣ
//    �޸�����:������Ϣ�ṹErrorStru���ӱ���iELevel
// 2. ...
//////////////////////////////////////////////////////////////////////
#ifndef _PUB_ERRDEFFUNC_H__
#define _PUB_ERRDEFFUNC_H__


// �ɴ�����������
#define MAXERRNUM 100
//----------ϵͳ�������Ŷ���---------------------------------------
// ϵͳ�������Ŷ���
// �ر��ļ������ڻ�򿪳����������󣬳����˳�
# define ERROR_FILEOPEN 1
// �ر��ļ����ڵ�Ϊ�գ��������󣬳����˳�
# define ERROR_FILEEMPTY 2
// Ԫ��������д����,�������ģ����ɳ�������
# define ERROR_PARA 3
// Ԫ�����������������󣬳����˳�
# define ERROR_PARAFATAL 4
// Ԫ�����Ͳ�����
# define ERROR_CTYP 5
// ��������Ų�����
# define ERROR_VARNO 6
// Ԫ���������������
# define ERROR_NOVAROUT 7
// �ڴ����ʧ��
# define ERROR_ALLOC 8
// �ڴ����·���ʧ��
# define ERROR_REALLOC 9
// �����ļ�ʧ��
# define ERROR_FILECREAT 10
// ��Ŀ������������ά�����������󣬳����˳�
# define ERROR_OUTOFMAX 11
// ��ģ����
# define ERROR_MODELING 12
// �ӳ������г���
# define ERROR_FUNC 13
// ����Ϊ0
# define ERROR_DIV 14
// ��������ʧ��
# define ERROR_INVMAT 15
// Ԫ�����˽ӵ�
# define ERROR_GNDBOTH 16
// �������죬LU�ֽ�ʧ��
# define ERROR_LUMAT 17
// ���ڳ���׼��������������������������������
# define ERROR_IOCALTYP 18
// �Զ������ʧ��
# define ERROR_UDMCAL 19
// ��ֵ�㷨��Ȼͨ�ϲ�ֵʱ�����ʧ��
# define ERROR_INTERPTIMECAL 20
// ���㲻����
# define ERROR_NOTCONVERGENT 21

//----------�洢������Ϣ����---------------------------------------
// ������Ϣ�ṹ
typedef struct
{
	int iENo;	// �����
	int iECTyp;	// ����Ԫ������
	int iECNo;	// ����Ԫ����
	int iELevel;// ����ȼ���
}ErrorStru;


//-----------�ⲿ��������--------------------------------------------
// ���������Ϣ����ֵ����
extern void clearSimError();
// ������ķ������ֵ��������Ϣ�����е�һ�������Ϊ0��λ��
extern void setSimError(ErrorStru *Err);
// ͳһ����������
extern int dealSimError();
// �ڴ����ɹ�����飬����ӡ��Ӧ��Ϣ
extern void checkPoint(void *p,char *str);

#endif