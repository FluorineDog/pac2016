����Intel�ں�ƽ̨�ĶԳ�ϡ�����Է�������Ⲣ���Ż�
1��������ܣ�
�����ṩ����Դ�ļ���LUSolve.cpp,LE_SymSprsMats.cpp�Լ�һ��Makefile�ļ���
���У�
a)LUSolve.cpp����������룬�ṩ���ݶ�д�����ü��㺯���������������ļ����벻���Ż����ݲ��֡�
b)LE_SymSprsMats.cpp��ϡ����������㺯�����ǻ����ں�ƽ̨�Ż����ݡ�
ֱ��make�������ɿ�ִ���ļ�test1��ֱ�����п�ִ���ļ��ó������ʱ��
���룺
$make
icc -c LUSolve.cpp
icc -c LE_SymSprsmats.cpp
icc LUSolve.o  LE_SymSprsmats.o -o test1
���У�
./test1
The program elapsed :   ����ִ��ʱ�䵥λ��

2����������:
 a) �Ż����ݲ�����������Ҷ���Ķ�д���̲�����
 b) �Ż������޶���ĳһ���Ҷ���B�����ķ�������ⲿ�֣�����չ��ѭ��������ѭ���ϲ�ͳһ�����Ż�;
 c) �����޸�LU_SymbolicSymG��LU_NumbericSymG��LE_FBackwardSym���������ӿڲ�������͸�ʽ��
 d) �Ż�ƽ̨�޶���Intel �ں�ƽ̨�ϡ�