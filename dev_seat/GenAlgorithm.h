//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  GenAlgorithm.h
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


#ifndef GenAlgorithmH
#define GenAlgorithmH
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "data.h"

#include <string>
#include <vector>
#include <bitset>
#include <algorithm>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//        NAME: func_fitness0(vChrom &chroms)
// DESCRIPTION: ������ѡֵ�ĺ���
//				��Ϊ�����������Ƚ��ɷ�������ǰ��ź���λ��Ȼ���Ż�ÿһ������λ
//   ARGUMENTS: int iter	- iteration times
//				vChrom &chroms	- population
//				int(*fitness)(Chrom& chrom)	- function to evaluate fitness of chrom, return the best
// 
// 
// 
//     RETURNS: Chrom	- the best Chrom find
//  
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class CTable{
public:
	CTable(){};
	~CTable(){};
	int label;
	std::vector<int> personID; // <person ID, seat>
	std::vector<int> seat;
	std::vector<std::vector<int> > Yij;
};


// Ⱦɫ�����ж���
class Chrom{
public:
	Chrom() : len(0), fit(2, 0), r_fit(2, 0), c_fit(2, 0), bit(""){};
	~Chrom(){};
	int len;	// length of Chrom, == bit.length();
	std::vector<int> fit;	// fitness val J1, J2
	std::vector<double> r_fit, c_fit; // rate of total fit

	std::string bit; // Ⱦɫ���ַ�
}; // �Ƿ���Ҫ����һ����Ⱥ��

typedef std::vector<Chrom>	vChrom;
typedef std::vector<int>	vInt;
bool CompareJ1(Chrom &a, Chrom &b);	// use to sort J1
bool CompareJ2(Chrom &a, Chrom &b);	// use to sort J2

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	���Ի�������ض���
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	�Ŵ��㷨��
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


class CGenAlgorithm
{
public:
	CGenAlgorithm(int iter, int random_loop, int pop_size);
	~CGenAlgorithm();

/****************************************
	����
*****************************************/
	// ѭ������
	int m_iter, m_loop, m_pop_size;
	int f_size;			// familie.size, total person = f_size * 2
	// person list like this; f1_w,f1_h ... fk_w, fk_h
	// each chrom has f_size * 2 * (8bits)
	// 4 bits for table, 4 bits for seat at that table
	int bits_one, bits_table, bits_seat;
	int chrom_len;		
	
	int check;

	// ����input
	CData data;		// input data, ID, K, families

	// member output
	vInt	best_fit, temp_best;	// (4)
	vChrom	best_chrom, temp_chrom;	// ��Ӧ��õ�Ⱦɫ��


	// member use
	vChrom chroms;	


/******************************************
	����
*****************************************/
	Chrom GenAlgorithm();
	// check valid & other function
	vInt RandVector(int num, int max);

	// �Ŵ���ز���
	// �����Ⱦɫ�������߼� 
	// wife_husband=0 ���޲�Ҫ��ͬ��
	// wife_husband=1 ����Ҫ��ͬ��
	bool CheckChrom(std::string s, int wife_husband);
	Chrom CreateOneChrom(int length, int wife_husband);	// ����һ������valid��Ⱦɫ��

	void InitPopulation(int length, int pop_size, int wife_husband); // �����ʼ��Ⱦɫ����Ⱥ

	
	bool encode(std::string &s, std::map<int, std::map<int, int> > &map_index);
	bool decode(std::string &s, std::map<int, std::map<int, int> > &map_index, int wife_husband);


	
	vInt func_fit(Chrom & chrom, int wife_husband);	// ��Ӧ���������� J1 , J2 ������ת������Ŵ����ֵ���� chroms
	vInt func_fitness(vChrom &chroms, int wife_husband, int J=1);	// ��Ӧ���������� J1 , J2 ������ת������Ŵ����ֵ���� chroms



	void PickChroms(vChrom &chroms, int J=1);	// �������̳�ȡ new Ⱦɫ��
	void CrossChrom(Chrom &chrom_a, Chrom &chrom_b, int wife_husband);	// ��������Ⱦɫ����������
	void CrossChroms(vChrom &chroms, int wife_husband); // ����Ⱦɫ���������
	void MutationChrom(vChrom &chroms, int wife_husband);		// �Խ�С�������б���


	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	//        NAME: GenAlgorithm(int iter)
	// DESCRIPTION: �Ŵ��㷨������������ںϺ���
	//				ÿ�θ�������Ⱥ��ģ���䣬ȡ���������е����ֵ
	//   ARGUMENTS: 
	//				int iter	- iteration times
	//				int length	- length of chromss
	//				int pop_size	- size of population
	//				int(*fitness)(Chrom& chrom)	- function to evaluate fitness of chrom, return the best
	//				CData data	- relationship between families
	//     RETURNS: 
	//				Chrom	- the best Chrom find
	//  
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	

};



























//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  end of file
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
