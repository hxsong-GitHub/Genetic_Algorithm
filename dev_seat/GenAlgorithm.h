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
// DESCRIPTION: 计算优选值的函数
//				分为两部分做，先将丈夫和妻子们安排好座位，然后优化每一桌的座位
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


// 染色体序列定义
class Chrom{
public:
	Chrom() : len(0), fit(2, 0), r_fit(2, 0), c_fit(2, 0), bit(""){};
	~Chrom(){};
	int len;	// length of Chrom, == bit.length();
	std::vector<int> fit;	// fitness val J1, J2
	std::vector<double> r_fit, c_fit; // rate of total fit

	std::string bit; // 染色体字符
}; // 是否需要定义一个种群？

typedef std::vector<Chrom>	vChrom;
typedef std::vector<int>	vInt;
bool CompareJ1(Chrom &a, Chrom &b);	// use to sort J1
bool CompareJ2(Chrom &a, Chrom &b);	// use to sort J2

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	个性化函数相关定义
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	遗传算法类
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


class CGenAlgorithm
{
public:
	CGenAlgorithm(int iter, int random_loop, int pop_size);
	~CGenAlgorithm();

/****************************************
	变量
*****************************************/
	// 循环控制
	int m_iter, m_loop, m_pop_size;
	int f_size;			// familie.size, total person = f_size * 2
	// person list like this; f1_w,f1_h ... fk_w, fk_h
	// each chrom has f_size * 2 * (8bits)
	// 4 bits for table, 4 bits for seat at that table
	int bits_one, bits_table, bits_seat;
	int chrom_len;		
	
	int check;

	// 数据input
	CData data;		// input data, ID, K, families

	// member output
	vInt	best_fit, temp_best;	// (4)
	vChrom	best_chrom, temp_chrom;	// 对应最好的染色体


	// member use
	vChrom chroms;	


/******************************************
	函数
*****************************************/
	Chrom GenAlgorithm();
	// check valid & other function
	vInt RandVector(int num, int max);

	// 遗传相关操作
	// 检查是染色体否符合逻辑 
	// wife_husband=0 夫妻不要求同桌
	// wife_husband=1 夫妻要求同桌
	bool CheckChrom(std::string s, int wife_husband);
	Chrom CreateOneChrom(int length, int wife_husband);	// 创建一条符合valid的染色体

	void InitPopulation(int length, int pop_size, int wife_husband); // 随机初始化染色体种群

	
	bool encode(std::string &s, std::map<int, std::map<int, int> > &map_index);
	bool decode(std::string &s, std::map<int, std::map<int, int> > &map_index, int wife_husband);


	
	vInt func_fit(Chrom & chrom, int wife_husband);	// 适应函数，返回 J1 , J2 ，并将转换后的遗传最大值赋予 chroms
	vInt func_fitness(vChrom &chroms, int wife_husband, int J=1);	// 适应函数，返回 J1 , J2 ，并将转换后的遗传最大值赋予 chroms



	void PickChroms(vChrom &chroms, int J=1);	// 按照轮盘抽取 new 染色体
	void CrossChrom(Chrom &chrom_a, Chrom &chrom_b, int wife_husband);	// 交叉两条染色体的随机部分
	void CrossChroms(vChrom &chroms, int wife_husband); // 所有染色体随机交叉
	void MutationChrom(vChrom &chroms, int wife_husband);		// 以较小比例进行变异


	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	//        NAME: GenAlgorithm(int iter)
	// DESCRIPTION: 遗传算法和随机搜索的融合函数
	//				每次更新中种群规模不变，取整个过程中的最佳值
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
