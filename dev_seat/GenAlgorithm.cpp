//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  GenAlgorithm.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "utility.h"
#include "GenAlgorithm.h"

#include <algorithm>
#include <fstream>
#include <time.h>
#include <cmath>

#define MAX(x,y) x>=y?x:y
#define MIN(x,y) x>=y?y:x


// 随机搜索下， max j1=40000， j2理论值 max=1944
#define TOP_J1 30000
#define TOP_J2 3000

#define best_J1 734
#define best_J2 128



//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	遗传算法类
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

// sort J1 & J2
bool CompareJ1(Chrom &a, Chrom &b)
{
	return a.fit[0] > b.fit[0];
}
bool CompareJ2(Chrom &a, Chrom &b)
{
	return a.fit[1] > b.fit[1];
}


CGenAlgorithm::CGenAlgorithm(int iter, int loop, int pop_size)
	: m_iter(iter), m_loop(loop), m_pop_size(pop_size),
	bits_one(8), bits_table(4), bits_seat(4), 
	best_fit(4, INT_MAX), temp_best(4), best_chrom(4), temp_chrom(4)
{
	data.ReadDataFromTXT("./input.txt");
	f_size		= data.families.size();
	chrom_len	= f_size * 2 * bits_one;

	check = 0;

}


Chrom CGenAlgorithm::GenAlgorithm()
{
	long start_t = clock(); 
	// 分别计算四类结果
	for(int J=0; J!=2;++J){
		for(int wife_husband=0; wife_husband!=2; ++wife_husband){
			// best J1, J2
			std::vector<int> best(2, INT_MAX);
			//init chroms
			InitPopulation(chrom_len, m_pop_size, wife_husband);
			int loop(0);
			for(int times=0; times!=m_iter; ++times){
				int m = J*2 + wife_husband;

				// 计算 fit value，并对求值之后的 chroms 对 fit value的降序排序
				std::vector<int> val = func_fitness(chroms, wife_husband, J);

				// 超过一定次数不迭代，则重置种群
				if( (MIN(best[0], val[0])) == best[0]){
					if(++loop >= m_loop){
						srand(clock());
						InitPopulation(chrom_len, m_pop_size, wife_husband);
						val = func_fitness(chroms, wife_husband, J);
						loop = 0;
						long now_t  = clock();
						long cost_t = now_t - start_t;
						start_t = clock();
						ut::OutputLine(times, "times best plan",m, "is", best_fit[m]);
						ut::OutputLine("time consume", cost_t, "seconds:", (double)cost_t/CLOCKS_PER_SEC);
					}
					best[0] = MIN(best[0], val[0]);
					best[1] = MIN(best[1], val[1]);
				}
				else{
					best[0] = MIN(best[0], val[0]);
					best[1] = MIN(best[1], val[1]);

				}

				PickChroms(chroms);		// 择优染色体
				CrossChroms(chroms, wife_husband);	// 交叉染色体

				//变异，设定极小的概率
				if(1){
					MutationChrom(chroms, wife_husband);
				}

				
				if(best_fit[m] > temp_best[m]){
					// 出现更好数据，更新值和染色体
					best_fit[m] = temp_best[m];
					best_chrom[m] = temp_chrom[m];

					// 输出值和染色体到txt
					std::map<int, std::map<int, int> > seat_list;
					for(int i=0; i!=data.animosity.size(); ++i){
						seat_list[i][0] = 0;
						seat_list[i][1] = 0;
					}
					int tt = decode(best_chrom[m].bit, seat_list, wife_husband);

					// 输出到屏幕
					ut::OutputLine("new best plan", m, "is", best_fit[m]);

					// 输出值和染色体到txt


					std::fstream fid, fbest;
					std::stringstream ss_temp;
					ss_temp << "./result_"<< J<< "_" << wife_husband << ".txt";
					std::string file_temp;
					file_temp = ss_temp.str();
					fbest.open(ss_temp.str().c_str(), std::fstream::in | std::fstream::out);
					int best_val(INT_MAX);
					fbest >> best_val;
					if( best_val > best_fit[m] ){
						fid.open(ss_temp.str().c_str(), std::fstream::out);
						// 最佳值
						fid << best_fit[m] << "\r\n\r\n";


						// 座位次序
						for(int i=0; i!=seat_list.size(); ++i){
							fid << seat_list[i][0] << "\t";
							fid << seat_list[i][1] << "\t";
							fid << "\r\n";
						}
						fid.close();
					}
					fbest.close();
				}
			

			}
		}
	}
	return chroms[0]; // 返回最后的优选值，或者返回运行期间的最佳值作为最终结果
}

// 产生 num 个随机数列， 0:max-1， num 可以 = max
vInt CGenAlgorithm::RandVector(int num, int max)
{
	//srand(clock());
	vInt ret, temp;
	for(int i=0; i!=max; ++i){
		temp.push_back(i);
	}
	for(int i=0; i!=num; ++i){
		int rand_int = rand()%(temp.size());
		ret.push_back(temp[rand_int]);
		temp.erase(temp.begin()+rand_int);
	}
	return ret;
}

bool CGenAlgorithm::CheckChrom(std::string s, int wife_husband)
{
	// seat_list : <table-seat, person_num>, if person_num>1, two person in a seat
	std::map<int, int> seat_list; 
	for(int i=0; i!=f_size*2; ++i){
		std::string table, seat;
		table = s.substr(i*bits_one, bits_table);
		seat  = s.substr(i*bits_one + bits_table, bits_seat);
		std::bitset<32> bit_t(table), bit_s(seat);

		// bit_t < K
		if(bit_t.to_ulong() > data.K){
			//ut::OutputLine("table number > K");
			return false;
		}
		if(bit_s.to_ulong() > 9){
			//ut::OutputLine("seat number > 10");
			return false;
		}
		seat_list[bit_t.to_ulong()*10 + bit_s.to_ulong()*1] += 1;
	}

	// map.size() == data.families.size()*2
	if(seat_list.size() != f_size*2){
		//ut::OutputLine("one seat > 2 person");
		return false;
	}

	for(int i=0; i!=data.K; ++i){
		if(!seat_list.count(i*10)){
			//ut::OutputLine("miss one head seat");
			return false;
		}
	}

	// check if wife and husband in one table
	if(wife_husband){
		for(int i=0; i!=data.families.size(); ++i){
			std::string w_table, h_table;
			w_table = s.substr(2*i*bits_one, bits_table);
			h_table = s.substr((2*i+1)*bits_one, bits_table);
			if(w_table != h_table){
				//ut::OutputLine("wife and husband isn't in the same table");
				return false;
			}
		}
	}

	return true;
}


Chrom CGenAlgorithm::CreateOneChrom(int length, int wife_husband)
{
	//srand(clock());
	Chrom new_chrom;
	new_chrom.len = length;
	while(1){
		new_chrom.bit.clear();
		int seats(f_size*2); //
		// 产生首席的 的可行解 K*10, 一串 w.h w.h
		vInt seat_list(seats, 0); // 1-12 人的座位号,eg: 02 - 0桌 2号

		if(wife_husband){
			/**
				some one seat in head, couples seat in same table
			*/
			// list of each table
			std::vector<vInt> each_table(data.K); 
			// random to pick table
			vInt table_rand = RandVector(f_size, f_size);	
			// random pick a family to each table
			for(int i=0; i!=data.K; ++i){
				each_table[i].push_back(2*table_rand[i]);
				each_table[i].push_back(2*table_rand[i]+1);
			}
			// each rest family pick a table
			for(int i=data.K; i!=table_rand.size(); ++i){
				int table = rand()%data.K;
				each_table[table].push_back(2*table_rand[i]);
				each_table[table].push_back(2*table_rand[i]+1);
			}

			// assign head and order for each table
			for(int i=0; i!=data.K; ++i){
				// assign head
				int seat_head = rand()%each_table[i].size();
				vInt seat_rand = RandVector(9, 9);
				for(int j=0; j!=each_table[i].size(); ++j){
					if(j==seat_head){
						// who seat in head
						seat_list[each_table[i][j]] = i*10;
					}
					else{
						// the rest seat by seat_rand
						seat_list[each_table[i][j]] = i*10 + (seat_rand[j] + 1);
					}
				}
			}


		}
		else{
			/**
				some one seat in head, couples don't need to seat in same table
			*/
			// random pick K person seat in head of each table
			vInt each_head = RandVector(data.K, seats); // 4/12
			// a random list for the rest
			vInt rest_list = RandVector(seats, data.K*10);
			for(int i=0; i!=seats; ++i){
				std::vector<int>::iterator iter = find(each_head.begin(), each_head.end(), i);
				if(iter!=each_head.end()){
					// assign head seat
					seat_list[i] = (iter-each_head.begin())*10;
				}
				else{
					int temp(0);
					while(temp%10 == 0){
						// remove head seat
						temp = rest_list.back();
						rest_list.pop_back();
					}
					seat_list[i] = temp;
				}
			}

		}

		// seat_list → one chrom
		for(int j=0; j!=seats; ++j){ // each chrom
			std::bitset<4> bit_t(seat_list[j]/10), bit_s(seat_list[j]%10);
			new_chrom.bit += bit_t.to_string() + bit_s.to_string();
		}

		if(CheckChrom(new_chrom.bit, wife_husband)){
			break;
		}
	}
	return new_chrom;
}


void  CGenAlgorithm::InitPopulation(int length, int pop_size, int wife_husband)
{
	chroms.clear();
	for(int i=0; i!=pop_size; ++i){ // each pop_size
		Chrom new_chrom;
		new_chrom = CreateOneChrom(length, wife_husband);
		chroms.push_back(new_chrom);
	}
}

bool CGenAlgorithm::encode(std::string &s, std::map<int, std::map<int, int> > &seat_list)
{
	// 从 families → string
	s.clear();
	for(int i=0; i!=seat_list.size(); ++i){
		std::bitset<4> bit_t(seat_list[i][0]), bit_s(seat_list[i][1]);
		s += bit_t.to_string() + bit_s.to_string();
	}
	return true;
}
bool CGenAlgorithm::decode(std::string &s, std::map<int,std::map<int, int> > &seat_list, int wife_husband)
{
	// map_index 0→f_size ， wife → husband
	// ( map_index[i][0], map_index[i][1] )   means the ith person's table-seat
	if(!CheckChrom(s, wife_husband)){
		return false;
	}
	// 从 string → families
	for(int i=0; i!=seat_list.size(); ++i){
		std::string table, seat;
		table = s.substr(i*bits_one, bits_table);
		seat  = s.substr(i*bits_one+bits_table, bits_seat);      

		std::bitset<32> bit_t(table), bit_s(seat);

		seat_list[i][0] = bit_t.to_ulong();
		seat_list[i][1] = bit_s.to_ulong();
	}
	return true;
}

vInt CGenAlgorithm::func_fit(Chrom & chrom, int wife_husband)
{
	vInt val_fit;

	std::vector<vInt> &ani_arr = data.animosity;

	// ( seat_list[i][0], seat_list[i][1] )   i-th person table-seat
	std::map<int, std::map<int, int> > seat_list;
	for(int i=0; i!=ani_arr.size(); ++i){
		seat_list[i][0] = 0;
		seat_list[i][1] = 0;
	}
	
	// chrom to seat_list
	//  error
	int tt = decode(chrom.bit, seat_list, wife_husband);

	tt = 1;

	// compute y_k[i][j] in each_table
	std::vector<CTable> each_table(data.K);

	for(int i=0; i!=seat_list.size(); ++i){
		each_table[seat_list[i][0]].personID.push_back(i);
		each_table[seat_list[i][0]].seat.push_back(seat_list[i][0]*10 + seat_list[i][1]);
	}
	int ttt = each_table[0].personID.size();


	ttt = 0;

	for(int k=0; k!=each_table.size(); ++k){
		each_table[k].label = k;
		int size = each_table[k].personID.size();
		each_table[k].Yij.resize(size, std::vector<int>(size, 0));
		for(int i=0; i!=size; ++i){
			int IDi, IDj;
			int Ri, Rj, dis;
			IDi = each_table[k].personID[i];
			Ri  = each_table[k].seat[i]%10==0? 5 : 1;
			for(int j=0; j!=size; ++j){
				IDj = each_table[k].personID[j];
				Rj = each_table[k].seat[j]%10==0? 5 : 1;
				int abs_dis = abs(each_table[k].seat[i] - each_table[k].seat[j]);
				dis = MAX(abs_dis, 10-abs_dis);

				each_table[k].Yij[i][j] = (Ri*ani_arr[IDi][IDj] + Rj*ani_arr[IDj][IDi])*pow(dis, 2);
			}
		}

	}

	// J1, J2
	int J2(0), J1(0);
	for(int k=0; k!=each_table.size(); ++k){
		for(int i=0; i!=each_table[k].Yij.size(); ++i){
			for(int j=0; j!=each_table[k].Yij.size(); ++j){
				J1 += each_table[k].Yij[i][j];
				J2 = MAX(J2, each_table[k].Yij[i][j]);
			}
		}
	}

	val_fit.push_back(J1);
	val_fit.push_back(J2);

	return val_fit;
}

vInt CGenAlgorithm::func_fitness(vChrom &chroms, int wife_husband, int J/* =1 */)
{
	// best fit for this term's chroms
	vInt best_fit(2, INT_MAX);
	for(vChrom::iterator iter=chroms.begin(); iter!=chroms.end(); ++iter){
		// get the minimize of J1 & J2
		std::vector<int> temp = func_fit(*iter, wife_husband);
		// convert to maximize for iter
		iter->fit[0] = TOP_J1 - temp[0];
		iter->fit[1] = TOP_J2 - temp[1];
		// if best return the minimize
		int temp0 = MIN(best_fit[0], temp[0]);
		int temp1 = MIN(best_fit[1], temp[1]);
		if(J==0){
			if(temp0 != best_fit[0]){
				temp_chrom[J*2+wife_husband] = *iter;
				temp_best[J*2+wife_husband]  = temp0;
			}
			best_fit[0] = temp0;
			best_fit[1] = temp1;
		}
		if(J==1){
			if(temp1 != best_fit[1]){
				temp_chrom[J*2+wife_husband] = *iter;
				temp_best[J*2+wife_husband]  = temp1;
			}
			best_fit[0] = temp0;
			best_fit[1] = temp1;
		}
	}

	// sort the max for Ji,  0->J1, 1->J2
	if(J == 0){
		std::sort(chroms.begin(), chroms.end(), CompareJ1);
	}
	else{
		std::sort(chroms.begin(), chroms.end(), CompareJ2);
	}	
	return best_fit;
}


void CGenAlgorithm::PickChroms(vChrom &chroms, int J /* =1 */)
{
	// int J = 0; // or 1

	// sum of fit(has converted to max)
	int sum(0);
	for(vChrom::iterator iter=chroms.begin(); iter!=chroms.end(); ++iter){
		sum += iter->fit[J];
	}
	// 确定比例
	for(vChrom::iterator iter=chroms.begin(); iter!=chroms.end(); ++iter){
		iter->r_fit[J] = (double)iter->fit[J]/(double)sum * 100;
		if(iter!=chroms.begin()){
			iter->c_fit[J] += (iter-1)->c_fit[J] + iter->r_fit[J];
		}
		else{
			iter->c_fit[J] = iter->r_fit[J];
		}
	}

	// 根据比例随机选取 size（） 个样本, 轮盘大的选取几率大
	vChrom new_chroms;
	//srand(clock());
	for(int i=0; i!=chroms.size(); ++i){
		int rate = rand()%100;
		int pick(0);
		while(1){
			if(rate <= chroms[pick].c_fit[J]){
				break;
			}
			++pick;
		}
		new_chroms.push_back(chroms[pick]);
	}
	chroms.swap(new_chroms);
}


void CGenAlgorithm::CrossChrom(Chrom &chrom_a, Chrom &chrom_b, int wife_husband)
{
	// random pick i-th person， 0:f_size*2
	// 随机选择交叉点 0:bits_one
	//srand(clock());
	vInt cross_list = RandVector(rand()%f_size*2+1, f_size*2);

	for(int i=0; i!=cross_list.size(); ++i){
		vInt points = RandVector(2, bits_one);
		std::string s_a, s_b;
		int s_start  = MIN(points[0], points[1]);
		int s_length = abs(points[0] - points[1]);
		s_a = chrom_a.bit.substr(cross_list[i]*bits_one+s_start, s_length);
		s_b = chrom_b.bit.substr(cross_list[i]*bits_one+s_start, s_length);
		chrom_b.bit.replace(cross_list[i]*bits_one+s_start, s_length, s_a);
		chrom_a.bit.replace(cross_list[i]*bits_one+s_start, s_length, s_b);
	}

	if(!CheckChrom(chrom_a.bit, wife_husband)){
		Chrom new_chrom = CreateOneChrom(chrom_a.len, wife_husband);
		chrom_a.bit = new_chrom.bit;
		++check;
	}
	if(!CheckChrom(chrom_b.bit, wife_husband)){
		Chrom new_chrom = CreateOneChrom(chrom_b.len, wife_husband);
		chrom_b.bit = new_chrom.bit;
		++check;
	}

}
void CGenAlgorithm::CrossChroms(vChrom &chroms, int wife_husband)
{
	// 产生不相同的 chroms.size()，两两组合
	vInt cross_list = RandVector(chroms.size(), chroms.size());

	for(int i=0; i!=chroms.size()/2; ++i){
		CrossChrom(chroms[i*2], chroms[i*2+1], wife_husband);
	}
}

void CGenAlgorithm::MutationChrom(vChrom &chroms, int wife_husband)
{
	//srand(clock());
	if(rand()%1000 == 1){
		Chrom *p = &chroms[rand()%chroms.size()];
		int mutation = rand()%chroms[0].len;
		p->bit[mutation] = p->bit[mutation]=='0'?'1':'0';

		if(!CheckChrom(p->bit, wife_husband)){
			Chrom new_chrom = CreateOneChrom(p->len, wife_husband);
			p->bit = new_chrom.bit;
		}
	}

	
}
























//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
