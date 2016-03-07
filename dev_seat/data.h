//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  input_data.h
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef inputH
#define inputH
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <string>
#include <vector>
#include <map>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  declarations
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

typedef std::vector<int> vInt;

class CFamily
{
public:
	CFamily(): w_table(0), w_seat(0), h_table(0), h_seat(0){};
	~CFamily(){};

	std::string name, short_name, wife, husband, others;

	int ani_w2h, ani_h2w; // animosity wife to husband, husband to wife, between families

	std::map<std::string, int> ani_families;

	// seat plan
	int w_table, w_seat, h_table, h_seat;
};


class CData
{
public:
	CData(): ID(0), K(0){};
	~CData(){};

	/**
		member
	*/
	unsigned int ID;
	unsigned int K;
	std::vector<CFamily> families;
	std::vector< vInt > animosity;
	/**
		function
	*/
	void ReadDataFromTXT(std::string path);
	std::vector<vInt> GetAnimosityArray(); // 得到全部人的憎恨矩阵

};




//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  end of file
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX








