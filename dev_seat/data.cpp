//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XX  input_data.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "data.h"
// #include "utility.h"

#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

// get key field from txt
template<typename T>
inline std::string GetWord(T &fid, std::string key)
{
	std::string temp;
	std::string next, start, end;
	start = "<" + key + ">";
	end	  = "</" + key + ">";
	// <key>
	while(next != start && !fid.eof()){
		fid >> next;
	}
// 	fid >> next;
// 	temp << next;
	fid >> next;
	while(next != end && !fid.eof()){
		next = temp.empty()? next : " "+next;
		temp += next;
		fid >> next;
	}
	return temp;
	// </key>
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//        NAME: CData::ReadDataFromTXT(std::string path)
// DESCRIPTION: Input data from txt. 
//				ID, relationship between families and people
//   ARGUMENTS: std::string path - txt path
//     RETURNS: none
//  
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void CData::ReadDataFromTXT(std::string path)
{
	std::fstream fid;
	fid.open(path.c_str(), std::fstream::in);
	if(!fid){
		//ut::OutputLine("Can not open input txt.");
		return;
	}

	std::string next;
	// <ID>
	fid >> next;
	if(next != "<ID>"){
		//ut::OutputLine("Data error.");
		return;
	}
	fid >> ID;
	int x = (ID%10)%3;
	int y = (ID/10)%10%3;
	fid >> next; 
	// </ID>

	// <K>
	fid >> next;
	fid >> K;
	fid >> next;
	// </K>

	while(!fid.eof()){
		// <family>
		fid >> next;
		if(next == "<family>"){
			CFamily new_family;

			// <name>
			new_family.name = GetWord(fid, "name");

			// <short>
			new_family.short_name = GetWord(fid, "short");

			// <wife>
			new_family.wife = GetWord(fid, "wife");

			// <husband>
			new_family.husband = GetWord(fid, "husband");

			// others
			new_family.others = GetWord(fid, "others");

			// animosity
			 // w & h
			new_family.ani_w2h = atoi(GetWord(fid, "w2h").c_str());

			new_family.ani_h2w = atoi(GetWord(fid, "h2w").c_str());

			 // Inter_families
			std::stringstream Inter_family;
			 Inter_family.str(GetWord(fid, "Inter_family"));
			if(!Inter_family.str().empty()){
				while(!Inter_family.eof()){
					//new_family.ani_families.insert()
					std::string name = GetWord(Inter_family, "name");
					if(name.empty()){
						continue;
					}
					int val;
					// value has x,y
					std::stringstream ss_val;
					ss_val.str(GetWord(Inter_family, "value"));
					ss_val >> val;
					if(ss_val.str().length() > 2){
						std::string s_temp;
						ss_val >> s_temp >> s_temp;
						val += s_temp=="x"? x : y;
					}

					// insert ani_families
					
					new_family.ani_families.insert(make_pair(name, val));
				}
			}
			families.push_back(new_family);
		}
		// </family>
		fid >> next >> next;	// </animosity> </family>
	}
	fid.close();


	animosity = GetAnimosityArray();

}

std::vector<vInt> CData::GetAnimosityArray()
{
	std::vector<CFamily> &f = families;
	int arr_size = f.size()*2;
	std::map<std::string, int> map_index;
	for(int i=0; i!=f.size(); ++i){
		map_index[f[i].name] = i;
		map_index[f[i].name + "s"] = i;
	}

	// create array
	std::vector< std::vector<int> > ret;
	ret.resize(arr_size, std::vector<int>(arr_size, 0));

	// ¥” data ∏≥÷µ
	for(int i=0; i!=f.size(); ++i){
		// between wife & husband
		ret[2*i][2*i+1] = f[i].ani_w2h;
		ret[2*i+1][2*i] = f[i].ani_h2w;

		// between families
		std::map<std::string, int> &ani_f = f[i].ani_families;
		for(std::map<std::string, int>::iterator iter=ani_f.begin(); iter!=ani_f.end(); ++iter){
			// 			auto col = map_index.size();
			// 			auto row = 2*i;
			ret[2*i][map_index[iter->first]*2] = iter->second;
			ret[2*i][map_index[iter->first]*2+1] = iter->second;
			ret[2*i+1][map_index[iter->first]*2] = iter->second;
			ret[2*i+1][map_index[iter->first]*2+1] = iter->second;
		}
	}

	// ¥Ú”°≤‚ ‘
	std::fstream fid;
	fid.open("ani_array.txt", std::fstream::out);

	for(int i=0; i!=arr_size; ++i){
		for(int j=0; j!=arr_size; ++j){
			fid << ret[i][j] << "\t";
		}
		fid << "\r\n";
	}
	fid.close();
	return ret;
}


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
