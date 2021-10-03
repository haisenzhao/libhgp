//#include "cgal.h"
#include "geom.h"

using namespace std;

void Test_CGAL_GEOM() 
{
#ifdef CGAL_ONCE
	//HMODULE hModule = Functs::LoadHMODULE("ppgl.dll");
	//auto read_mesh = (CGAL_Vector_Base)GetProcAddress(hModule, "CGAL_Vector_Base");
	//Vector3d result;
	//read_mesh(Vector3d(1.0, 323.0, 0.0), result);
	//std::cerr << result[0] << " " << result[1] << " " << result[2] << std::endl;
#endif

#ifdef geom_hpp
	Vector3d result;
	CGAL_Vector_Base(Vector3d(1.0, 323.0, 0.0), result);
	std::cerr << result[0] << " " << result[1] << " " << result[2] << std::endl;

#endif

};
/*
void PostProcess() 
{
	system("copy E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");
	system("copy E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");

	//read
	VectorStr1 funct_titles;
	VectorStr1 funct_paras;
	std::ifstream file("E:\\Task2\\personal-pack-geom-lib\\ppgl\\geom.h");

	for (std::string line; std::getline(file, line); )
	{
		if (Functs::StringContain(line, "CGAL"))
		{
			auto title = line.substr(line.find("CGAL"), line.find("(") - line.find("CGAL"));
			auto para = line.substr(line.find("("));
			funct_titles.push_back(title);
			funct_paras.push_back(para);
			std::cerr << "Function: " << title <<" para: " <<para<< std::endl;
		}
	}
	file.close();


	std::ofstream cgal_file("E:\\Task2\\personal-pack-geom-lib\\ppgl\\cgal.h");
	cgal_file << "#ifndef CGAL_ONCE" << std::endl;
	cgal_file << "#define CGAL_ONCE" << std::endl;
	cgal_file << "#pragma once" << std::endl;
	cgal_file << "#include <pgl_functs.hpp>" << std::endl;
	cgal_file << "using namespace std;" << std::endl;
	cgal_file << "using namespace PGL;" << std::endl;

	//typedef void (*CGAL_Vector_Base)(Vector3d n, Vector3d&);
	for (int i = 0; i < funct_titles.size(); i++)
	{
		auto title = funct_titles[i];
		auto para = funct_titles[i];
		cgal_file << "typedef void (*" << title << ")" << para << std::endl;
	}


	cgal_file << "#endif" << std::endl;


	cgal_file.close();

};

*/
int main()
{
	//PostProcess();

	
	system("pause");
	return 0;
}

