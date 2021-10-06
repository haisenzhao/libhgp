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

void Generate_CGAL_H() 
{
	//parse geom.h
	VectorStr1 funct_values, funct_titles, funct_paras;
	std::map<int, VectorStr1> funct_notations;
	{
		std::ifstream file("E:\\Task2\\personal-pack-geom-lib\\ppgl\\geom.h");
		for (std::string line; std::getline(file, line); )
		{
			if (Functs::StringContain(line, "CGAL"))
			{
				auto value = line.substr(line.find("PPGL_EXPORT") + 11, line.find("CGAL") - line.find("PPGL_EXPORT") - 12);
				auto title = line.substr(line.find("CGAL"), line.find("(") - line.find("CGAL"));
				auto para = line.substr(line.find("("));
				funct_values.push_back(value);
				funct_titles.push_back(title);
				funct_paras.push_back(para);
				std::cerr << "Function: " << title << " para: " << para << std::endl;
			}
			else
			{
				if (Functs::StringContain(line, "//"))
				{
					funct_notations[funct_values.size()-1].push_back(line);
				}
			}

		}
		file.close();
	}

	//if (mt.find(t1) == mt.end())

	//output cgal.h
	std::ofstream cgal_file("E:\\Task2\\personal-pack-geom-lib\\ppgl\\cgal.h");
	cgal_file << "#ifndef CGAL_ONCE" << std::endl;
	cgal_file << "#define CGAL_ONCE" << std::endl;
	cgal_file << "#pragma once" << std::endl;
	cgal_file << "#include <pgl_functs.hpp>" << std::endl;
	cgal_file << "using namespace std;" << std::endl;
	cgal_file << "using namespace PGL;" << std::endl;

	//define functions
	if (funct_notations.find(-1) != funct_notations.end())
		for (auto& line : funct_notations[-1])
			cgal_file << line << std::endl;
	for (int i = 0; i < funct_titles.size(); i++)
	{
		auto value = funct_values[i];
		auto title = funct_titles[i];
		auto para = funct_paras[i];
		cgal_file << "typedef "+value+" (*" << title << ")" << para << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				cgal_file << line << std::endl;
	}

	//define class
	cgal_file << std::endl;
	cgal_file << std::endl;
	cgal_file << "class PL" << std::endl;
	cgal_file << "{" << std::endl;
	cgal_file << "	public:" << std::endl;
	cgal_file << "	PL()" << std::endl;
	cgal_file << "	{" << std::endl;
	cgal_file << "		hModule = Functs::LoadHMODULE(\"ppgl.dll\");" << std::endl;
	for (int i = 0; i < funct_titles.size(); i++)
	{
		auto title = funct_titles[i];
		cgal_file << "		" << title << "_C = ("<<title<<")GetProcAddress(hModule, \""<<title<<"\");" << std::endl;
	}
	cgal_file << "	};" << std::endl << std::endl;
	cgal_file << "	static PL& Instance()" << std::endl;
	cgal_file << "	{" << std::endl;
	cgal_file << "		static PL instance;" << std::endl;
	cgal_file << "		return instance;" << std::endl;
	cgal_file << "	};" << std::endl << std::endl;


	cgal_file << "	HMODULE hModule;" << std::endl;

	for (int i = 0; i < funct_titles.size(); i++)
	{
		auto title = funct_titles[i];
		cgal_file <<"	" << title << " " << title << "_C;" << std::endl;
	}

	cgal_file << "};" << std::endl;


	cgal_file << "#endif" << std::endl;


	cgal_file.close();

};


int main(int argc, char* argv[])
{
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\Release\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\Release\\");

	Generate_CGAL_H();
	system("pause");
	return 0;
}

