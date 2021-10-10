#include "pgl_functs.hpp"

using namespace std;
using namespace PGL;

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
				std::cerr << "Function: " << title << " para: " << para.substr(0,6) << std::endl;
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
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				cgal_file << line << std::endl;
	}
	cgal_file << "	};" << std::endl << std::endl;
	cgal_file << "	static PL& Inst()" << std::endl;
	cgal_file << "	{" << std::endl;
	cgal_file << "		static PL instance;" << std::endl;
	cgal_file << "		return instance;" << std::endl;
	cgal_file << "	};" << std::endl << std::endl;


	cgal_file << "	HMODULE hModule;" << std::endl;

	for (int i = 0; i < funct_titles.size(); i++)
	{
		auto title = funct_titles[i];
		cgal_file <<"	" << title << " " << title << "_C;" << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				cgal_file << line << std::endl;
	}

	cgal_file << "};" << std::endl;
	cgal_file << "#endif" << std::endl;
	cgal_file.close();
};


int main(int argc, char* argv[])
{
	//copy dll
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\Release\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\Release\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\post\\Release\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\post\\Release\\");

	//generate cgal head file
	Generate_CGAL_H();

	Functs::MAssert("Successfully finishing the post processing...");
	return 0;
}

