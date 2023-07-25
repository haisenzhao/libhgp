#include "pgl_functs.hpp"

#include <iostream>

using namespace std;
using namespace PGL;

void Generate_CGAL_H(std::string input_path, std::string output_path)
{
	//parse geom.h
	VectorStr1 funct_values, funct_titles, funct_paras;
	std::map<int, VectorStr1> funct_notations;
	{
		std::ifstream file(input_path);
		for (std::string line; std::getline(file, line); )
		{
			if (Functs::StringContain(line, "CGAL")&& Functs::StringContain(line, "extern \"C\" PPGL_EXPORT"))
			{
				if (!Functs::StringContain(line, "(") || !Functs::StringContain(line, ")"))
				{
					Functs::MAssert("Does not include both ( and )");
				}

				auto value = line.substr(line.find("PPGL_EXPORT") + 11, line.find("CGAL") - line.find("PPGL_EXPORT") - 12);
				auto title = line.substr(line.find("CGAL"), line.find("(") - line.find("CGAL"));
				auto para = line.substr(line.find("("));
				funct_values.push_back(value);
				funct_titles.push_back(title);
				funct_paras.push_back(para);
				std::cerr << "Function: " << title  << std::endl;
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
	std::ofstream cgal_file(output_path);
	cgal_file << "#ifndef CGAL_ONCE" << std::endl;
	cgal_file << "#define CGAL_ONCE" << std::endl;
	cgal_file << "#pragma once" << std::endl;
	cgal_file << "#include <pgl_functs.hpp>" << std::endl;
	cgal_file << "using namespace std;" << std::endl;
	
	cgal_file << "namespace PPGL {" << std::endl;
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
	cgal_file << "class CGALPL" << std::endl;
	cgal_file << "{" << std::endl;
	cgal_file << "	public:" << std::endl;
	cgal_file << "	CGALPL()" << std::endl;
	cgal_file << "	{" << std::endl;
	cgal_file << "		hModule = Functs::LoadHMODULE(\"ppgl.dll\");" << std::endl;

	for (int i = 0; i < funct_titles.size(); i++)
	{
		std::string pre_str = "		";
		auto title = funct_titles[i];
		cgal_file << pre_str << title << "_C = ("<<title<<")GetProcAddress(hModule, \""<<title<<"\");" << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				cgal_file << pre_str << line << std::endl;
	}
	cgal_file << "	};" << std::endl << std::endl;

	cgal_file << "	static CGALPL& Inst()" << std::endl;
	cgal_file << "	{" << std::endl;
	cgal_file << "		static CGALPL instance;" << std::endl;
	cgal_file << "		return instance;" << std::endl;
	cgal_file << "	};" << std::endl << std::endl;


	cgal_file << "	HMODULE hModule;" << std::endl;

	for (int i = 0; i < funct_titles.size(); i++)
	{
		std::string pre_str = "	";
		auto title = funct_titles[i];
		cgal_file << pre_str << title << " " << title << "_C;" << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				cgal_file << pre_str << line << std::endl;
	}

	cgal_file << "};" << std::endl;

	cgal_file << "#define PL() CGALPL::Inst()" << std::endl;

	cgal_file << "}" << std::endl;

	cgal_file << "#endif" << std::endl;
	cgal_file.close();
};

int main(int argc, char* argv[])
{
	//get root path
	std::string root_path = Functs::WinGetCurDirectory().substr(0, Functs::WinGetCurDirectory().find_last_of("\\"));
	root_path = root_path.substr(0, root_path.find_last_of("\\")+1);

	//generate cgal.h
	Generate_CGAL_H(root_path + "ppgl\\geom.h", root_path + "port\\cgal.h");

	//copy dll - Release
	VectorStr1 pds = {"Debug", "Release", "RelWithDebInfo"};
	for (auto pd : pds)
	{
		if (Functs::DetectExisting(root_path + "build\\ppgl\\" + pd) && 
			Functs::DetectExisting(root_path + "port\\" + pd))
		{
			Functs::WinCopy(root_path + "build\\ppgl\\" + pd, root_path + "port\\"+pd);

			if (Functs::DetectExisting(root_path + "port\\" + pd + "\\ppgl.pdb"))
				Functs::WinDel(root_path + "port\\" + pd + "\\ppgl.pdb");

			if (Functs::DetectExisting(root_path + "port\\" + pd + "\\ppgl.exp"))
				Functs::WinDel(root_path + "port\\" + pd + "\\ppgl.exp");
		}
	}

	Functs::MAssert("Successfully finishing the post processing...", 1.5);
	return 0;
}

