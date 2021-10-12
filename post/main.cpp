#include "pgl_functs.hpp"

using namespace std;
using namespace PGL;

void Check_Not_Included() 
{

	auto Get_CGAL_Functions = [](const std::string& path, 
		VectorStr1& funct_titles, VectorStr1& funct_lines)
	{
		std::ifstream file(path);
		for (std::string line; std::getline(file, line); )
		{
			if (Functs::StringContain(line, "CGAL"))
			{
				auto title = line.substr(line.find("CGAL"), line.find("(") - line.find("CGAL"));
				funct_titles.push_back(title);
				funct_lines.push_back(line);
				//std::cerr << title << std::endl;
			}
		}
		file.close();

	};

	//parse geom.h
	VectorStr1 funct_titles;
	VectorStr1 funct_lines;

	Get_CGAL_Functions("E:\\Task2\\personal-pack-geom-lib\\ppgl\\geom.h", funct_titles, funct_lines);


	VectorStr1 sfc_titles;
	VectorStr1 sfc_lines;
	Get_CGAL_Functions("E:\\Task2\\SmartCFS\\cgal\\cgalpackage.h",sfc_titles, sfc_lines);

	int debug_int = 0;

	ofstream file("E:\\out.txt");
	int i = 0;
	for (auto& sfc_title : sfc_titles)
	{
		bool b = true;
		for (auto& title : funct_titles)
		{
			if (Functs::StringContain(title, sfc_title))
			{
				b = false;
				break;
			}
		}

		std::cerr << i << std::endl;

		if (b)
		{
			//void CGAL_2D_Convex_Hulls(Vector2d1 & vec, Vector2d1 & hull_points);
			auto line = sfc_lines[i];

			if (!Functs::StringContain(line, "(")|| !Functs::StringContain(line, ")"))
			{
				std::cerr << line << std::endl;
				Functs::MAssert(line);
			}

			if (Functs::SplitStr(line, "(").size() != 2 || Functs::SplitStr(line, ")").size() != 2 || Functs::StringContain(line, "const"))
			{
				file << "extern \"C\" PPGL_EXPORT " << line << std::endl;
			}
			else
			{
				auto str_0 = Functs::SplitStr(line, "(")[0];
				auto str_1 = Functs::SplitStr(line, "(")[1];
				auto str_2 = Functs::SplitStr(str_1, ")")[0];
				auto str_3 = Functs::SplitStr(str_1, ")")[1];

				//std::string path, Vector3d1 &ves, Vector3d1 &ners

				str_2 = Functs::StringReplace(str_2, "&", "");
				str_2 = Functs::StringReplace(str_2, ",   ", ",");
				str_2 = Functs::StringReplace(str_2, ",  ", ",");
				str_2 = Functs::StringReplace(str_2, ", ", ",");
				
				auto paras = Functs::SplitStr(str_2,",");
				str_2 = "";
				int ii = 0;
				for (auto para : paras)
				{
					auto defs=Functs::SplitStr(para, " ");

					if (defs.size() == 2)
					{
						if (ii == paras.size() - 1)
							str_2 += "const " + defs[0] + "& " + defs[1];
						else
							str_2 += "const " + defs[0] + "& " + defs[1] + ", ";
					}
					else
					{
						if (defs.size() == 1)
						{
							if (ii == paras.size() - 1)
								str_2 += "const " + defs[0] + "& ";
							else
								str_2 += "const " + defs[0] + "& " + ", ";
						}
						else
						{
							if (ii == paras.size() - 1)
								str_2 += "const " + para;
							else
								str_2 += "const " + para + ", ";
						}
					}

					ii++;
				}

				//str_2 = Functs::StringReplace(str_2, ",", ", const");

				auto str = str_0 + "(" + str_2 + ")" + str_3;
				//const
				//&
				//line = Functs::StringReplace(line, "(", "(const ");
				//line = Functs::StringReplace(line, ",", ", const");
				file << "extern \"C\" PPGL_EXPORT " << str << std::endl;
			}
		}
		i++;
	}
	
	//Functs::CerrLine();
	file.close();
	system("pause");

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
		std::string pre_str = "		";
		auto title = funct_titles[i];
		cgal_file << pre_str << title << "_C = ("<<title<<")GetProcAddress(hModule, \""<<title<<"\");" << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				cgal_file << pre_str << line << std::endl;
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
		std::string pre_str = "	";
		auto title = funct_titles[i];
		cgal_file << pre_str << title << " " << title << "_C;" << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				cgal_file << pre_str << line << std::endl;
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

	//RelWithDebInfo
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\RelWithDebInfo\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\RelWithDebInfo\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\RelWithDebInfo\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\RelWithDebInfo\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\RelWithDebInfo\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\RelWithDebInfo\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\test\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\RelWithDebInfo\\ppgl.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\post\\RelWithDebInfo\\");
	Functs::WinCopy("E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\RelWithDebInfo\\gmp.dll", "E:\\Task2\\personal-pack-geom-lib\\build\\post\\RelWithDebInfo\\");

	if(argc==1)
		Generate_CGAL_H();
	else
		Check_Not_Included();

	//Functs::MAssert("Successfully finishing the post processing...");
	return 0;
}

