#include "liblgp/liblgp.hpp"

#include <iostream>

using namespace std;
using namespace liblgp;

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
	cgal_file << "#include <iostream>" << std::endl;
	cgal_file << "#include <vector>" << std::endl;
	cgal_file << "#include <windows.h>" << std::endl;
	cgal_file << "#include <direct.h>" << std::endl;
	cgal_file << "#include <tchar.h>" << std::endl;
	cgal_file << "#include <string>" << std::endl;
	cgal_file << "#include \"glm/glm.hpp\"" << std::endl;

	cgal_file << "using namespace std;" << std::endl;
	
	cgal_file << "namespace PPGL {" << std::endl;

	cgal_file << "" << std::endl;
	cgal_file << "template <typename datum>" << std::endl;
	cgal_file << "using Vector1 = std::vector<datum>;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "template <typename datum>" << std::endl;
	cgal_file << "using Vector2 = std::vector<std::vector<datum>>;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "template <typename datum>" << std::endl;
	cgal_file << "using Vector3 = std::vector<std::vector<std::vector<datum>>>;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef glm::highp_dvec2 Vector2d;" << std::endl;
	cgal_file << "typedef glm::highp_dvec3 Vector3d;" << std::endl;
	cgal_file << "typedef glm::highp_ivec2 Vector2i;" << std::endl;
	cgal_file << "typedef glm::highp_ivec3 Vector3i;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<Vector2d> Vector2d1;" << std::endl;
	cgal_file << "typedef Vector2<Vector2d> Vector2d2;" << std::endl;
	cgal_file << "typedef Vector3<Vector2d> Vector2d3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<Vector3d> Vector3d1;" << std::endl;
	cgal_file << "typedef Vector2<Vector3d> Vector3d2;" << std::endl;
	cgal_file << "typedef Vector3<Vector3d> Vector3d3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<bool> Vector1b1;" << std::endl;
	cgal_file << "typedef Vector2<bool> Vector1b2;" << std::endl;
	cgal_file << "typedef Vector3<bool> Vector1b3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<int> Vector1i1;" << std::endl;
	cgal_file << "typedef Vector2<int> Vector1i2;" << std::endl;
	cgal_file << "typedef Vector3<int> Vector1i3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<double> Vector1d1;" << std::endl;
	cgal_file << "typedef Vector2<double> Vector1d2;" << std::endl;
	cgal_file << "typedef Vector3<double> Vector1d3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<std::string> VectorStr1;" << std::endl;
	cgal_file << "typedef Vector2<std::string> VectorStr2;" << std::endl;
	cgal_file << "typedef Vector3<std::string> VectorStr3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<Vector2i> Vector2i1;" << std::endl;
	cgal_file << "typedef Vector2<Vector2i> Vector2i2;" << std::endl;
	cgal_file << "typedef Vector3<Vector2i> Vector2i3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<Vector3i> Vector3i1;" << std::endl;
	cgal_file << "typedef Vector2<Vector3i> Vector3i2;" << std::endl;
	cgal_file << "typedef Vector3<Vector3i> Vector3i3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<std::pair<int, int>> VectorPI1;" << std::endl;
	cgal_file << "typedef Vector2<std::pair<int, int>> VectorPI2;" << std::endl;
	cgal_file << "typedef Vector3<std::pair<int, int>> VectorPI3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef Vector1<std::pair<bool, bool>> VectorPB1;" << std::endl;
	cgal_file << "typedef Vector2<std::pair<bool, bool>> VectorPB2;" << std::endl;
	cgal_file << "typedef Vector3<std::pair<bool, bool>> VectorPB3;" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "typedef std::tuple<int, int, int> TI3;" << std::endl;
	cgal_file << "typedef Vector1<std::tuple<int, int, int>> VectorTI3;" << std::endl;
	cgal_file << "" << std::endl;

	cgal_file << "static HMODULE LoadHMODULE(const string& dll_path)" << std::endl;
	cgal_file << "{" << std::endl;
	cgal_file << "	struct stat buffer;" << std::endl;
	cgal_file << "	if (!(stat(dll_path.c_str(), &buffer) == 0))" << std::endl;
	cgal_file << "	{" << std::endl;
	cgal_file << "		char tmp[256];" << std::endl;
	cgal_file << "		if (_getcwd(tmp, 256)) {};" << std::endl;
	cgal_file << "		std::string root_path = std::string(tmp);" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "		std::string str;" << std::endl;
	cgal_file << "		str += \"The dll does not exist: \" + dll_path + \"; \\n\";" << std::endl;
	cgal_file << "		str += \"The current running directory : \" + root_path + \"; \\n\";" << std::endl;
	cgal_file << "		str += \"Please gurrentee the dll is in the right place; \\n\";" << std::endl;
	cgal_file << "		std::cerr << str << std::endl;" << std::endl;
	cgal_file << "	}" << std::endl;
	cgal_file << "" << std::endl;
	cgal_file << "	HMODULE hModule = LoadLibraryA(LPCSTR(dll_path.c_str()));" << std::endl;
	cgal_file << "	if (!hModule)" << std::endl;
	cgal_file << "	{" << std::endl;
	cgal_file << "		DWORD dw = GetLastError(); // returns 0xc1 (193)" << std::endl;
	cgal_file << "		std::cerr << \"LoadLibrary failed with error code \" + std::to_string(dw) << std::endl;" << std::endl;
	cgal_file << "	}" << std::endl;
	cgal_file << "	else" << std::endl;
	cgal_file << "		std::cerr << \"LoadLibrary success\\n\";" << std::endl;
	cgal_file << "	return hModule;" << std::endl;
	cgal_file << "};" << std::endl;


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
	cgal_file << "		hModule = LoadHMODULE(\"ppgl.dll\");" << std::endl;

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
