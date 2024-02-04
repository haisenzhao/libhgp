#include "liblgp/liblgp.hpp"

#include <iostream>

using namespace std;
using namespace liblgp;

void Generate_LIBHGP_H(const std::string & input_path, const std::string & output_path)
{
	if (!Functs::DetectExisting(input_path))
	{
		Functs::MAssert(Functs::DetectExisting(input_path), "File is not existing: "+ input_path);
		return;
	}

	//parse geom.h
	VectorStr1 funct_values, funct_titles, funct_paras;
	std::map<int, VectorStr1> funct_notations;
	{
		std::ifstream file(input_path);
		for (std::string line; std::getline(file, line); )
		{
			if (Functs::StringContain(line, "CGAL")&& Functs::StringContain(line, "extern \"C\" LIBHGP_EXPORT"))
			{
				if (!Functs::StringContain(line, "(") || !Functs::StringContain(line, ")"))
				{
					Functs::MAssert("Does not include both ( and )");
				}

				auto value = line.substr(line.find("LIBHGP_EXPORT") + 14, line.find("CGAL") - line.find("LIBHGP_EXPORT") - 15);
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

	//output libhgp.h
	std::ofstream libhgp_file(output_path);
	libhgp_file << "#ifndef LIBHGP_ONCE" << std::endl;
	libhgp_file << "#define LIBHGP_ONCE" << std::endl;
	libhgp_file << "#pragma once" << std::endl;
	libhgp_file << "#include <iostream>" << std::endl;
	libhgp_file << "#include <vector>" << std::endl;
	libhgp_file << "#include <windows.h>" << std::endl;
	libhgp_file << "#include <direct.h>" << std::endl;
	libhgp_file << "#include <tchar.h>" << std::endl;
	libhgp_file << "#include <string>" << std::endl;
	libhgp_file << "#include \"glm/glm.hpp\"" << std::endl;

	libhgp_file << "using namespace std;" << std::endl;
	
	libhgp_file << "namespace libhgp {" << std::endl;

	libhgp_file << "" << std::endl;
	libhgp_file << "template <typename datum>" << std::endl;
	libhgp_file << "using Vector1 = std::vector<datum>;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "template <typename datum>" << std::endl;
	libhgp_file << "using Vector2 = std::vector<std::vector<datum>>;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "template <typename datum>" << std::endl;
	libhgp_file << "using Vector3 = std::vector<std::vector<std::vector<datum>>>;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef glm::highp_dvec2 Vector2d;" << std::endl;
	libhgp_file << "typedef glm::highp_dvec3 Vector3d;" << std::endl;
	libhgp_file << "typedef glm::highp_ivec2 Vector2i;" << std::endl;
	libhgp_file << "typedef glm::highp_ivec3 Vector3i;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<Vector2d> Vector2d1;" << std::endl;
	libhgp_file << "typedef Vector2<Vector2d> Vector2d2;" << std::endl;
	libhgp_file << "typedef Vector3<Vector2d> Vector2d3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<Vector3d> Vector3d1;" << std::endl;
	libhgp_file << "typedef Vector2<Vector3d> Vector3d2;" << std::endl;
	libhgp_file << "typedef Vector3<Vector3d> Vector3d3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<bool> Vector1b1;" << std::endl;
	libhgp_file << "typedef Vector2<bool> Vector1b2;" << std::endl;
	libhgp_file << "typedef Vector3<bool> Vector1b3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<int> Vector1i1;" << std::endl;
	libhgp_file << "typedef Vector2<int> Vector1i2;" << std::endl;
	libhgp_file << "typedef Vector3<int> Vector1i3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<double> Vector1d1;" << std::endl;
	libhgp_file << "typedef Vector2<double> Vector1d2;" << std::endl;
	libhgp_file << "typedef Vector3<double> Vector1d3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<std::string> VectorStr1;" << std::endl;
	libhgp_file << "typedef Vector2<std::string> VectorStr2;" << std::endl;
	libhgp_file << "typedef Vector3<std::string> VectorStr3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<Vector2i> Vector2i1;" << std::endl;
	libhgp_file << "typedef Vector2<Vector2i> Vector2i2;" << std::endl;
	libhgp_file << "typedef Vector3<Vector2i> Vector2i3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<Vector3i> Vector3i1;" << std::endl;
	libhgp_file << "typedef Vector2<Vector3i> Vector3i2;" << std::endl;
	libhgp_file << "typedef Vector3<Vector3i> Vector3i3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<std::pair<int, int>> VectorPI1;" << std::endl;
	libhgp_file << "typedef Vector2<std::pair<int, int>> VectorPI2;" << std::endl;
	libhgp_file << "typedef Vector3<std::pair<int, int>> VectorPI3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef Vector1<std::pair<bool, bool>> VectorPB1;" << std::endl;
	libhgp_file << "typedef Vector2<std::pair<bool, bool>> VectorPB2;" << std::endl;
	libhgp_file << "typedef Vector3<std::pair<bool, bool>> VectorPB3;" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "typedef std::tuple<int, int, int> TI3;" << std::endl;
	libhgp_file << "typedef Vector1<std::tuple<int, int, int>> VectorTI3;" << std::endl;
	libhgp_file << "" << std::endl;

	libhgp_file << "static HMODULE LoadHMODULE(const string& dll_path)" << std::endl;
	libhgp_file << "{" << std::endl;
	libhgp_file << "	struct stat buffer;" << std::endl;
	libhgp_file << "	if (!(stat(dll_path.c_str(), &buffer) == 0))" << std::endl;
	libhgp_file << "	{" << std::endl;
	libhgp_file << "		char tmp[256];" << std::endl;
	libhgp_file << "		if (_getcwd(tmp, 256)) {};" << std::endl;
	libhgp_file << "		std::string root_path = std::string(tmp);" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "		std::string str;" << std::endl;
	libhgp_file << "		str += \"The dll does not exist: \" + dll_path + \"; \\n\";" << std::endl;
	libhgp_file << "		str += \"The current running directory : \" + root_path + \"; \\n\";" << std::endl;
	libhgp_file << "		str += \"Please gurrentee the dll is in the right place; \\n\";" << std::endl;
	libhgp_file << "		std::cerr << str << std::endl;" << std::endl;
	libhgp_file << "	}" << std::endl;
	libhgp_file << "" << std::endl;
	libhgp_file << "	HMODULE hModule = LoadLibraryA(LPCSTR(dll_path.c_str()));" << std::endl;
	libhgp_file << "	if (!hModule)" << std::endl;
	libhgp_file << "	{" << std::endl;
	libhgp_file << "		DWORD dw = GetLastError(); // returns 0xc1 (193)" << std::endl;
	libhgp_file << "		std::cerr << \"LoadLibrary failed with error code \" + std::to_string(dw) << std::endl;" << std::endl;
	libhgp_file << "	}" << std::endl;
	libhgp_file << "	else" << std::endl;
	libhgp_file << "		std::cerr << \"LoadLibrary success\\n\";" << std::endl;
	libhgp_file << "	return hModule;" << std::endl;
	libhgp_file << "};" << std::endl;


	//define functions
	if (funct_notations.find(-1) != funct_notations.end())
		for (auto& line : funct_notations[-1])
			libhgp_file << line << std::endl;
	for (int i = 0; i < funct_titles.size(); i++)
	{
		auto value = funct_values[i];
		auto title = funct_titles[i];
		auto para = funct_paras[i];
		libhgp_file << "typedef "+value+" (*" << title << ")" << para << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				libhgp_file << line << std::endl;
	}

	//define class
	libhgp_file << std::endl;
	libhgp_file << std::endl;
	libhgp_file << "class CGALPL" << std::endl;
	libhgp_file << "{" << std::endl;
	libhgp_file << "	public:" << std::endl;
	libhgp_file << "	CGALPL()" << std::endl;
	libhgp_file << "	{" << std::endl;
	libhgp_file << "		hModule = LoadHMODULE(\"libhgp.dll\");" << std::endl;

	for (int i = 0; i < funct_titles.size(); i++)
	{
		std::string pre_str = "		";
		auto title = funct_titles[i];
		libhgp_file << pre_str << title << "_C = ("<<title<<")GetProcAddress(hModule, \""<<title<<"\");" << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				libhgp_file << pre_str << line << std::endl;
	}
	libhgp_file << "	};" << std::endl << std::endl;

	libhgp_file << "	static CGALPL& Inst()" << std::endl;
	libhgp_file << "	{" << std::endl;
	libhgp_file << "		static CGALPL instance;" << std::endl;
	libhgp_file << "		return instance;" << std::endl;
	libhgp_file << "	};" << std::endl << std::endl;


	libhgp_file << "	HMODULE hModule;" << std::endl;

	for (int i = 0; i < funct_titles.size(); i++)
	{
		std::string pre_str = "	";
		auto title = funct_titles[i];
		libhgp_file << pre_str << title << " " << title << "_C;" << std::endl;
		if (funct_notations.find(i) != funct_notations.end())
			for (auto& line : funct_notations[i])
				libhgp_file << pre_str << line << std::endl;
	}

	libhgp_file << "};" << std::endl;

	libhgp_file << "#define PL() CGALPL::Inst()" << std::endl;

	libhgp_file << "}" << std::endl;

	libhgp_file << "#endif" << std::endl;
	libhgp_file.close();
};

int main(int argc, char* argv[])
{
	//get root path
	std::string build_path = Functs::WinGetCurDirectory().substr(0, Functs::WinGetCurDirectory().find_last_of("\\"));
	std::string root_path = build_path.substr(0, build_path.find_last_of("\\")+1);

	//generate libhgp.h
	Generate_LIBHGP_H(root_path + "dev\\geom.h", root_path + "libhgp\\libhgp.h");

	//copy dll - Release
	VectorStr1 pds = {"Debug", "Release", "RelWithDebInfo"};
	for (auto pd : pds)
	{
		std::string dev_path = build_path + "\\dev\\" + pd;
		std::string lib_path = root_path + "libhgp\\" + pd;

		auto bool_dev = Functs::DetectExisting(dev_path);
		auto bool_lib = Functs::DetectExisting(lib_path);

		if (bool_dev && bool_lib)
		{
			Functs::ClearFolder(lib_path);
			Functs::WinCopy(dev_path, lib_path);
			if (Functs::DetectExisting(lib_path + "\\libhgp.pdb")) Functs::WinDel(lib_path + "\\libhgp.pdb");
			if (Functs::DetectExisting(lib_path + "\\libhgp.exp")) Functs::WinDel(lib_path + "\\libhgp.exp");
		}
	}

	Functs::MAssert("Successfully finishing the post processing...", 1.5);
	system("pause");
	return 0;
}
