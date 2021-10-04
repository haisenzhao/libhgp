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

void PostProcess() 
{


	//read
	VectorStr1 funct_values;
	VectorStr1 funct_titles;
	VectorStr1 funct_paras;
	std::ifstream file("E:\\Task2\\personal-pack-geom-lib\\ppgl\\geom.h");

	for (std::string line; std::getline(file, line); )
	{
		if (Functs::StringContain(line, "CGAL"))
		{
			auto value = line.substr(line.find("PPGL_EXPORT")+11, line.find("CGAL") - line.find("PPGL_EXPORT")-12);
			auto title = line.substr(line.find("CGAL"), line.find("(") - line.find("CGAL"));
			auto para = line.substr(line.find("("));
			funct_values.push_back(value);
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
		auto value = funct_values[i];
		auto title = funct_titles[i];
		auto para = funct_paras[i];
		cgal_file << "typedef "+value+" (*" << title << ")" << para << std::endl;
	}

	//	auto ReadMesh = (CGAL_3D_Read_Triangle_Mesh)GetProcAddress(CompilerConfig::Instance().GetHModule(), "CGAL_3D_Read_Triangle_Mesh");

	cgal_file << std::endl;
	cgal_file << std::endl;
	cgal_file << "class PL" << std::endl;
	cgal_file << "{" << std::endl;
	cgal_file << "	public:" << std::endl;
	cgal_file << "	PL()" << std::endl;
	cgal_file << "	{" << std::endl;
	cgal_file << "		hModule = Functs::LoadHMODULE(\"ppgl.dll\");" << std::endl;

	//ReadMesh = (CGAL_2D_Distance_Point_Point)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Point");
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



	//static CompilerConfig& Instance();


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

void PreProcess()
{
	std::string folder("E:\\Dropbox\\Mold\\microstructures");
	auto objs = Functs::GetFilesInDirectory(folder);

	//std::string path, Vector3d1 &vecs,Vector1i1 &face_id_0, Vector1i1 &face_id_1, Vector1i1 &face_id_2

	for (int i = 0; i < objs.size(); i++)
	{
		//input
		Vector3d1 in_vecs;
		Vector1i1 in_face_id_0, in_face_id_1, in_face_id_2;

		auto in_path = folder + "\\" + objs[i];
		CGAL_3D_Read_Triangle_Mesh(in_path, in_vecs, in_face_id_0, in_face_id_1, in_face_id_2);

		Vector3d1 vecms;
		{
			Vector3d minc, maxc;
			Functs::GetBoundingBox(in_vecs, minc, maxc);
			auto bb_size = Functs::GetMax(maxc - minc);
			auto center = (minc + maxc) / 2.0;

			std::cerr << i << " : " << objs[i] << " bb_size: " << bb_size << " center: " << Functs::VectorString(center, " ") << " minc: " << Functs::VectorString(minc, " ") << std::endl;


			//transformation
			auto tm = Functs::TranslationMatrix(-center);
			auto sm = Functs::ScaleMatrix(Vector3d(100.0 / bb_size, 100.0 / bb_size, 100.0 / bb_size));
			vecms = Functs::PosApplyM(in_vecs, sm * tm);
		}

		{
			Vector3d minc, maxc;
			Functs::GetBoundingBox(vecms, minc, maxc);
			auto bb_size = Functs::GetMax(maxc - minc);
			auto center = (minc + maxc) / 2.0;

			std::cerr << i << " : " << objs[i] << " bb_size: " << bb_size << " center: " << Functs::VectorString(center, " ") << " minc: " << Functs::VectorString(minc, " ") << std::endl;

		}


		//output
		auto out_path = folder + "\\temp\\ts_" + objs[i];



		CGAL_3D_Output_Triangle_Mesh(out_path, vecms, in_face_id_0, in_face_id_1, in_face_id_2);



		//Functs::OutputObj3d();
	}
};
int main()
{

	system("copy E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");
	system("copy E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");

	system("copy E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll E:\\Task2\\personal-pack-geom-lib\\build\\test\\Release\\");
	system("copy E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll E:\\Task2\\personal-pack-geom-lib\\build\\test\\Release\\");

	PreProcess();
	//PostProcess();

	
	system("pause");
	return 0;
}

