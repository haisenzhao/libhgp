#include "cgal.h"
//#include "geom.h"

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
	system("copy E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\ppgl.dll E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");
	system("copy E:\\Task2\\personal-pack-geom-lib\\build\\ppgl\\Release\\gmp.dll E:\\Task2\\personal-pack-geom-lib\\ppgl\\dll\\");

	std::ifstream file("E:\\Task2\\personal-pack-geom-lib\\ppgl\\geom.h");

	for (std::string line; std::getline(file, line); )
	{
		std::cout << line << std::endl;
	}

	file.close();
};

int main()
{
	PostProcess();


	system("pause");
	return 0;
}

