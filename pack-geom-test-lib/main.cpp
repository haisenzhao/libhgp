#include "cgal.h"

using namespace std;

int main()
{
	HMODULE hModule = Functs::LoadHMODULE("carpentry_geom.dll");
	auto read_mesh = (CGAL_Vector_Base)GetProcAddress(hModule, "CGAL_Vector_Base");
	Vector3d result;
	read_mesh(Vector3d(1.0, 323.0, 0.0), result);
	std::cerr << result[0] << " " << result[1] << " " << result[2] << std::endl;
	system("pause");
	return 0;
}

