#include "cgal.h"
//#include "geom.h"

using namespace std;

int main()
{
#ifdef CGAL_ONCE

	HMODULE hModule = Functs::LoadHMODULE("carpentry_geom.dll");
	auto read_mesh = (CGAL_Vector_Base)GetProcAddress(hModule, "CGAL_Vector_Base");
	Vector3d result;
	read_mesh(Vector3d(1.0, 323.0, 0.0), result);
	std::cerr << result[0] << " " << result[1] << " " << result[2] << std::endl;
#endif

#ifdef geom_hpp
	Vector3d result;
	CGAL_Vector_Base(Vector3d(1.0, 323.0, 0.0), result);
	std::cerr << result[0] << " " << result[1] << " " << result[2] << std::endl;

#endif

	system("pause");
	return 0;
}

