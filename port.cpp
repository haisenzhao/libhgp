//#include "../port/cgal.h" //call functions from a dll
//#include "geom.h" //call functions directly
#include"port/cgal.h"

using namespace std;

using namespace PGL;
using namespace PPGL;


int main(int argc, char* argv[])
{
	auto cur_path = Functs::WinGetCurDirectory();

	//PGL functions
	auto d = Functs::DoubleString(123.456);

	// PPGL functions
	auto inter = PL().CGAL_3D_Projection_Point_Segment_C(Vector3d(0, 0, 0), Vector3d(17, 2, 2), Vector3d{ 32,23,234 });

	system("pause");
	return 0;
}

