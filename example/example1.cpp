#include"libhgp.h"

using namespace libhgp;

int main(int argc, char* argv[])
{
	Vector3d n;
	PL().HGP_Vector_Base_C(Vector3d(1, 0, 0), n);
	std::cerr << "Vector base of (1,0,0): " << n[0] << " " << n[1] << " " << n[2] << std::endl;

	// PPGL functions
	auto b = PL().HGP_2D_Distance_Point_Point_C(Vector2d(0, 0), Vector2d(1, 1));
	std::cerr << "Distance from point (0,0) to point (1,1): " << b << std::endl;

	// PPGL functions
	auto p = PL().HGP_3D_Projection_Point_Segment_C(Vector3d(0, 0, 0), Vector3d(17, 2, 2), Vector3d{ 32,23,234 });
	std::cerr << "Project point (0,0,0) to segment (17,2,2)-(32,23,234): " << p[0] << " " << p[1] << " " << p[2] << std::endl;

	system("pause");
	return 0;
}