#include "cgal.h" //call functions from a dll
//#include "geom.h" //call functions directly

using namespace std;
#ifdef CGAL_ONCE
using namespace PPGL;
#endif

int main(int argc, char* argv[])
{
#ifdef CGAL_ONCE
	//HMODULE hModule = Functs::LoadHMODULE("ppgl.dll");
	//auto read_mesh = (CGAL_Vector_Base)GetProcAddress(hModule, "CGAL_Vector_Base");
	//Vector3d result;
	//read_mesh(Vector3d(1.0, 323.0, 0.0), result);
	//std::cerr << result[0] << " " << result[1] << " " << result[2] << std::endl;

	//Vector3d p, Vector3d s_s, Vector3d s_e
	auto inter=PL().CGAL_3D_Projection_Point_Segment_C(Vector3d(0, 0, 0), Vector3d(17, 2, 2), Vector3d{32,23,234});

	Functs::CerrLine(Functs::VectorString(inter," "));

#endif

#ifdef geom_hpp
	Vector3d result;
	CGAL_Vector_Base(Vector3d(1.0, 323.0, 0.0), result);
	std::cerr << result[0] << " " << result[1] << " " << result[2] << std::endl;

	std::string root = "E:\\Dropbox\\Mold\\microstructures\\demo";
	std::string cube_file = "1_cube.obj";
	std::string micro_file = "1_micro.obj";
	Vector3d1 sampling_points;
	

	Vector3d1 origin_vecs;
	Vector1i1 origin_face_id_0, origin_face_id_1, origin_face_id_2;
	if (Functs::CerrLine("Read Input File.", 1))
	{
		std::cerr << "Path: " << root + "\\0_origin.obj" << std::endl;
		//CGAL_3D_Read_Triangle_Mesh(root + "\\0_origin.obj", origin_vecs, origin_face_id_0, origin_face_id_1, origin_face_id_2);
		FF::LoadObj3d(std::string(root + "\\0_origin.obj").c_str(), origin_vecs, origin_face_id_0, origin_face_id_1, origin_face_id_2);
	}

	std::cerr << origin_vecs.size() << std::endl;

#endif

	system("pause");
	return 0;
}

