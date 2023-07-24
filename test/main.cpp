#include "../port/cgal.h" //call functions from a dll
//#include "geom.h" //call functions directly

using namespace std;

#ifdef CGAL_ONCE
using namespace PPGL;
#endif



int Check2DIntersection()
{
	auto Load = [](Vector2d2& contours)
	{
		string address = "D:/test.txt";
		ifstream infile;
		infile.open(address);
		if (!infile.is_open()) {
			std::cout << "文件打开失败" << endl;
			return;
		}

		string line;

		while (getline(infile, line))
		{//每次从文件读取一行
			istringstream iss(line);
			Vector2d1 points;
			int n;
			iss >> n;
			double x, y;
			for (int i = 0; i < n; i++)
			{
				iss >> x >> y;
				points.push_back(Vector2d(x, y));
			}
			contours.push_back(points);
		}
	};

	Vector2d2 contours;
	Load(contours);


	for (int i = 0; i < contours.size(); i++)
	{
		if (PL().CGAL_2D_Polygon_Is_Clockwise_Oriented_C(contours[i]))
		{
			std::reverse(contours[i].begin(), contours[i].end());
		}
		//Functs::Export_Segment(contours[i],);

		Vector2d center = Functs::GetCenter(contours[i]);

		for (int j = 0; j < contours[i].size(); j++)
		{
			contours[i][j] = contours[i][j] - center;
		}

		auto c3d = Functs::Vector2d3d(contours[i], 0.0);
		Functs::OutputObj3d("D:\\debug_" + Functs::IntString(i) + ".obj", c3d);

		bool b = PL().CGAL_2D_Polygon_Simple_C(contours[i]);

		std::cerr << "Simple test: " << i << " " << b << std::endl;
	}

	Vector1d1 des;
	for (int i = 0; i < contours.size(); i++)
	{
		for (int j = i + 1; j < contours.size(); j++)
		{
			double ds = PL().CGAL_2D_Distance_Polygon_Polygon_C(contours[i], contours[j]);
			//double cs = CGAL_2D_Intersection_Polygon_Polygon(contours[i], contours[j]);

			double cs = PL().CGAL_2D_Two_Polygons_Intersection_C(contours[i], contours[j]);

			std::cerr << i << " " << j << " dis: " << ds << " collision: " << cs << std::endl;
			des.push_back(ds);
		}
	}

	return 0;
}

int main(int argc, char* argv[])
{
#ifdef CGAL_ONCE

	//HMODULE hModule = Functs::LoadHMODULE("ppgl.dll");
	//auto read_mesh = (CGAL_Vector_Base)GetProcAddress(hModule, "CGAL_Vector_Base");
	//Vector3d result;
	//read_mesh(Vector3d(1.0, 323.0, 0.0), result);
	//std::cerr << result[0] << " " << result[1] << " " << result[2] << std::endl;

	auto cur_path = Functs::WinGetCurDirectory();

	//Vector3d p, Vector3d s_s, Vector3d s_e
	auto inter=PL().CGAL_3D_Projection_Point_Segment_C(Vector3d(0, 0, 0), Vector3d(17, 2, 2), Vector3d{32,23,234});

	Check2DIntersection();

	Functs::CerrLine(Functs::VectorString(inter," "));

	CGALPL cgal;
	double abcdd = cgal.CGAL_2D_Distance_Point_Point_C(Vector2d(0, 0), Vector2d(1, 1));
	cout << abcdd << endl;
	const char* path = "D:\\cube.off";
	auto aaa = Functs::DetectExisting(path);
	// C:\Users\86139\Desktop\Graduation_Design\fileTransform\obj
	Vector3d plane_normal(1.0, 0.0, 0.0);
	std::vector<double> plane_d;
	plane_d.push_back(1.0);
	Vector3d3 offsetses;
	Vector3d2 offsets;
	cgal.CGAL_Slicer_Mesh_C(path, plane_normal, plane_d, offsetses, offsets);

#endif

#ifdef geom_hpp


	

	/*
	Vector3d plane_normal(1.0, 0.0, 0.0);
	std::vector<double> plane_d;
	plane_d.push_back(1.0);
	Vector3d3 offsetses;
	Vector3d2 offsets;
	//(const char* path, const Vector3d & plane_normal, const std::vector<double> &plane_d, Vector3d3 & offsetses, Vector3d2 & offsets)
	CGAL_Slicer_Mesh(path, plane_normal, plane_d, offsetses, offsets);
	*/
	
	/*
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
	*/
#endif

	system("pause");
	return 0;
}

