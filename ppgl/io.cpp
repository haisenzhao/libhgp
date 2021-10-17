#include "geom.h"

Vector3d PointVector3d(Point_3 p) {
    return Vector3d(p[0], p[1], p[2]);
}

Point_3 VectorPoint3d(Vector3d p) {
    return Point_3(p[0], p[1], p[2]);
}

Vector2d PointVector2d(Point_2 p) {
    return Vector2d(p[0], p[1]);
}

Point_2 VectorPoint2d(Vector2d p) {
    return Point_2(p[0], p[1]);
}

extern "C" PPGL_EXPORT void CGAL_Test_PGL(const Vector3d& n) 
{
    Functs::MAssert("Test PGL ...");
}

extern "C" PPGL_EXPORT void CGAL_Vector_Base(const Vector3d& n, Vector3d &result) {
    Plane_3 plane(Point_3(0.0, 0.0, 0.0), Vector_3(n[0], n[1], n[2]));
    Vector_3 v = plane.base1();
    result = Vector3d(v[0], v[1], v[2]);
}

//3d plane related
Vector3d CGAL_3D_Plane_Base_1(Vector3d plane_p, Vector3d plane_n) {
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    Vector_3 v = plane.base1();
    return Vector3d(v[0], v[1], v[2]);
}

extern "C" PPGL_EXPORT void CGAL_Export_Path_Segment(std::ofstream &export_file_output, int &export_index,
                                                               const std::string s_name, const double r, const double g, const double b,
                                                               const Vector3d &start, const Vector3d &end, const double radius) {
    Vector3d normal = end - start;
    Vector3d base_1;
    CGAL_Vector_Base(normal, base_1);
    double length_base_1 = glm::length(base_1);

    base_1[0] = base_1[0] / length_base_1 * radius;
    base_1[1] = base_1[1] / length_base_1 * radius;
    base_1[2] = base_1[2] / length_base_1 * radius;

    Vector3d1 vecs;
    for (int i = 0; i < 4; i++) {
        double angle = (double)(i) * Math_PI / 2.0;
        Vector3d v = Functs::RotationAxis(normal + base_1, angle, normal);
        vecs.push_back(v + start);
    }
    for (int i = 0; i < 4; i++) {
        vecs.push_back(vecs[i] - normal);
    }

    Vector1i2 faces;

    int face_index_0[4] = {0, 1, 2, 3};
    int face_index_1[4] = {5, 1, 0, 4};
    int face_index_2[4] = {4, 0, 3, 7};
    int face_index_3[4] = {5, 4, 7, 6};
    int face_index_4[4] = {7, 3, 2, 6};
    int face_index_5[4] = {6, 2, 1, 5};

    faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
    faces.push_back(std::vector<int>(face_index_1, face_index_1 + 4));
    faces.push_back(std::vector<int>(face_index_2, face_index_2 + 4));
    faces.push_back(std::vector<int>(face_index_3, face_index_3 + 4));
    faces.push_back(std::vector<int>(face_index_4, face_index_4 + 4));
    faces.push_back(std::vector<int>(face_index_5, face_index_5 + 4));

    export_file_output << "g " + s_name << std::endl;

    for (int i = 0; i < vecs.size(); i++) {
        export_file_output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << r << " " << g
                           << " " << b << std::endl;
    }

    for (int i = 0; i < faces.size(); i++) {
        export_file_output << "f ";

        for (int j = 0; j < faces[i].size(); j++) {
            export_file_output << faces[i][j] + export_index << " ";
        }
        export_file_output << "" << std::endl;
    }

    export_index += 8;
}

extern "C" PPGL_EXPORT void CGAL_Export_Path_Point
(std::ofstream &export_file_output, int &export_index,const std::string s_name, const double r, const double g, const double b,const Vector3d point, const double radius)
{
    Vector3d1 vecs;
    vecs.push_back(Vector3d(0.5, 0.5, 0.5));
    vecs.push_back(Vector3d(-0.5, 0.5, 0.5));
    vecs.push_back(Vector3d(-0.5, 0.5, -0.5));
    vecs.push_back(Vector3d(0.5, 0.5, -0.5));

    vecs.push_back(Vector3d(0.5, -0.5, 0.5));
    vecs.push_back(Vector3d(-0.5, -0.5, 0.5));
    vecs.push_back(Vector3d(-0.5, -0.5, -0.5));
    vecs.push_back(Vector3d(0.5, -0.5, -0.5));

    Vector1i2 faces;

    int face_index_0[4] = {0, 1, 2, 3};
    int face_index_1[4] = {5, 1, 0, 4};
    int face_index_2[4] = {4, 0, 3, 7};
    int face_index_3[4] = {5, 4, 7, 6};
    int face_index_4[4] = {7, 3, 2, 6};
    int face_index_5[4] = {6, 2, 1, 5};

    faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
    faces.push_back(std::vector<int>(face_index_1, face_index_1 + 4));
    faces.push_back(std::vector<int>(face_index_2, face_index_2 + 4));
    faces.push_back(std::vector<int>(face_index_3, face_index_3 + 4));
    faces.push_back(std::vector<int>(face_index_4, face_index_4 + 4));
    faces.push_back(std::vector<int>(face_index_5, face_index_5 + 4));

    export_file_output << "g " + s_name << std::endl;

    for (int i = 0; i < vecs.size(); i++) {
        vecs[i][0] = vecs[i][0] * radius;
        vecs[i][1] = vecs[i][1] * radius;
        vecs[i][2] = vecs[i][2] * radius;

        vecs[i][0] += point[0];
        vecs[i][1] += point[1];
        vecs[i][2] += point[2];

        export_file_output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << r << " " << g
                           << " " << b << std::endl;
    }
    for (int i = 0; i < faces.size(); i++) {
        export_file_output << "f ";

        for (int j = (int)faces[i].size() - 1; j >= 0; j--) {
            export_file_output << faces[i][j] + export_index << " ";
        }
        export_file_output << "" << std::endl;
    }
    export_index += 8;
}





extern "C" PPGL_EXPORT void CGAL_Output_Obj_C1(const std::string & path, const Vector3d1 & vecs)
{
	if (vecs.size() < 3)
	{
		std::cout << "CGAL_Output_Obj error: vecs.size() < 3 " << std::endl;
		return;
	}

	std::ofstream file(path);
	for (int i = 0; i < vecs.size(); i++)
	{
		file << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
	}

	file.close();
}
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C2(const std::string & path, const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2)
{

	if (vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
	{
		std::cout << "CGAL_Output_Obj error: vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
		return;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
		{
			std::cout << "CGAL_Output_Obj error: index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
			return;
		}
	}

	std::ofstream file(path);
	for (int i = 0; i < vecs.size(); i++)
	{
		Vector3d v = vecs[i];
		//file << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
		file << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 != index_1 && index_0 != index_2 && index_1 != index_2)
			file << "f " << index_0 + 1 << " " << index_1 + 1 << " " << index_2 + 1 << std::endl;
	}
	file.close();
}
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C3(const std::string & path, const Vector3d1 & vecs, const std::vector<std::vector<int>>&face_ids)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}

	CGAL_Output_Obj_C2(path, vecs, face_id_0, face_id_1, face_id_2);
}
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C4(const std::string & path, const Vector3d1 & vecs, const std::vector<std::vector<int>>&face_ids, const std::vector<int>&triangles_lables, const int& index)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	std::vector<int> lables(vecs.size(), -1);
	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
		if (triangles_lables[i] == index)
		{
			lables[face_ids[i][0]] = 0;
			lables[face_ids[i][1]] = 0;
			lables[face_ids[i][2]] = 0;
		}
	}
	Vector3d1 new_vecs;
	std::vector<int> new_face_id_0;
	std::vector<int> new_face_id_1;
	std::vector<int> new_face_id_2;

	int vertices_nb = 0;
	for (int i = 0; i < vecs.size(); i++)
	{
		if (lables[i] == 0)
		{
			Vector3d v = vecs[i];
			new_vecs.push_back(v);
			lables[i] = vertices_nb;
			vertices_nb++;
		}
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		if (triangles_lables[i] == index)
		{
			new_face_id_0.push_back(lables[face_id_0[i]]);
			new_face_id_1.push_back(lables[face_id_1[i]]);
			new_face_id_2.push_back(lables[face_id_2[i]]);
		}
	}

	CGAL_Output_Obj_C2(path, new_vecs, new_face_id_0, new_face_id_1, new_face_id_2);
}
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C5(const std::string & path, const Vector3d1 & vecs, const Vector3d1 & colors, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2)
{

	if (vecs.size() < 3 || colors.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
	{
		std::cout << "CGAL_Output_Obj error: vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
		return;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
		{
			std::cout << "CGAL_Output_Obj error: index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
			return;
		}
	}

	std::ofstream file(path);
	for (int i = 0; i < vecs.size(); i++)
	{
		file << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << std::endl;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 != index_1 && index_0 != index_2 && index_1 != index_2)
			file << "f " << index_0 + 1 << " " << index_1 + 1 << " " << index_2 + 1 << std::endl;
	}
	file.close();
}

extern "C" PPGL_EXPORT void CGAL_Output_Obj_C6(const std::string & path, const Vector3d1 & vecs, const Vector3d1 & colors, const std::vector<std::vector<int>>&face_ids)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}
	CGAL_Output_Obj_C5(path, vecs, colors, face_id_0, face_id_1, face_id_2);
}


extern "C" PPGL_EXPORT void CGAL_Output_Off(const std::string & path, const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2)
{
	if (vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
	{
		std::cout << "CGAL_Output_Off error: vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
		return;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
		{
			std::cout << "CGAL_Output_Off error: index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
			return;
		}
	}
	std::ofstream file(path);
	file << "OFF" << std::endl;
	file << vecs.size() << " " << face_id_0.size() << " 0" << std::endl;
	for (int i = 0; i < vecs.size(); i++)
		file << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
	for (int i = 0; i < face_id_0.size(); i++)
		file << "3 " << face_id_0[i] << " " << face_id_1[i] << " " << face_id_2[i] << " " << std::endl;
	file.close();
}

extern "C" PPGL_EXPORT void CGAL_Load_Obj(const std::string & path, std::vector<double>&coords, std::vector<int>&tris)
{
	std::cerr << 4 << std::endl;

	auto get_first_integer = [](const char* v)
	{
		int ival;
		std::string s(v);
		std::replace(s.begin(), s.end(), '/', ' ');
		sscanf(s.c_str(), "%d", &ival);
		return ival;
	};

	double x, y, z;
	char line[1024], v0[1024], v1[1024], v2[1024];
	std::cerr << 5 << std::endl;

	std::cerr << path << std::endl;

	// open the file, return if open fails
	FILE* fp = fopen(path.c_str(), "r");
	if (!fp)
	{
		Functs::MAssert("This file does not exist: "+path);
		return;
	};

	std::cerr << 6 << std::endl;


	int i = 0;
	while (fgets(line, 1024, fp)) {


		std::cerr << line << std::endl;
		std::cerr << coords.size() << std::endl;
		if (line[0] == 'v') {
			sscanf(line, "%*s%lf%lf%lf", &x, &y, &z);
			coords.push_back(x);
			coords.push_back(y);
			coords.push_back(z);
		}
		else if (line[0] == 'f') {
			sscanf(line, "%*s%s%s%s", v0, v1, v2);
			tris.push_back(get_first_integer(v0) - 1);
			tris.push_back(get_first_integer(v1) - 1);
			tris.push_back(get_first_integer(v2) - 1);
		}
	}
	fclose(fp);
}
