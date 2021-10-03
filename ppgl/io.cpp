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

extern "C" PPGL_EXPORT void Test_PGL(Vector3d n) 
{
    Functs::MAssert("Test PGL ...");
}

extern "C" PPGL_EXPORT void CGAL_Vector_Base(Vector3d n, Vector3d &result) {
    Plane_3 plane(Point_3(0.0, 0.0, 0.0), Vector_3(n[0], n[1], n[2]));
    Vector_3 v = plane.base1();
    result = Vector3d(v[0], v[1], v[2]);
}

//3d plane relaed
Vector3d CGAL_3D_Plane_Base_1(Vector3d plane_p, Vector3d plane_n) {
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    Vector_3 v = plane.base1();
    return Vector3d(v[0], v[1], v[2]);
}

extern "C" PPGL_EXPORT void CGAL_Export_Path_Segment(std::ofstream &export_file_output, int &export_index,
                                                               std::string s_name, double r, double g, double b,
                                                               Vector3d &start, Vector3d &end, double radius) {
    Vector3d normal = end - start;
    Vector3d base_1;
    CGAL_Vector_Base(normal, base_1);
    double length_base_1 = glm::length(base_1);

    base_1[0] = base_1[0] / length_base_1 * radius;
    base_1[1] = base_1[1] / length_base_1 * radius;
    base_1[2] = base_1[2] / length_base_1 * radius;

    std::vector<Vector3d> vecs;
    for (int i = 0; i < 4; i++) {
        double angle = i * 2 * Math_PI / 4;
        Vector3d v = Functs::RotationAxis(normal + base_1, angle, normal);
        vecs.push_back(v + start);
    }
    for (int i = 0; i < 4; i++) {
        vecs.push_back(vecs[i] - normal);
    }

    std::vector<std::vector<int>> faces;

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

extern "C" PPGL_EXPORT void CGAL_Export_Path_Point(std::ofstream &export_file_output, int &export_index,
                                                             std::string s_name, double r, double g, double b,
                                                             Vector3d point, double radius) {
    std::vector<Vector3d> vecs;
    vecs.push_back(Vector3d(0.5, 0.5, 0.5));
    vecs.push_back(Vector3d(-0.5, 0.5, 0.5));
    vecs.push_back(Vector3d(-0.5, 0.5, -0.5));
    vecs.push_back(Vector3d(0.5, 0.5, -0.5));

    vecs.push_back(Vector3d(0.5, -0.5, 0.5));
    vecs.push_back(Vector3d(-0.5, -0.5, 0.5));
    vecs.push_back(Vector3d(-0.5, -0.5, -0.5));
    vecs.push_back(Vector3d(0.5, -0.5, -0.5));

    std::vector<std::vector<int>> faces;

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

        for (int j = faces[i].size() - 1; j >= 0; j--) {
            export_file_output << faces[i][j] + export_index << " ";
        }
        export_file_output << "" << std::endl;
    }
    export_index += 8;
}
