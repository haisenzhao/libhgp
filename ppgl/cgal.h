#ifndef CGAL_ONCE
#define CGAL_ONCE
#pragma once

#include <stdio.h>
#include <tchar.h>
#include <stdio.h>
#include <tchar.h>
#include "iostream"
#include <windows.h>
#include <vector>
#include <set>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <functional>
#include <queue>
#include <sstream>
#include <math.h>
#include <fstream>

#include <glm/glm.hpp>
#include <random>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/transform.hpp>

#include <pgl_functs.hpp>


using namespace std;
using namespace PGL;


typedef void (*Test_PGL)(Vector3d n);


//implementation in "io.cpp"
//####################################################################################
typedef void (*CGAL_Vector_Base)(Vector3d n, Vector3d&);
typedef void (*CGAL_Export_Path_Segment)(std::ofstream& export_file_output, int& export_index,
	std::string s_name, double r, double g, double b,
	Vector3d& start,
	Vector3d& end, double radius);
typedef void (*CGAL_Export_Path_Point)(std::ofstream& export_file_output, int& export_index,
	std::string s_name, double r, double g, double b,
	Vector3d point,
	double radius);

//implementation in "twoD.cpp"
//####################################################################################
typedef double (*CGAL_2D_Distance_Point_Point)(Vector2d p_0, Vector2d p_1);
typedef double (*CGAL_2D_Distance_Point_Line)(Vector2d v, Vector2d l_0, Vector2d l_1);
typedef double (*CGAL_2D_Distance_Point_Segment)(Vector2d v, Vector2d s_0, Vector2d s_1);
typedef double (*CGAL_2D_Distance_Segment_Segment)(Vector2d s_0, Vector2d s_1, Vector2d e_0, Vector2d e_1);
typedef bool (*CGAL_2D_Location_Point_Polygon)(Vector2d p, Vector2d1 py);
typedef bool (*CGAL_2D_Location_Points_Polygon)(const Vector2d1& ps, const Vector2d1& py);
typedef double (*CGAL_2D_Distance_Point_Polygon)(Vector2d p, Vector2d1 py);
typedef bool (*CGAL_2D_Intersection_Segment_Segment)(Vector2d s_0_s, Vector2d s_0_e, Vector2d s_1_s, Vector2d s_1_e, Vector2d& inter);
typedef bool (*CGAL_2D_Intersection_Line_Line)(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
typedef bool (*CGAL_2D_Intersection_Segment_Polygon)(Vector2d s_s, Vector2d s_e, Vector2d1& p);
typedef bool (*CGAL_2D_Polygon_Is_Clockwise_Oriented)(Vector2d1& ps);
typedef double (*CGAL_2D_Two_Polygons_Union)(Vector2d1 poly_0, Vector2d1 poly_1, Vector2d2& inter_polygons);
typedef double (*CGAL_2D_Two_Polygons_Intersection)(const Vector2d1& poly_0, const Vector2d1& poly_1);
typedef Vector1i1(*CGAL_Decompose_Polyline)(Vector2d1& polyline, double threshold);
typedef bool (*CGAL_Identify_Polycut_Extend)(const Vector2d1& polygon, const Vector2d& s, const Vector2d& e, Vector2d& ns, Vector2d& ne);
typedef bool (*CGAL_Identify_Polycut_NotExtend)(const Vector2d1& polygon, const Vector2d& s, const Vector2d& e);
typedef bool (*CGAL_Identify_Polycut)(const Vector2d1& polygon, const Vector2d1& cutLine, VectorPB1& result);
typedef void (*CGAL_2D_Polygon_One_Offsets)(Vector2d1& poly, double d, Vector2d2& offset_polys);
typedef bool (*CGAL_Construct_InOutSide_Polygon)(const Vector2d1& py, const Vector2d& p, const Vector2d& q, bool& isPInside, bool& isQInside);
typedef bool (*CGAL_2D_Intersection_Ray_Segment)(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
typedef double (*GetAngleKerfOffsetTan)(const Vector2d& a, const Vector2d& b);


//implementation in "threeD.cpp"
//####################################################################################
typedef double (*CGAL_3D_Distance_Point_Segment)(Vector3d p, Vector3d s_s, Vector3d s_e);
typedef void (*CGAL_3D_Plane_Fitting)(Vector3d1& points, Vector3d& plane_p, Vector3d& plane_n);
typedef void (*CGAL_3D_Plane_Point_Projection)(Vector3d& plane_p, Vector3d& plane_n, Vector3d& p, Vector3d& result);
typedef void (*CGAL_3D_Plane_Points_Projection)(Vector3d& plane_p, Vector3d& plane_n, Vector3d1& points, Vector3d1& project_points);
typedef void (*CGAL_3D_Plane_3D_to_2D_Point)(Vector3d& plane_p, Vector3d& plane_n, Vector3d& point_3d, Vector2d& result);
typedef void (*CGAL_3D_Plane_2D_to_3D_Point)(Vector3d& plane_p, Vector3d& plane_n, Vector2d& points_2d, Vector3d& result);
typedef void (*CGAL_3D_Plane_3D_to_2D_Points)(Vector3d& plane_p, Vector3d& plane_n, Vector3d1& points_3d, Vector2d1& points_2d);
typedef void (*CGAL_3D_Plane_2D_to_3D_Points)(Vector3d& plane_p, Vector3d& plane_n, Vector2d1& points_2d, Vector3d1& points_3d);
typedef Vector3d(*CGAL_3D_Projection_Point_Segment)(Vector3d p, Vector3d s_s, Vector3d s_e);
typedef double (*CGAL_3D_Distance_Point_Point)(const Vector3d& v0, const Vector3d& v1);
typedef double (*CGAL_3D_Distance_Point_Polygon)(const Vector3d1& py, const Vector3d& p);
typedef void (*CGAL_2D_Polygon_Triangulation)(const Vector2d2& polys, Vector1i2& faces);


//implementation in "mesh.cpp"
//####################################################################################

typedef void (*CGAL_Remesh_Surface_by_Adding_Feature)(const Vector3d1& feature,
	const Vector1i1& face_ids,
	const Vector3d1& vecs,
	const Vector1i1& face_id_0,
	const Vector1i1& face_id_1,
	const Vector1i1& face_id_2,
	Vector1i1& igl_cutting_0_edges,
	Vector1i1& igl_cutting_1_edges,
	Vector3d1& igl_cutting_points,
	Vector1i2& cutting_faces);

typedef void (*CGAL_3D_Read_Triangle_Mesh)(std::string path, Vector3d1& vecs,
	Vector1i1& face_id_0,
	Vector1i1& face_id_1,
	Vector1i1& face_id_2);
typedef void (*CGAL_Mesh_Edges)(std::string path);

typedef bool (*CGAL_3D_Intersection_Sphere_Ray)(double, double, double, double,
	double, double, double, double, double, double, std::vector<double>&, std::vector<double>&, std::vector<double>&);
typedef bool (*CGAL_3D_Intersection_Ray_Triangle)(Vector3d p, Vector3d n, Vector3d p0, Vector3d p1, Vector3d p2);
typedef bool (*CGAL_3D_Intersection_Ray_Mesh)(Vector3d p, Vector3d n, std::string path);
typedef void (*CGAL_3D_Intersection_Rays_Mesh)(Vector3d1 ps, Vector3d1 ns, std::string path, Vector3d1& inters);

class CGAL
{
public:

};

#endif
