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
typedef double
(*CGAL_2D_Distance_Segment_Segment)(Vector2d s_0, Vector2d s_1, Vector2d e_0, Vector2d e_1);
typedef bool (*CGAL_2D_Location_Point_Polygon)(Vector2d p, Vector2d1 py);
typedef bool (*CGAL_2D_Location_Points_Polygon)(const Vector2d1& ps,
	const Vector2d1& py);
typedef double (*CGAL_2D_Distance_Point_Polygon)(Vector2d p, Vector2d1 py);

typedef bool (*CGAL_2D_Intersection_Segment_Segment
	)(Vector2d s_0_s, Vector2d s_0_e, Vector2d s_1_s, Vector2d s_1_e, Vector2d& inter);

typedef bool (*CGAL_2D_Intersection_Line_Line
	)(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);

typedef bool
(*CGAL_2D_Intersection_Segment_Polygon)(Vector2d s_s, Vector2d s_e, Vector2d1& p);
typedef bool (*CGAL_2D_Polygon_Is_Clockwise_Oriented)(Vector2d1& ps);
typedef double
(*CGAL_2D_Two_Polygons_Union)(Vector2d1 poly_0, Vector2d1 poly_1,
	std::vector<Vector2d1 >& inter_polygons);

typedef double (*CGAL_2D_Two_Polygons_Intersection)(const Vector2d1& poly_0,
	const Vector2d1& poly_1);

typedef std::vector<int>
(*CGAL_Decompose_Polyline)(Vector2d1& polyline, double threshold);
typedef bool
(*CGAL_Identify_Polycut_Extend)(const Vector2d1& polygon, const Vector2d& s,
	const Vector2d& e, Vector2d& ns, Vector2d& ne);
typedef bool
(*CGAL_Identify_Polycut_NotExtend)(const Vector2d1& polygon, const Vector2d& s,
	const Vector2d& e);
typedef bool (*CGAL_Identify_Polycut)(const Vector2d1& polygon,
	const Vector2d1& cutLine,
	std::vector<std::pair<bool, bool> >& result);

typedef void (*CGAL_2D_Polygon_One_Offsets)(Vector2d1& poly,
	double d,
	std::vector<Vector2d1 >& offset_polys);

typedef bool
(*CGAL_Construct_InOutSide_Polygon)(const Vector2d1& py, const Vector2d& p, const Vector2d& q, bool& isPInside,
	bool& isQInside);
typedef bool
(*CGAL_2D_Intersection_Ray_Segment)(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s,
	const Vector2d& s_1_e, Vector2d& inter);

typedef double (*GetAngleKerfOffsetTan)(const Vector2d& a, const Vector2d& b);


//implementation in "threeD.cpp"
//####################################################################################
typedef double (*CGAL_3D_Distance_Point_Segment)(Vector3d p, Vector3d s_s, Vector3d s_e);
typedef void
(*CGAL_3D_Plane_Fitting)(std::vector<Vector3d>& points, Vector3d& plane_p, Vector3d& plane_n);
typedef void (*CGAL_3D_Plane_Point_Projection)(Vector3d& plane_p, Vector3d& plane_n, Vector3d& p,
	Vector3d& result);
typedef void (*CGAL_3D_Plane_Points_Projection)(Vector3d& plane_p, Vector3d& plane_n,
	std::vector<Vector3d>& points,
	std::vector<Vector3d>& project_points);

typedef void
(*CGAL_3D_Plane_3D_to_2D_Point)(Vector3d& plane_p, Vector3d& plane_n, Vector3d& point_3d,
	Vector2d& result);
typedef void
(*CGAL_3D_Plane_2D_to_3D_Point)(Vector3d& plane_p, Vector3d& plane_n, Vector2d& points_2d,
	Vector3d& result);

typedef void (*CGAL_3D_Plane_3D_to_2D_Points)(Vector3d& plane_p, Vector3d& plane_n,
	std::vector<Vector3d>& points_3d,
	Vector2d1& points_2d);
typedef void (*CGAL_3D_Plane_2D_to_3D_Points)(Vector3d& plane_p, Vector3d& plane_n,
	Vector2d1& points_2d,
	std::vector<Vector3d>& points_3d);

typedef Vector3d(*CGAL_3D_Projection_Point_Segment)(Vector3d p, Vector3d s_s, Vector3d s_e);
typedef double
(*CGAL_3D_Distance_Point_Point)(double p_0_x, double p_0_y, double p_0_z, double p_1_x,
	double p_1_y, double p_1_z);
typedef double
(*CGAL_3D_Distance_Point_Polygon)(const std::vector<Vector3d>& py, const Vector3d& p);

typedef void (*CGAL_2D_Polygon_Triangulation)(const Vector2d2& polys, std::vector<std::vector<int>>& faces);

//implementation in "mesh.cpp"
//####################################################################################

typedef void (*CGAL_Remesh_Surface_by_Adding_Feature)(const std::vector<Vector3d>& feature,
	const std::vector<int>& face_ids,
	const std::vector<Vector3d>& vecs,
	const std::vector<int>& face_id_0,
	const std::vector<int>& face_id_1,
	const std::vector<int>& face_id_2,
	std::vector<int>& igl_cutting_0_edges,
	std::vector<int>& igl_cutting_1_edges,
	std::vector<Vector3d>& igl_cutting_points,
	std::vector<std::vector<int> >& cutting_faces);

typedef void (*CGAL_3D_Read_Triangle_Mesh)(std::string path, std::vector<Vector3d>& vecs,
	std::vector<int>& face_id_0,
	std::vector<int>& face_id_1,
	std::vector<int>& face_id_2);
typedef void (*CGAL_Mesh_Edges)(std::string path);


#endif
