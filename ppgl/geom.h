#ifndef geom_hpp
#define geom_hpp
#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <algorithm>
#include <list>
#include <ppgl_export.h>

#include "pgl_functs.hpp"
#include "include_cgal.h"

using namespace PGL;

Vector3d PointVector3d(Point_3 p);
Point_3 VectorPoint3d(Vector3d p);
Vector2d PointVector2d(Point_2 p);
Point_2 VectorPoint2d(Vector2d p);

void  Construct_Polyhedron(Polyhedron_3& polyhedron, const Vector3d1& vecs, const Vector1i1& face_id_0, const Vector1i1& face_id_1, const Vector1i1& face_id_2);
void  Construct_Polyhedron(Polyhedron_3& polyhedron, const std::string& path);
void  Construct_Polyhedron(Polyhedron_3& polyhedron, const std::string& path, Vector3d1& vecs, Vector1i1& face_id_0, Vector1i1& face_id_1, Vector1i1& face_id_2);

extern "C" PPGL_EXPORT void CGAL_Test_PGL(const Vector3d& n);

//implementation in "io.cpp"
//####################################################################################
extern "C" PPGL_EXPORT void CGAL_Vector_Base(const Vector3d& n, Vector3d &);
extern "C" PPGL_EXPORT void CGAL_Export_Path_Segment(std::ofstream & export_file_output, int& export_index,const std::string s_name, const double r, const double g, const double b,const Vector3d & start, const Vector3d & end, const double radius);
extern "C" PPGL_EXPORT void CGAL_Export_Path_Point(std::ofstream & export_file_output, int& export_index, const std::string s_name, const double r, const double g, const double b, const Vector3d point, const double radius);


extern "C" PPGL_EXPORT void CGAL_Output_Obj_C1(const std::string & path, const Vector3d1 & vecs);
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C2(const std::string & path, const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2);
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C3(const std::string & path, const Vector3d1 & vecs, const std::vector<std::vector<int>>&face_ids);
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C4(const std::string & path, const Vector3d1 & vecs, const std::vector<std::vector<int>>&face_ids, const std::vector<int>&triangles_lables, const int& index);
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C5(const std::string & path, const Vector3d1 & vecs, const Vector3d1 & colors, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2);
extern "C" PPGL_EXPORT void CGAL_Output_Obj_C6(const std::string & path, const Vector3d1 & vecs, const Vector3d1 & colors, const std::vector<std::vector<int>>&face_ids);

//implementation in "twoD.cpp"
//####################################################################################
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Point(const Vector2d & p_0, const Vector2d& p_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Line(const Vector2d & v, const Vector2d & l_0, const Vector2d & l_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Segment(const Vector2d & v, const Vector2d & s_0, const Vector2d & s_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Segment_Segment(const Vector2d & s_0, const Vector2d & s_1, const Vector2d & e_0, const Vector2d & e_1);
extern "C" PPGL_EXPORT bool CGAL_2D_Location_Point_Polygon(const Vector2d & p, const Vector2d1 & py);
extern "C" PPGL_EXPORT bool CGAL_2D_Location_Points_Polygon(const Vector2d1 &ps, const Vector2d1 &py);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_2D_Polygon_Dart_Sampling(const Vector2d1& py, const double& d, Vector2d1& sampling_points, const int& total_iter);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Polygon(const Vector2d & p, const Vector2d1 & py);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Polygons(const Vector2d & p, const Vector2d2 & pys);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Segment(const Vector2d & s_0_s, const Vector2d & s_0_e, const Vector2d & s_1_s, const Vector2d & s_1_e, Vector2d &inter);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Line_Line(const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s, const Vector2d &s_1_e, Vector2d &inter);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Line(const Vector2d& s_s, const Vector2d & s_e, const Vector2d & l_s, const Vector2d & l_e, Vector2d& inter);

extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Polygon(const Vector2d & s_s, const Vector2d & s_e, Vector2d1 &p);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Is_Clockwise_Oriented(const Vector2d1 &ps);
extern "C" PPGL_EXPORT double CGAL_2D_Two_Polygons_Union(const Vector2d1 & poly_0, const Vector2d1 & poly_1, Vector2d2 & inter_polygons);
extern "C" PPGL_EXPORT double CGAL_2D_Two_Polygons_Intersection(const Vector2d1 &poly_0, const Vector2d1 &poly_1);
extern "C" PPGL_EXPORT void CGAL_Decompose_Polyline(const Vector2d1 & polyline, const double& threshold, Vector1i1 & result);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_Extend(const Vector2d1 &polygon, const Vector2d &s,const Vector2d &e, Vector2d &ns, Vector2d &ne);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_NotExtend(const Vector2d1 &polygon, const Vector2d &s,const Vector2d &e);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut(const Vector2d1 &polygon, const Vector2d1 &cutLine, VectorPB1 &result);

extern "C" PPGL_EXPORT bool CGAL_Construct_InOutSide_Polygon(const Vector2d1 &py, const Vector2d &p, const Vector2d &q, bool &isPInside,bool &isQInside);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Ray_Segment(const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s,const Vector2d &s_1_e, Vector2d &inter);
extern "C" PPGL_EXPORT double CGAL_Get_Angle_Kerf_Offset_Tan(const Vector2d &a, const Vector2d &b);
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Projection_Point_Segment(const Vector2d& p, const Vector2d& s, const Vector2d& e);

extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C1(const Vector2d1& outside_py,  const Vector2d& p);
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C2(const Vector2d1& outside_py,  const Vector2d1 & inside_py);
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C3(const Vector2d2& outside_pys, const Vector2d& p);
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C4(const Vector2d2& outside_pys, const Vector2d1& inside_py);
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C5(const Vector2d2& outside_pys, const Vector2d2& inside_pys);

extern "C" PPGL_EXPORT double CGAL_2D_Distance_Polygon_Polygon(const Vector2d1& poly_0, const Vector2d1& poly_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Polygons_Polygons(const Vector2d2& poly_0, const Vector2d2& poly_1);

extern "C" PPGL_EXPORT Vector2d CGAL_2D_Nearest_Point_Polygon_C1(const Vector2d& v, const Vector2d1& poly);
extern "C" PPGL_EXPORT void CGAL_2D_Nearest_Point_Polygon_C2(const Vector2d & v, const Vector2d1 & poly, Vector2d& p, double& min_d);
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Nearest_Point_Polygons(const Vector2d & v, const Vector2d2& polys);

extern "C" PPGL_EXPORT void CGAL_2d_Polygon_Boundingbox(const Vector2d1& ps, Vector2d& min_corner, Vector2d& max_corner);
extern "C" PPGL_EXPORT double CGAL_2D_Polygon_Area(const Vector2d1& py);
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Polygon_Inside_Point_C1(const Vector2d1& poly);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Inside_Point_C2(const Vector2d2& polys, Vector2d& inner_vec);

extern "C" PPGL_EXPORT void CGAL_2D_Polygon_One_Offsets(const Vector2d1 & poly, const double& d, Vector2d2 & offset_polys);
extern "C" PPGL_EXPORT void CGAL_2D_Polygons_One_Offsets(const Vector2d2 & polys, const double& d, Vector2d2 & offset_polys);

extern "C" PPGL_EXPORT bool CGAL_2D_Polygons_Simple(const Vector2d2& poly);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Simple(const Vector2d1 & poly);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Simple_Inter(const Vector2d1& poly);

extern "C" PPGL_EXPORT void CGAL_2D_Convex_Hulls(const Vector2d1 & vec, Vector2d1 & hull_points);

extern "C" PPGL_EXPORT void CGAL_2D_OBB_Box(const Vector2d1 & vec, Vector2d & center, Vector2d & axis_0, Vector2d & axis_1, double& entent_0, double& entent_1);



//implementation in "threeD.cpp"
//####################################################################################
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Segment(const Vector3d & p, const Vector3d & s_s, const Vector3d & s_e);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Fitting(const Vector3d1 & points, Vector3d &plane_p, Vector3d &plane_n);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Point_Projection(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d & p, Vector3d & result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Points_Projection(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d1 & points, Vector3d1 & project_points);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_3D_to_2D_Point(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d & point_3d, Vector2d & result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_2D_to_3D_Point(const Vector3d & plane_p, const Vector3d & plane_n, const Vector2d & points_2d, Vector3d & result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_3D_to_2D_Points(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d1 & points_3d, Vector2d1 & points_2d);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_2D_to_3D_Points(const Vector3d & plane_p, const Vector3d & plane_n, const Vector2d1 & points_2d, Vector3d1 & points_3d);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Segment(const Vector3d & p, const Vector3d & s_s, const Vector3d & s_e);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Point(const Vector3d & v0, const Vector3d & v1);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Polygon(const Vector3d1 &py, const Vector3d &p);
extern "C" PPGL_EXPORT void CGAL_2D_Polygon_Triangulation(const Vector2d2 &polys, Vector1i2 &faces);

extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Line(const Vector3d & p, const Vector3d & l_s, const Vector3d & l_e);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Line(const Vector3d & p, const Vector3d & l_s, const Vector3d & l_e);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Segment_Segment(const Vector3d & s_0_s, const Vector3d & s_0_e, const Vector3d & s_1_s, const Vector3d & s_1_e);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Plane(const Vector3d & v, const Vector3d & plane_p, const Vector3d & plane_n);

extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Line(const Vector3d & s_s, const Vector3d & s_e, const Vector3d & l_s, const Vector3d & l_e, Vector3d & inter);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Segment(const Vector3d & s_0_s, const Vector3d & s_0_e, const Vector3d & s_1_s, const Vector3d & s_1_e, Vector3d & iter);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Plane(const Vector3d & s_s, const Vector3d & s_e, const Vector3d & plane_p, const Vector3d & plane_n, Vector3d & inter);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Line_Plane(const Vector3d & l_s, const Vector3d & l_e, const Vector3d & plane_p, const Vector3d & plane_n, Vector3d & inter);

//implementation in "mesh.cpp"
//####################################################################################
extern "C" PPGL_EXPORT void CGAL_Remesh_Surface_by_Adding_Feature(const Vector3d1 &feature,const Vector1i1 &face_ids, const Vector3d1 &vecs, const Vector1i1 &face_id_0,const Vector1i1 &face_id_1,const Vector1i1 &face_id_2, Vector1i1 &igl_cutting_0_edges,Vector1i1 &igl_cutting_1_edges, Vector3d1 &igl_cutting_points,Vector1i2 &cutting_faces);
extern "C" PPGL_EXPORT void CGAL_3D_Output_Triangle_Mesh(const std::string & path, const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2);
extern "C" PPGL_EXPORT void CGAL_3D_Read_Triangle_Mesh(const std::string& path, Vector3d1 &vecs,Vector1i1 &face_id_0, Vector1i1 &face_id_1, Vector1i1 &face_id_2);
extern "C" PPGL_EXPORT void CGAL_Mesh_Edges(const std::string & path);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Sphere_Ray(const double& center_x, const double& center_y, const double& center_z, const double& radius,const double& ray_origin_x, const double& ray_origin_y, const double& ray_origin_z, const double& ray_direction_x, const double& ray_direction_y, const double& ray_direction_z,std::vector<double>&i_x, std::vector<double>&i_y, std::vector<double>&i_z);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Triangle(const Vector3d & p, const Vector3d & n, const Vector3d & p0, const Vector3d & p1, const Vector3d & p2);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Mesh(const Vector3d & p, const Vector3d & n, const std::string & path);
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_Vector3d(const Vector3d1 & ps, const Vector3d1 & ns, const std::string & path, Vector3d1 & inters);
//test each group directions (nes[i]) for each point in ps
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_C1_Bool(const Vector3d1 & ps, const Vector3d2 & nes, const std::string & path, Vector1b2 & inters);
//test all directions (ns) for each point in ps
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_C2_Bool(const Vector3d1 & ps, const Vector3d1 & ns, const std::string & path, Vector1b2 & inters);
extern "C" PPGL_EXPORT void CGAL_3D_Points_Inside_Triangles_C1_Bool(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points, std::vector<bool>& insides);
extern "C" PPGL_EXPORT void CGAL_3D_Points_Inside_Triangles_C2_Bool(const std::string& path, const Vector3d1& points, std::vector<bool>& insides);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Dart_Sampling_C1(const std::string & outside_path, const double& d, Vector3d1 & sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Dart_Sampling_C2(const std::string & outside_path, const std::string & inside_path, const double& d, Vector3d1 & sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C1(const std::string & outside_path, const double& d, Vector3d1 & sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C2(const std::string & outside_path, const std::string & inside_path, const double& d, Vector3d1 & sampling_points);

extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Triangle(const Vector3d & p, const Vector3d & t_0, const Vector3d & t_1, const Vector3d & t_2);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Triangles(const Vector3d & p, const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Nearest_Point_Triangles(const Vector3d & p, const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2);
extern "C" PPGL_EXPORT void CGAL_3D_Distance_Point_Mesh(const std::string & path, const Vector3d1 & query_points, std::vector<double>&distances);
extern "C" PPGL_EXPORT void CGAL_3D_Neareast_Point_Mesh(const std::string & path, const Vector3d1 & ves, Vector3d1 & ners);
extern "C" PPGL_EXPORT void  CGAL_3D_Mesh_Near_Triangles(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const Vector3d1 & points, const double& d, std::vector<std::vector<int>>&triangles);

extern "C" PPGL_EXPORT void CGAL_3D_Points_inside_Triangles_C1(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const Vector3d1 & points, std::vector<bool>&insides);
extern "C" PPGL_EXPORT void CGAL_3D_Points_inside_Triangles_C2(const std::string & path, const Vector3d1 & points, std::vector<bool>&insides);

extern "C" PPGL_EXPORT void CGAL_Mesh_Subdivision(const std::string & in_path, const std::string & sub, const int& step, const std::string & out_path);
extern "C" PPGL_EXPORT void CGAL_Mesh_Loop_Subdivision_One_Step(Vector3d1 & vecs, std::vector<int>&face_id_0, std::vector<int>&face_id_1, std::vector<int>&face_id_2);

/////////////////////////////////////////////////////////////


extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const std::vector<bool>&bools)
{
	std::vector<bool>().swap(bools);

	std::vector<std::vector<int>> vecs_neigbor(vecs.size(), std::vector<int>());
	std::vector<std::vector<int>> vecs_neigbor_lable(vecs.size(), std::vector<int>());
	std::vector<int> edges;
	for (int i = 0; i < face_id_0.size(); i++) {
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		edges.push_back(index_0);
		edges.push_back(index_1);
		edges.push_back(index_1);
		edges.push_back(index_2);
		edges.push_back(index_2);
		edges.push_back(index_0);

		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
	}

	for (int i = 0; i < edges.size(); i = i + 2) {
		int index_0 = edges[i];
		int index_1 = edges[i + 1];

		int search_0 = -1;
		for (int j = 0; j < vecs_neigbor[index_0].size() && search_0 < 0; j++) {
			if (vecs_neigbor[index_0][j] == index_1) {
				search_0 = j;
				vecs_neigbor_lable[index_0][j]++;
			}
		}

		if (search_0 < 0) {
			vecs_neigbor[index_0].push_back(index_1);
			vecs_neigbor_lable[index_0].push_back(1);
		}
	}

	for (int i = 0; i < vecs.size(); i++) {
		bool b = false;
		for (int j = 0; j < vecs_neigbor_lable[i].size() & !b; j++) {
			if (vecs_neigbor_lable[i][j] == 1) {
				b = true;
			}
		}
		bools.push_back(b);
	}

	std::vector<std::vector<int>>().swap(vecs_neigbor);
	std::vector<std::vector<int>>().swap(vecs_neigbor_lable);
	std::vector<int>().swap(edges);
}
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary(const std::string & path, const std::vector<bool>&bools)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, bools);
}

extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary(Vector3d1 & vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, Vector3d2 & boundaries, Vector3d & inside = Vector3d(0.0, 0.0, 0.0))
{
	Vector3d2 segments;

	std::vector<std::vector<int>> vecs_neigbor(vecs.size(), std::vector<int>());
	std::vector<std::vector<int>> vecs_neigbor_lable(vecs.size(), std::vector<int>());
	std::vector<int> edges;
	for (int i = 0; i < face_id_0.size(); i++) {
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		edges.push_back(index_0);
		edges.push_back(index_1);
		edges.push_back(index_1);
		edges.push_back(index_2);
		edges.push_back(index_2);
		edges.push_back(index_0);

		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
	}

	for (int i = 0; i < edges.size(); i = i + 2) {
		int index_0 = edges[i];
		int index_1 = edges[i + 1];

		int search_0 = -1;
		for (int j = 0; j < vecs_neigbor[index_0].size() && search_0 < 0; j++) {
			if (vecs_neigbor[index_0][j] == index_1) {
				search_0 = j;
				vecs_neigbor_lable[index_0][j]++;
			}
		}

		if (search_0 < 0) {
			vecs_neigbor[index_0].push_back(index_1);
			vecs_neigbor_lable[index_0].push_back(1);
		}
	}

	for (int i = 0; i < vecs.size(); i++)
	{
		int index_0 = i;

		for (int j = 0; j < vecs_neigbor_lable[i].size(); j++)
		{
			if (vecs_neigbor_lable[i][j] == 1)
			{
				int index_1 = vecs_neigbor[i][j];

				Vector3d1 segment;
				segment.push_back(vecs[index_0]);
				segment.push_back(vecs[index_1]);

				segments.push_back(segment);

				//delete
				for (int k = 0; k < vecs_neigbor[index_1].size(); k++)
				{
					if (vecs_neigbor[index_1][k] == index_0)
					{
						vecs_neigbor_lable[index_1][k] = 0;
					}
				}

			}
		}
	}

	CGAL_3D_Connecting_Segments(segments, boundaries);

	//CGAL_3D_Triangel_Mesh_Most_Inside_Point(vecs,face_id_0,face_id_1,face_id_2,inside);

	for (int i = 0; i < segments.size(); i++)
		Vector3d1().swap(segments[i]);
	Vector3d2().swap(segments);
	std::vector<std::vector<int>>().swap(vecs_neigbor);
	std::vector<std::vector<int>>().swap(vecs_neigbor_lable);
	std::vector<int>().swap(edges);
}
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary(const Vector3d1 & vecs, const std::vector<std::vector<int>>&face_ids, const Vector3d2 & boundaries)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (auto& face : face_ids)
	{
		face_id_0.emplace_back(face[0]);
		face_id_1.emplace_back(face[1]);
		face_id_2.emplace_back(face[2]);
	}

	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, boundaries);
}
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary(std::string path, Vector3d2 & boundaries, Vector3d & inside = Vector3d(0.0, 0.0, 0.0))
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, boundaries, inside);
}


extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth(const std::string & in_path, const std::string & out_path, const int laplace_nb)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(in_path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_Mesh_Laplace_Smooth(vecs, face_id_0, face_id_1, face_id_2, laplace_nb);
	CGAL_Output_Obj_C2(out_path, vecs, face_id_0, face_id_1, face_id_2);
}

extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth(Vector3d1 & vecs, std::vector<int>&face_id_0, std::vector<int>&face_id_1, std::vector<int>&face_id_2, const int laplace_nb)
{
	std::vector<bool> vertices_boundary;
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, vertices_boundary);
	std::vector<std::vector<int>> neighs;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(vecs, face_id_0, face_id_1, face_id_2, neighs);

	for (int iter = 0; iter < laplace_nb; iter++)
	{
		Vector3d1 new_vecs;
		for (int i = 0; i < vecs.size(); i++)
		{
			if (!vertices_boundary[i])
			{
				Vector3d v(0.0, 0.0, 0.0);
				double w = 0.0;
				for (int j = 0; j < neighs[i].size(); j++)
				{
					double d = CGAL_3D_Distance_Point_Point(vecs[neighs[i][j]], vecs[i]);
					v += vecs[neighs[i][j]] * (double)(1.0 / d);
					w += (1.0 / d);
				}
				v = vecs[i] * (double)0.5 + (double)0.5 * v / (double)w;

				new_vecs.push_back(v);
			}
			else
			{
				new_vecs.push_back(vecs[i]);
			}

		}
		vecs = new_vecs;
	}
}

extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth_by_Curvature(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const double& low_curvature)
{
	std::vector<double> max_curvature;
	std::vector<double> min_curvature;
	Vector3d1 max_curvature_direction;
	Vector3d1 min_curvature_direction;

	std::vector<bool> boundary;
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, boundary);

	Vector3d1 vecs_normals;

	int last_nb = 0;

	std::vector<std::vector<int>> vecs_neighbors;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(vecs, face_id_0, face_id_1, face_id_2, vecs_neighbors);

	std::vector<std::vector<std::vector<int>>> surface_vectices_to_neighbor_edges;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges(vecs, face_id_0, face_id_1, face_id_2, surface_vectices_to_neighbor_edges);

	int stop = 0;

	double par_0 = 0.1;
	double par_1 = 0.6;
	double par_2 = 0.3;

	for (int iter = 0; iter < 500; iter++)
		//while (true)
	{
		//compute vertices curvature
		std::vector<double>().swap(max_curvature);
		std::vector<double>().swap(min_curvature);
		Vector3d1().swap(max_curvature_direction);
		Vector3d1().swap(min_curvature_direction);
		Vector3d1().swap(vecs_normals);

		CGAL_3D_Mesh_Curvature(vecs, face_id_0, face_id_1, face_id_2, max_curvature, min_curvature, max_curvature_direction, min_curvature_direction, vecs_normals);

		int nb = 0;
		double minimal_cur = 100000.0;
		for (int i = 0; i < vecs.size(); i++)
		{
			if (min_curvature[i] < low_curvature && !boundary[i])
			{
				nb++;
				minimal_cur = std::min(minimal_cur, min_curvature[i]);
			}
		}


		//terminal condition

		if (nb < 5)break;

		if (abs(last_nb - nb) < 2)
			stop++;
		else
			stop = 0;

		if (stop == 60)
		{
			//break;
			par_0 = 0.10;
			par_1 = 0.65;
			par_2 = 0.25;
		}

		if (stop == 100)
		{
			par_0 = 0.10;
			par_1 = 0.70;
			par_2 = 0.20;
		}

		if (stop == 300)
		{
			break;
		}

		std::cout << "Current low curvature points number: " << stop << " " << nb << " " << minimal_cur << std::endl;

		last_nb = nb;

		//one iteration
		Vector3d1 iteration_vecs;

		for (int i = 0; i < vecs.size(); i++)
		{
			if (boundary[i])
			{
				iteration_vecs.push_back(vecs[i]);
			}
			else
			{
				bool run = min_curvature[i] < low_curvature + 0.1;

				for (int j = 0; j < vecs_neighbors[i].size() && !run; j++)
				{
					if (min_curvature[vecs_neighbors[i][j]] < low_curvature)
					{
						run = true;
					}

					if (true)
					{
						int vertice_id = vecs_neighbors[i][j];

						for (int k = 0; k < vecs_neighbors[vertice_id].size() && !run; k++)
						{
							if (min_curvature[vecs_neighbors[vertice_id][k]] < low_curvature)
							{
								run = true;
							}
						}
					}
				}

				if (run)
				{
					Vector3d cur_v = vecs[i] + Math::SetVectorLength(vecs_normals[i], 0.001 * min_curvature[i] / low_curvature);

					Vector3d smooth_v(0.0, 0.0, 0.0);

					double weight = 0.0;
					for (int j = 0; j < vecs_neighbors[i].size(); j++)
					{
						double l = CGAL_3D_Distance_Point_Point(vecs[vecs_neighbors[i][j]], vecs[i]);
						smooth_v += vecs[vecs_neighbors[i][j]] * (float)l;
						weight += l;
					}
					smooth_v = smooth_v / (float)weight;

					if (min_curvature[i] < 0.0 && max_curvature[i]>0.0)
					{
						//smooth_v = vecs[i];
						iteration_vecs.push_back((float)par_0 * vecs[i] + (float)par_1 * cur_v + (float)par_2 * smooth_v);
					}
					else
						iteration_vecs.push_back((float)par_0 * vecs[i] + (float)par_1 * cur_v + (float)par_2 * smooth_v);

					//iteration_vecs.push_back(LaplaceMeshSmoothForOnePoint(i, vecs, min_curvature_direction, vecs_normals, surface_vectices_to_neighbor_edges));
				}
				else
				{
					iteration_vecs.push_back(vecs[i]);
				}
			}
		}

		Vector3d1().swap(vecs);
		vecs = iteration_vecs;
		Vector3d1().swap(iteration_vecs);
	}
}

extern "C" PPGL_EXPORT void CGAL_Mesh_Loop_Subdivision_Own_Version(const std::string & in_path, const int& step, const std::string & out_path, const int& laplace_nb = 0)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	CGAL_3D_Read_Triangle_Mesh(in_path, vecs, face_id_0, face_id_1, face_id_2);

	for (int i = 0; i < step; i++)
	{
		CGAL_Mesh_Loop_Subdivision_One_Step(vecs, face_id_0, face_id_1, face_id_2);
		CGAL_Mesh_Laplace_Smooth(vecs, face_id_0, face_id_1, face_id_2, laplace_nb);
	}

	CGAL_Output_Obj(out_path, vecs, face_id_0, face_id_1, face_id_2);

	Vector3d1().swap(vecs);
	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);
}


#endif
