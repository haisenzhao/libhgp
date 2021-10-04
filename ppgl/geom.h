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
void  Construct_Polyhedron(Polyhedron_3& polyhedron, std::string path);
void  Construct_Polyhedron(Polyhedron_3& polyhedron, std::string path, Vector3d1& vecs, Vector1i1& face_id_0, Vector1i1& face_id_1, Vector1i1& face_id_2);

extern "C" PPGL_EXPORT void CGAL_Test_PGL(Vector3d n);


//implementation in "io.cpp"
//####################################################################################
extern "C" PPGL_EXPORT void CGAL_Vector_Base(Vector3d n, Vector3d &);
extern "C" PPGL_EXPORT void CGAL_Export_Path_Segment(std::ofstream &export_file_output, int &export_index,std::string s_name, double r, double g, double b,Vector3d &start,Vector3d &end, double radius);
extern "C" PPGL_EXPORT void CGAL_Export_Path_Point(std::ofstream &export_file_output, int &export_index,std::string s_name, double r, double g, double b, Vector3d point, double radius);

//implementation in "twoD.cpp"
//####################################################################################
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Point(Vector2d p_0, Vector2d p_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Line(Vector2d v, Vector2d l_0, Vector2d l_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Segment(Vector2d v, Vector2d s_0, Vector2d s_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Segment_Segment(Vector2d s_0, Vector2d s_1, Vector2d e_0, Vector2d e_1);
extern "C" PPGL_EXPORT bool CGAL_2D_Location_Point_Polygon(Vector2d p, Vector2d1 py);
extern "C" PPGL_EXPORT bool CGAL_2D_Location_Points_Polygon(const Vector2d1 &ps, const Vector2d1 &py);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Polygon(Vector2d p, Vector2d1 py);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Segment(Vector2d s_0_s, Vector2d s_0_e, Vector2d s_1_s, Vector2d s_1_e, Vector2d &inter);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Line_Line(const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s, const Vector2d &s_1_e, Vector2d &inter);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Polygon(Vector2d s_s, Vector2d s_e, Vector2d1 &p);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Is_Clockwise_Oriented(Vector2d1 &ps);
extern "C" PPGL_EXPORT double CGAL_2D_Two_Polygons_Union(Vector2d1 poly_0, Vector2d1 poly_1,Vector2d2 &inter_polygons);
extern "C" PPGL_EXPORT double CGAL_2D_Two_Polygons_Intersection(const Vector2d1 &poly_0, const Vector2d1 &poly_1);
extern "C" PPGL_EXPORT Vector1i1 CGAL_Decompose_Polyline(Vector2d1 &polyline, double threshold);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_Extend(const Vector2d1 &polygon, const Vector2d &s,const Vector2d &e, Vector2d &ns, Vector2d &ne);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_NotExtend(const Vector2d1 &polygon, const Vector2d &s,const Vector2d &e);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut(const Vector2d1 &polygon, const Vector2d1 &cutLine, VectorPB1 &result);
extern "C" PPGL_EXPORT void CGAL_2D_Polygon_One_Offsets(Vector2d1 &poly, double d,Vector2d2  &offset_polys);
extern "C" PPGL_EXPORT bool CGAL_Construct_InOutSide_Polygon(const Vector2d1 &py, const Vector2d &p, const Vector2d &q, bool &isPInside,bool &isQInside);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Ray_Segment(const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s,const Vector2d &s_1_e, Vector2d &inter);
extern "C" PPGL_EXPORT double CGAL_Get_Angle_Kerf_Offset_Tan(const Vector2d &a, const Vector2d &b);


//implementation in "threeD.cpp"
//####################################################################################
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Segment(Vector3d p, Vector3d s_s, Vector3d s_e);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Fitting(Vector3d1 &points, Vector3d &plane_p, Vector3d &plane_n);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Point_Projection(Vector3d &plane_p, Vector3d &plane_n, Vector3d &p,Vector3d &result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Points_Projection(Vector3d &plane_p, Vector3d &plane_n, Vector3d1 &points, Vector3d1 &project_points);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_3D_to_2D_Point(Vector3d &plane_p, Vector3d &plane_n, Vector3d &point_3d,Vector2d &result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_2D_to_3D_Point(Vector3d &plane_p, Vector3d &plane_n, Vector2d &points_2d, Vector3d &result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_3D_to_2D_Points(Vector3d &plane_p, Vector3d &plane_n, Vector3d1 &points_3d,Vector2d1 &points_2d);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_2D_to_3D_Points(Vector3d &plane_p, Vector3d &plane_n, Vector2d1 &points_2d,Vector3d1 &points_3d);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Segment(Vector3d p, Vector3d s_s, Vector3d s_e);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Point(const Vector3d & v0, const Vector3d & v1);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Polygon(const Vector3d1 &py, const Vector3d &p);
extern "C" PPGL_EXPORT void CGAL_2D_Polygon_Triangulation(const Vector2d2 &polys, Vector1i2 &faces);


//implementation in "mesh.cpp"
//####################################################################################

extern "C" PPGL_EXPORT void CGAL_Remesh_Surface_by_Adding_Feature(const Vector3d1 &feature,const Vector1i1 &face_ids, const Vector3d1 &vecs, const Vector1i1 &face_id_0,const Vector1i1 &face_id_1,const Vector1i1 &face_id_2, Vector1i1 &igl_cutting_0_edges,Vector1i1 &igl_cutting_1_edges, Vector3d1 &igl_cutting_points,Vector1i2 &cutting_faces);
extern "C" PPGL_EXPORT void CGAL_3D_Read_Triangle_Mesh(std::string path, Vector3d1 &vecs,Vector1i1 &face_id_0, Vector1i1 &face_id_1, Vector1i1 &face_id_2);
extern "C" PPGL_EXPORT void CGAL_Mesh_Edges(std::string path);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Sphere_Ray(double, double, double, double, double, double, double, double, double, double, std::vector<double>&, std::vector<double>&, std::vector<double>&);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Triangle(Vector3d p, Vector3d n, Vector3d p0, Vector3d p1, Vector3d p2);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Mesh(Vector3d p, Vector3d n, std::string path);
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_Vector3d(Vector3d1 ps, Vector3d1 ns, std::string path, Vector3d1 & inters);
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_Bool(Vector3d1 ps, Vector3d2 nes, std::string path, Vector1b2 & inters);


#endif
