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

struct GridEdge {
	int s_i, s_j;
	int e_i, e_j;
	int type;
	GridEdge(int s_i_0, int s_j_0, int e_i_0, int e_j_0, int t) {
		s_i = s_i_0;
		s_j = s_j_0;
		e_i = e_i_0;
		e_j = e_j_0;
		type = t;
	}
};

Vector3d PointVector3d(Point_3 p);
Point_3 VectorPoint3d(Vector3d p);
Vector2d PointVector2d(Point_2 p);
Point_2 VectorPoint2d(Vector2d p);

void  Construct_Polyhedron(Polyhedron_3& polyhedron, const Vector3d1& vecs, const Vector1i1& face_id_0, const Vector1i1& face_id_1, const Vector1i1& face_id_2);
void  Construct_Polyhedron(Polyhedron_3& polyhedron, const char* path);
void  Construct_Polyhedron(Polyhedron_3& polyhedron, const char* path, Vector3d1& vecs, Vector1i1& face_id_0, Vector1i1& face_id_1, Vector1i1& face_id_2);


Point_3 point_to_3d(const Point_2& p, Plane_3& pl);
Point_2 point_to_2d(const Point_3& p, Plane_3& pl);

//Project p onto the planar surface of 3d triangle
//Checking the position relationship between the p and 3d triangle
//face: 3d triangle
//p: 3d point
//return true: inside
//return false: outside
bool OutsidePointInsideTriangle(Poly_facet_iterator& face, Vector3d p);
bool Intersection(Halfedge_handle& hh, int nb, Vector3d inside, Vector3d outside, Halfedge_handle& handle, Vector3d& intersection);
Vector3d RelatedFaceNormal(Polyhedron_3& polyhedron, Tree& tree, Vector3d1& normals, Vector3d source);
void RelatedFaceAndBarycentric(const Polyhedron_3& polyhedron, const Tree& tree, const Vector3d& source, double& u, double& v, double& w, Point_3& nearest_point, face_iterator& face_it);
Vector3d RelatedFaceNormal(Polyhedron_3& polyhedron, Tree& tree, Vector3d1& normals, Vector3d source);

bool DebugInformation();
void ComputeEdgeLables(const int size_of_vertice, Halfedge_handle& start_hh, std::vector<std::pair<int, int>>& edges, std::vector<int>& lables);
void CGAL_2D_Polygon_Simple_0(Vector2d1 points_2d);
void ComputeRemeshTriangles(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2,
	const std::vector<std::pair<int, int>>& edges, const Vector3d1& cutting_points, const std::vector<int>& multi_cps, const std::vector<int>& two_class_lables, const std::string output_path);
void insert_polygon(CDT& cdt, const Polygon_2& polygon, Vector1i1& indexInt);
void mark_domains(CDT& ct, CDT::Face_handle start, int index, std::list<CDT::Edge>& border);
void mark_domains(CDT& cdt);

Polygon_2 Polygon2D(const Vector2d1& py);

int GetIndex(std::vector<GridEdge>& grid_edges, std::pair<int, int> i_0, std::pair<int, int> i_1);

extern "C" PPGL_EXPORT void CGAL_Test_PGL(const Vector3d& n, const char* str, const char* char_);

//implementation in "io.cpp"
//####################################################################################
extern "C" PPGL_EXPORT void CGAL_Vector_Base(const Vector3d& n, Vector3d&);

//implementation in "twoD.cpp"
//####################################################################################
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Point(const Vector2d& p_0, const Vector2d& p_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Line(const Vector2d& v, const Vector2d& l_0, const Vector2d& l_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Segment(const Vector2d& v, const Vector2d& s_0, const Vector2d& s_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Segment_Segment(const Vector2d& s_0, const Vector2d& s_1, const Vector2d& e_0, const Vector2d& e_1);
extern "C" PPGL_EXPORT bool CGAL_2D_Location_Point_Polygon(const Vector2d& p, const Vector2d1& py);
extern "C" PPGL_EXPORT bool CGAL_2D_Location_Points_Polygon(const Vector2d1& ps, const Vector2d1& py);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_2D_Polygon_Dart_Sampling(const Vector2d1& py, const double& d, Vector2d1& sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Polygon_Regular_Sampling_C1(const Vector2d1& py, const double& d);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Polygon_Regular_Sampling_C2(const Vector2d1& py, const double& d, VectorPI1& neighbors);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Polygon_Regular_Sampling_C3(const Vector2d1& py, const double& d, VectorPI1& neighbors, const bool& compute_neighbors);

extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Square_Regular_Sampling_C1(const double& d);
extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Square_Regular_Sampling_C2(const double& d, VectorPI1& neighbors);
extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Square_Regular_Sampling_C3(const double& d, VectorPI1& neighbors, const bool& compute_neighbors);



extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Polygon(const Vector2d& p, const Vector2d1& py);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Polygons(const Vector2d& p, const Vector2d2& pys);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Segment(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Line_Line(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Line(const Vector2d& s_s, const Vector2d& s_e, const Vector2d& l_s, const Vector2d& l_e, Vector2d& inter);

extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Polygon(const Vector2d& s_s, const Vector2d& s_e, const Vector2d1& p);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Polygon_Polygon(const Vector2d1 & p1, const Vector2d1 & p2);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Is_Clockwise_Oriented(const Vector2d1& ps);
extern "C" PPGL_EXPORT double CGAL_2D_Two_Polygons_Union(const Vector2d1& poly_0, const Vector2d1& poly_1, Vector2d2& inter_polygons);
extern "C" PPGL_EXPORT double CGAL_2D_Two_Polygons_Intersection(const Vector2d1& poly_0, const Vector2d1& poly_1);
extern "C" PPGL_EXPORT void CGAL_Decompose_Polyline(const Vector2d1& polyline, const double& threshold, Vector1i1& result);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_Extend(const Vector2d1& polygon, const Vector2d& s, const Vector2d& e, Vector2d& ns, Vector2d& ne);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_NotExtend(const Vector2d1& polygon, const Vector2d& s, const Vector2d& e);
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut(const Vector2d1& polygon, const Vector2d1& cutLine, VectorPB1& result);

extern "C" PPGL_EXPORT bool CGAL_Construct_InOutSide_Polygon(const Vector2d1& py, const Vector2d& p, const Vector2d& q, bool& isPInside, bool& isQInside);
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Ray_Segment(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
extern "C" PPGL_EXPORT double CGAL_Get_Angle_Kerf_Offset_Tan(const Vector2d& a, const Vector2d& b);
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Projection_Point_Segment(const Vector2d& p, const Vector2d& s, const Vector2d& e);

extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C1(const Vector2d1& outside_py, const Vector2d& p);
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C2(const Vector2d1& outside_py, const Vector2d1& inside_py);
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C3(const Vector2d2& outside_pys, const Vector2d& p);
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C4(const Vector2d2& outside_pys, const Vector2d1& inside_py);
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C5(const Vector2d2& outside_pys, const Vector2d2& inside_pys);

extern "C" PPGL_EXPORT double CGAL_2D_Distance_Polygon_Polygon(const Vector2d1& poly_0, const Vector2d1& poly_1);
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Polygons_Polygons(const Vector2d2& poly_0, const Vector2d2& poly_1);

extern "C" PPGL_EXPORT Vector2d CGAL_2D_Nearest_Point_Polygon_C1(const Vector2d& v, const Vector2d1& poly);
extern "C" PPGL_EXPORT void CGAL_2D_Nearest_Point_Polygon_C2(const Vector2d& v, const Vector2d1& poly, Vector2d& p, double& min_d);
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Nearest_Point_Polygons(const Vector2d& v, const Vector2d2& polys);

extern "C" PPGL_EXPORT void CGAL_2d_Polygon_Boundingbox(const Vector2d1& ps, Vector2d& min_corner, Vector2d& max_corner);
extern "C" PPGL_EXPORT double CGAL_2D_Polygon_Area(const Vector2d1& py);
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Polygon_Inside_Point_C1(const Vector2d1& poly);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Inside_Point_C2(const Vector2d2& polys, Vector2d& inner_vec);

extern "C" PPGL_EXPORT void CGAL_2D_Polygon_One_Offsets(const Vector2d1& poly, const double& d, Vector2d2& offset_polys);
extern "C" PPGL_EXPORT void CGAL_2D_Polygons_One_Offsets(const Vector2d2& polys, const double& d, Vector2d2& offset_polys);

extern "C" PPGL_EXPORT bool CGAL_2D_Polygons_Simple(const Vector2d2& poly);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Simple(const Vector2d1& poly);
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Simple_Inter(const Vector2d1& poly);

extern "C" PPGL_EXPORT void CGAL_2D_Convex_Hulls(const Vector2d1& vec, Vector2d1& hull_points);

extern "C" PPGL_EXPORT void CGAL_2D_OBB_Box(const Vector2d1& vec, Vector2d& center, Vector2d& axis_0, Vector2d& axis_1, double& entent_0, double& entent_1);

extern "C" PPGL_EXPORT void CGAL_Image_Grid_Decomposition_C1(std::vector<std::vector<int>>& image, std::vector<std::vector<double>>& boundary_xs, std::vector<std::vector<double>>& boundary_ys);
extern "C" PPGL_EXPORT void CGAL_Image_Grid_Decomposition_Conservative_C1(std::vector<std::vector<int>>& image, std::vector<std::vector<double>>& boundary_xs, std::vector<std::vector<double>>& boundary_ys);
extern "C" PPGL_EXPORT void CGAL_Image_Grid_Decomposition_C2(std::vector<std::vector<int>>& image, Vector2d2& boundaries);
extern "C" PPGL_EXPORT void CGAL_Image_Grid_Decomposition_Conservative_C2(std::vector<std::vector<int>>& image, Vector2d2& boundaries);



//implementation in "threeD.cpp"
//####################################################################################
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Segment(const Vector3d& p, const Vector3d& s_s, const Vector3d& s_e);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Fitting(const Vector3d1& points, Vector3d& plane_p, Vector3d& plane_n);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Point_Projection(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d& p, Vector3d& result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_Points_Projection(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d1& points, Vector3d1& project_points);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_3D_to_2D_Point(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d& point_3d, Vector2d& result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_2D_to_3D_Point(const Vector3d& plane_p, const Vector3d& plane_n, const Vector2d& points_2d, Vector3d& result);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_3D_to_2D_Points(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d1& points_3d, Vector2d1& points_2d);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_3Ds_to_2Ds_Points(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d2& points_3d, Vector2d2& points_2d);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_2D_to_3D_Points(const Vector3d& plane_p, const Vector3d& plane_n, const Vector2d1& points_2d, Vector3d1& points_3d);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Segment(const Vector3d& p, const Vector3d& s_s, const Vector3d& s_e);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Point(const Vector3d& v0, const Vector3d& v1);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Polygon(const Vector3d1& py, const Vector3d& p);
extern "C" PPGL_EXPORT void CGAL_2D_Polygon_Triangulation_C1(const Vector2d2& polys, Vector1i2& faces);

extern "C" PPGL_EXPORT std::vector<std::vector<int>> CGAL_2D_Polygon_Triangulation_C2(const Vector2d2& polys);
extern "C" PPGL_EXPORT std::vector<std::vector<int>> CGAL_2D_Polygon_Triangulation_C3(const Vector2d1& poly);

extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Line(const Vector3d& p, const Vector3d& l_s, const Vector3d& l_e);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Line(const Vector3d& p, const Vector3d& l_s, const Vector3d& l_e);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Segment_Segment(const Vector3d& s_0_s, const Vector3d& s_0_e, const Vector3d& s_1_s, const Vector3d& s_1_e);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Plane(const Vector3d& v, const Vector3d& plane_p, const Vector3d& plane_n);

extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Line(const Vector3d& s_s, const Vector3d& s_e, const Vector3d& l_s, const Vector3d& l_e, Vector3d& inter);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Segment(const Vector3d& s_0_s, const Vector3d& s_0_e, const Vector3d& s_1_s, const Vector3d& s_1_e, Vector3d& iter);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Plane(const Vector3d& s_s, const Vector3d& s_e, const Vector3d& plane_p, const Vector3d& plane_n, Vector3d& inter);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Line_Plane(const Vector3d& l_s, const Vector3d& l_e, const Vector3d& plane_p, const Vector3d& plane_n, Vector3d& inter);

extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Plane_C1(const Vector3d& p, const Vector3d& plane_p, const Vector3d& plane_n);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Plane_C2(const Vector3d& p, const Vector3d& plane_p_0, const Vector3d& plane_p_1, const Vector3d& plane_p_2);
extern "C" PPGL_EXPORT Vector2d CGAL_3D_Projection_3D_Point_Plane_2D_C1(const Vector3d& p, const Vector3d& plane_p, const Vector3d& plane_n);
extern "C" PPGL_EXPORT Vector2d CGAL_3D_Projection_3D_Point_Plane_2D_C2(const Vector3d& p, const Vector3d& plane_p_0, const Vector3d& plane_p_1, const Vector3d& plane_p_2);
extern "C" PPGL_EXPORT void CGAL_3D_Plane_ABCD(const Vector3d& plane_p, const Vector3d& plane_n, double& a, double& b, double& c, double& d);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Plane_Base_1(const Vector3d& plane_p, const Vector3d& plane_n);


extern "C" PPGL_EXPORT Vector3d CGAL_Face_Normal(const Vector3d& source, const Vector3d& tri_0, const Vector3d& tri_1, const Vector3d& tri_2, Vector3d& normal_0, Vector3d& normal_1, Vector3d& normal_2);


//implementation in "mesh.cpp"
//####################################################################################
extern "C" PPGL_EXPORT void CGAL_Remesh_Surface_by_Adding_Feature(const Vector3d1& feature, const Vector1i1& face_ids, const Vector3d1& vecs, const Vector1i1& face_id_0, const Vector1i1& face_id_1, const Vector1i1& face_id_2, Vector1i1& igl_cutting_0_edges, Vector1i1& igl_cutting_1_edges, Vector3d1& igl_cutting_points, Vector1i2& cutting_faces);
extern "C" PPGL_EXPORT void CGAL_Mesh_Edges(const char* path);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Sphere_Ray(const double& center_x, const double& center_y, const double& center_z, const double& radius, const double& ray_origin_x, const double& ray_origin_y, const double& ray_origin_z, const double& ray_direction_x, const double& ray_direction_y, const double& ray_direction_z, std::vector<double>& i_x, std::vector<double>& i_y, std::vector<double>& i_z);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Triangle(const Vector3d& p, const Vector3d& n, const Vector3d& p0, const Vector3d& p1, const Vector3d& p2);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Mesh(const Vector3d& p, const Vector3d& n, const char* path);
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Mesh(const Vector3d& s, const Vector3d& e, const char* path);
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Segments_Mesh(const Vector3d1& ss, const Vector3d1& ee, const char* path, Vector1b1& inters);
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Polygons_Mesh(const Vector3d2& polygons, const char* path, Vector1b1& inters);
//check whether there is a polygon intersected with the input mesh
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Polygons_Mesh_Bool(const Vector3d2& polygons, const char* path);
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_Vector3d(const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector3d1& inters);
//test each group directions (nes[i]) for each point in ps
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_C1_Bool(const Vector3d1& ps, const Vector3d2& nes, const char* path, Vector1b2& inters);
//test all directions (ns) for each point in ps
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_C2_Bool(const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector1b2& inters);
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d(const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector1d2& inters);
extern "C" PPGL_EXPORT void CGAL_3D_Points_Inside_Triangles_C1_Bool(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points, std::vector<bool>& insides);
extern "C" PPGL_EXPORT void CGAL_3D_Points_Inside_Triangles_C2_Bool(const char* path, const Vector3d1& points, std::vector<bool>& insides);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Dart_Sampling_C1(const char* outside_path, const double& d, Vector3d1& sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Dart_Sampling_C2(const char* outside_path, const char* inside_path, const double& d, Vector3d1& sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C1(const char* outside_path, const double& d, Vector3d1& sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C2(const char* outside_path, const char* inside_path, const double& d, Vector3d1& sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Cube_Surface_Sampling_C1(const double& cube_size, const double& d, Vector3d2& sampling_points, VectorPI2& neighbors, const bool& compute_neighbors);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Cube_Surface_Sampling_C2(const double& cube_size, const double& d, Vector3d2& sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Cube_Surface_Sampling_C3(const double& cube_size, const double& d, Vector3d2& sampling_points, VectorPI2& neighbors);

//with neighboring
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C3(const char* outside_path, const char* inside_path, const double& d, Vector3d1& sampling_points, VectorPI1& neighbors);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Triangle(const Vector3d& p, const Vector3d& t_0, const Vector3d& t_1, const Vector3d& t_2);
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Triangles(const Vector3d& p, const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Nearest_Point_Triangles(const Vector3d& p, const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2);
extern "C" PPGL_EXPORT void CGAL_3D_Distance_Point_Mesh(const char* path, const Vector3d1& query_points, std::vector<double>& distances);
extern "C" PPGL_EXPORT void CGAL_3D_Neareast_Point_Mesh(const char* path, const Vector3d1& ves, Vector3d1& ners);
extern "C" PPGL_EXPORT void  CGAL_3D_Mesh_Near_Triangles(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points, const double& d, std::vector<std::vector<int>>& triangles);

extern "C" PPGL_EXPORT void CGAL_3D_Points_inside_Triangles_C1(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points, std::vector<bool>& insides);
extern "C" PPGL_EXPORT void CGAL_3D_Points_inside_Triangles_C2(const char* path, const Vector3d1& points, std::vector<bool>& insides);

extern "C" PPGL_EXPORT void CGAL_Mesh_Subdivision(const char* in_path, const char* sub, const int& step, const char* out_path);
extern "C" PPGL_EXPORT void CGAL_Mesh_Loop_Subdivision_One_Step(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2);

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C1(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, std::vector<double>& max_curs, std::vector<double>& min_curs);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C2(const Vector3d1& vecs, const std::vector<std::vector<int>>& face_ids, std::vector<double>& max_curs, std::vector<double>& min_curs);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C3(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C4(const Vector3d1& vecs, const std::vector<std::vector<int>>& face_ids, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C5(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions, Vector3d1& normals);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C6(const Vector3d1& vecs, const std::vector<std::vector<int>>& face_ids, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions, Vector3d1& normals);
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C1(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, std::vector<bool>& bools);
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C2(const char* path, std::vector<bool>& bools);
extern "C" PPGL_EXPORT void CGAL_3D_Connecting_Segments_C1(Vector2d2& segments, Vector2d2& lines);
extern "C" PPGL_EXPORT void CGAL_3D_Connecting_Segments_C2(Vector3d2& segments, Vector3d2& lines);
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C3(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, Vector3d2& boundaries);
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C4(Vector3d1& vecs, std::vector<std::vector<int>>& face_ids, Vector3d2& boundaries);
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C5(const char* path, Vector3d2& boundaries);
extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth_C1(const char* in_path, const char* out_path, const int laplace_nb);
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Vecs_Neighbors(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<std::vector<int>>& neighs);
extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth_C2(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, const int laplace_nb);
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Vecs_Faces(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<std::vector<int>>& surface_vectices_to_face);
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<std::vector<std::vector<int>>>& surface_vectices_to_neighbor_edges);
extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth_by_Curvature(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, double& low_curvature);
extern "C" PPGL_EXPORT void CGAL_Mesh_Loop_Subdivision_Own_Version(const char* in_path, const int& step, const char* out_path, const int& laplace_nb);
extern "C" PPGL_EXPORT void CGAL_Rotation_Obj(const char* path, const double& angle, const Vector3d& axis);
extern "C" PPGL_EXPORT void CGAL_Slicer_Mesh(const char* path, const Vector3d& plane_normal, const std::vector<double>& plane_d, Vector3d3& offsetses, Vector3d2& offsets);
extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C1(const char* path, Vector3d1& xyzs);
extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C3(const char* path, Vector3d source, Vector3d target, Vector3d1& output);
extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C4(const char* path, Vector3d1 sources, Vector3d1 targets, Vector3d2& xyzes);
extern "C" PPGL_EXPORT double CGAL_Geodesic_Distance(const char* path, const Vector3d& source, const Vector3d& target);
extern "C" PPGL_EXPORT Vector3d1 CGAL_Project_Points_Onto_Surface_C1(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points);
extern "C" PPGL_EXPORT Vector3d1 CGAL_Project_Points_Onto_Surface_C2(const char* path, const Vector3d1& points);
extern "C" PPGL_EXPORT void CGAL_3D_Triangel_Mesh_Most_Inside_Point(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, Vector3d& inside);
extern "C" PPGL_EXPORT double CGAL_3D_One_Triangle_Area(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2);
extern "C" PPGL_EXPORT double CGAL_3D_Triangle_Mesh_Area(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2);
extern "C" PPGL_EXPORT void CGAL_3D_Convex_Hulls_C1(const Vector3d1& vec, Vector3d1& hull_points);
extern "C" PPGL_EXPORT void CGAL_3D_Convex_Hulls_C2(const Vector3d1& vec, Vector3d1& hull_points, std::vector<int>& hulls_surface_0, std::vector<int>& hulls_surface_1, std::vector<int>& hulls_surface_2);
extern "C" PPGL_EXPORT void CGAL_3D_Convex_Hulls_C3(const Vector3d1& vec, Vector3d1& hull_points, Vector3d1& plane_p, Vector3d1& plane_n);
extern "C" PPGL_EXPORT void CGAL_3D_Convex_Hulls_C4(const Vector3d1& vec, Vector3d1& hull_points, std::vector<int>& hulls_surface_0, std::vector<int>& hulls_surface_1, std::vector<int>& hulls_surface_2, Vector3d1& plane_p, Vector3d1& plane_n);
extern "C" PPGL_EXPORT void CGAL_Mesh_Field_Query_C1(const char* path, const Vector3d1& gradients, const Vector3d1& input_points, Vector3d1& points_gradients);
extern "C" PPGL_EXPORT void CGAL_Mesh_Field_Query_C2(const char* path, const std::vector<double>& gradient_values, const Vector3d1& input_points, std::vector<double>& points_gradient_values);
extern "C" PPGL_EXPORT void CGAL_Mesh_Field_Query_C3(const char* path, const std::vector<double>& gradient_values, const Vector3d2& input_point_es, std::vector<std::vector<double>>& points_gradient_value_es);
extern "C" PPGL_EXPORT void CGAL_Curvature_Mesh(const char* path, const Vector3d1& input_points, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions);

extern "C" PPGL_EXPORT void CGAL_Normal_Mesh_C1(const char* path, const Vector3d1& mesh_points, Vector3d1& mesh_normals);
extern "C" PPGL_EXPORT void CGAL_Normal_Mesh_C2(const char* path, const Vector3d2& mesh_pointses, Vector3d2& mesh_normalses);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Normal_C1(const Vector3d1& ps, const std::vector<std::vector<int>>& face_ids, Vector3d1& normals);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Normal_C2(const Vector3d1& ps, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, Vector3d1& normals);

extern "C" PPGL_EXPORT Vector3d CGAL_3D_Mesh_Center_C1(const Vector3d2& ps);
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Mesh_Center_C2(const Vector3d1& ps);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Boundingbox_C1(const Vector3d2& ps, Vector3d& min_corner, Vector3d& max_corner);
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Boundingbox_C2(const Vector3d1& ps, Vector3d& min_corner, Vector3d& max_corner);

extern "C" PPGL_EXPORT void CGAL_Surface_Decomposition(const char* path, std::vector<double>& face_sdf, int& regions_nb, std::vector<int>& face_segments);

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Gradient(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const std::vector<double>& psd, Vector3d1& vecs_gradients, Vector3d1& faces_gradients);

extern "C" PPGL_EXPORT void CGAL_Intergral_Curvature(const Vector2d1& input_points, const int& sampling_points_nb, const double& radius, const double& thresholder, Vector2d1& output_points, std::vector<double>& output_rates);
extern "C" PPGL_EXPORT bool CGAL_3D_Mesh_Extract_Isoline(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const std::vector<double>& psd, const double& d, Vector3d2& isolines);
extern "C" PPGL_EXPORT void CGAL_BSplineCurveFit(const Vector3d1& samples, Vector3d1& output);

extern "C" PPGL_EXPORT void CGAL_Cut_Surface(const Vector3d1& boundary, const Vector3d& inside_point, const char* full_path, char* output_path);
extern "C" PPGL_EXPORT void CGAL_Cut_Surface_by_Multi_Boundaries(const Vector3d2& multi_boundary, const Vector3d& inside_point, const char* full_path, char* output_path);



/////////////////////////////////////////////////////////////
//

#endif
