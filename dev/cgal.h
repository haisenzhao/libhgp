#ifndef CGAL_ONCE
#define CGAL_ONCE
#pragma once
#include <pgl_functs.hpp>
using namespace std;
using namespace PGL;
namespace PPGL {
//Project p onto the planar surface of 3d triangle
//Checking the position relationship between the p and 3d triangle
//face: 3d triangle
//p: 3d point
//return true: inside
//return false: outside
typedef  void (*CGAL_Test_PGL)(const Vector3d& n, const char* str, const char* char_);
//implementation in "io.cpp"
//####################################################################################
typedef  void (*CGAL_Vector_Base)(const Vector3d& n, Vector3d&);
//implementation in "twoD.cpp"
//####################################################################################
typedef  double (*CGAL_2D_Distance_Point_Point)(const Vector2d& p_0, const Vector2d& p_1);
typedef  double (*CGAL_2D_Distance_Point_Line)(const Vector2d& v, const Vector2d& l_0, const Vector2d& l_1);
typedef  double (*CGAL_2D_Distance_Point_Segment)(const Vector2d& v, const Vector2d& s_0, const Vector2d& s_1);
typedef  double (*CGAL_2D_Distance_Segment_Segment)(const Vector2d& s_0, const Vector2d& s_1, const Vector2d& e_0, const Vector2d& e_1);
typedef  bool (*CGAL_2D_Location_Point_Polygon)(const Vector2d& p, const Vector2d1& py);
typedef  bool (*CGAL_2D_Location_Points_Polygon)(const Vector2d1& ps, const Vector2d1& py);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  void (*CGAL_2D_Polygon_Dart_Sampling)(const Vector2d1& py, const double& d, Vector2d1& sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  Vector2d1 (*CGAL_2D_Polygon_Regular_Sampling_C1)(const Vector2d1& py, const double& d);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  Vector2d1 (*CGAL_2D_Polygon_Regular_Sampling_C2)(const Vector2d1& py, const double& d, VectorPI1& neighbors);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  Vector2d1 (*CGAL_2D_Polygon_Regular_Sampling_C3)(const Vector2d1& py, const double& d, VectorPI1& neighbors, const bool& compute_neighbors);
typedef  Vector2d1 (*CGAL_2D_Square_Regular_Sampling_C1)(const double& d);
typedef  Vector2d1 (*CGAL_2D_Square_Regular_Sampling_C2)(const double& d, VectorPI1& neighbors);
typedef  Vector2d1 (*CGAL_2D_Square_Regular_Sampling_C3)(const double& d, VectorPI1& neighbors, const bool& compute_neighbors);
typedef  double (*CGAL_2D_Distance_Point_Polygon)(const Vector2d& p, const Vector2d1& py);
typedef  double (*CGAL_2D_Distance_Point_Polygons)(const Vector2d& p, const Vector2d2& pys);
typedef  bool (*CGAL_2D_Intersection_Segment_Segment)(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
typedef  bool (*CGAL_2D_Intersection_Line_Line)(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
typedef  bool (*CGAL_2D_Intersection_Segment_Line)(const Vector2d& s_s, const Vector2d& s_e, const Vector2d& l_s, const Vector2d& l_e, Vector2d& inter);
typedef  bool (*CGAL_2D_Intersection_Segment_Polygon)(const Vector2d& s_s, const Vector2d& s_e, const Vector2d1& p);
typedef  bool (*CGAL_2D_Polygon_Is_Clockwise_Oriented)(const Vector2d1& ps);
typedef  double (*CGAL_2D_Two_Polygons_Union)(const Vector2d1& poly_0, const Vector2d1& poly_1, Vector2d2& inter_polygons);
typedef  double (*CGAL_2D_Two_Polygons_Intersection)(const Vector2d1& poly_0, const Vector2d1& poly_1);
typedef  void (*CGAL_Decompose_Polyline)(const Vector2d1& polyline, const double& threshold, Vector1i1& result);
typedef  bool (*CGAL_Identify_Polycut_Extend)(const Vector2d1& polygon, const Vector2d& s, const Vector2d& e, Vector2d& ns, Vector2d& ne);
typedef  bool (*CGAL_Identify_Polycut_NotExtend)(const Vector2d1& polygon, const Vector2d& s, const Vector2d& e);
typedef  bool (*CGAL_Identify_Polycut)(const Vector2d1& polygon, const Vector2d1& cutLine, VectorPB1& result);
typedef  bool (*CGAL_Construct_InOutSide_Polygon)(const Vector2d1& py, const Vector2d& p, const Vector2d& q, bool& isPInside, bool& isQInside);
typedef  bool (*CGAL_2D_Intersection_Ray_Segment)(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
typedef  double (*CGAL_Get_Angle_Kerf_Offset_Tan)(const Vector2d& a, const Vector2d& b);
typedef  Vector2d (*CGAL_2D_Projection_Point_Segment)(const Vector2d& p, const Vector2d& s, const Vector2d& e);
typedef  bool (*CGAL_2D_Detect_Polygon_Inside_C1)(const Vector2d1& outside_py, const Vector2d& p);
typedef  bool (*CGAL_2D_Detect_Polygon_Inside_C2)(const Vector2d1& outside_py, const Vector2d1& inside_py);
typedef  bool (*CGAL_2D_Detect_Polygon_Inside_C3)(const Vector2d2& outside_pys, const Vector2d& p);
typedef  bool (*CGAL_2D_Detect_Polygon_Inside_C4)(const Vector2d2& outside_pys, const Vector2d1& inside_py);
typedef  bool (*CGAL_2D_Detect_Polygon_Inside_C5)(const Vector2d2& outside_pys, const Vector2d2& inside_pys);
typedef  double (*CGAL_2D_Distance_Polygon_Polygon)(const Vector2d1& poly_0, const Vector2d1& poly_1);
typedef  double (*CGAL_2D_Distance_Polygons_Polygons)(const Vector2d2& poly_0, const Vector2d2& poly_1);
typedef  Vector2d (*CGAL_2D_Nearest_Point_Polygon_C1)(const Vector2d& v, const Vector2d1& poly);
typedef  void (*CGAL_2D_Nearest_Point_Polygon_C2)(const Vector2d& v, const Vector2d1& poly, Vector2d& p, double& min_d);
typedef  Vector2d (*CGAL_2D_Nearest_Point_Polygons)(const Vector2d& v, const Vector2d2& polys);
typedef  void (*CGAL_2d_Polygon_Boundingbox)(const Vector2d1& ps, Vector2d& min_corner, Vector2d& max_corner);
typedef  double (*CGAL_2D_Polygon_Area)(const Vector2d1& py);
typedef  Vector2d (*CGAL_2D_Polygon_Inside_Point_C1)(const Vector2d1& poly);
typedef  bool (*CGAL_2D_Polygon_Inside_Point_C2)(const Vector2d2& polys, Vector2d& inner_vec);
typedef  void (*CGAL_2D_Polygon_One_Offsets)(const Vector2d1& poly, const double& d, Vector2d2& offset_polys);
typedef  void (*CGAL_2D_Polygons_One_Offsets)(const Vector2d2& polys, const double& d, Vector2d2& offset_polys);
typedef  bool (*CGAL_2D_Polygons_Simple)(const Vector2d2& poly);
typedef  bool (*CGAL_2D_Polygon_Simple)(const Vector2d1& poly);
typedef  bool (*CGAL_2D_Polygon_Simple_Inter)(const Vector2d1& poly);
typedef  void (*CGAL_2D_Convex_Hulls)(const Vector2d1& vec, Vector2d1& hull_points);
typedef  void (*CGAL_2D_OBB_Box)(const Vector2d1& vec, Vector2d& center, Vector2d& axis_0, Vector2d& axis_1, double& entent_0, double& entent_1);
typedef  void (*CGAL_Image_Grid_Decomposition_C1)(std::vector<std::vector<int>>& image, std::vector<std::vector<double>>& boundary_xs, std::vector<std::vector<double>>& boundary_ys);
typedef  void (*CGAL_Image_Grid_Decomposition_Conservative_C1)(std::vector<std::vector<int>>& image, std::vector<std::vector<double>>& boundary_xs, std::vector<std::vector<double>>& boundary_ys);
typedef  void (*CGAL_Image_Grid_Decomposition_C2)(std::vector<std::vector<int>>& image, Vector2d2& boundaries);
typedef  void (*CGAL_Image_Grid_Decomposition_Conservative_C2)(std::vector<std::vector<int>>& image, Vector2d2& boundaries);
//implementation in "threeD.cpp"
//####################################################################################
typedef  double (*CGAL_3D_Distance_Point_Segment)(const Vector3d& p, const Vector3d& s_s, const Vector3d& s_e);
typedef  void (*CGAL_3D_Plane_Fitting)(const Vector3d1& points, Vector3d& plane_p, Vector3d& plane_n);
typedef  void (*CGAL_3D_Plane_Point_Projection)(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d& p, Vector3d& result);
typedef  void (*CGAL_3D_Plane_Points_Projection)(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d1& points, Vector3d1& project_points);
typedef  void (*CGAL_3D_Plane_3D_to_2D_Point)(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d& point_3d, Vector2d& result);
typedef  void (*CGAL_3D_Plane_2D_to_3D_Point)(const Vector3d& plane_p, const Vector3d& plane_n, const Vector2d& points_2d, Vector3d& result);
typedef  void (*CGAL_3D_Plane_3D_to_2D_Points)(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d1& points_3d, Vector2d1& points_2d);
typedef  void (*CGAL_3D_Plane_3Ds_to_2Ds_Points)(const Vector3d& plane_p, const Vector3d& plane_n, const Vector3d2& points_3d, Vector2d2& points_2d);
typedef  void (*CGAL_3D_Plane_2D_to_3D_Points)(const Vector3d& plane_p, const Vector3d& plane_n, const Vector2d1& points_2d, Vector3d1& points_3d);
typedef  Vector3d (*CGAL_3D_Projection_Point_Segment)(const Vector3d& p, const Vector3d& s_s, const Vector3d& s_e);
typedef  double (*CGAL_3D_Distance_Point_Point)(const Vector3d& v0, const Vector3d& v1);
typedef  double (*CGAL_3D_Distance_Point_Polygon)(const Vector3d1& py, const Vector3d& p);
typedef  void (*CGAL_2D_Polygon_Triangulation_C1)(const Vector2d2& polys, Vector1i2& faces);
typedef  std::vector<std::vector<int>> (*CGAL_2D_Polygon_Triangulation_C2)(const Vector2d2& polys);
typedef  std::vector<std::vector<int>> (*CGAL_2D_Polygon_Triangulation_C3)(const Vector2d1& poly);
typedef  double (*CGAL_3D_Distance_Point_Line)(const Vector3d& p, const Vector3d& l_s, const Vector3d& l_e);
typedef  Vector3d (*CGAL_3D_Projection_Point_Line)(const Vector3d& p, const Vector3d& l_s, const Vector3d& l_e);
typedef  double (*CGAL_3D_Distance_Segment_Segment)(const Vector3d& s_0_s, const Vector3d& s_0_e, const Vector3d& s_1_s, const Vector3d& s_1_e);
typedef  double (*CGAL_3D_Distance_Point_Plane)(const Vector3d& v, const Vector3d& plane_p, const Vector3d& plane_n);
typedef  bool (*CGAL_3D_Intersection_Segment_Line)(const Vector3d& s_s, const Vector3d& s_e, const Vector3d& l_s, const Vector3d& l_e, Vector3d& inter);
typedef  bool (*CGAL_3D_Intersection_Segment_Segment)(const Vector3d& s_0_s, const Vector3d& s_0_e, const Vector3d& s_1_s, const Vector3d& s_1_e, Vector3d& iter);
typedef  bool (*CGAL_3D_Intersection_Segment_Plane)(const Vector3d& s_s, const Vector3d& s_e, const Vector3d& plane_p, const Vector3d& plane_n, Vector3d& inter);
typedef  bool (*CGAL_3D_Intersection_Line_Plane)(const Vector3d& l_s, const Vector3d& l_e, const Vector3d& plane_p, const Vector3d& plane_n, Vector3d& inter);
typedef  Vector3d (*CGAL_3D_Projection_Point_Plane_C1)(const Vector3d& p, const Vector3d& plane_p, const Vector3d& plane_n);
typedef  Vector3d (*CGAL_3D_Projection_Point_Plane_C2)(const Vector3d& p, const Vector3d& plane_p_0, const Vector3d& plane_p_1, const Vector3d& plane_p_2);
typedef  Vector2d (*CGAL_3D_Projection_3D_Point_Plane_2D_C1)(const Vector3d& p, const Vector3d& plane_p, const Vector3d& plane_n);
typedef  Vector2d (*CGAL_3D_Projection_3D_Point_Plane_2D_C2)(const Vector3d& p, const Vector3d& plane_p_0, const Vector3d& plane_p_1, const Vector3d& plane_p_2);
typedef  void (*CGAL_3D_Plane_ABCD)(const Vector3d& plane_p, const Vector3d& plane_n, double& a, double& b, double& c, double& d);
typedef  Vector3d (*CGAL_3D_Plane_Base_1)(const Vector3d& plane_p, const Vector3d& plane_n);
typedef  Vector3d (*CGAL_Face_Normal)(const Vector3d& source, const Vector3d& tri_0, const Vector3d& tri_1, const Vector3d& tri_2, Vector3d& normal_0, Vector3d& normal_1, Vector3d& normal_2);
//implementation in "mesh.cpp"
//####################################################################################
typedef  void (*CGAL_Remesh_Surface_by_Adding_Feature)(const Vector3d1& feature, const Vector1i1& face_ids, const Vector3d1& vecs, const Vector1i1& face_id_0, const Vector1i1& face_id_1, const Vector1i1& face_id_2, Vector1i1& igl_cutting_0_edges, Vector1i1& igl_cutting_1_edges, Vector3d1& igl_cutting_points, Vector1i2& cutting_faces);
typedef  void (*CGAL_Mesh_Edges)(const char* path);
typedef  bool (*CGAL_3D_Intersection_Sphere_Ray)(const double& center_x, const double& center_y, const double& center_z, const double& radius, const double& ray_origin_x, const double& ray_origin_y, const double& ray_origin_z, const double& ray_direction_x, const double& ray_direction_y, const double& ray_direction_z, std::vector<double>& i_x, std::vector<double>& i_y, std::vector<double>& i_z);
typedef  bool (*CGAL_3D_Intersection_Ray_Triangle)(const Vector3d& p, const Vector3d& n, const Vector3d& p0, const Vector3d& p1, const Vector3d& p2);
typedef  bool (*CGAL_3D_Intersection_Ray_Mesh)(const Vector3d& p, const Vector3d& n, const char* path);
typedef  bool (*CGAL_3D_Intersection_Segment_Mesh)(const Vector3d& s, const Vector3d& e, const char* path);
typedef  void (*CGAL_3D_Intersection_Segments_Mesh)(const Vector3d1& ss, const Vector3d1& ee, const char* path, Vector1b1& inters);
typedef  void (*CGAL_3D_Intersection_Polygons_Mesh)(const Vector3d2& polygons, const char* path, Vector1b1& inters);
//check whether there is a polygon intersected with the input mesh
typedef  bool (*CGAL_3D_Intersection_Polygons_Mesh_Bool)(const Vector3d2& polygons, const char* path);
typedef  void (*CGAL_3D_Intersection_Rays_Mesh_Vector3d)(const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector3d1& inters);
//test each group directions (nes[i]) for each point in ps
typedef  void (*CGAL_3D_Intersection_Rays_Mesh_C1_Bool)(const Vector3d1& ps, const Vector3d2& nes, const char* path, Vector1b2& inters);
//test all directions (ns) for each point in ps
typedef  void (*CGAL_3D_Intersection_Rays_Mesh_C2_Bool)(const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector1b2& inters);
typedef  void (*CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d)(const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector1d2& inters);
typedef  void (*CGAL_3D_Points_Inside_Triangles_C1_Bool)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points, std::vector<bool>& insides);
typedef  void (*CGAL_3D_Points_Inside_Triangles_C2_Bool)(const char* path, const Vector3d1& points, std::vector<bool>& insides);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  void (*CGAL_3D_Mesh_Dart_Sampling_C1)(const char* outside_path, const double& d, Vector3d1& sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  void (*CGAL_3D_Mesh_Dart_Sampling_C2)(const char* outside_path, const char* inside_path, const double& d, Vector3d1& sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  void (*CGAL_3D_Mesh_Regular_Sampling_C1)(const char* outside_path, const double& d, Vector3d1& sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  void (*CGAL_3D_Mesh_Regular_Sampling_C2)(const char* outside_path, const char* inside_path, const double& d, Vector3d1& sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  void (*CGAL_3D_Cube_Surface_Sampling_C1)(const double& cube_size, const double& d, Vector3d2& sampling_points, VectorPI2& neighbors, const bool& compute_neighbors);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  void (*CGAL_3D_Cube_Surface_Sampling_C2)(const double& cube_size, const double& d, Vector3d2& sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
typedef  void (*CGAL_3D_Cube_Surface_Sampling_C3)(const double& cube_size, const double& d, Vector3d2& sampling_points, VectorPI2& neighbors);
//with neighboring
typedef  void (*CGAL_3D_Mesh_Regular_Sampling_C3)(const char* outside_path, const char* inside_path, const double& d, Vector3d1& sampling_points, VectorPI1& neighbors);
typedef  double (*CGAL_3D_Distance_Point_Triangle)(const Vector3d& p, const Vector3d& t_0, const Vector3d& t_1, const Vector3d& t_2);
typedef  double (*CGAL_3D_Distance_Point_Triangles)(const Vector3d& p, const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2);
typedef  Vector3d (*CGAL_3D_Nearest_Point_Triangles)(const Vector3d& p, const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2);
typedef  void (*CGAL_3D_Distance_Point_Mesh)(const char* path, const Vector3d1& query_points, std::vector<double>& distances);
typedef  void (*CGAL_3D_Neareast_Point_Mesh)(const char* path, const Vector3d1& ves, Vector3d1& ners);
typedef  void  (*CGAL_3D_Mesh_Near_Triangles)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points, const double& d, std::vector<std::vector<int>>& triangles);
typedef  void (*CGAL_3D_Points_inside_Triangles_C1)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points, std::vector<bool>& insides);
typedef  void (*CGAL_3D_Points_inside_Triangles_C2)(const char* path, const Vector3d1& points, std::vector<bool>& insides);
typedef  void (*CGAL_Mesh_Subdivision)(const char* in_path, const char* sub, const int& step, const char* out_path);
typedef  void (*CGAL_Mesh_Loop_Subdivision_One_Step)(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2);
typedef  void (*CGAL_3D_Mesh_Curvature_C1)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, std::vector<double>& max_curs, std::vector<double>& min_curs);
typedef  void (*CGAL_3D_Mesh_Curvature_C2)(const Vector3d1& vecs, const std::vector<std::vector<int>>& face_ids, std::vector<double>& max_curs, std::vector<double>& min_curs);
typedef  void (*CGAL_3D_Mesh_Curvature_C3)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions);
typedef  void (*CGAL_3D_Mesh_Curvature_C4)(const Vector3d1& vecs, const std::vector<std::vector<int>>& face_ids, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions);
typedef  void (*CGAL_3D_Mesh_Curvature_C5)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions, Vector3d1& normals);
typedef  void (*CGAL_3D_Mesh_Curvature_C6)(const Vector3d1& vecs, const std::vector<std::vector<int>>& face_ids, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions, Vector3d1& normals);
typedef  void (*CGAL_3D_Triangle_Mesh_Boundary_C1)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, std::vector<bool>& bools);
typedef  void (*CGAL_3D_Triangle_Mesh_Boundary_C2)(const char* path, std::vector<bool>& bools);
typedef  void (*CGAL_3D_Connecting_Segments_C1)(Vector2d2& segments, Vector2d2& lines);
typedef  void (*CGAL_3D_Connecting_Segments_C2)(Vector3d2& segments, Vector3d2& lines);
typedef  void (*CGAL_3D_Triangle_Mesh_Boundary_C3)(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, Vector3d2& boundaries);
typedef  void (*CGAL_3D_Triangle_Mesh_Boundary_C4)(Vector3d1& vecs, std::vector<std::vector<int>>& face_ids, Vector3d2& boundaries);
typedef  void (*CGAL_3D_Triangle_Mesh_Boundary_C5)(const char* path, Vector3d2& boundaries);
typedef  void (*CGAL_Mesh_Laplace_Smooth_C1)(const char* in_path, const char* out_path, const int laplace_nb);
typedef  void (*CGAL_3D_Triangle_Mesh_Vecs_Neighbors)(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<std::vector<int>>& neighs);
typedef  void (*CGAL_Mesh_Laplace_Smooth_C2)(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, const int laplace_nb);
typedef  void (*CGAL_3D_Triangle_Mesh_Vecs_Faces)(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<std::vector<int>>& surface_vectices_to_face);
typedef  void (*CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges)(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<std::vector<std::vector<int>>>& surface_vectices_to_neighbor_edges);
typedef  void (*CGAL_Mesh_Laplace_Smooth_by_Curvature)(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, double& low_curvature);
typedef  void (*CGAL_Mesh_Loop_Subdivision_Own_Version)(const char* in_path, const int& step, const char* out_path, const int& laplace_nb);
typedef  void (*CGAL_Rotation_Obj)(const char* path, const double& angle, const Vector3d& axis);
typedef  void (*CGAL_Slicer_Mesh)(const char* path, const Vector3d& plane_normal, const std::vector<double>& plane_d, Vector3d3& offsetses, Vector3d2& offsets);
typedef  void (*CGAL_Shortest_Geodesic_Path_C1)(const char* path, Vector3d1& xyzs);
typedef  void (*CGAL_Shortest_Geodesic_Path_C3)(const char* path, Vector3d source, Vector3d target, Vector3d1& output);
typedef  void (*CGAL_Shortest_Geodesic_Path_C4)(const char* path, Vector3d1 sources, Vector3d1 targets, Vector3d2& xyzes);
typedef  double (*CGAL_Geodesic_Distance)(const char* path, const Vector3d& source, const Vector3d& target);
typedef  Vector3d1 (*CGAL_Project_Points_Onto_Surface_C1)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points);
typedef  Vector3d1 (*CGAL_Project_Points_Onto_Surface_C2)(const char* path, const Vector3d1& points);
typedef  void (*CGAL_3D_Triangel_Mesh_Most_Inside_Point)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, Vector3d& inside);
typedef  double (*CGAL_3D_One_Triangle_Area)(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2);
typedef  double (*CGAL_3D_Triangle_Mesh_Area)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2);
typedef  void (*CGAL_3D_Convex_Hulls_C1)(const Vector3d1& vec, Vector3d1& hull_points);
typedef  void (*CGAL_3D_Convex_Hulls_C2)(const Vector3d1& vec, Vector3d1& hull_points, std::vector<int>& hulls_surface_0, std::vector<int>& hulls_surface_1, std::vector<int>& hulls_surface_2);
typedef  void (*CGAL_3D_Convex_Hulls_C3)(const Vector3d1& vec, Vector3d1& hull_points, Vector3d1& plane_p, Vector3d1& plane_n);
typedef  void (*CGAL_3D_Convex_Hulls_C4)(const Vector3d1& vec, Vector3d1& hull_points, std::vector<int>& hulls_surface_0, std::vector<int>& hulls_surface_1, std::vector<int>& hulls_surface_2, Vector3d1& plane_p, Vector3d1& plane_n);
typedef  void (*CGAL_Mesh_Field_Query_C1)(const char* path, const Vector3d1& gradients, const Vector3d1& input_points, Vector3d1& points_gradients);
typedef  void (*CGAL_Mesh_Field_Query_C2)(const char* path, const std::vector<double>& gradient_values, const Vector3d1& input_points, std::vector<double>& points_gradient_values);
typedef  void (*CGAL_Mesh_Field_Query_C3)(const char* path, const std::vector<double>& gradient_values, const Vector3d2& input_point_es, std::vector<std::vector<double>>& points_gradient_value_es);
typedef  void (*CGAL_Curvature_Mesh)(const char* path, const Vector3d1& input_points, std::vector<double>& max_curs, std::vector<double>& min_curs, Vector3d1& max_curs_directions, Vector3d1& min_curs_directions);
typedef  void (*CGAL_Normal_Mesh_C1)(const char* path, const Vector3d1& mesh_points, Vector3d1& mesh_normals);
typedef  void (*CGAL_Normal_Mesh_C2)(const char* path, const Vector3d2& mesh_pointses, Vector3d2& mesh_normalses);
typedef  void (*CGAL_3D_Mesh_Normal_C1)(const Vector3d1& ps, const std::vector<std::vector<int>>& face_ids, Vector3d1& normals);
typedef  void (*CGAL_3D_Mesh_Normal_C2)(const Vector3d1& ps, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, Vector3d1& normals);
typedef  Vector3d (*CGAL_3D_Mesh_Center_C1)(const Vector3d2& ps);
typedef  Vector3d (*CGAL_3D_Mesh_Center_C2)(const Vector3d1& ps);
typedef  void (*CGAL_3D_Mesh_Boundingbox_C1)(const Vector3d2& ps, Vector3d& min_corner, Vector3d& max_corner);
typedef  void (*CGAL_3D_Mesh_Boundingbox_C2)(const Vector3d1& ps, Vector3d& min_corner, Vector3d& max_corner);
typedef  void (*CGAL_Surface_Decomposition)(const char* path, std::vector<double>& face_sdf, int& regions_nb, std::vector<int>& face_segments);
typedef  void (*CGAL_3D_Mesh_Gradient)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const std::vector<double>& psd, Vector3d1& vecs_gradients, Vector3d1& faces_gradients);
typedef  void (*CGAL_Intergral_Curvature)(const Vector2d1& input_points, const int& sampling_points_nb, const double& radius, const double& thresholder, Vector2d1& output_points, std::vector<double>& output_rates);
typedef  bool (*CGAL_3D_Mesh_Extract_Isoline)(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const std::vector<double>& psd, const double& d, Vector3d2& isolines);
typedef  void (*CGAL_BSplineCurveFit)(const Vector3d1& samples, Vector3d1& output);
typedef  void (*CGAL_Cut_Surface)(const Vector3d1& boundary, const Vector3d& inside_point, const char* full_path, char* output_path);
typedef  void (*CGAL_Cut_Surface_by_Multi_Boundaries)(const Vector3d2& multi_boundary, const Vector3d& inside_point, const char* full_path, char* output_path);
/////////////////////////////////////////////////////////////
//


class CGALPL
{
	public:
	CGALPL()
	{
		hModule = Functs::LoadHMODULE("ppgl.dll");
		CGAL_Test_PGL_C = (CGAL_Test_PGL)GetProcAddress(hModule, "CGAL_Test_PGL");
		//implementation in "io.cpp"
		//####################################################################################
		CGAL_Vector_Base_C = (CGAL_Vector_Base)GetProcAddress(hModule, "CGAL_Vector_Base");
		//implementation in "twoD.cpp"
		//####################################################################################
		CGAL_2D_Distance_Point_Point_C = (CGAL_2D_Distance_Point_Point)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Point");
		CGAL_2D_Distance_Point_Line_C = (CGAL_2D_Distance_Point_Line)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Line");
		CGAL_2D_Distance_Point_Segment_C = (CGAL_2D_Distance_Point_Segment)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Segment");
		CGAL_2D_Distance_Segment_Segment_C = (CGAL_2D_Distance_Segment_Segment)GetProcAddress(hModule, "CGAL_2D_Distance_Segment_Segment");
		CGAL_2D_Location_Point_Polygon_C = (CGAL_2D_Location_Point_Polygon)GetProcAddress(hModule, "CGAL_2D_Location_Point_Polygon");
		CGAL_2D_Location_Points_Polygon_C = (CGAL_2D_Location_Points_Polygon)GetProcAddress(hModule, "CGAL_2D_Location_Points_Polygon");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_2D_Polygon_Dart_Sampling_C = (CGAL_2D_Polygon_Dart_Sampling)GetProcAddress(hModule, "CGAL_2D_Polygon_Dart_Sampling");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_2D_Polygon_Regular_Sampling_C1_C = (CGAL_2D_Polygon_Regular_Sampling_C1)GetProcAddress(hModule, "CGAL_2D_Polygon_Regular_Sampling_C1");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_2D_Polygon_Regular_Sampling_C2_C = (CGAL_2D_Polygon_Regular_Sampling_C2)GetProcAddress(hModule, "CGAL_2D_Polygon_Regular_Sampling_C2");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_2D_Polygon_Regular_Sampling_C3_C = (CGAL_2D_Polygon_Regular_Sampling_C3)GetProcAddress(hModule, "CGAL_2D_Polygon_Regular_Sampling_C3");
		CGAL_2D_Square_Regular_Sampling_C1_C = (CGAL_2D_Square_Regular_Sampling_C1)GetProcAddress(hModule, "CGAL_2D_Square_Regular_Sampling_C1");
		CGAL_2D_Square_Regular_Sampling_C2_C = (CGAL_2D_Square_Regular_Sampling_C2)GetProcAddress(hModule, "CGAL_2D_Square_Regular_Sampling_C2");
		CGAL_2D_Square_Regular_Sampling_C3_C = (CGAL_2D_Square_Regular_Sampling_C3)GetProcAddress(hModule, "CGAL_2D_Square_Regular_Sampling_C3");
		CGAL_2D_Distance_Point_Polygon_C = (CGAL_2D_Distance_Point_Polygon)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Polygon");
		CGAL_2D_Distance_Point_Polygons_C = (CGAL_2D_Distance_Point_Polygons)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Polygons");
		CGAL_2D_Intersection_Segment_Segment_C = (CGAL_2D_Intersection_Segment_Segment)GetProcAddress(hModule, "CGAL_2D_Intersection_Segment_Segment");
		CGAL_2D_Intersection_Line_Line_C = (CGAL_2D_Intersection_Line_Line)GetProcAddress(hModule, "CGAL_2D_Intersection_Line_Line");
		CGAL_2D_Intersection_Segment_Line_C = (CGAL_2D_Intersection_Segment_Line)GetProcAddress(hModule, "CGAL_2D_Intersection_Segment_Line");
		CGAL_2D_Intersection_Segment_Polygon_C = (CGAL_2D_Intersection_Segment_Polygon)GetProcAddress(hModule, "CGAL_2D_Intersection_Segment_Polygon");
		CGAL_2D_Polygon_Is_Clockwise_Oriented_C = (CGAL_2D_Polygon_Is_Clockwise_Oriented)GetProcAddress(hModule, "CGAL_2D_Polygon_Is_Clockwise_Oriented");
		CGAL_2D_Two_Polygons_Union_C = (CGAL_2D_Two_Polygons_Union)GetProcAddress(hModule, "CGAL_2D_Two_Polygons_Union");
		CGAL_2D_Two_Polygons_Intersection_C = (CGAL_2D_Two_Polygons_Intersection)GetProcAddress(hModule, "CGAL_2D_Two_Polygons_Intersection");
		CGAL_Decompose_Polyline_C = (CGAL_Decompose_Polyline)GetProcAddress(hModule, "CGAL_Decompose_Polyline");
		CGAL_Identify_Polycut_Extend_C = (CGAL_Identify_Polycut_Extend)GetProcAddress(hModule, "CGAL_Identify_Polycut_Extend");
		CGAL_Identify_Polycut_NotExtend_C = (CGAL_Identify_Polycut_NotExtend)GetProcAddress(hModule, "CGAL_Identify_Polycut_NotExtend");
		CGAL_Identify_Polycut_C = (CGAL_Identify_Polycut)GetProcAddress(hModule, "CGAL_Identify_Polycut");
		CGAL_Construct_InOutSide_Polygon_C = (CGAL_Construct_InOutSide_Polygon)GetProcAddress(hModule, "CGAL_Construct_InOutSide_Polygon");
		CGAL_2D_Intersection_Ray_Segment_C = (CGAL_2D_Intersection_Ray_Segment)GetProcAddress(hModule, "CGAL_2D_Intersection_Ray_Segment");
		CGAL_Get_Angle_Kerf_Offset_Tan_C = (CGAL_Get_Angle_Kerf_Offset_Tan)GetProcAddress(hModule, "CGAL_Get_Angle_Kerf_Offset_Tan");
		CGAL_2D_Projection_Point_Segment_C = (CGAL_2D_Projection_Point_Segment)GetProcAddress(hModule, "CGAL_2D_Projection_Point_Segment");
		CGAL_2D_Detect_Polygon_Inside_C1_C = (CGAL_2D_Detect_Polygon_Inside_C1)GetProcAddress(hModule, "CGAL_2D_Detect_Polygon_Inside_C1");
		CGAL_2D_Detect_Polygon_Inside_C2_C = (CGAL_2D_Detect_Polygon_Inside_C2)GetProcAddress(hModule, "CGAL_2D_Detect_Polygon_Inside_C2");
		CGAL_2D_Detect_Polygon_Inside_C3_C = (CGAL_2D_Detect_Polygon_Inside_C3)GetProcAddress(hModule, "CGAL_2D_Detect_Polygon_Inside_C3");
		CGAL_2D_Detect_Polygon_Inside_C4_C = (CGAL_2D_Detect_Polygon_Inside_C4)GetProcAddress(hModule, "CGAL_2D_Detect_Polygon_Inside_C4");
		CGAL_2D_Detect_Polygon_Inside_C5_C = (CGAL_2D_Detect_Polygon_Inside_C5)GetProcAddress(hModule, "CGAL_2D_Detect_Polygon_Inside_C5");
		CGAL_2D_Distance_Polygon_Polygon_C = (CGAL_2D_Distance_Polygon_Polygon)GetProcAddress(hModule, "CGAL_2D_Distance_Polygon_Polygon");
		CGAL_2D_Distance_Polygons_Polygons_C = (CGAL_2D_Distance_Polygons_Polygons)GetProcAddress(hModule, "CGAL_2D_Distance_Polygons_Polygons");
		CGAL_2D_Nearest_Point_Polygon_C1_C = (CGAL_2D_Nearest_Point_Polygon_C1)GetProcAddress(hModule, "CGAL_2D_Nearest_Point_Polygon_C1");
		CGAL_2D_Nearest_Point_Polygon_C2_C = (CGAL_2D_Nearest_Point_Polygon_C2)GetProcAddress(hModule, "CGAL_2D_Nearest_Point_Polygon_C2");
		CGAL_2D_Nearest_Point_Polygons_C = (CGAL_2D_Nearest_Point_Polygons)GetProcAddress(hModule, "CGAL_2D_Nearest_Point_Polygons");
		CGAL_2d_Polygon_Boundingbox_C = (CGAL_2d_Polygon_Boundingbox)GetProcAddress(hModule, "CGAL_2d_Polygon_Boundingbox");
		CGAL_2D_Polygon_Area_C = (CGAL_2D_Polygon_Area)GetProcAddress(hModule, "CGAL_2D_Polygon_Area");
		CGAL_2D_Polygon_Inside_Point_C1_C = (CGAL_2D_Polygon_Inside_Point_C1)GetProcAddress(hModule, "CGAL_2D_Polygon_Inside_Point_C1");
		CGAL_2D_Polygon_Inside_Point_C2_C = (CGAL_2D_Polygon_Inside_Point_C2)GetProcAddress(hModule, "CGAL_2D_Polygon_Inside_Point_C2");
		CGAL_2D_Polygon_One_Offsets_C = (CGAL_2D_Polygon_One_Offsets)GetProcAddress(hModule, "CGAL_2D_Polygon_One_Offsets");
		CGAL_2D_Polygons_One_Offsets_C = (CGAL_2D_Polygons_One_Offsets)GetProcAddress(hModule, "CGAL_2D_Polygons_One_Offsets");
		CGAL_2D_Polygons_Simple_C = (CGAL_2D_Polygons_Simple)GetProcAddress(hModule, "CGAL_2D_Polygons_Simple");
		CGAL_2D_Polygon_Simple_C = (CGAL_2D_Polygon_Simple)GetProcAddress(hModule, "CGAL_2D_Polygon_Simple");
		CGAL_2D_Polygon_Simple_Inter_C = (CGAL_2D_Polygon_Simple_Inter)GetProcAddress(hModule, "CGAL_2D_Polygon_Simple_Inter");
		CGAL_2D_Convex_Hulls_C = (CGAL_2D_Convex_Hulls)GetProcAddress(hModule, "CGAL_2D_Convex_Hulls");
		CGAL_2D_OBB_Box_C = (CGAL_2D_OBB_Box)GetProcAddress(hModule, "CGAL_2D_OBB_Box");
		CGAL_Image_Grid_Decomposition_C1_C = (CGAL_Image_Grid_Decomposition_C1)GetProcAddress(hModule, "CGAL_Image_Grid_Decomposition_C1");
		CGAL_Image_Grid_Decomposition_Conservative_C1_C = (CGAL_Image_Grid_Decomposition_Conservative_C1)GetProcAddress(hModule, "CGAL_Image_Grid_Decomposition_Conservative_C1");
		CGAL_Image_Grid_Decomposition_C2_C = (CGAL_Image_Grid_Decomposition_C2)GetProcAddress(hModule, "CGAL_Image_Grid_Decomposition_C2");
		CGAL_Image_Grid_Decomposition_Conservative_C2_C = (CGAL_Image_Grid_Decomposition_Conservative_C2)GetProcAddress(hModule, "CGAL_Image_Grid_Decomposition_Conservative_C2");
		//implementation in "threeD.cpp"
		//####################################################################################
		CGAL_3D_Distance_Point_Segment_C = (CGAL_3D_Distance_Point_Segment)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Segment");
		CGAL_3D_Plane_Fitting_C = (CGAL_3D_Plane_Fitting)GetProcAddress(hModule, "CGAL_3D_Plane_Fitting");
		CGAL_3D_Plane_Point_Projection_C = (CGAL_3D_Plane_Point_Projection)GetProcAddress(hModule, "CGAL_3D_Plane_Point_Projection");
		CGAL_3D_Plane_Points_Projection_C = (CGAL_3D_Plane_Points_Projection)GetProcAddress(hModule, "CGAL_3D_Plane_Points_Projection");
		CGAL_3D_Plane_3D_to_2D_Point_C = (CGAL_3D_Plane_3D_to_2D_Point)GetProcAddress(hModule, "CGAL_3D_Plane_3D_to_2D_Point");
		CGAL_3D_Plane_2D_to_3D_Point_C = (CGAL_3D_Plane_2D_to_3D_Point)GetProcAddress(hModule, "CGAL_3D_Plane_2D_to_3D_Point");
		CGAL_3D_Plane_3D_to_2D_Points_C = (CGAL_3D_Plane_3D_to_2D_Points)GetProcAddress(hModule, "CGAL_3D_Plane_3D_to_2D_Points");
		CGAL_3D_Plane_3Ds_to_2Ds_Points_C = (CGAL_3D_Plane_3Ds_to_2Ds_Points)GetProcAddress(hModule, "CGAL_3D_Plane_3Ds_to_2Ds_Points");
		CGAL_3D_Plane_2D_to_3D_Points_C = (CGAL_3D_Plane_2D_to_3D_Points)GetProcAddress(hModule, "CGAL_3D_Plane_2D_to_3D_Points");
		CGAL_3D_Projection_Point_Segment_C = (CGAL_3D_Projection_Point_Segment)GetProcAddress(hModule, "CGAL_3D_Projection_Point_Segment");
		CGAL_3D_Distance_Point_Point_C = (CGAL_3D_Distance_Point_Point)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Point");
		CGAL_3D_Distance_Point_Polygon_C = (CGAL_3D_Distance_Point_Polygon)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Polygon");
		CGAL_2D_Polygon_Triangulation_C1_C = (CGAL_2D_Polygon_Triangulation_C1)GetProcAddress(hModule, "CGAL_2D_Polygon_Triangulation_C1");
		CGAL_2D_Polygon_Triangulation_C2_C = (CGAL_2D_Polygon_Triangulation_C2)GetProcAddress(hModule, "CGAL_2D_Polygon_Triangulation_C2");
		CGAL_2D_Polygon_Triangulation_C3_C = (CGAL_2D_Polygon_Triangulation_C3)GetProcAddress(hModule, "CGAL_2D_Polygon_Triangulation_C3");
		CGAL_3D_Distance_Point_Line_C = (CGAL_3D_Distance_Point_Line)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Line");
		CGAL_3D_Projection_Point_Line_C = (CGAL_3D_Projection_Point_Line)GetProcAddress(hModule, "CGAL_3D_Projection_Point_Line");
		CGAL_3D_Distance_Segment_Segment_C = (CGAL_3D_Distance_Segment_Segment)GetProcAddress(hModule, "CGAL_3D_Distance_Segment_Segment");
		CGAL_3D_Distance_Point_Plane_C = (CGAL_3D_Distance_Point_Plane)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Plane");
		CGAL_3D_Intersection_Segment_Line_C = (CGAL_3D_Intersection_Segment_Line)GetProcAddress(hModule, "CGAL_3D_Intersection_Segment_Line");
		CGAL_3D_Intersection_Segment_Segment_C = (CGAL_3D_Intersection_Segment_Segment)GetProcAddress(hModule, "CGAL_3D_Intersection_Segment_Segment");
		CGAL_3D_Intersection_Segment_Plane_C = (CGAL_3D_Intersection_Segment_Plane)GetProcAddress(hModule, "CGAL_3D_Intersection_Segment_Plane");
		CGAL_3D_Intersection_Line_Plane_C = (CGAL_3D_Intersection_Line_Plane)GetProcAddress(hModule, "CGAL_3D_Intersection_Line_Plane");
		CGAL_3D_Projection_Point_Plane_C1_C = (CGAL_3D_Projection_Point_Plane_C1)GetProcAddress(hModule, "CGAL_3D_Projection_Point_Plane_C1");
		CGAL_3D_Projection_Point_Plane_C2_C = (CGAL_3D_Projection_Point_Plane_C2)GetProcAddress(hModule, "CGAL_3D_Projection_Point_Plane_C2");
		CGAL_3D_Projection_3D_Point_Plane_2D_C1_C = (CGAL_3D_Projection_3D_Point_Plane_2D_C1)GetProcAddress(hModule, "CGAL_3D_Projection_3D_Point_Plane_2D_C1");
		CGAL_3D_Projection_3D_Point_Plane_2D_C2_C = (CGAL_3D_Projection_3D_Point_Plane_2D_C2)GetProcAddress(hModule, "CGAL_3D_Projection_3D_Point_Plane_2D_C2");
		CGAL_3D_Plane_ABCD_C = (CGAL_3D_Plane_ABCD)GetProcAddress(hModule, "CGAL_3D_Plane_ABCD");
		CGAL_3D_Plane_Base_1_C = (CGAL_3D_Plane_Base_1)GetProcAddress(hModule, "CGAL_3D_Plane_Base_1");
		CGAL_Face_Normal_C = (CGAL_Face_Normal)GetProcAddress(hModule, "CGAL_Face_Normal");
		//implementation in "mesh.cpp"
		//####################################################################################
		CGAL_Remesh_Surface_by_Adding_Feature_C = (CGAL_Remesh_Surface_by_Adding_Feature)GetProcAddress(hModule, "CGAL_Remesh_Surface_by_Adding_Feature");
		CGAL_Mesh_Edges_C = (CGAL_Mesh_Edges)GetProcAddress(hModule, "CGAL_Mesh_Edges");
		CGAL_3D_Intersection_Sphere_Ray_C = (CGAL_3D_Intersection_Sphere_Ray)GetProcAddress(hModule, "CGAL_3D_Intersection_Sphere_Ray");
		CGAL_3D_Intersection_Ray_Triangle_C = (CGAL_3D_Intersection_Ray_Triangle)GetProcAddress(hModule, "CGAL_3D_Intersection_Ray_Triangle");
		CGAL_3D_Intersection_Ray_Mesh_C = (CGAL_3D_Intersection_Ray_Mesh)GetProcAddress(hModule, "CGAL_3D_Intersection_Ray_Mesh");
		CGAL_3D_Intersection_Segment_Mesh_C = (CGAL_3D_Intersection_Segment_Mesh)GetProcAddress(hModule, "CGAL_3D_Intersection_Segment_Mesh");
		CGAL_3D_Intersection_Segments_Mesh_C = (CGAL_3D_Intersection_Segments_Mesh)GetProcAddress(hModule, "CGAL_3D_Intersection_Segments_Mesh");
		CGAL_3D_Intersection_Polygons_Mesh_C = (CGAL_3D_Intersection_Polygons_Mesh)GetProcAddress(hModule, "CGAL_3D_Intersection_Polygons_Mesh");
		//check whether there is a polygon intersected with the input mesh
		CGAL_3D_Intersection_Polygons_Mesh_Bool_C = (CGAL_3D_Intersection_Polygons_Mesh_Bool)GetProcAddress(hModule, "CGAL_3D_Intersection_Polygons_Mesh_Bool");
		CGAL_3D_Intersection_Rays_Mesh_Vector3d_C = (CGAL_3D_Intersection_Rays_Mesh_Vector3d)GetProcAddress(hModule, "CGAL_3D_Intersection_Rays_Mesh_Vector3d");
		//test each group directions (nes[i]) for each point in ps
		CGAL_3D_Intersection_Rays_Mesh_C1_Bool_C = (CGAL_3D_Intersection_Rays_Mesh_C1_Bool)GetProcAddress(hModule, "CGAL_3D_Intersection_Rays_Mesh_C1_Bool");
		//test all directions (ns) for each point in ps
		CGAL_3D_Intersection_Rays_Mesh_C2_Bool_C = (CGAL_3D_Intersection_Rays_Mesh_C2_Bool)GetProcAddress(hModule, "CGAL_3D_Intersection_Rays_Mesh_C2_Bool");
		CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d_C = (CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d)GetProcAddress(hModule, "CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d");
		CGAL_3D_Points_Inside_Triangles_C1_Bool_C = (CGAL_3D_Points_Inside_Triangles_C1_Bool)GetProcAddress(hModule, "CGAL_3D_Points_Inside_Triangles_C1_Bool");
		CGAL_3D_Points_Inside_Triangles_C2_Bool_C = (CGAL_3D_Points_Inside_Triangles_C2_Bool)GetProcAddress(hModule, "CGAL_3D_Points_Inside_Triangles_C2_Bool");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_3D_Mesh_Dart_Sampling_C1_C = (CGAL_3D_Mesh_Dart_Sampling_C1)GetProcAddress(hModule, "CGAL_3D_Mesh_Dart_Sampling_C1");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_3D_Mesh_Dart_Sampling_C2_C = (CGAL_3D_Mesh_Dart_Sampling_C2)GetProcAddress(hModule, "CGAL_3D_Mesh_Dart_Sampling_C2");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_3D_Mesh_Regular_Sampling_C1_C = (CGAL_3D_Mesh_Regular_Sampling_C1)GetProcAddress(hModule, "CGAL_3D_Mesh_Regular_Sampling_C1");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_3D_Mesh_Regular_Sampling_C2_C = (CGAL_3D_Mesh_Regular_Sampling_C2)GetProcAddress(hModule, "CGAL_3D_Mesh_Regular_Sampling_C2");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_3D_Cube_Surface_Sampling_C1_C = (CGAL_3D_Cube_Surface_Sampling_C1)GetProcAddress(hModule, "CGAL_3D_Cube_Surface_Sampling_C1");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_3D_Cube_Surface_Sampling_C2_C = (CGAL_3D_Cube_Surface_Sampling_C2)GetProcAddress(hModule, "CGAL_3D_Cube_Surface_Sampling_C2");
		//d: percentage value of the length of the diagonal of the bounding box.
		CGAL_3D_Cube_Surface_Sampling_C3_C = (CGAL_3D_Cube_Surface_Sampling_C3)GetProcAddress(hModule, "CGAL_3D_Cube_Surface_Sampling_C3");
		//with neighboring
		CGAL_3D_Mesh_Regular_Sampling_C3_C = (CGAL_3D_Mesh_Regular_Sampling_C3)GetProcAddress(hModule, "CGAL_3D_Mesh_Regular_Sampling_C3");
		CGAL_3D_Distance_Point_Triangle_C = (CGAL_3D_Distance_Point_Triangle)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Triangle");
		CGAL_3D_Distance_Point_Triangles_C = (CGAL_3D_Distance_Point_Triangles)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Triangles");
		CGAL_3D_Nearest_Point_Triangles_C = (CGAL_3D_Nearest_Point_Triangles)GetProcAddress(hModule, "CGAL_3D_Nearest_Point_Triangles");
		CGAL_3D_Distance_Point_Mesh_C = (CGAL_3D_Distance_Point_Mesh)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Mesh");
		CGAL_3D_Neareast_Point_Mesh_C = (CGAL_3D_Neareast_Point_Mesh)GetProcAddress(hModule, "CGAL_3D_Neareast_Point_Mesh");
		CGAL_3D_Mesh_Near_Triangles_C = (CGAL_3D_Mesh_Near_Triangles)GetProcAddress(hModule, "CGAL_3D_Mesh_Near_Triangles");
		CGAL_3D_Points_inside_Triangles_C1_C = (CGAL_3D_Points_inside_Triangles_C1)GetProcAddress(hModule, "CGAL_3D_Points_inside_Triangles_C1");
		CGAL_3D_Points_inside_Triangles_C2_C = (CGAL_3D_Points_inside_Triangles_C2)GetProcAddress(hModule, "CGAL_3D_Points_inside_Triangles_C2");
		CGAL_Mesh_Subdivision_C = (CGAL_Mesh_Subdivision)GetProcAddress(hModule, "CGAL_Mesh_Subdivision");
		CGAL_Mesh_Loop_Subdivision_One_Step_C = (CGAL_Mesh_Loop_Subdivision_One_Step)GetProcAddress(hModule, "CGAL_Mesh_Loop_Subdivision_One_Step");
		CGAL_3D_Mesh_Curvature_C1_C = (CGAL_3D_Mesh_Curvature_C1)GetProcAddress(hModule, "CGAL_3D_Mesh_Curvature_C1");
		CGAL_3D_Mesh_Curvature_C2_C = (CGAL_3D_Mesh_Curvature_C2)GetProcAddress(hModule, "CGAL_3D_Mesh_Curvature_C2");
		CGAL_3D_Mesh_Curvature_C3_C = (CGAL_3D_Mesh_Curvature_C3)GetProcAddress(hModule, "CGAL_3D_Mesh_Curvature_C3");
		CGAL_3D_Mesh_Curvature_C4_C = (CGAL_3D_Mesh_Curvature_C4)GetProcAddress(hModule, "CGAL_3D_Mesh_Curvature_C4");
		CGAL_3D_Mesh_Curvature_C5_C = (CGAL_3D_Mesh_Curvature_C5)GetProcAddress(hModule, "CGAL_3D_Mesh_Curvature_C5");
		CGAL_3D_Mesh_Curvature_C6_C = (CGAL_3D_Mesh_Curvature_C6)GetProcAddress(hModule, "CGAL_3D_Mesh_Curvature_C6");
		CGAL_3D_Triangle_Mesh_Boundary_C1_C = (CGAL_3D_Triangle_Mesh_Boundary_C1)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Boundary_C1");
		CGAL_3D_Triangle_Mesh_Boundary_C2_C = (CGAL_3D_Triangle_Mesh_Boundary_C2)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Boundary_C2");
		CGAL_3D_Connecting_Segments_C1_C = (CGAL_3D_Connecting_Segments_C1)GetProcAddress(hModule, "CGAL_3D_Connecting_Segments_C1");
		CGAL_3D_Connecting_Segments_C2_C = (CGAL_3D_Connecting_Segments_C2)GetProcAddress(hModule, "CGAL_3D_Connecting_Segments_C2");
		CGAL_3D_Triangle_Mesh_Boundary_C3_C = (CGAL_3D_Triangle_Mesh_Boundary_C3)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Boundary_C3");
		CGAL_3D_Triangle_Mesh_Boundary_C4_C = (CGAL_3D_Triangle_Mesh_Boundary_C4)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Boundary_C4");
		CGAL_3D_Triangle_Mesh_Boundary_C5_C = (CGAL_3D_Triangle_Mesh_Boundary_C5)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Boundary_C5");
		CGAL_Mesh_Laplace_Smooth_C1_C = (CGAL_Mesh_Laplace_Smooth_C1)GetProcAddress(hModule, "CGAL_Mesh_Laplace_Smooth_C1");
		CGAL_3D_Triangle_Mesh_Vecs_Neighbors_C = (CGAL_3D_Triangle_Mesh_Vecs_Neighbors)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Vecs_Neighbors");
		CGAL_Mesh_Laplace_Smooth_C2_C = (CGAL_Mesh_Laplace_Smooth_C2)GetProcAddress(hModule, "CGAL_Mesh_Laplace_Smooth_C2");
		CGAL_3D_Triangle_Mesh_Vecs_Faces_C = (CGAL_3D_Triangle_Mesh_Vecs_Faces)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Vecs_Faces");
		CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges_C = (CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges");
		CGAL_Mesh_Laplace_Smooth_by_Curvature_C = (CGAL_Mesh_Laplace_Smooth_by_Curvature)GetProcAddress(hModule, "CGAL_Mesh_Laplace_Smooth_by_Curvature");
		CGAL_Mesh_Loop_Subdivision_Own_Version_C = (CGAL_Mesh_Loop_Subdivision_Own_Version)GetProcAddress(hModule, "CGAL_Mesh_Loop_Subdivision_Own_Version");
		CGAL_Rotation_Obj_C = (CGAL_Rotation_Obj)GetProcAddress(hModule, "CGAL_Rotation_Obj");
		CGAL_Slicer_Mesh_C = (CGAL_Slicer_Mesh)GetProcAddress(hModule, "CGAL_Slicer_Mesh");
		CGAL_Shortest_Geodesic_Path_C1_C = (CGAL_Shortest_Geodesic_Path_C1)GetProcAddress(hModule, "CGAL_Shortest_Geodesic_Path_C1");
		CGAL_Shortest_Geodesic_Path_C3_C = (CGAL_Shortest_Geodesic_Path_C3)GetProcAddress(hModule, "CGAL_Shortest_Geodesic_Path_C3");
		CGAL_Shortest_Geodesic_Path_C4_C = (CGAL_Shortest_Geodesic_Path_C4)GetProcAddress(hModule, "CGAL_Shortest_Geodesic_Path_C4");
		CGAL_Geodesic_Distance_C = (CGAL_Geodesic_Distance)GetProcAddress(hModule, "CGAL_Geodesic_Distance");
		CGAL_Project_Points_Onto_Surface_C1_C = (CGAL_Project_Points_Onto_Surface_C1)GetProcAddress(hModule, "CGAL_Project_Points_Onto_Surface_C1");
		CGAL_Project_Points_Onto_Surface_C2_C = (CGAL_Project_Points_Onto_Surface_C2)GetProcAddress(hModule, "CGAL_Project_Points_Onto_Surface_C2");
		CGAL_3D_Triangel_Mesh_Most_Inside_Point_C = (CGAL_3D_Triangel_Mesh_Most_Inside_Point)GetProcAddress(hModule, "CGAL_3D_Triangel_Mesh_Most_Inside_Point");
		CGAL_3D_One_Triangle_Area_C = (CGAL_3D_One_Triangle_Area)GetProcAddress(hModule, "CGAL_3D_One_Triangle_Area");
		CGAL_3D_Triangle_Mesh_Area_C = (CGAL_3D_Triangle_Mesh_Area)GetProcAddress(hModule, "CGAL_3D_Triangle_Mesh_Area");
		CGAL_3D_Convex_Hulls_C1_C = (CGAL_3D_Convex_Hulls_C1)GetProcAddress(hModule, "CGAL_3D_Convex_Hulls_C1");
		CGAL_3D_Convex_Hulls_C2_C = (CGAL_3D_Convex_Hulls_C2)GetProcAddress(hModule, "CGAL_3D_Convex_Hulls_C2");
		CGAL_3D_Convex_Hulls_C3_C = (CGAL_3D_Convex_Hulls_C3)GetProcAddress(hModule, "CGAL_3D_Convex_Hulls_C3");
		CGAL_3D_Convex_Hulls_C4_C = (CGAL_3D_Convex_Hulls_C4)GetProcAddress(hModule, "CGAL_3D_Convex_Hulls_C4");
		CGAL_Mesh_Field_Query_C1_C = (CGAL_Mesh_Field_Query_C1)GetProcAddress(hModule, "CGAL_Mesh_Field_Query_C1");
		CGAL_Mesh_Field_Query_C2_C = (CGAL_Mesh_Field_Query_C2)GetProcAddress(hModule, "CGAL_Mesh_Field_Query_C2");
		CGAL_Mesh_Field_Query_C3_C = (CGAL_Mesh_Field_Query_C3)GetProcAddress(hModule, "CGAL_Mesh_Field_Query_C3");
		CGAL_Curvature_Mesh_C = (CGAL_Curvature_Mesh)GetProcAddress(hModule, "CGAL_Curvature_Mesh");
		CGAL_Normal_Mesh_C1_C = (CGAL_Normal_Mesh_C1)GetProcAddress(hModule, "CGAL_Normal_Mesh_C1");
		CGAL_Normal_Mesh_C2_C = (CGAL_Normal_Mesh_C2)GetProcAddress(hModule, "CGAL_Normal_Mesh_C2");
		CGAL_3D_Mesh_Normal_C1_C = (CGAL_3D_Mesh_Normal_C1)GetProcAddress(hModule, "CGAL_3D_Mesh_Normal_C1");
		CGAL_3D_Mesh_Normal_C2_C = (CGAL_3D_Mesh_Normal_C2)GetProcAddress(hModule, "CGAL_3D_Mesh_Normal_C2");
		CGAL_3D_Mesh_Center_C1_C = (CGAL_3D_Mesh_Center_C1)GetProcAddress(hModule, "CGAL_3D_Mesh_Center_C1");
		CGAL_3D_Mesh_Center_C2_C = (CGAL_3D_Mesh_Center_C2)GetProcAddress(hModule, "CGAL_3D_Mesh_Center_C2");
		CGAL_3D_Mesh_Boundingbox_C1_C = (CGAL_3D_Mesh_Boundingbox_C1)GetProcAddress(hModule, "CGAL_3D_Mesh_Boundingbox_C1");
		CGAL_3D_Mesh_Boundingbox_C2_C = (CGAL_3D_Mesh_Boundingbox_C2)GetProcAddress(hModule, "CGAL_3D_Mesh_Boundingbox_C2");
		CGAL_Surface_Decomposition_C = (CGAL_Surface_Decomposition)GetProcAddress(hModule, "CGAL_Surface_Decomposition");
		CGAL_3D_Mesh_Gradient_C = (CGAL_3D_Mesh_Gradient)GetProcAddress(hModule, "CGAL_3D_Mesh_Gradient");
		CGAL_Intergral_Curvature_C = (CGAL_Intergral_Curvature)GetProcAddress(hModule, "CGAL_Intergral_Curvature");
		CGAL_3D_Mesh_Extract_Isoline_C = (CGAL_3D_Mesh_Extract_Isoline)GetProcAddress(hModule, "CGAL_3D_Mesh_Extract_Isoline");
		CGAL_BSplineCurveFit_C = (CGAL_BSplineCurveFit)GetProcAddress(hModule, "CGAL_BSplineCurveFit");
		CGAL_Cut_Surface_C = (CGAL_Cut_Surface)GetProcAddress(hModule, "CGAL_Cut_Surface");
		CGAL_Cut_Surface_by_Multi_Boundaries_C = (CGAL_Cut_Surface_by_Multi_Boundaries)GetProcAddress(hModule, "CGAL_Cut_Surface_by_Multi_Boundaries");
		/////////////////////////////////////////////////////////////
		//
	};

	static CGALPL& Inst()
	{
		static CGALPL instance;
		return instance;
	};

	HMODULE hModule;
	CGAL_Test_PGL CGAL_Test_PGL_C;
	//implementation in "io.cpp"
	//####################################################################################
	CGAL_Vector_Base CGAL_Vector_Base_C;
	//implementation in "twoD.cpp"
	//####################################################################################
	CGAL_2D_Distance_Point_Point CGAL_2D_Distance_Point_Point_C;
	CGAL_2D_Distance_Point_Line CGAL_2D_Distance_Point_Line_C;
	CGAL_2D_Distance_Point_Segment CGAL_2D_Distance_Point_Segment_C;
	CGAL_2D_Distance_Segment_Segment CGAL_2D_Distance_Segment_Segment_C;
	CGAL_2D_Location_Point_Polygon CGAL_2D_Location_Point_Polygon_C;
	CGAL_2D_Location_Points_Polygon CGAL_2D_Location_Points_Polygon_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_2D_Polygon_Dart_Sampling CGAL_2D_Polygon_Dart_Sampling_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_2D_Polygon_Regular_Sampling_C1 CGAL_2D_Polygon_Regular_Sampling_C1_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_2D_Polygon_Regular_Sampling_C2 CGAL_2D_Polygon_Regular_Sampling_C2_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_2D_Polygon_Regular_Sampling_C3 CGAL_2D_Polygon_Regular_Sampling_C3_C;
	CGAL_2D_Square_Regular_Sampling_C1 CGAL_2D_Square_Regular_Sampling_C1_C;
	CGAL_2D_Square_Regular_Sampling_C2 CGAL_2D_Square_Regular_Sampling_C2_C;
	CGAL_2D_Square_Regular_Sampling_C3 CGAL_2D_Square_Regular_Sampling_C3_C;
	CGAL_2D_Distance_Point_Polygon CGAL_2D_Distance_Point_Polygon_C;
	CGAL_2D_Distance_Point_Polygons CGAL_2D_Distance_Point_Polygons_C;
	CGAL_2D_Intersection_Segment_Segment CGAL_2D_Intersection_Segment_Segment_C;
	CGAL_2D_Intersection_Line_Line CGAL_2D_Intersection_Line_Line_C;
	CGAL_2D_Intersection_Segment_Line CGAL_2D_Intersection_Segment_Line_C;
	CGAL_2D_Intersection_Segment_Polygon CGAL_2D_Intersection_Segment_Polygon_C;
	CGAL_2D_Polygon_Is_Clockwise_Oriented CGAL_2D_Polygon_Is_Clockwise_Oriented_C;
	CGAL_2D_Two_Polygons_Union CGAL_2D_Two_Polygons_Union_C;
	CGAL_2D_Two_Polygons_Intersection CGAL_2D_Two_Polygons_Intersection_C;
	CGAL_Decompose_Polyline CGAL_Decompose_Polyline_C;
	CGAL_Identify_Polycut_Extend CGAL_Identify_Polycut_Extend_C;
	CGAL_Identify_Polycut_NotExtend CGAL_Identify_Polycut_NotExtend_C;
	CGAL_Identify_Polycut CGAL_Identify_Polycut_C;
	CGAL_Construct_InOutSide_Polygon CGAL_Construct_InOutSide_Polygon_C;
	CGAL_2D_Intersection_Ray_Segment CGAL_2D_Intersection_Ray_Segment_C;
	CGAL_Get_Angle_Kerf_Offset_Tan CGAL_Get_Angle_Kerf_Offset_Tan_C;
	CGAL_2D_Projection_Point_Segment CGAL_2D_Projection_Point_Segment_C;
	CGAL_2D_Detect_Polygon_Inside_C1 CGAL_2D_Detect_Polygon_Inside_C1_C;
	CGAL_2D_Detect_Polygon_Inside_C2 CGAL_2D_Detect_Polygon_Inside_C2_C;
	CGAL_2D_Detect_Polygon_Inside_C3 CGAL_2D_Detect_Polygon_Inside_C3_C;
	CGAL_2D_Detect_Polygon_Inside_C4 CGAL_2D_Detect_Polygon_Inside_C4_C;
	CGAL_2D_Detect_Polygon_Inside_C5 CGAL_2D_Detect_Polygon_Inside_C5_C;
	CGAL_2D_Distance_Polygon_Polygon CGAL_2D_Distance_Polygon_Polygon_C;
	CGAL_2D_Distance_Polygons_Polygons CGAL_2D_Distance_Polygons_Polygons_C;
	CGAL_2D_Nearest_Point_Polygon_C1 CGAL_2D_Nearest_Point_Polygon_C1_C;
	CGAL_2D_Nearest_Point_Polygon_C2 CGAL_2D_Nearest_Point_Polygon_C2_C;
	CGAL_2D_Nearest_Point_Polygons CGAL_2D_Nearest_Point_Polygons_C;
	CGAL_2d_Polygon_Boundingbox CGAL_2d_Polygon_Boundingbox_C;
	CGAL_2D_Polygon_Area CGAL_2D_Polygon_Area_C;
	CGAL_2D_Polygon_Inside_Point_C1 CGAL_2D_Polygon_Inside_Point_C1_C;
	CGAL_2D_Polygon_Inside_Point_C2 CGAL_2D_Polygon_Inside_Point_C2_C;
	CGAL_2D_Polygon_One_Offsets CGAL_2D_Polygon_One_Offsets_C;
	CGAL_2D_Polygons_One_Offsets CGAL_2D_Polygons_One_Offsets_C;
	CGAL_2D_Polygons_Simple CGAL_2D_Polygons_Simple_C;
	CGAL_2D_Polygon_Simple CGAL_2D_Polygon_Simple_C;
	CGAL_2D_Polygon_Simple_Inter CGAL_2D_Polygon_Simple_Inter_C;
	CGAL_2D_Convex_Hulls CGAL_2D_Convex_Hulls_C;
	CGAL_2D_OBB_Box CGAL_2D_OBB_Box_C;
	CGAL_Image_Grid_Decomposition_C1 CGAL_Image_Grid_Decomposition_C1_C;
	CGAL_Image_Grid_Decomposition_Conservative_C1 CGAL_Image_Grid_Decomposition_Conservative_C1_C;
	CGAL_Image_Grid_Decomposition_C2 CGAL_Image_Grid_Decomposition_C2_C;
	CGAL_Image_Grid_Decomposition_Conservative_C2 CGAL_Image_Grid_Decomposition_Conservative_C2_C;
	//implementation in "threeD.cpp"
	//####################################################################################
	CGAL_3D_Distance_Point_Segment CGAL_3D_Distance_Point_Segment_C;
	CGAL_3D_Plane_Fitting CGAL_3D_Plane_Fitting_C;
	CGAL_3D_Plane_Point_Projection CGAL_3D_Plane_Point_Projection_C;
	CGAL_3D_Plane_Points_Projection CGAL_3D_Plane_Points_Projection_C;
	CGAL_3D_Plane_3D_to_2D_Point CGAL_3D_Plane_3D_to_2D_Point_C;
	CGAL_3D_Plane_2D_to_3D_Point CGAL_3D_Plane_2D_to_3D_Point_C;
	CGAL_3D_Plane_3D_to_2D_Points CGAL_3D_Plane_3D_to_2D_Points_C;
	CGAL_3D_Plane_3Ds_to_2Ds_Points CGAL_3D_Plane_3Ds_to_2Ds_Points_C;
	CGAL_3D_Plane_2D_to_3D_Points CGAL_3D_Plane_2D_to_3D_Points_C;
	CGAL_3D_Projection_Point_Segment CGAL_3D_Projection_Point_Segment_C;
	CGAL_3D_Distance_Point_Point CGAL_3D_Distance_Point_Point_C;
	CGAL_3D_Distance_Point_Polygon CGAL_3D_Distance_Point_Polygon_C;
	CGAL_2D_Polygon_Triangulation_C1 CGAL_2D_Polygon_Triangulation_C1_C;
	CGAL_2D_Polygon_Triangulation_C2 CGAL_2D_Polygon_Triangulation_C2_C;
	CGAL_2D_Polygon_Triangulation_C3 CGAL_2D_Polygon_Triangulation_C3_C;
	CGAL_3D_Distance_Point_Line CGAL_3D_Distance_Point_Line_C;
	CGAL_3D_Projection_Point_Line CGAL_3D_Projection_Point_Line_C;
	CGAL_3D_Distance_Segment_Segment CGAL_3D_Distance_Segment_Segment_C;
	CGAL_3D_Distance_Point_Plane CGAL_3D_Distance_Point_Plane_C;
	CGAL_3D_Intersection_Segment_Line CGAL_3D_Intersection_Segment_Line_C;
	CGAL_3D_Intersection_Segment_Segment CGAL_3D_Intersection_Segment_Segment_C;
	CGAL_3D_Intersection_Segment_Plane CGAL_3D_Intersection_Segment_Plane_C;
	CGAL_3D_Intersection_Line_Plane CGAL_3D_Intersection_Line_Plane_C;
	CGAL_3D_Projection_Point_Plane_C1 CGAL_3D_Projection_Point_Plane_C1_C;
	CGAL_3D_Projection_Point_Plane_C2 CGAL_3D_Projection_Point_Plane_C2_C;
	CGAL_3D_Projection_3D_Point_Plane_2D_C1 CGAL_3D_Projection_3D_Point_Plane_2D_C1_C;
	CGAL_3D_Projection_3D_Point_Plane_2D_C2 CGAL_3D_Projection_3D_Point_Plane_2D_C2_C;
	CGAL_3D_Plane_ABCD CGAL_3D_Plane_ABCD_C;
	CGAL_3D_Plane_Base_1 CGAL_3D_Plane_Base_1_C;
	CGAL_Face_Normal CGAL_Face_Normal_C;
	//implementation in "mesh.cpp"
	//####################################################################################
	CGAL_Remesh_Surface_by_Adding_Feature CGAL_Remesh_Surface_by_Adding_Feature_C;
	CGAL_Mesh_Edges CGAL_Mesh_Edges_C;
	CGAL_3D_Intersection_Sphere_Ray CGAL_3D_Intersection_Sphere_Ray_C;
	CGAL_3D_Intersection_Ray_Triangle CGAL_3D_Intersection_Ray_Triangle_C;
	CGAL_3D_Intersection_Ray_Mesh CGAL_3D_Intersection_Ray_Mesh_C;
	CGAL_3D_Intersection_Segment_Mesh CGAL_3D_Intersection_Segment_Mesh_C;
	CGAL_3D_Intersection_Segments_Mesh CGAL_3D_Intersection_Segments_Mesh_C;
	CGAL_3D_Intersection_Polygons_Mesh CGAL_3D_Intersection_Polygons_Mesh_C;
	//check whether there is a polygon intersected with the input mesh
	CGAL_3D_Intersection_Polygons_Mesh_Bool CGAL_3D_Intersection_Polygons_Mesh_Bool_C;
	CGAL_3D_Intersection_Rays_Mesh_Vector3d CGAL_3D_Intersection_Rays_Mesh_Vector3d_C;
	//test each group directions (nes[i]) for each point in ps
	CGAL_3D_Intersection_Rays_Mesh_C1_Bool CGAL_3D_Intersection_Rays_Mesh_C1_Bool_C;
	//test all directions (ns) for each point in ps
	CGAL_3D_Intersection_Rays_Mesh_C2_Bool CGAL_3D_Intersection_Rays_Mesh_C2_Bool_C;
	CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d_C;
	CGAL_3D_Points_Inside_Triangles_C1_Bool CGAL_3D_Points_Inside_Triangles_C1_Bool_C;
	CGAL_3D_Points_Inside_Triangles_C2_Bool CGAL_3D_Points_Inside_Triangles_C2_Bool_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_3D_Mesh_Dart_Sampling_C1 CGAL_3D_Mesh_Dart_Sampling_C1_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_3D_Mesh_Dart_Sampling_C2 CGAL_3D_Mesh_Dart_Sampling_C2_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_3D_Mesh_Regular_Sampling_C1 CGAL_3D_Mesh_Regular_Sampling_C1_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_3D_Mesh_Regular_Sampling_C2 CGAL_3D_Mesh_Regular_Sampling_C2_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_3D_Cube_Surface_Sampling_C1 CGAL_3D_Cube_Surface_Sampling_C1_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_3D_Cube_Surface_Sampling_C2 CGAL_3D_Cube_Surface_Sampling_C2_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	CGAL_3D_Cube_Surface_Sampling_C3 CGAL_3D_Cube_Surface_Sampling_C3_C;
	//with neighboring
	CGAL_3D_Mesh_Regular_Sampling_C3 CGAL_3D_Mesh_Regular_Sampling_C3_C;
	CGAL_3D_Distance_Point_Triangle CGAL_3D_Distance_Point_Triangle_C;
	CGAL_3D_Distance_Point_Triangles CGAL_3D_Distance_Point_Triangles_C;
	CGAL_3D_Nearest_Point_Triangles CGAL_3D_Nearest_Point_Triangles_C;
	CGAL_3D_Distance_Point_Mesh CGAL_3D_Distance_Point_Mesh_C;
	CGAL_3D_Neareast_Point_Mesh CGAL_3D_Neareast_Point_Mesh_C;
	CGAL_3D_Mesh_Near_Triangles CGAL_3D_Mesh_Near_Triangles_C;
	CGAL_3D_Points_inside_Triangles_C1 CGAL_3D_Points_inside_Triangles_C1_C;
	CGAL_3D_Points_inside_Triangles_C2 CGAL_3D_Points_inside_Triangles_C2_C;
	CGAL_Mesh_Subdivision CGAL_Mesh_Subdivision_C;
	CGAL_Mesh_Loop_Subdivision_One_Step CGAL_Mesh_Loop_Subdivision_One_Step_C;
	CGAL_3D_Mesh_Curvature_C1 CGAL_3D_Mesh_Curvature_C1_C;
	CGAL_3D_Mesh_Curvature_C2 CGAL_3D_Mesh_Curvature_C2_C;
	CGAL_3D_Mesh_Curvature_C3 CGAL_3D_Mesh_Curvature_C3_C;
	CGAL_3D_Mesh_Curvature_C4 CGAL_3D_Mesh_Curvature_C4_C;
	CGAL_3D_Mesh_Curvature_C5 CGAL_3D_Mesh_Curvature_C5_C;
	CGAL_3D_Mesh_Curvature_C6 CGAL_3D_Mesh_Curvature_C6_C;
	CGAL_3D_Triangle_Mesh_Boundary_C1 CGAL_3D_Triangle_Mesh_Boundary_C1_C;
	CGAL_3D_Triangle_Mesh_Boundary_C2 CGAL_3D_Triangle_Mesh_Boundary_C2_C;
	CGAL_3D_Connecting_Segments_C1 CGAL_3D_Connecting_Segments_C1_C;
	CGAL_3D_Connecting_Segments_C2 CGAL_3D_Connecting_Segments_C2_C;
	CGAL_3D_Triangle_Mesh_Boundary_C3 CGAL_3D_Triangle_Mesh_Boundary_C3_C;
	CGAL_3D_Triangle_Mesh_Boundary_C4 CGAL_3D_Triangle_Mesh_Boundary_C4_C;
	CGAL_3D_Triangle_Mesh_Boundary_C5 CGAL_3D_Triangle_Mesh_Boundary_C5_C;
	CGAL_Mesh_Laplace_Smooth_C1 CGAL_Mesh_Laplace_Smooth_C1_C;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors CGAL_3D_Triangle_Mesh_Vecs_Neighbors_C;
	CGAL_Mesh_Laplace_Smooth_C2 CGAL_Mesh_Laplace_Smooth_C2_C;
	CGAL_3D_Triangle_Mesh_Vecs_Faces CGAL_3D_Triangle_Mesh_Vecs_Faces_C;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges_C;
	CGAL_Mesh_Laplace_Smooth_by_Curvature CGAL_Mesh_Laplace_Smooth_by_Curvature_C;
	CGAL_Mesh_Loop_Subdivision_Own_Version CGAL_Mesh_Loop_Subdivision_Own_Version_C;
	CGAL_Rotation_Obj CGAL_Rotation_Obj_C;
	CGAL_Slicer_Mesh CGAL_Slicer_Mesh_C;
	CGAL_Shortest_Geodesic_Path_C1 CGAL_Shortest_Geodesic_Path_C1_C;
	CGAL_Shortest_Geodesic_Path_C3 CGAL_Shortest_Geodesic_Path_C3_C;
	CGAL_Shortest_Geodesic_Path_C4 CGAL_Shortest_Geodesic_Path_C4_C;
	CGAL_Geodesic_Distance CGAL_Geodesic_Distance_C;
	CGAL_Project_Points_Onto_Surface_C1 CGAL_Project_Points_Onto_Surface_C1_C;
	CGAL_Project_Points_Onto_Surface_C2 CGAL_Project_Points_Onto_Surface_C2_C;
	CGAL_3D_Triangel_Mesh_Most_Inside_Point CGAL_3D_Triangel_Mesh_Most_Inside_Point_C;
	CGAL_3D_One_Triangle_Area CGAL_3D_One_Triangle_Area_C;
	CGAL_3D_Triangle_Mesh_Area CGAL_3D_Triangle_Mesh_Area_C;
	CGAL_3D_Convex_Hulls_C1 CGAL_3D_Convex_Hulls_C1_C;
	CGAL_3D_Convex_Hulls_C2 CGAL_3D_Convex_Hulls_C2_C;
	CGAL_3D_Convex_Hulls_C3 CGAL_3D_Convex_Hulls_C3_C;
	CGAL_3D_Convex_Hulls_C4 CGAL_3D_Convex_Hulls_C4_C;
	CGAL_Mesh_Field_Query_C1 CGAL_Mesh_Field_Query_C1_C;
	CGAL_Mesh_Field_Query_C2 CGAL_Mesh_Field_Query_C2_C;
	CGAL_Mesh_Field_Query_C3 CGAL_Mesh_Field_Query_C3_C;
	CGAL_Curvature_Mesh CGAL_Curvature_Mesh_C;
	CGAL_Normal_Mesh_C1 CGAL_Normal_Mesh_C1_C;
	CGAL_Normal_Mesh_C2 CGAL_Normal_Mesh_C2_C;
	CGAL_3D_Mesh_Normal_C1 CGAL_3D_Mesh_Normal_C1_C;
	CGAL_3D_Mesh_Normal_C2 CGAL_3D_Mesh_Normal_C2_C;
	CGAL_3D_Mesh_Center_C1 CGAL_3D_Mesh_Center_C1_C;
	CGAL_3D_Mesh_Center_C2 CGAL_3D_Mesh_Center_C2_C;
	CGAL_3D_Mesh_Boundingbox_C1 CGAL_3D_Mesh_Boundingbox_C1_C;
	CGAL_3D_Mesh_Boundingbox_C2 CGAL_3D_Mesh_Boundingbox_C2_C;
	CGAL_Surface_Decomposition CGAL_Surface_Decomposition_C;
	CGAL_3D_Mesh_Gradient CGAL_3D_Mesh_Gradient_C;
	CGAL_Intergral_Curvature CGAL_Intergral_Curvature_C;
	CGAL_3D_Mesh_Extract_Isoline CGAL_3D_Mesh_Extract_Isoline_C;
	CGAL_BSplineCurveFit CGAL_BSplineCurveFit_C;
	CGAL_Cut_Surface CGAL_Cut_Surface_C;
	CGAL_Cut_Surface_by_Multi_Boundaries CGAL_Cut_Surface_by_Multi_Boundaries_C;
	/////////////////////////////////////////////////////////////
	//
};
#define PL() CGALPL::Inst()
}
#endif
