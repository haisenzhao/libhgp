#ifndef LIBHGP_ONCE
#define LIBHGP_ONCE
#pragma once
#include <iostream>
#include <vector>
#include <windows.h>
#include <direct.h>
#include <tchar.h>
#include <string>
#include "glm/glm.hpp"
using namespace std;
namespace libhgp {

template <typename datum>
using Vector1 = std::vector<datum>;

template <typename datum>
using Vector2 = std::vector<std::vector<datum>>;

template <typename datum>
using Vector3 = std::vector<std::vector<std::vector<datum>>>;

typedef glm::highp_dvec2 Vector2d;
typedef glm::highp_dvec3 Vector3d;
typedef glm::highp_ivec2 Vector2i;
typedef glm::highp_ivec3 Vector3i;

typedef Vector1<Vector2d> Vector2d1;
typedef Vector2<Vector2d> Vector2d2;
typedef Vector3<Vector2d> Vector2d3;

typedef Vector1<Vector3d> Vector3d1;
typedef Vector2<Vector3d> Vector3d2;
typedef Vector3<Vector3d> Vector3d3;

typedef Vector1<bool> Vector1b1;
typedef Vector2<bool> Vector1b2;
typedef Vector3<bool> Vector1b3;

typedef Vector1<int> Vector1i1;
typedef Vector2<int> Vector1i2;
typedef Vector3<int> Vector1i3;

typedef Vector1<double> Vector1d1;
typedef Vector2<double> Vector1d2;
typedef Vector3<double> Vector1d3;

typedef Vector1<std::string> VectorStr1;
typedef Vector2<std::string> VectorStr2;
typedef Vector3<std::string> VectorStr3;

typedef Vector1<Vector2i> Vector2i1;
typedef Vector2<Vector2i> Vector2i2;
typedef Vector3<Vector2i> Vector2i3;

typedef Vector1<Vector3i> Vector3i1;
typedef Vector2<Vector3i> Vector3i2;
typedef Vector3<Vector3i> Vector3i3;

typedef Vector1<std::pair<int, int>> VectorPI1;
typedef Vector2<std::pair<int, int>> VectorPI2;
typedef Vector3<std::pair<int, int>> VectorPI3;

typedef Vector1<std::pair<bool, bool>> VectorPB1;
typedef Vector2<std::pair<bool, bool>> VectorPB2;
typedef Vector3<std::pair<bool, bool>> VectorPB3;

typedef std::tuple<int, int, int> TI3;
typedef Vector1<std::tuple<int, int, int>> VectorTI3;

static HMODULE LoadHMODULE(const string& dll_path)
{
	struct stat buffer;
	if (!(stat(dll_path.c_str(), &buffer) == 0))
	{
		char tmp[256];
		if (_getcwd(tmp, 256)) {};
		std::string root_path = std::string(tmp);

		std::string str;
		str += "The dll does not exist: " + dll_path + "; \n";
		str += "The current running directory : " + root_path + "; \n";
		str += "Please gurrentee the dll is in the right place; \n";
		std::cerr << str << std::endl;
	}

	HMODULE hModule = LoadLibraryA(LPCSTR(dll_path.c_str()));
	if (!hModule)
	{
		DWORD dw = GetLastError(); // returns 0xc1 (193)
		std::cerr << "LoadLibrary failed with error code " + std::to_string(dw) << std::endl;
	}
	else
		std::cerr << "LoadLibrary success\n";
	return hModule;
};
//implementation in "io.cpp"
//####################################################################################
//Usage: Get a vector which is perpenticular to the input vector;
//Input: input vector "n"; 
//Output: output vector "r";
typedef void (*HGP_Vector_Base)(const Vector3d & n, Vector3d & r);
typedef void (*HGP_Test_PGL)(const Vector3d & n, const char* str, const char* char_);
//implementation in "hgp2d.cpp"
//####################################################################################
typedef double (*HGP_2D_Distance_Point_Point)(const Vector2d & p_0, const Vector2d & p_1);
typedef double (*HGP_2D_Distance_Point_Line)(const Vector2d & v, const Vector2d & l_0, const Vector2d & l_1);
typedef double (*HGP_2D_Distance_Point_Segment)(const Vector2d & v, const Vector2d & s_0, const Vector2d & s_1);
typedef double (*HGP_2D_Distance_Segment_Segment)(const Vector2d & s_0, const Vector2d & s_1, const Vector2d & e_0, const Vector2d & e_1);
typedef bool (*HGP_2D_Location_Point_Polygon)(const Vector2d & p, const Vector2d1 & py);
typedef bool (*HGP_2D_Location_Points_Polygon)(const Vector2d1 & ps, const Vector2d1 & py);
//d: percentage value of the length of the diagonal of the bounding box.
typedef void (*HGP_2D_Polygon_Dart_Sampling)(const Vector2d1 & py, const double& d, Vector2d1 & sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
typedef Vector2d1 (*HGP_2D_Polygon_Regular_Sampling_C1)(const Vector2d1 & py, const double& d);
//d: percentage value of the length of the diagonal of the bounding box.
typedef Vector2d1 (*HGP_2D_Polygon_Regular_Sampling_C2)(const Vector2d1 & py, const double& d, VectorPI1 & neighbors);
//d: percentage value of the length of the diagonal of the bounding box.
typedef Vector2d1 (*HGP_2D_Polygon_Regular_Sampling_C3)(const Vector2d1 & py, const double& d, VectorPI1 & neighbors, const bool& compute_neighbors);
typedef Vector2d1 (*HGP_2D_Square_Regular_Sampling_C1)(const double& d);
typedef Vector2d1 (*HGP_2D_Square_Regular_Sampling_C2)(const double& d, VectorPI1 & neighbors);
typedef Vector2d1 (*HGP_2D_Square_Regular_Sampling_C3)(const double& d, VectorPI1 & neighbors, const bool& compute_neighbors);
typedef double (*HGP_2D_Distance_Point_Polygon)(const Vector2d & p, const Vector2d1 & py);
typedef double (*HGP_2D_Distance_Point_Polygons)(const Vector2d & p, const Vector2d2 & pys);
typedef bool (*HGP_2D_Intersection_Segment_Segment)(const Vector2d & s_0_s, const Vector2d & s_0_e, const Vector2d & s_1_s, const Vector2d & s_1_e, Vector2d & inter);
typedef bool (*HGP_2D_Intersection_Line_Line)(const Vector2d & s_0_s, const Vector2d & s_0_e, const Vector2d & s_1_s, const Vector2d & s_1_e, Vector2d & inter);
typedef bool (*HGP_2D_Intersection_Segment_Line)(const Vector2d & s_s, const Vector2d & s_e, const Vector2d & l_s, const Vector2d & l_e, Vector2d & inter);
typedef bool (*HGP_2D_Intersection_Segment_Polygon)(const Vector2d & s_s, const Vector2d & s_e, const Vector2d1 & p);
typedef bool (*HGP_2D_Intersection_Polygon_Polygon)(const Vector2d1 & p1, const Vector2d1 & p2);
typedef bool (*HGP_2D_Polygon_Is_Clockwise_Oriented)(const Vector2d1 & ps);
typedef double (*HGP_2D_Two_Polygons_Union)(const Vector2d1 & poly_0, const Vector2d1 & poly_1, Vector2d2 & inter_polygons);
typedef double (*HGP_2D_Two_Polygons_Intersection)(const Vector2d1 & poly_0, const Vector2d1 & poly_1);
typedef void (*HGP_Decompose_Polyline)(const Vector2d1 & polyline, const double& threshold, Vector1i1 & result);
typedef bool (*HGP_Identify_Polycut_Extend)(const Vector2d1 & polygon, const Vector2d & s, const Vector2d & e, Vector2d & ns, Vector2d & ne);
typedef bool (*HGP_Identify_Polycut_NotExtend)(const Vector2d1 & polygon, const Vector2d & s, const Vector2d & e);
typedef bool (*HGP_Identify_Polycut)(const Vector2d1 & polygon, const Vector2d1 & cutLine, VectorPB1 & result);
typedef bool (*HGP_Construct_InOutSide_Polygon)(const Vector2d1 & py, const Vector2d & p, const Vector2d & q, bool& isPInside, bool& isQInside);
typedef bool (*HGP_2D_Intersection_Ray_Segment)(const Vector2d & s_0_s, const Vector2d & s_0_e, const Vector2d & s_1_s, const Vector2d & s_1_e, Vector2d & inter);
typedef double (*HGP_Get_Angle_Kerf_Offset_Tan)(const Vector2d & a, const Vector2d & b);
typedef Vector2d (*HGP_2D_Projection_Point_Segment)(const Vector2d & p, const Vector2d & s, const Vector2d & e);
typedef bool (*HGP_2D_Detect_Polygon_Inside_C1)(const Vector2d1 & outside_py, const Vector2d & p);
typedef bool (*HGP_2D_Detect_Polygon_Inside_C2)(const Vector2d1 & outside_py, const Vector2d1 & inside_py);
typedef bool (*HGP_2D_Detect_Polygon_Inside_C3)(const Vector2d2 & outside_pys, const Vector2d & p);
typedef bool (*HGP_2D_Detect_Polygon_Inside_C4)(const Vector2d2 & outside_pys, const Vector2d1 & inside_py);
typedef bool (*HGP_2D_Detect_Polygon_Inside_C5)(const Vector2d2 & outside_pys, const Vector2d2 & inside_pys);
typedef double (*HGP_2D_Distance_Polygon_Polygon)(const Vector2d1 & poly_0, const Vector2d1 & poly_1);
typedef double (*HGP_2D_Distance_Polygons_Polygons)(const Vector2d2 & poly_0, const Vector2d2 & poly_1);
typedef Vector2d (*HGP_2D_Nearest_Point_Polygon_C1)(const Vector2d & v, const Vector2d1 & poly);
typedef void (*HGP_2D_Nearest_Point_Polygon_C2)(const Vector2d & v, const Vector2d1 & poly, Vector2d & p, double& min_d);
typedef Vector2d (*HGP_2D_Nearest_Point_Polygons)(const Vector2d & v, const Vector2d2 & polys);
typedef void (*HGP_2d_Polygon_Boundingbox)(const Vector2d1 & ps, Vector2d & min_corner, Vector2d & max_corner);
typedef double (*HGP_2D_Polygon_Area)(const Vector2d1 & py);
typedef Vector2d (*HGP_2D_Polygon_Inside_Point_C1)(const Vector2d1 & poly);
typedef bool (*HGP_2D_Polygon_Inside_Point_C2)(const Vector2d2 & polys, Vector2d & inner_vec);
typedef void (*HGP_2D_Polygon_One_Offsets)(const Vector2d1 & poly, const double& d, Vector2d2 & offset_polys);
typedef void (*HGP_2D_Polygons_One_Offsets)(const Vector2d2 & polys, const double& d, Vector2d2 & offset_polys);
typedef bool (*HGP_2D_Polygons_Simple)(const Vector2d2 & poly);
typedef bool (*HGP_2D_Polygon_Simple)(const Vector2d1 & poly);
typedef bool (*HGP_2D_Polygon_Simple_Inter)(const Vector2d1 & poly);
typedef void (*HGP_2D_Convex_Hulls)(const Vector2d1 & vec, Vector2d1 & hull_points);
typedef void (*HGP_2D_OBB_Box)(const Vector2d1 & vec, Vector2d & center, Vector2d & axis_0, Vector2d & axis_1, double& entent_0, double& entent_1);
typedef void (*HGP_Image_Grid_Decomposition_C1)(Vector1i2 & image, Vector1d2 & boundary_xs, Vector1d2 & boundary_ys);
typedef void (*HGP_Image_Grid_Decomposition_Conservative_C1)(Vector1i2 & image, Vector1d2 & boundary_xs, Vector1d2 & boundary_ys);
typedef void (*HGP_Image_Grid_Decomposition_C2)(Vector1i2 & image, Vector2d2 & boundaries);
typedef void (*HGP_Image_Grid_Decomposition_Conservative_C2)(Vector1i2 & image, Vector2d2 & boundaries);
//implementation in "hgp3d.cpp"
//####################################################################################
typedef double (*HGP_3D_Distance_Point_Segment)(const Vector3d & p, const Vector3d & s_s, const Vector3d & s_e);
typedef void (*HGP_3D_Plane_Fitting)(const Vector3d1 & points, Vector3d & plane_p, Vector3d & plane_n);
typedef void (*HGP_3D_Plane_Point_Projection)(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d & p, Vector3d & result);
typedef void (*HGP_3D_Plane_Points_Projection)(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d1 & points, Vector3d1 & project_points);
typedef void (*HGP_3D_Plane_3D_to_2D_Point)(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d & point_3d, Vector2d & result);
typedef void (*HGP_3D_Plane_2D_to_3D_Point)(const Vector3d & plane_p, const Vector3d & plane_n, const Vector2d & points_2d, Vector3d & result);
typedef void (*HGP_3D_Plane_3D_to_2D_Points)(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d1 & points_3d, Vector2d1 & points_2d);
typedef void (*HGP_3D_Plane_3Ds_to_2Ds_Points)(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d2 & points_3d, Vector2d2 & points_2d);
typedef void (*HGP_3D_Plane_2D_to_3D_Points)(const Vector3d & plane_p, const Vector3d & plane_n, const Vector2d1 & points_2d, Vector3d1 & points_3d);
typedef Vector3d (*HGP_3D_Projection_Point_Segment)(const Vector3d & p, const Vector3d & s_s, const Vector3d & s_e);
typedef double (*HGP_3D_Distance_Point_Point)(const Vector3d & v0, const Vector3d & v1);
typedef double (*HGP_3D_Distance_Point_Polygon)(const Vector3d1 & py, const Vector3d & p);
typedef void (*HGP_2D_Polygon_Triangulation_C1)(const Vector2d2 & polys, Vector1i2 & faces);
typedef Vector1i2 (*HGP_2D_Polygon_Triangulation_C2)(const Vector2d2 & polys);
typedef Vector1i2 (*HGP_2D_Polygon_Triangulation_C3)(const Vector2d1 & poly);
typedef double (*HGP_3D_Distance_Point_Line)(const Vector3d & p, const Vector3d & l_s, const Vector3d & l_e);
typedef Vector3d (*HGP_3D_Projection_Point_Line)(const Vector3d & p, const Vector3d & l_s, const Vector3d & l_e);
typedef double (*HGP_3D_Distance_Segment_Segment)(const Vector3d & s_0_s, const Vector3d & s_0_e, const Vector3d & s_1_s, const Vector3d & s_1_e);
typedef double (*HGP_3D_Distance_Point_Plane)(const Vector3d & v, const Vector3d & plane_p, const Vector3d & plane_n);
typedef bool (*HGP_3D_Intersection_Segment_Line)(const Vector3d & s_s, const Vector3d & s_e, const Vector3d & l_s, const Vector3d & l_e, Vector3d & inter);
typedef bool (*HGP_3D_Intersection_Segment_Segment)(const Vector3d & s_0_s, const Vector3d & s_0_e, const Vector3d & s_1_s, const Vector3d & s_1_e, Vector3d & iter);
typedef bool (*HGP_3D_Intersection_Segment_Plane)(const Vector3d & s_s, const Vector3d & s_e, const Vector3d & plane_p, const Vector3d & plane_n, Vector3d & inter);
typedef bool (*HGP_3D_Intersection_Line_Plane)(const Vector3d & l_s, const Vector3d & l_e, const Vector3d & plane_p, const Vector3d & plane_n, Vector3d & inter);
typedef Vector3d (*HGP_3D_Projection_Point_Plane_C1)(const Vector3d & p, const Vector3d & plane_p, const Vector3d & plane_n);
typedef Vector3d (*HGP_3D_Projection_Point_Plane_C2)(const Vector3d & p, const Vector3d & plane_p_0, const Vector3d & plane_p_1, const Vector3d & plane_p_2);
typedef Vector2d (*HGP_3D_Projection_3D_Point_Plane_2D_C1)(const Vector3d & p, const Vector3d & plane_p, const Vector3d & plane_n);
typedef Vector2d (*HGP_3D_Projection_3D_Point_Plane_2D_C2)(const Vector3d & p, const Vector3d & plane_p_0, const Vector3d & plane_p_1, const Vector3d & plane_p_2);
typedef void (*HGP_3D_Plane_ABCD)(const Vector3d & plane_p, const Vector3d & plane_n, double& a, double& b, double& c, double& d);
typedef Vector3d (*HGP_3D_Plane_Base_1)(const Vector3d & plane_p, const Vector3d & plane_n);
typedef Vector3d (*HGP_Face_Normal)(const Vector3d & source, const Vector3d & tri_0, const Vector3d & tri_1, const Vector3d & tri_2, Vector3d & normal_0, Vector3d & normal_1, Vector3d & normal_2);
//implementation in "hgpmesh.cpp"
//####################################################################################
typedef void (*HGP_Remesh_Surface_by_Adding_Feature)(const Vector3d1 & feature, const Vector1i1 & face_ids, const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, Vector1i1 & igl_cutting_0_edges, Vector1i1 & igl_cutting_1_edges, Vector3d1 & igl_cutting_points, Vector1i2 & cutting_faces);
typedef void (*HGP_Mesh_Edges)(const char* path);
typedef bool (*HGP_3D_Intersection_Sphere_Ray)(const double& center_x, const double& center_y, const double& center_z, const double& radius, const double& ray_origin_x, const double& ray_origin_y, const double& ray_origin_z, const double& ray_direction_x, const double& ray_direction_y, const double& ray_direction_z, Vector1d1 & i_x, Vector1d1 & i_y, Vector1d1 & i_z);
typedef bool (*HGP_3D_Intersection_Ray_Triangle)(const Vector3d & p, const Vector3d & n, const Vector3d & p0, const Vector3d & p1, const Vector3d & p2);
typedef bool (*HGP_3D_Intersection_Ray_Mesh)(const Vector3d & p, const Vector3d & n, const char* path);
typedef bool (*HGP_3D_Intersection_Segment_Mesh)(const Vector3d & s, const Vector3d & e, const char* path);
typedef void (*HGP_3D_Intersection_Segments_Mesh)(const Vector3d1 & ss, const Vector3d1 & ee, const char* path, Vector1b1 & inters);
typedef void (*HGP_3D_Intersection_Polygons_Mesh)(const Vector3d2 & polygons, const char* path, Vector1b1 & inters);
//check whether there is a polygon intersected with the input mesh
typedef bool (*HGP_3D_Intersection_Polygons_Mesh_Bool)(const Vector3d2 & polygons, const char* path);
typedef void (*HGP_3D_Intersection_Rays_Mesh_Vector3d)(const Vector3d1 & ps, const Vector3d1 & ns, const char* path, Vector3d1 & inters);
//test each group directions (nes[i]) for each point in ps
typedef void (*HGP_3D_Intersection_Rays_Mesh_C1_Bool)(const Vector3d1 & ps, const Vector3d2 & nes, const char* path, Vector1b2 & inters);
//test all directions (ns) for each point in ps
typedef void (*HGP_3D_Intersection_Rays_Mesh_C2_Bool)(const Vector3d1 & ps, const Vector3d1 & ns, const char* path, Vector1b2 & inters);
typedef void (*HGP_3D_Intersection_Rays_Mesh_C2_Vector3d)(const Vector3d1 & ps, const Vector3d1 & ns, const char* path, Vector1d2 & inters);
typedef void (*HGP_3D_Points_Inside_Triangles_C1_Bool)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, const Vector3d1 & points, Vector1b1 & insides);
typedef void (*HGP_3D_Points_Inside_Triangles_C2_Bool)(const char* path, const Vector3d1 & points, Vector1b1 & insides);
//d: percentage value of the length of the diagonal of the bounding box.
typedef void (*HGP_3D_Mesh_Dart_Sampling_C1)(const char* outside_path, const double& d, Vector3d1 & sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
typedef void (*HGP_3D_Mesh_Dart_Sampling_C2)(const char* outside_path, const char* inside_path, const double& d, Vector3d1 & sampling_points, const int& total_iter);
//d: percentage value of the length of the diagonal of the bounding box.
typedef void (*HGP_3D_Mesh_Regular_Sampling_C1)(const char* outside_path, const double& d, Vector3d1 & sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
typedef void (*HGP_3D_Mesh_Regular_Sampling_C2)(const char* outside_path, const char* inside_path, const double& d, Vector3d1 & sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
typedef void (*HGP_3D_Cube_Surface_Sampling_C1)(const double& cube_size, const double& d, Vector3d2 & sampling_points, VectorPI2 & neighbors, const bool& compute_neighbors);
//d: percentage value of the length of the diagonal of the bounding box.
typedef void (*HGP_3D_Cube_Surface_Sampling_C2)(const double& cube_size, const double& d, Vector3d2 & sampling_points);
//d: percentage value of the length of the diagonal of the bounding box.
typedef void (*HGP_3D_Cube_Surface_Sampling_C3)(const double& cube_size, const double& d, Vector3d2 & sampling_points, VectorPI2 & neighbors);
//with neighboring
typedef void (*HGP_3D_Mesh_Regular_Sampling_C3)(const char* outside_path, const char* inside_path, const double& d, Vector3d1 & sampling_points, VectorPI1 & neighbors);
typedef double (*HGP_3D_Distance_Point_Triangle)(const Vector3d & p, const Vector3d & t_0, const Vector3d & t_1, const Vector3d & t_2);
typedef double (*HGP_3D_Distance_Point_Triangles)(const Vector3d & p, const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2);
typedef Vector3d (*HGP_3D_Nearest_Point_Triangles)(const Vector3d & p, const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2);
typedef void (*HGP_3D_Distance_Point_Mesh)(const char* path, const Vector3d1 & query_points, Vector1d1 & distances);
typedef void (*HGP_3D_Neareast_Point_Mesh)(const char* path, const Vector3d1 & ves, Vector3d1 & ners);
typedef void  (*HGP_3D_Mesh_Near_Triangles)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, const Vector3d1 & points, const double& d, Vector1i2 & triangles);
typedef void (*HGP_3D_Points_inside_Triangles_C1)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, const Vector3d1 & points, Vector1b1 & insides);
typedef void (*HGP_3D_Points_inside_Triangles_C2)(const char* path, const Vector3d1 & points, Vector1b1 & insides);
typedef void (*HGP_Mesh_Subdivision)(const char* in_path, const char* sub, const int& step, const char* out_path);
typedef void (*HGP_Mesh_Loop_Subdivision_One_Step)(Vector3d1 & vecs, Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2);
typedef void (*HGP_3D_Mesh_Curvature_C1)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, Vector1d1 & max_curs, Vector1d1 & min_curs);
typedef void (*HGP_3D_Mesh_Curvature_C2)(const Vector3d1 & vecs, const Vector1i2 & face_ids, Vector1d1 & max_curs, Vector1d1 & min_curs);
typedef void (*HGP_3D_Mesh_Curvature_C3)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, Vector1d1 & max_curs, Vector1d1 & min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions);
typedef void (*HGP_3D_Mesh_Curvature_C4)(const Vector3d1 & vecs, const Vector1i2 & face_ids, Vector1d1 & max_curs, Vector1d1 & min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions);
typedef void (*HGP_3D_Mesh_Curvature_C5)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, Vector1d1 & max_curs, Vector1d1 & min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions, Vector3d1 & normals);
typedef void (*HGP_3D_Mesh_Curvature_C6)(const Vector3d1 & vecs, const Vector1i2 & face_ids, Vector1d1 & max_curs, Vector1d1 & min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions, Vector3d1 & normals);
typedef void (*HGP_3D_Triangle_Mesh_Boundary_C1)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, Vector1b1 & bools);
typedef void (*HGP_3D_Triangle_Mesh_Boundary_C2)(const char* path, Vector1b1 & bools);
typedef void (*HGP_3D_Connecting_Segments_C1)(Vector2d2 & segments, Vector2d2 & lines);
typedef void (*HGP_3D_Connecting_Segments_C2)(Vector3d2 & segments, Vector3d2 & lines);
typedef void (*HGP_3D_Triangle_Mesh_Boundary_C3)(Vector3d1 & vecs, Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2, Vector3d2 & boundaries);
typedef void (*HGP_3D_Triangle_Mesh_Boundary_C4)(Vector3d1 & vecs, Vector1i2 & face_ids, Vector3d2 & boundaries);
typedef void (*HGP_3D_Triangle_Mesh_Boundary_C5)(const char* path, Vector3d2 & boundaries);
typedef void (*HGP_Mesh_Laplace_Smooth_C1)(const char* in_path, const char* out_path, const int laplace_nb);
typedef void (*HGP_3D_Triangle_Mesh_Vecs_Neighbors)(Vector3d1 & vecs, Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2, Vector1i2 & neighs);
typedef void (*HGP_Mesh_Laplace_Smooth_C2)(Vector3d1 & vecs, Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2, const int laplace_nb);
typedef void (*HGP_3D_Triangle_Mesh_Vecs_Faces)(Vector3d1 & vecs, Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2, Vector1i2 & surface_vectices_to_face);
typedef void (*HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges)(Vector3d1 & vecs, Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2, Vector1i3 & surface_vectices_to_neighbor_edges);
typedef void (*HGP_Mesh_Laplace_Smooth_by_Curvature)(Vector3d1 & vecs, Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2, double& low_curvature);
typedef void (*HGP_Mesh_Loop_Subdivision_Own_Version)(const char* in_path, const int& step, const char* out_path, const int& laplace_nb);
typedef void (*HGP_Rotation_Obj)(const char* path, const double& angle, const Vector3d & axis);
typedef void (*HGP_Slicer_Mesh)(const char* path, const Vector3d & plane_normal, const Vector1d1 & plane_d, Vector3d3 & offsetses, Vector3d2 & offsets);
typedef void (*HGP_Shortest_Geodesic_Path_C1)(const char* path, Vector3d1 & xyzs);
typedef void (*HGP_Shortest_Geodesic_Path_C3)(const char* path, Vector3d source, Vector3d target, Vector3d1 & output);
typedef void (*HGP_Shortest_Geodesic_Path_C4)(const char* path, Vector3d1 sources, Vector3d1 targets, Vector3d2 & xyzes);
typedef double (*HGP_Geodesic_Distance)(const char* path, const Vector3d & source, const Vector3d & target);
typedef Vector3d1 (*HGP_Project_Points_Onto_Surface_C1)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, const Vector3d1 & points);
typedef Vector3d1 (*HGP_Project_Points_Onto_Surface_C2)(const char* path, const Vector3d1 & points);
typedef void (*HGP_3D_Triangel_Mesh_Most_Inside_Point)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, Vector3d & inside);
typedef double (*HGP_3D_One_Triangle_Area)(const Vector3d & v0, const Vector3d & v1, const Vector3d & v2);
typedef double (*HGP_3D_Triangle_Mesh_Area)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2);
typedef void (*HGP_3D_Convex_Hulls_C1)(const Vector3d1 & vec, Vector3d1 & hull_points);
typedef void (*HGP_3D_Convex_Hulls_C2)(const Vector3d1 & vec, Vector3d1 & hull_points, Vector1i1 & hulls_surface_0, Vector1i1 & hulls_surface_1, Vector1i1 & hulls_surface_2);
typedef void (*HGP_3D_Convex_Hulls_C3)(const Vector3d1 & vec, Vector3d1 & hull_points, Vector3d1 & plane_p, Vector3d1 & plane_n);
typedef void (*HGP_3D_Convex_Hulls_C4)(const Vector3d1 & vec, Vector3d1 & hull_points, Vector1i1 & hulls_surface_0, Vector1i1 & hulls_surface_1, Vector1i1 & hulls_surface_2, Vector3d1 & plane_p, Vector3d1 & plane_n);
typedef void (*HGP_Mesh_Field_Query_C1)(const char* path, const Vector3d1 & gradients, const Vector3d1 & input_points, Vector3d1 & points_gradients);
typedef void (*HGP_Mesh_Field_Query_C2)(const char* path, const Vector1d1 & gradient_values, const Vector3d1 & input_points, Vector1d1 & points_gradient_values);
typedef void (*HGP_Mesh_Field_Query_C3)(const char* path, const Vector1d1 & gradient_values, const Vector3d2 & input_point_es, Vector1d2 & points_gradient_value_es);
typedef void (*HGP_Curvature_Mesh)(const char* path, const Vector3d1 & input_points, Vector1d1 & max_curs, Vector1d1 & min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions);
typedef void (*HGP_Normal_Mesh_C1)(const char* path, const Vector3d1 & mesh_points, Vector3d1 & mesh_normals);
typedef void (*HGP_Normal_Mesh_C2)(const char* path, const Vector3d2 & mesh_pointses, Vector3d2 & mesh_normalses);
typedef void (*HGP_3D_Mesh_Normal_C1)(const Vector3d1 & ps, const Vector1i2 & face_ids, Vector3d1 & normals);
typedef void (*HGP_3D_Mesh_Normal_C2)(const Vector3d1 & ps, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, Vector3d1 & normals);
typedef Vector3d (*HGP_3D_Mesh_Center_C1)(const Vector3d2 & ps);
typedef Vector3d (*HGP_3D_Mesh_Center_C2)(const Vector3d1 & ps);
typedef void (*HGP_3D_Mesh_Boundingbox_C1)(const Vector3d2 & ps, Vector3d & min_corner, Vector3d & max_corner);
typedef void (*HGP_3D_Mesh_Boundingbox_C2)(const Vector3d1 & ps, Vector3d & min_corner, Vector3d & max_corner);
typedef void (*HGP_Surface_Decomposition)(const char* path, Vector1d1 & face_sdf, int& regions_nb, Vector1i1 & face_segments);
typedef void (*HGP_3D_Mesh_Gradient)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, const Vector1d1 & psd, Vector3d1 & vecs_gradients, Vector3d1 & faces_gradients);
typedef void (*HGP_Intergral_Curvature)(const Vector2d1 & input_points, const int& sampling_points_nb, const double& radius, const double& thresholder, Vector2d1 & output_points, Vector1d1 & output_rates);
typedef bool (*HGP_3D_Mesh_Extract_Isoline)(const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2, const Vector1d1 & psd, const double& d, Vector3d2 & isolines);
typedef void (*HGP_BSplineCurveFit)(const Vector3d1 & samples, Vector3d1 & output);
typedef void (*HGP_Cut_Surface)(const Vector3d1 & boundary, const Vector3d & inside_point, const char* full_path, char* output_path);
typedef void (*HGP_Cut_Surface_by_Multi_Boundaries)(const Vector3d2 & multi_boundary, const Vector3d & inside_point, const char* full_path, char* output_path);
/////////////////////////////////////////////////////////////
//


class HGPPL
{
	public:
	HGPPL()
	{
		hModule = LoadHMODULE("libhgp.dll");
		HGP_Vector_Base_C = (HGP_Vector_Base)GetProcAddress(hModule, "HGP_Vector_Base");
		HGP_Test_PGL_C = (HGP_Test_PGL)GetProcAddress(hModule, "HGP_Test_PGL");
		//implementation in "hgp2d.cpp"
		//####################################################################################
		HGP_2D_Distance_Point_Point_C = (HGP_2D_Distance_Point_Point)GetProcAddress(hModule, "HGP_2D_Distance_Point_Point");
		HGP_2D_Distance_Point_Line_C = (HGP_2D_Distance_Point_Line)GetProcAddress(hModule, "HGP_2D_Distance_Point_Line");
		HGP_2D_Distance_Point_Segment_C = (HGP_2D_Distance_Point_Segment)GetProcAddress(hModule, "HGP_2D_Distance_Point_Segment");
		HGP_2D_Distance_Segment_Segment_C = (HGP_2D_Distance_Segment_Segment)GetProcAddress(hModule, "HGP_2D_Distance_Segment_Segment");
		HGP_2D_Location_Point_Polygon_C = (HGP_2D_Location_Point_Polygon)GetProcAddress(hModule, "HGP_2D_Location_Point_Polygon");
		HGP_2D_Location_Points_Polygon_C = (HGP_2D_Location_Points_Polygon)GetProcAddress(hModule, "HGP_2D_Location_Points_Polygon");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_2D_Polygon_Dart_Sampling_C = (HGP_2D_Polygon_Dart_Sampling)GetProcAddress(hModule, "HGP_2D_Polygon_Dart_Sampling");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_2D_Polygon_Regular_Sampling_C1_C = (HGP_2D_Polygon_Regular_Sampling_C1)GetProcAddress(hModule, "HGP_2D_Polygon_Regular_Sampling_C1");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_2D_Polygon_Regular_Sampling_C2_C = (HGP_2D_Polygon_Regular_Sampling_C2)GetProcAddress(hModule, "HGP_2D_Polygon_Regular_Sampling_C2");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_2D_Polygon_Regular_Sampling_C3_C = (HGP_2D_Polygon_Regular_Sampling_C3)GetProcAddress(hModule, "HGP_2D_Polygon_Regular_Sampling_C3");
		HGP_2D_Square_Regular_Sampling_C1_C = (HGP_2D_Square_Regular_Sampling_C1)GetProcAddress(hModule, "HGP_2D_Square_Regular_Sampling_C1");
		HGP_2D_Square_Regular_Sampling_C2_C = (HGP_2D_Square_Regular_Sampling_C2)GetProcAddress(hModule, "HGP_2D_Square_Regular_Sampling_C2");
		HGP_2D_Square_Regular_Sampling_C3_C = (HGP_2D_Square_Regular_Sampling_C3)GetProcAddress(hModule, "HGP_2D_Square_Regular_Sampling_C3");
		HGP_2D_Distance_Point_Polygon_C = (HGP_2D_Distance_Point_Polygon)GetProcAddress(hModule, "HGP_2D_Distance_Point_Polygon");
		HGP_2D_Distance_Point_Polygons_C = (HGP_2D_Distance_Point_Polygons)GetProcAddress(hModule, "HGP_2D_Distance_Point_Polygons");
		HGP_2D_Intersection_Segment_Segment_C = (HGP_2D_Intersection_Segment_Segment)GetProcAddress(hModule, "HGP_2D_Intersection_Segment_Segment");
		HGP_2D_Intersection_Line_Line_C = (HGP_2D_Intersection_Line_Line)GetProcAddress(hModule, "HGP_2D_Intersection_Line_Line");
		HGP_2D_Intersection_Segment_Line_C = (HGP_2D_Intersection_Segment_Line)GetProcAddress(hModule, "HGP_2D_Intersection_Segment_Line");
		HGP_2D_Intersection_Segment_Polygon_C = (HGP_2D_Intersection_Segment_Polygon)GetProcAddress(hModule, "HGP_2D_Intersection_Segment_Polygon");
		HGP_2D_Intersection_Polygon_Polygon_C = (HGP_2D_Intersection_Polygon_Polygon)GetProcAddress(hModule, "HGP_2D_Intersection_Polygon_Polygon");
		HGP_2D_Polygon_Is_Clockwise_Oriented_C = (HGP_2D_Polygon_Is_Clockwise_Oriented)GetProcAddress(hModule, "HGP_2D_Polygon_Is_Clockwise_Oriented");
		HGP_2D_Two_Polygons_Union_C = (HGP_2D_Two_Polygons_Union)GetProcAddress(hModule, "HGP_2D_Two_Polygons_Union");
		HGP_2D_Two_Polygons_Intersection_C = (HGP_2D_Two_Polygons_Intersection)GetProcAddress(hModule, "HGP_2D_Two_Polygons_Intersection");
		HGP_Decompose_Polyline_C = (HGP_Decompose_Polyline)GetProcAddress(hModule, "HGP_Decompose_Polyline");
		HGP_Identify_Polycut_Extend_C = (HGP_Identify_Polycut_Extend)GetProcAddress(hModule, "HGP_Identify_Polycut_Extend");
		HGP_Identify_Polycut_NotExtend_C = (HGP_Identify_Polycut_NotExtend)GetProcAddress(hModule, "HGP_Identify_Polycut_NotExtend");
		HGP_Identify_Polycut_C = (HGP_Identify_Polycut)GetProcAddress(hModule, "HGP_Identify_Polycut");
		HGP_Construct_InOutSide_Polygon_C = (HGP_Construct_InOutSide_Polygon)GetProcAddress(hModule, "HGP_Construct_InOutSide_Polygon");
		HGP_2D_Intersection_Ray_Segment_C = (HGP_2D_Intersection_Ray_Segment)GetProcAddress(hModule, "HGP_2D_Intersection_Ray_Segment");
		HGP_Get_Angle_Kerf_Offset_Tan_C = (HGP_Get_Angle_Kerf_Offset_Tan)GetProcAddress(hModule, "HGP_Get_Angle_Kerf_Offset_Tan");
		HGP_2D_Projection_Point_Segment_C = (HGP_2D_Projection_Point_Segment)GetProcAddress(hModule, "HGP_2D_Projection_Point_Segment");
		HGP_2D_Detect_Polygon_Inside_C1_C = (HGP_2D_Detect_Polygon_Inside_C1)GetProcAddress(hModule, "HGP_2D_Detect_Polygon_Inside_C1");
		HGP_2D_Detect_Polygon_Inside_C2_C = (HGP_2D_Detect_Polygon_Inside_C2)GetProcAddress(hModule, "HGP_2D_Detect_Polygon_Inside_C2");
		HGP_2D_Detect_Polygon_Inside_C3_C = (HGP_2D_Detect_Polygon_Inside_C3)GetProcAddress(hModule, "HGP_2D_Detect_Polygon_Inside_C3");
		HGP_2D_Detect_Polygon_Inside_C4_C = (HGP_2D_Detect_Polygon_Inside_C4)GetProcAddress(hModule, "HGP_2D_Detect_Polygon_Inside_C4");
		HGP_2D_Detect_Polygon_Inside_C5_C = (HGP_2D_Detect_Polygon_Inside_C5)GetProcAddress(hModule, "HGP_2D_Detect_Polygon_Inside_C5");
		HGP_2D_Distance_Polygon_Polygon_C = (HGP_2D_Distance_Polygon_Polygon)GetProcAddress(hModule, "HGP_2D_Distance_Polygon_Polygon");
		HGP_2D_Distance_Polygons_Polygons_C = (HGP_2D_Distance_Polygons_Polygons)GetProcAddress(hModule, "HGP_2D_Distance_Polygons_Polygons");
		HGP_2D_Nearest_Point_Polygon_C1_C = (HGP_2D_Nearest_Point_Polygon_C1)GetProcAddress(hModule, "HGP_2D_Nearest_Point_Polygon_C1");
		HGP_2D_Nearest_Point_Polygon_C2_C = (HGP_2D_Nearest_Point_Polygon_C2)GetProcAddress(hModule, "HGP_2D_Nearest_Point_Polygon_C2");
		HGP_2D_Nearest_Point_Polygons_C = (HGP_2D_Nearest_Point_Polygons)GetProcAddress(hModule, "HGP_2D_Nearest_Point_Polygons");
		HGP_2d_Polygon_Boundingbox_C = (HGP_2d_Polygon_Boundingbox)GetProcAddress(hModule, "HGP_2d_Polygon_Boundingbox");
		HGP_2D_Polygon_Area_C = (HGP_2D_Polygon_Area)GetProcAddress(hModule, "HGP_2D_Polygon_Area");
		HGP_2D_Polygon_Inside_Point_C1_C = (HGP_2D_Polygon_Inside_Point_C1)GetProcAddress(hModule, "HGP_2D_Polygon_Inside_Point_C1");
		HGP_2D_Polygon_Inside_Point_C2_C = (HGP_2D_Polygon_Inside_Point_C2)GetProcAddress(hModule, "HGP_2D_Polygon_Inside_Point_C2");
		HGP_2D_Polygon_One_Offsets_C = (HGP_2D_Polygon_One_Offsets)GetProcAddress(hModule, "HGP_2D_Polygon_One_Offsets");
		HGP_2D_Polygons_One_Offsets_C = (HGP_2D_Polygons_One_Offsets)GetProcAddress(hModule, "HGP_2D_Polygons_One_Offsets");
		HGP_2D_Polygons_Simple_C = (HGP_2D_Polygons_Simple)GetProcAddress(hModule, "HGP_2D_Polygons_Simple");
		HGP_2D_Polygon_Simple_C = (HGP_2D_Polygon_Simple)GetProcAddress(hModule, "HGP_2D_Polygon_Simple");
		HGP_2D_Polygon_Simple_Inter_C = (HGP_2D_Polygon_Simple_Inter)GetProcAddress(hModule, "HGP_2D_Polygon_Simple_Inter");
		HGP_2D_Convex_Hulls_C = (HGP_2D_Convex_Hulls)GetProcAddress(hModule, "HGP_2D_Convex_Hulls");
		HGP_2D_OBB_Box_C = (HGP_2D_OBB_Box)GetProcAddress(hModule, "HGP_2D_OBB_Box");
		HGP_Image_Grid_Decomposition_C1_C = (HGP_Image_Grid_Decomposition_C1)GetProcAddress(hModule, "HGP_Image_Grid_Decomposition_C1");
		HGP_Image_Grid_Decomposition_Conservative_C1_C = (HGP_Image_Grid_Decomposition_Conservative_C1)GetProcAddress(hModule, "HGP_Image_Grid_Decomposition_Conservative_C1");
		HGP_Image_Grid_Decomposition_C2_C = (HGP_Image_Grid_Decomposition_C2)GetProcAddress(hModule, "HGP_Image_Grid_Decomposition_C2");
		HGP_Image_Grid_Decomposition_Conservative_C2_C = (HGP_Image_Grid_Decomposition_Conservative_C2)GetProcAddress(hModule, "HGP_Image_Grid_Decomposition_Conservative_C2");
		//implementation in "hgp3d.cpp"
		//####################################################################################
		HGP_3D_Distance_Point_Segment_C = (HGP_3D_Distance_Point_Segment)GetProcAddress(hModule, "HGP_3D_Distance_Point_Segment");
		HGP_3D_Plane_Fitting_C = (HGP_3D_Plane_Fitting)GetProcAddress(hModule, "HGP_3D_Plane_Fitting");
		HGP_3D_Plane_Point_Projection_C = (HGP_3D_Plane_Point_Projection)GetProcAddress(hModule, "HGP_3D_Plane_Point_Projection");
		HGP_3D_Plane_Points_Projection_C = (HGP_3D_Plane_Points_Projection)GetProcAddress(hModule, "HGP_3D_Plane_Points_Projection");
		HGP_3D_Plane_3D_to_2D_Point_C = (HGP_3D_Plane_3D_to_2D_Point)GetProcAddress(hModule, "HGP_3D_Plane_3D_to_2D_Point");
		HGP_3D_Plane_2D_to_3D_Point_C = (HGP_3D_Plane_2D_to_3D_Point)GetProcAddress(hModule, "HGP_3D_Plane_2D_to_3D_Point");
		HGP_3D_Plane_3D_to_2D_Points_C = (HGP_3D_Plane_3D_to_2D_Points)GetProcAddress(hModule, "HGP_3D_Plane_3D_to_2D_Points");
		HGP_3D_Plane_3Ds_to_2Ds_Points_C = (HGP_3D_Plane_3Ds_to_2Ds_Points)GetProcAddress(hModule, "HGP_3D_Plane_3Ds_to_2Ds_Points");
		HGP_3D_Plane_2D_to_3D_Points_C = (HGP_3D_Plane_2D_to_3D_Points)GetProcAddress(hModule, "HGP_3D_Plane_2D_to_3D_Points");
		HGP_3D_Projection_Point_Segment_C = (HGP_3D_Projection_Point_Segment)GetProcAddress(hModule, "HGP_3D_Projection_Point_Segment");
		HGP_3D_Distance_Point_Point_C = (HGP_3D_Distance_Point_Point)GetProcAddress(hModule, "HGP_3D_Distance_Point_Point");
		HGP_3D_Distance_Point_Polygon_C = (HGP_3D_Distance_Point_Polygon)GetProcAddress(hModule, "HGP_3D_Distance_Point_Polygon");
		HGP_2D_Polygon_Triangulation_C1_C = (HGP_2D_Polygon_Triangulation_C1)GetProcAddress(hModule, "HGP_2D_Polygon_Triangulation_C1");
		HGP_2D_Polygon_Triangulation_C2_C = (HGP_2D_Polygon_Triangulation_C2)GetProcAddress(hModule, "HGP_2D_Polygon_Triangulation_C2");
		HGP_2D_Polygon_Triangulation_C3_C = (HGP_2D_Polygon_Triangulation_C3)GetProcAddress(hModule, "HGP_2D_Polygon_Triangulation_C3");
		HGP_3D_Distance_Point_Line_C = (HGP_3D_Distance_Point_Line)GetProcAddress(hModule, "HGP_3D_Distance_Point_Line");
		HGP_3D_Projection_Point_Line_C = (HGP_3D_Projection_Point_Line)GetProcAddress(hModule, "HGP_3D_Projection_Point_Line");
		HGP_3D_Distance_Segment_Segment_C = (HGP_3D_Distance_Segment_Segment)GetProcAddress(hModule, "HGP_3D_Distance_Segment_Segment");
		HGP_3D_Distance_Point_Plane_C = (HGP_3D_Distance_Point_Plane)GetProcAddress(hModule, "HGP_3D_Distance_Point_Plane");
		HGP_3D_Intersection_Segment_Line_C = (HGP_3D_Intersection_Segment_Line)GetProcAddress(hModule, "HGP_3D_Intersection_Segment_Line");
		HGP_3D_Intersection_Segment_Segment_C = (HGP_3D_Intersection_Segment_Segment)GetProcAddress(hModule, "HGP_3D_Intersection_Segment_Segment");
		HGP_3D_Intersection_Segment_Plane_C = (HGP_3D_Intersection_Segment_Plane)GetProcAddress(hModule, "HGP_3D_Intersection_Segment_Plane");
		HGP_3D_Intersection_Line_Plane_C = (HGP_3D_Intersection_Line_Plane)GetProcAddress(hModule, "HGP_3D_Intersection_Line_Plane");
		HGP_3D_Projection_Point_Plane_C1_C = (HGP_3D_Projection_Point_Plane_C1)GetProcAddress(hModule, "HGP_3D_Projection_Point_Plane_C1");
		HGP_3D_Projection_Point_Plane_C2_C = (HGP_3D_Projection_Point_Plane_C2)GetProcAddress(hModule, "HGP_3D_Projection_Point_Plane_C2");
		HGP_3D_Projection_3D_Point_Plane_2D_C1_C = (HGP_3D_Projection_3D_Point_Plane_2D_C1)GetProcAddress(hModule, "HGP_3D_Projection_3D_Point_Plane_2D_C1");
		HGP_3D_Projection_3D_Point_Plane_2D_C2_C = (HGP_3D_Projection_3D_Point_Plane_2D_C2)GetProcAddress(hModule, "HGP_3D_Projection_3D_Point_Plane_2D_C2");
		HGP_3D_Plane_ABCD_C = (HGP_3D_Plane_ABCD)GetProcAddress(hModule, "HGP_3D_Plane_ABCD");
		HGP_3D_Plane_Base_1_C = (HGP_3D_Plane_Base_1)GetProcAddress(hModule, "HGP_3D_Plane_Base_1");
		HGP_Face_Normal_C = (HGP_Face_Normal)GetProcAddress(hModule, "HGP_Face_Normal");
		//implementation in "hgpmesh.cpp"
		//####################################################################################
		HGP_Remesh_Surface_by_Adding_Feature_C = (HGP_Remesh_Surface_by_Adding_Feature)GetProcAddress(hModule, "HGP_Remesh_Surface_by_Adding_Feature");
		HGP_Mesh_Edges_C = (HGP_Mesh_Edges)GetProcAddress(hModule, "HGP_Mesh_Edges");
		HGP_3D_Intersection_Sphere_Ray_C = (HGP_3D_Intersection_Sphere_Ray)GetProcAddress(hModule, "HGP_3D_Intersection_Sphere_Ray");
		HGP_3D_Intersection_Ray_Triangle_C = (HGP_3D_Intersection_Ray_Triangle)GetProcAddress(hModule, "HGP_3D_Intersection_Ray_Triangle");
		HGP_3D_Intersection_Ray_Mesh_C = (HGP_3D_Intersection_Ray_Mesh)GetProcAddress(hModule, "HGP_3D_Intersection_Ray_Mesh");
		HGP_3D_Intersection_Segment_Mesh_C = (HGP_3D_Intersection_Segment_Mesh)GetProcAddress(hModule, "HGP_3D_Intersection_Segment_Mesh");
		HGP_3D_Intersection_Segments_Mesh_C = (HGP_3D_Intersection_Segments_Mesh)GetProcAddress(hModule, "HGP_3D_Intersection_Segments_Mesh");
		HGP_3D_Intersection_Polygons_Mesh_C = (HGP_3D_Intersection_Polygons_Mesh)GetProcAddress(hModule, "HGP_3D_Intersection_Polygons_Mesh");
		//check whether there is a polygon intersected with the input mesh
		HGP_3D_Intersection_Polygons_Mesh_Bool_C = (HGP_3D_Intersection_Polygons_Mesh_Bool)GetProcAddress(hModule, "HGP_3D_Intersection_Polygons_Mesh_Bool");
		HGP_3D_Intersection_Rays_Mesh_Vector3d_C = (HGP_3D_Intersection_Rays_Mesh_Vector3d)GetProcAddress(hModule, "HGP_3D_Intersection_Rays_Mesh_Vector3d");
		//test each group directions (nes[i]) for each point in ps
		HGP_3D_Intersection_Rays_Mesh_C1_Bool_C = (HGP_3D_Intersection_Rays_Mesh_C1_Bool)GetProcAddress(hModule, "HGP_3D_Intersection_Rays_Mesh_C1_Bool");
		//test all directions (ns) for each point in ps
		HGP_3D_Intersection_Rays_Mesh_C2_Bool_C = (HGP_3D_Intersection_Rays_Mesh_C2_Bool)GetProcAddress(hModule, "HGP_3D_Intersection_Rays_Mesh_C2_Bool");
		HGP_3D_Intersection_Rays_Mesh_C2_Vector3d_C = (HGP_3D_Intersection_Rays_Mesh_C2_Vector3d)GetProcAddress(hModule, "HGP_3D_Intersection_Rays_Mesh_C2_Vector3d");
		HGP_3D_Points_Inside_Triangles_C1_Bool_C = (HGP_3D_Points_Inside_Triangles_C1_Bool)GetProcAddress(hModule, "HGP_3D_Points_Inside_Triangles_C1_Bool");
		HGP_3D_Points_Inside_Triangles_C2_Bool_C = (HGP_3D_Points_Inside_Triangles_C2_Bool)GetProcAddress(hModule, "HGP_3D_Points_Inside_Triangles_C2_Bool");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_3D_Mesh_Dart_Sampling_C1_C = (HGP_3D_Mesh_Dart_Sampling_C1)GetProcAddress(hModule, "HGP_3D_Mesh_Dart_Sampling_C1");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_3D_Mesh_Dart_Sampling_C2_C = (HGP_3D_Mesh_Dart_Sampling_C2)GetProcAddress(hModule, "HGP_3D_Mesh_Dart_Sampling_C2");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_3D_Mesh_Regular_Sampling_C1_C = (HGP_3D_Mesh_Regular_Sampling_C1)GetProcAddress(hModule, "HGP_3D_Mesh_Regular_Sampling_C1");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_3D_Mesh_Regular_Sampling_C2_C = (HGP_3D_Mesh_Regular_Sampling_C2)GetProcAddress(hModule, "HGP_3D_Mesh_Regular_Sampling_C2");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_3D_Cube_Surface_Sampling_C1_C = (HGP_3D_Cube_Surface_Sampling_C1)GetProcAddress(hModule, "HGP_3D_Cube_Surface_Sampling_C1");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_3D_Cube_Surface_Sampling_C2_C = (HGP_3D_Cube_Surface_Sampling_C2)GetProcAddress(hModule, "HGP_3D_Cube_Surface_Sampling_C2");
		//d: percentage value of the length of the diagonal of the bounding box.
		HGP_3D_Cube_Surface_Sampling_C3_C = (HGP_3D_Cube_Surface_Sampling_C3)GetProcAddress(hModule, "HGP_3D_Cube_Surface_Sampling_C3");
		//with neighboring
		HGP_3D_Mesh_Regular_Sampling_C3_C = (HGP_3D_Mesh_Regular_Sampling_C3)GetProcAddress(hModule, "HGP_3D_Mesh_Regular_Sampling_C3");
		HGP_3D_Distance_Point_Triangle_C = (HGP_3D_Distance_Point_Triangle)GetProcAddress(hModule, "HGP_3D_Distance_Point_Triangle");
		HGP_3D_Distance_Point_Triangles_C = (HGP_3D_Distance_Point_Triangles)GetProcAddress(hModule, "HGP_3D_Distance_Point_Triangles");
		HGP_3D_Nearest_Point_Triangles_C = (HGP_3D_Nearest_Point_Triangles)GetProcAddress(hModule, "HGP_3D_Nearest_Point_Triangles");
		HGP_3D_Distance_Point_Mesh_C = (HGP_3D_Distance_Point_Mesh)GetProcAddress(hModule, "HGP_3D_Distance_Point_Mesh");
		HGP_3D_Neareast_Point_Mesh_C = (HGP_3D_Neareast_Point_Mesh)GetProcAddress(hModule, "HGP_3D_Neareast_Point_Mesh");
		HGP_3D_Mesh_Near_Triangles_C = (HGP_3D_Mesh_Near_Triangles)GetProcAddress(hModule, "HGP_3D_Mesh_Near_Triangles");
		HGP_3D_Points_inside_Triangles_C1_C = (HGP_3D_Points_inside_Triangles_C1)GetProcAddress(hModule, "HGP_3D_Points_inside_Triangles_C1");
		HGP_3D_Points_inside_Triangles_C2_C = (HGP_3D_Points_inside_Triangles_C2)GetProcAddress(hModule, "HGP_3D_Points_inside_Triangles_C2");
		HGP_Mesh_Subdivision_C = (HGP_Mesh_Subdivision)GetProcAddress(hModule, "HGP_Mesh_Subdivision");
		HGP_Mesh_Loop_Subdivision_One_Step_C = (HGP_Mesh_Loop_Subdivision_One_Step)GetProcAddress(hModule, "HGP_Mesh_Loop_Subdivision_One_Step");
		HGP_3D_Mesh_Curvature_C1_C = (HGP_3D_Mesh_Curvature_C1)GetProcAddress(hModule, "HGP_3D_Mesh_Curvature_C1");
		HGP_3D_Mesh_Curvature_C2_C = (HGP_3D_Mesh_Curvature_C2)GetProcAddress(hModule, "HGP_3D_Mesh_Curvature_C2");
		HGP_3D_Mesh_Curvature_C3_C = (HGP_3D_Mesh_Curvature_C3)GetProcAddress(hModule, "HGP_3D_Mesh_Curvature_C3");
		HGP_3D_Mesh_Curvature_C4_C = (HGP_3D_Mesh_Curvature_C4)GetProcAddress(hModule, "HGP_3D_Mesh_Curvature_C4");
		HGP_3D_Mesh_Curvature_C5_C = (HGP_3D_Mesh_Curvature_C5)GetProcAddress(hModule, "HGP_3D_Mesh_Curvature_C5");
		HGP_3D_Mesh_Curvature_C6_C = (HGP_3D_Mesh_Curvature_C6)GetProcAddress(hModule, "HGP_3D_Mesh_Curvature_C6");
		HGP_3D_Triangle_Mesh_Boundary_C1_C = (HGP_3D_Triangle_Mesh_Boundary_C1)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Boundary_C1");
		HGP_3D_Triangle_Mesh_Boundary_C2_C = (HGP_3D_Triangle_Mesh_Boundary_C2)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Boundary_C2");
		HGP_3D_Connecting_Segments_C1_C = (HGP_3D_Connecting_Segments_C1)GetProcAddress(hModule, "HGP_3D_Connecting_Segments_C1");
		HGP_3D_Connecting_Segments_C2_C = (HGP_3D_Connecting_Segments_C2)GetProcAddress(hModule, "HGP_3D_Connecting_Segments_C2");
		HGP_3D_Triangle_Mesh_Boundary_C3_C = (HGP_3D_Triangle_Mesh_Boundary_C3)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Boundary_C3");
		HGP_3D_Triangle_Mesh_Boundary_C4_C = (HGP_3D_Triangle_Mesh_Boundary_C4)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Boundary_C4");
		HGP_3D_Triangle_Mesh_Boundary_C5_C = (HGP_3D_Triangle_Mesh_Boundary_C5)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Boundary_C5");
		HGP_Mesh_Laplace_Smooth_C1_C = (HGP_Mesh_Laplace_Smooth_C1)GetProcAddress(hModule, "HGP_Mesh_Laplace_Smooth_C1");
		HGP_3D_Triangle_Mesh_Vecs_Neighbors_C = (HGP_3D_Triangle_Mesh_Vecs_Neighbors)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Vecs_Neighbors");
		HGP_Mesh_Laplace_Smooth_C2_C = (HGP_Mesh_Laplace_Smooth_C2)GetProcAddress(hModule, "HGP_Mesh_Laplace_Smooth_C2");
		HGP_3D_Triangle_Mesh_Vecs_Faces_C = (HGP_3D_Triangle_Mesh_Vecs_Faces)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Vecs_Faces");
		HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges_C = (HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges");
		HGP_Mesh_Laplace_Smooth_by_Curvature_C = (HGP_Mesh_Laplace_Smooth_by_Curvature)GetProcAddress(hModule, "HGP_Mesh_Laplace_Smooth_by_Curvature");
		HGP_Mesh_Loop_Subdivision_Own_Version_C = (HGP_Mesh_Loop_Subdivision_Own_Version)GetProcAddress(hModule, "HGP_Mesh_Loop_Subdivision_Own_Version");
		HGP_Rotation_Obj_C = (HGP_Rotation_Obj)GetProcAddress(hModule, "HGP_Rotation_Obj");
		HGP_Slicer_Mesh_C = (HGP_Slicer_Mesh)GetProcAddress(hModule, "HGP_Slicer_Mesh");
		HGP_Shortest_Geodesic_Path_C1_C = (HGP_Shortest_Geodesic_Path_C1)GetProcAddress(hModule, "HGP_Shortest_Geodesic_Path_C1");
		HGP_Shortest_Geodesic_Path_C3_C = (HGP_Shortest_Geodesic_Path_C3)GetProcAddress(hModule, "HGP_Shortest_Geodesic_Path_C3");
		HGP_Shortest_Geodesic_Path_C4_C = (HGP_Shortest_Geodesic_Path_C4)GetProcAddress(hModule, "HGP_Shortest_Geodesic_Path_C4");
		HGP_Geodesic_Distance_C = (HGP_Geodesic_Distance)GetProcAddress(hModule, "HGP_Geodesic_Distance");
		HGP_Project_Points_Onto_Surface_C1_C = (HGP_Project_Points_Onto_Surface_C1)GetProcAddress(hModule, "HGP_Project_Points_Onto_Surface_C1");
		HGP_Project_Points_Onto_Surface_C2_C = (HGP_Project_Points_Onto_Surface_C2)GetProcAddress(hModule, "HGP_Project_Points_Onto_Surface_C2");
		HGP_3D_Triangel_Mesh_Most_Inside_Point_C = (HGP_3D_Triangel_Mesh_Most_Inside_Point)GetProcAddress(hModule, "HGP_3D_Triangel_Mesh_Most_Inside_Point");
		HGP_3D_One_Triangle_Area_C = (HGP_3D_One_Triangle_Area)GetProcAddress(hModule, "HGP_3D_One_Triangle_Area");
		HGP_3D_Triangle_Mesh_Area_C = (HGP_3D_Triangle_Mesh_Area)GetProcAddress(hModule, "HGP_3D_Triangle_Mesh_Area");
		HGP_3D_Convex_Hulls_C1_C = (HGP_3D_Convex_Hulls_C1)GetProcAddress(hModule, "HGP_3D_Convex_Hulls_C1");
		HGP_3D_Convex_Hulls_C2_C = (HGP_3D_Convex_Hulls_C2)GetProcAddress(hModule, "HGP_3D_Convex_Hulls_C2");
		HGP_3D_Convex_Hulls_C3_C = (HGP_3D_Convex_Hulls_C3)GetProcAddress(hModule, "HGP_3D_Convex_Hulls_C3");
		HGP_3D_Convex_Hulls_C4_C = (HGP_3D_Convex_Hulls_C4)GetProcAddress(hModule, "HGP_3D_Convex_Hulls_C4");
		HGP_Mesh_Field_Query_C1_C = (HGP_Mesh_Field_Query_C1)GetProcAddress(hModule, "HGP_Mesh_Field_Query_C1");
		HGP_Mesh_Field_Query_C2_C = (HGP_Mesh_Field_Query_C2)GetProcAddress(hModule, "HGP_Mesh_Field_Query_C2");
		HGP_Mesh_Field_Query_C3_C = (HGP_Mesh_Field_Query_C3)GetProcAddress(hModule, "HGP_Mesh_Field_Query_C3");
		HGP_Curvature_Mesh_C = (HGP_Curvature_Mesh)GetProcAddress(hModule, "HGP_Curvature_Mesh");
		HGP_Normal_Mesh_C1_C = (HGP_Normal_Mesh_C1)GetProcAddress(hModule, "HGP_Normal_Mesh_C1");
		HGP_Normal_Mesh_C2_C = (HGP_Normal_Mesh_C2)GetProcAddress(hModule, "HGP_Normal_Mesh_C2");
		HGP_3D_Mesh_Normal_C1_C = (HGP_3D_Mesh_Normal_C1)GetProcAddress(hModule, "HGP_3D_Mesh_Normal_C1");
		HGP_3D_Mesh_Normal_C2_C = (HGP_3D_Mesh_Normal_C2)GetProcAddress(hModule, "HGP_3D_Mesh_Normal_C2");
		HGP_3D_Mesh_Center_C1_C = (HGP_3D_Mesh_Center_C1)GetProcAddress(hModule, "HGP_3D_Mesh_Center_C1");
		HGP_3D_Mesh_Center_C2_C = (HGP_3D_Mesh_Center_C2)GetProcAddress(hModule, "HGP_3D_Mesh_Center_C2");
		HGP_3D_Mesh_Boundingbox_C1_C = (HGP_3D_Mesh_Boundingbox_C1)GetProcAddress(hModule, "HGP_3D_Mesh_Boundingbox_C1");
		HGP_3D_Mesh_Boundingbox_C2_C = (HGP_3D_Mesh_Boundingbox_C2)GetProcAddress(hModule, "HGP_3D_Mesh_Boundingbox_C2");
		HGP_Surface_Decomposition_C = (HGP_Surface_Decomposition)GetProcAddress(hModule, "HGP_Surface_Decomposition");
		HGP_3D_Mesh_Gradient_C = (HGP_3D_Mesh_Gradient)GetProcAddress(hModule, "HGP_3D_Mesh_Gradient");
		HGP_Intergral_Curvature_C = (HGP_Intergral_Curvature)GetProcAddress(hModule, "HGP_Intergral_Curvature");
		HGP_3D_Mesh_Extract_Isoline_C = (HGP_3D_Mesh_Extract_Isoline)GetProcAddress(hModule, "HGP_3D_Mesh_Extract_Isoline");
		HGP_BSplineCurveFit_C = (HGP_BSplineCurveFit)GetProcAddress(hModule, "HGP_BSplineCurveFit");
		HGP_Cut_Surface_C = (HGP_Cut_Surface)GetProcAddress(hModule, "HGP_Cut_Surface");
		HGP_Cut_Surface_by_Multi_Boundaries_C = (HGP_Cut_Surface_by_Multi_Boundaries)GetProcAddress(hModule, "HGP_Cut_Surface_by_Multi_Boundaries");
		/////////////////////////////////////////////////////////////
		//
	};

	static HGPPL& Inst()
	{
		static HGPPL instance;
		return instance;
	};

	HMODULE hModule;
	HGP_Vector_Base HGP_Vector_Base_C;
	HGP_Test_PGL HGP_Test_PGL_C;
	//implementation in "hgp2d.cpp"
	//####################################################################################
	HGP_2D_Distance_Point_Point HGP_2D_Distance_Point_Point_C;
	HGP_2D_Distance_Point_Line HGP_2D_Distance_Point_Line_C;
	HGP_2D_Distance_Point_Segment HGP_2D_Distance_Point_Segment_C;
	HGP_2D_Distance_Segment_Segment HGP_2D_Distance_Segment_Segment_C;
	HGP_2D_Location_Point_Polygon HGP_2D_Location_Point_Polygon_C;
	HGP_2D_Location_Points_Polygon HGP_2D_Location_Points_Polygon_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_2D_Polygon_Dart_Sampling HGP_2D_Polygon_Dart_Sampling_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_2D_Polygon_Regular_Sampling_C1 HGP_2D_Polygon_Regular_Sampling_C1_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_2D_Polygon_Regular_Sampling_C2 HGP_2D_Polygon_Regular_Sampling_C2_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_2D_Polygon_Regular_Sampling_C3 HGP_2D_Polygon_Regular_Sampling_C3_C;
	HGP_2D_Square_Regular_Sampling_C1 HGP_2D_Square_Regular_Sampling_C1_C;
	HGP_2D_Square_Regular_Sampling_C2 HGP_2D_Square_Regular_Sampling_C2_C;
	HGP_2D_Square_Regular_Sampling_C3 HGP_2D_Square_Regular_Sampling_C3_C;
	HGP_2D_Distance_Point_Polygon HGP_2D_Distance_Point_Polygon_C;
	HGP_2D_Distance_Point_Polygons HGP_2D_Distance_Point_Polygons_C;
	HGP_2D_Intersection_Segment_Segment HGP_2D_Intersection_Segment_Segment_C;
	HGP_2D_Intersection_Line_Line HGP_2D_Intersection_Line_Line_C;
	HGP_2D_Intersection_Segment_Line HGP_2D_Intersection_Segment_Line_C;
	HGP_2D_Intersection_Segment_Polygon HGP_2D_Intersection_Segment_Polygon_C;
	HGP_2D_Intersection_Polygon_Polygon HGP_2D_Intersection_Polygon_Polygon_C;
	HGP_2D_Polygon_Is_Clockwise_Oriented HGP_2D_Polygon_Is_Clockwise_Oriented_C;
	HGP_2D_Two_Polygons_Union HGP_2D_Two_Polygons_Union_C;
	HGP_2D_Two_Polygons_Intersection HGP_2D_Two_Polygons_Intersection_C;
	HGP_Decompose_Polyline HGP_Decompose_Polyline_C;
	HGP_Identify_Polycut_Extend HGP_Identify_Polycut_Extend_C;
	HGP_Identify_Polycut_NotExtend HGP_Identify_Polycut_NotExtend_C;
	HGP_Identify_Polycut HGP_Identify_Polycut_C;
	HGP_Construct_InOutSide_Polygon HGP_Construct_InOutSide_Polygon_C;
	HGP_2D_Intersection_Ray_Segment HGP_2D_Intersection_Ray_Segment_C;
	HGP_Get_Angle_Kerf_Offset_Tan HGP_Get_Angle_Kerf_Offset_Tan_C;
	HGP_2D_Projection_Point_Segment HGP_2D_Projection_Point_Segment_C;
	HGP_2D_Detect_Polygon_Inside_C1 HGP_2D_Detect_Polygon_Inside_C1_C;
	HGP_2D_Detect_Polygon_Inside_C2 HGP_2D_Detect_Polygon_Inside_C2_C;
	HGP_2D_Detect_Polygon_Inside_C3 HGP_2D_Detect_Polygon_Inside_C3_C;
	HGP_2D_Detect_Polygon_Inside_C4 HGP_2D_Detect_Polygon_Inside_C4_C;
	HGP_2D_Detect_Polygon_Inside_C5 HGP_2D_Detect_Polygon_Inside_C5_C;
	HGP_2D_Distance_Polygon_Polygon HGP_2D_Distance_Polygon_Polygon_C;
	HGP_2D_Distance_Polygons_Polygons HGP_2D_Distance_Polygons_Polygons_C;
	HGP_2D_Nearest_Point_Polygon_C1 HGP_2D_Nearest_Point_Polygon_C1_C;
	HGP_2D_Nearest_Point_Polygon_C2 HGP_2D_Nearest_Point_Polygon_C2_C;
	HGP_2D_Nearest_Point_Polygons HGP_2D_Nearest_Point_Polygons_C;
	HGP_2d_Polygon_Boundingbox HGP_2d_Polygon_Boundingbox_C;
	HGP_2D_Polygon_Area HGP_2D_Polygon_Area_C;
	HGP_2D_Polygon_Inside_Point_C1 HGP_2D_Polygon_Inside_Point_C1_C;
	HGP_2D_Polygon_Inside_Point_C2 HGP_2D_Polygon_Inside_Point_C2_C;
	HGP_2D_Polygon_One_Offsets HGP_2D_Polygon_One_Offsets_C;
	HGP_2D_Polygons_One_Offsets HGP_2D_Polygons_One_Offsets_C;
	HGP_2D_Polygons_Simple HGP_2D_Polygons_Simple_C;
	HGP_2D_Polygon_Simple HGP_2D_Polygon_Simple_C;
	HGP_2D_Polygon_Simple_Inter HGP_2D_Polygon_Simple_Inter_C;
	HGP_2D_Convex_Hulls HGP_2D_Convex_Hulls_C;
	HGP_2D_OBB_Box HGP_2D_OBB_Box_C;
	HGP_Image_Grid_Decomposition_C1 HGP_Image_Grid_Decomposition_C1_C;
	HGP_Image_Grid_Decomposition_Conservative_C1 HGP_Image_Grid_Decomposition_Conservative_C1_C;
	HGP_Image_Grid_Decomposition_C2 HGP_Image_Grid_Decomposition_C2_C;
	HGP_Image_Grid_Decomposition_Conservative_C2 HGP_Image_Grid_Decomposition_Conservative_C2_C;
	//implementation in "hgp3d.cpp"
	//####################################################################################
	HGP_3D_Distance_Point_Segment HGP_3D_Distance_Point_Segment_C;
	HGP_3D_Plane_Fitting HGP_3D_Plane_Fitting_C;
	HGP_3D_Plane_Point_Projection HGP_3D_Plane_Point_Projection_C;
	HGP_3D_Plane_Points_Projection HGP_3D_Plane_Points_Projection_C;
	HGP_3D_Plane_3D_to_2D_Point HGP_3D_Plane_3D_to_2D_Point_C;
	HGP_3D_Plane_2D_to_3D_Point HGP_3D_Plane_2D_to_3D_Point_C;
	HGP_3D_Plane_3D_to_2D_Points HGP_3D_Plane_3D_to_2D_Points_C;
	HGP_3D_Plane_3Ds_to_2Ds_Points HGP_3D_Plane_3Ds_to_2Ds_Points_C;
	HGP_3D_Plane_2D_to_3D_Points HGP_3D_Plane_2D_to_3D_Points_C;
	HGP_3D_Projection_Point_Segment HGP_3D_Projection_Point_Segment_C;
	HGP_3D_Distance_Point_Point HGP_3D_Distance_Point_Point_C;
	HGP_3D_Distance_Point_Polygon HGP_3D_Distance_Point_Polygon_C;
	HGP_2D_Polygon_Triangulation_C1 HGP_2D_Polygon_Triangulation_C1_C;
	HGP_2D_Polygon_Triangulation_C2 HGP_2D_Polygon_Triangulation_C2_C;
	HGP_2D_Polygon_Triangulation_C3 HGP_2D_Polygon_Triangulation_C3_C;
	HGP_3D_Distance_Point_Line HGP_3D_Distance_Point_Line_C;
	HGP_3D_Projection_Point_Line HGP_3D_Projection_Point_Line_C;
	HGP_3D_Distance_Segment_Segment HGP_3D_Distance_Segment_Segment_C;
	HGP_3D_Distance_Point_Plane HGP_3D_Distance_Point_Plane_C;
	HGP_3D_Intersection_Segment_Line HGP_3D_Intersection_Segment_Line_C;
	HGP_3D_Intersection_Segment_Segment HGP_3D_Intersection_Segment_Segment_C;
	HGP_3D_Intersection_Segment_Plane HGP_3D_Intersection_Segment_Plane_C;
	HGP_3D_Intersection_Line_Plane HGP_3D_Intersection_Line_Plane_C;
	HGP_3D_Projection_Point_Plane_C1 HGP_3D_Projection_Point_Plane_C1_C;
	HGP_3D_Projection_Point_Plane_C2 HGP_3D_Projection_Point_Plane_C2_C;
	HGP_3D_Projection_3D_Point_Plane_2D_C1 HGP_3D_Projection_3D_Point_Plane_2D_C1_C;
	HGP_3D_Projection_3D_Point_Plane_2D_C2 HGP_3D_Projection_3D_Point_Plane_2D_C2_C;
	HGP_3D_Plane_ABCD HGP_3D_Plane_ABCD_C;
	HGP_3D_Plane_Base_1 HGP_3D_Plane_Base_1_C;
	HGP_Face_Normal HGP_Face_Normal_C;
	//implementation in "hgpmesh.cpp"
	//####################################################################################
	HGP_Remesh_Surface_by_Adding_Feature HGP_Remesh_Surface_by_Adding_Feature_C;
	HGP_Mesh_Edges HGP_Mesh_Edges_C;
	HGP_3D_Intersection_Sphere_Ray HGP_3D_Intersection_Sphere_Ray_C;
	HGP_3D_Intersection_Ray_Triangle HGP_3D_Intersection_Ray_Triangle_C;
	HGP_3D_Intersection_Ray_Mesh HGP_3D_Intersection_Ray_Mesh_C;
	HGP_3D_Intersection_Segment_Mesh HGP_3D_Intersection_Segment_Mesh_C;
	HGP_3D_Intersection_Segments_Mesh HGP_3D_Intersection_Segments_Mesh_C;
	HGP_3D_Intersection_Polygons_Mesh HGP_3D_Intersection_Polygons_Mesh_C;
	//check whether there is a polygon intersected with the input mesh
	HGP_3D_Intersection_Polygons_Mesh_Bool HGP_3D_Intersection_Polygons_Mesh_Bool_C;
	HGP_3D_Intersection_Rays_Mesh_Vector3d HGP_3D_Intersection_Rays_Mesh_Vector3d_C;
	//test each group directions (nes[i]) for each point in ps
	HGP_3D_Intersection_Rays_Mesh_C1_Bool HGP_3D_Intersection_Rays_Mesh_C1_Bool_C;
	//test all directions (ns) for each point in ps
	HGP_3D_Intersection_Rays_Mesh_C2_Bool HGP_3D_Intersection_Rays_Mesh_C2_Bool_C;
	HGP_3D_Intersection_Rays_Mesh_C2_Vector3d HGP_3D_Intersection_Rays_Mesh_C2_Vector3d_C;
	HGP_3D_Points_Inside_Triangles_C1_Bool HGP_3D_Points_Inside_Triangles_C1_Bool_C;
	HGP_3D_Points_Inside_Triangles_C2_Bool HGP_3D_Points_Inside_Triangles_C2_Bool_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_3D_Mesh_Dart_Sampling_C1 HGP_3D_Mesh_Dart_Sampling_C1_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_3D_Mesh_Dart_Sampling_C2 HGP_3D_Mesh_Dart_Sampling_C2_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_3D_Mesh_Regular_Sampling_C1 HGP_3D_Mesh_Regular_Sampling_C1_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_3D_Mesh_Regular_Sampling_C2 HGP_3D_Mesh_Regular_Sampling_C2_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_3D_Cube_Surface_Sampling_C1 HGP_3D_Cube_Surface_Sampling_C1_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_3D_Cube_Surface_Sampling_C2 HGP_3D_Cube_Surface_Sampling_C2_C;
	//d: percentage value of the length of the diagonal of the bounding box.
	HGP_3D_Cube_Surface_Sampling_C3 HGP_3D_Cube_Surface_Sampling_C3_C;
	//with neighboring
	HGP_3D_Mesh_Regular_Sampling_C3 HGP_3D_Mesh_Regular_Sampling_C3_C;
	HGP_3D_Distance_Point_Triangle HGP_3D_Distance_Point_Triangle_C;
	HGP_3D_Distance_Point_Triangles HGP_3D_Distance_Point_Triangles_C;
	HGP_3D_Nearest_Point_Triangles HGP_3D_Nearest_Point_Triangles_C;
	HGP_3D_Distance_Point_Mesh HGP_3D_Distance_Point_Mesh_C;
	HGP_3D_Neareast_Point_Mesh HGP_3D_Neareast_Point_Mesh_C;
	HGP_3D_Mesh_Near_Triangles HGP_3D_Mesh_Near_Triangles_C;
	HGP_3D_Points_inside_Triangles_C1 HGP_3D_Points_inside_Triangles_C1_C;
	HGP_3D_Points_inside_Triangles_C2 HGP_3D_Points_inside_Triangles_C2_C;
	HGP_Mesh_Subdivision HGP_Mesh_Subdivision_C;
	HGP_Mesh_Loop_Subdivision_One_Step HGP_Mesh_Loop_Subdivision_One_Step_C;
	HGP_3D_Mesh_Curvature_C1 HGP_3D_Mesh_Curvature_C1_C;
	HGP_3D_Mesh_Curvature_C2 HGP_3D_Mesh_Curvature_C2_C;
	HGP_3D_Mesh_Curvature_C3 HGP_3D_Mesh_Curvature_C3_C;
	HGP_3D_Mesh_Curvature_C4 HGP_3D_Mesh_Curvature_C4_C;
	HGP_3D_Mesh_Curvature_C5 HGP_3D_Mesh_Curvature_C5_C;
	HGP_3D_Mesh_Curvature_C6 HGP_3D_Mesh_Curvature_C6_C;
	HGP_3D_Triangle_Mesh_Boundary_C1 HGP_3D_Triangle_Mesh_Boundary_C1_C;
	HGP_3D_Triangle_Mesh_Boundary_C2 HGP_3D_Triangle_Mesh_Boundary_C2_C;
	HGP_3D_Connecting_Segments_C1 HGP_3D_Connecting_Segments_C1_C;
	HGP_3D_Connecting_Segments_C2 HGP_3D_Connecting_Segments_C2_C;
	HGP_3D_Triangle_Mesh_Boundary_C3 HGP_3D_Triangle_Mesh_Boundary_C3_C;
	HGP_3D_Triangle_Mesh_Boundary_C4 HGP_3D_Triangle_Mesh_Boundary_C4_C;
	HGP_3D_Triangle_Mesh_Boundary_C5 HGP_3D_Triangle_Mesh_Boundary_C5_C;
	HGP_Mesh_Laplace_Smooth_C1 HGP_Mesh_Laplace_Smooth_C1_C;
	HGP_3D_Triangle_Mesh_Vecs_Neighbors HGP_3D_Triangle_Mesh_Vecs_Neighbors_C;
	HGP_Mesh_Laplace_Smooth_C2 HGP_Mesh_Laplace_Smooth_C2_C;
	HGP_3D_Triangle_Mesh_Vecs_Faces HGP_3D_Triangle_Mesh_Vecs_Faces_C;
	HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges_C;
	HGP_Mesh_Laplace_Smooth_by_Curvature HGP_Mesh_Laplace_Smooth_by_Curvature_C;
	HGP_Mesh_Loop_Subdivision_Own_Version HGP_Mesh_Loop_Subdivision_Own_Version_C;
	HGP_Rotation_Obj HGP_Rotation_Obj_C;
	HGP_Slicer_Mesh HGP_Slicer_Mesh_C;
	HGP_Shortest_Geodesic_Path_C1 HGP_Shortest_Geodesic_Path_C1_C;
	HGP_Shortest_Geodesic_Path_C3 HGP_Shortest_Geodesic_Path_C3_C;
	HGP_Shortest_Geodesic_Path_C4 HGP_Shortest_Geodesic_Path_C4_C;
	HGP_Geodesic_Distance HGP_Geodesic_Distance_C;
	HGP_Project_Points_Onto_Surface_C1 HGP_Project_Points_Onto_Surface_C1_C;
	HGP_Project_Points_Onto_Surface_C2 HGP_Project_Points_Onto_Surface_C2_C;
	HGP_3D_Triangel_Mesh_Most_Inside_Point HGP_3D_Triangel_Mesh_Most_Inside_Point_C;
	HGP_3D_One_Triangle_Area HGP_3D_One_Triangle_Area_C;
	HGP_3D_Triangle_Mesh_Area HGP_3D_Triangle_Mesh_Area_C;
	HGP_3D_Convex_Hulls_C1 HGP_3D_Convex_Hulls_C1_C;
	HGP_3D_Convex_Hulls_C2 HGP_3D_Convex_Hulls_C2_C;
	HGP_3D_Convex_Hulls_C3 HGP_3D_Convex_Hulls_C3_C;
	HGP_3D_Convex_Hulls_C4 HGP_3D_Convex_Hulls_C4_C;
	HGP_Mesh_Field_Query_C1 HGP_Mesh_Field_Query_C1_C;
	HGP_Mesh_Field_Query_C2 HGP_Mesh_Field_Query_C2_C;
	HGP_Mesh_Field_Query_C3 HGP_Mesh_Field_Query_C3_C;
	HGP_Curvature_Mesh HGP_Curvature_Mesh_C;
	HGP_Normal_Mesh_C1 HGP_Normal_Mesh_C1_C;
	HGP_Normal_Mesh_C2 HGP_Normal_Mesh_C2_C;
	HGP_3D_Mesh_Normal_C1 HGP_3D_Mesh_Normal_C1_C;
	HGP_3D_Mesh_Normal_C2 HGP_3D_Mesh_Normal_C2_C;
	HGP_3D_Mesh_Center_C1 HGP_3D_Mesh_Center_C1_C;
	HGP_3D_Mesh_Center_C2 HGP_3D_Mesh_Center_C2_C;
	HGP_3D_Mesh_Boundingbox_C1 HGP_3D_Mesh_Boundingbox_C1_C;
	HGP_3D_Mesh_Boundingbox_C2 HGP_3D_Mesh_Boundingbox_C2_C;
	HGP_Surface_Decomposition HGP_Surface_Decomposition_C;
	HGP_3D_Mesh_Gradient HGP_3D_Mesh_Gradient_C;
	HGP_Intergral_Curvature HGP_Intergral_Curvature_C;
	HGP_3D_Mesh_Extract_Isoline HGP_3D_Mesh_Extract_Isoline_C;
	HGP_BSplineCurveFit HGP_BSplineCurveFit_C;
	HGP_Cut_Surface HGP_Cut_Surface_C;
	HGP_Cut_Surface_by_Multi_Boundaries HGP_Cut_Surface_by_Multi_Boundaries_C;
	/////////////////////////////////////////////////////////////
	//
};
#define PL() HGPPL::Inst()
}
#endif
