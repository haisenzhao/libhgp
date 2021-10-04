#ifndef CGAL_ONCE
#define CGAL_ONCE
#pragma once
#include <pgl_functs.hpp>
using namespace std;
using namespace PGL;
typedef  void (*CGAL_Test_PGL)(Vector3d n);
typedef  void (*CGAL_Vector_Base)(Vector3d n, Vector3d &);
typedef  void (*CGAL_Export_Path_Segment)(std::ofstream &export_file_output, int &export_index,std::string s_name, double r, double g, double b,Vector3d &start,Vector3d &end, double radius);
typedef  void (*CGAL_Export_Path_Point)(std::ofstream &export_file_output, int &export_index,std::string s_name, double r, double g, double b, Vector3d point, double radius);
typedef  double (*CGAL_2D_Distance_Point_Point)(Vector2d p_0, Vector2d p_1);
typedef  double (*CGAL_2D_Distance_Point_Line)(Vector2d v, Vector2d l_0, Vector2d l_1);
typedef  double (*CGAL_2D_Distance_Point_Segment)(Vector2d v, Vector2d s_0, Vector2d s_1);
typedef  double (*CGAL_2D_Distance_Segment_Segment)(Vector2d s_0, Vector2d s_1, Vector2d e_0, Vector2d e_1);
typedef  bool (*CGAL_2D_Location_Point_Polygon)(Vector2d p, Vector2d1 py);
typedef  bool (*CGAL_2D_Location_Points_Polygon)(const Vector2d1 &ps, const Vector2d1 &py);
typedef  double (*CGAL_2D_Distance_Point_Polygon)(Vector2d p, Vector2d1 py);
typedef  bool (*CGAL_2D_Intersection_Segment_Segment)(Vector2d s_0_s, Vector2d s_0_e, Vector2d s_1_s, Vector2d s_1_e, Vector2d &inter);
typedef  bool (*CGAL_2D_Intersection_Line_Line)(const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s, const Vector2d &s_1_e, Vector2d &inter);
typedef  bool (*CGAL_2D_Intersection_Segment_Polygon)(Vector2d s_s, Vector2d s_e, Vector2d1 &p);
typedef  bool (*CGAL_2D_Polygon_Is_Clockwise_Oriented)(Vector2d1 &ps);
typedef  double (*CGAL_2D_Two_Polygons_Union)(Vector2d1 poly_0, Vector2d1 poly_1,Vector2d2 &inter_polygons);
typedef  double (*CGAL_2D_Two_Polygons_Intersection)(const Vector2d1 &poly_0, const Vector2d1 &poly_1);
typedef  Vector1i1 (*CGAL_Decompose_Polyline)(Vector2d1 &polyline, double threshold);
typedef  bool (*CGAL_Identify_Polycut_Extend)(const Vector2d1 &polygon, const Vector2d &s,const Vector2d &e, Vector2d &ns, Vector2d &ne);
typedef  bool (*CGAL_Identify_Polycut_NotExtend)(const Vector2d1 &polygon, const Vector2d &s,const Vector2d &e);
typedef  bool (*CGAL_Identify_Polycut)(const Vector2d1 &polygon, const Vector2d1 &cutLine, VectorPB1 &result);
typedef  void (*CGAL_2D_Polygon_One_Offsets)(Vector2d1 &poly, double d,Vector2d2  &offset_polys);
typedef  bool (*CGAL_Construct_InOutSide_Polygon)(const Vector2d1 &py, const Vector2d &p, const Vector2d &q, bool &isPInside,bool &isQInside);
typedef  bool (*CGAL_2D_Intersection_Ray_Segment)(const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s,const Vector2d &s_1_e, Vector2d &inter);
typedef  double (*CGAL_Get_Angle_Kerf_Offset_Tan)(const Vector2d &a, const Vector2d &b);
typedef  double (*CGAL_3D_Distance_Point_Segment)(Vector3d p, Vector3d s_s, Vector3d s_e);
typedef  void (*CGAL_3D_Plane_Fitting)(Vector3d1 &points, Vector3d &plane_p, Vector3d &plane_n);
typedef  void (*CGAL_3D_Plane_Point_Projection)(Vector3d &plane_p, Vector3d &plane_n, Vector3d &p,Vector3d &result);
typedef  void (*CGAL_3D_Plane_Points_Projection)(Vector3d &plane_p, Vector3d &plane_n, Vector3d1 &points, Vector3d1 &project_points);
typedef  void (*CGAL_3D_Plane_3D_to_2D_Point)(Vector3d &plane_p, Vector3d &plane_n, Vector3d &point_3d,Vector2d &result);
typedef  void (*CGAL_3D_Plane_2D_to_3D_Point)(Vector3d &plane_p, Vector3d &plane_n, Vector2d &points_2d, Vector3d &result);
typedef  void (*CGAL_3D_Plane_3D_to_2D_Points)(Vector3d &plane_p, Vector3d &plane_n, Vector3d1 &points_3d,Vector2d1 &points_2d);
typedef  void (*CGAL_3D_Plane_2D_to_3D_Points)(Vector3d &plane_p, Vector3d &plane_n, Vector2d1 &points_2d,Vector3d1 &points_3d);
typedef  Vector3d (*CGAL_3D_Projection_Point_Segment)(Vector3d p, Vector3d s_s, Vector3d s_e);
typedef  double (*CGAL_3D_Distance_Point_Point)(const Vector3d & v0, const Vector3d & v1);
typedef  double (*CGAL_3D_Distance_Point_Polygon)(const Vector3d1 &py, const Vector3d &p);
typedef  void (*CGAL_2D_Polygon_Triangulation)(const Vector2d2 &polys, Vector1i2 &faces);
typedef  void (*CGAL_Remesh_Surface_by_Adding_Feature)(const Vector3d1 &feature,const Vector1i1 &face_ids, const Vector3d1 &vecs, const Vector1i1 &face_id_0,const Vector1i1 &face_id_1,const Vector1i1 &face_id_2, Vector1i1 &igl_cutting_0_edges,Vector1i1 &igl_cutting_1_edges, Vector3d1 &igl_cutting_points,Vector1i2 &cutting_faces);
typedef  void (*CGAL_3D_Read_Triangle_Mesh)(std::string path, Vector3d1 &vecs,Vector1i1 &face_id_0, Vector1i1 &face_id_1, Vector1i1 &face_id_2);
typedef  void (*CGAL_Mesh_Edges)(std::string path);
typedef  bool (*CGAL_3D_Intersection_Sphere_Ray)(double, double, double, double, double, double, double, double, double, double, std::vector<double>&, std::vector<double>&, std::vector<double>&);
typedef  bool (*CGAL_3D_Intersection_Ray_Triangle)(Vector3d p, Vector3d n, Vector3d p0, Vector3d p1, Vector3d p2);
typedef  bool (*CGAL_3D_Intersection_Ray_Mesh)(Vector3d p, Vector3d n, std::string path);
typedef  void (*CGAL_3D_Intersection_Rays_Mesh_Vector3d)(Vector3d1 ps, Vector3d1 ns, std::string path, Vector3d1 & inters);
typedef  void (*CGAL_3D_Intersection_Rays_Mesh_C1_Bool)(Vector3d1 ps, Vector3d2 nes, std::string path, Vector1b2 & inters);
typedef  void (*CGAL_3D_Intersection_Rays_Mesh_C2_Bool)(Vector3d1 ps, Vector3d1 ns, std::string path, Vector1b2 & inters);


class PL
{
	public:
	PL()
	{
		hModule = Functs::LoadHMODULE("ppgl.dll");
		CGAL_Test_PGL_C = (CGAL_Test_PGL)GetProcAddress(hModule, "CGAL_Test_PGL");
		CGAL_Vector_Base_C = (CGAL_Vector_Base)GetProcAddress(hModule, "CGAL_Vector_Base");
		CGAL_Export_Path_Segment_C = (CGAL_Export_Path_Segment)GetProcAddress(hModule, "CGAL_Export_Path_Segment");
		CGAL_Export_Path_Point_C = (CGAL_Export_Path_Point)GetProcAddress(hModule, "CGAL_Export_Path_Point");
		CGAL_2D_Distance_Point_Point_C = (CGAL_2D_Distance_Point_Point)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Point");
		CGAL_2D_Distance_Point_Line_C = (CGAL_2D_Distance_Point_Line)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Line");
		CGAL_2D_Distance_Point_Segment_C = (CGAL_2D_Distance_Point_Segment)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Segment");
		CGAL_2D_Distance_Segment_Segment_C = (CGAL_2D_Distance_Segment_Segment)GetProcAddress(hModule, "CGAL_2D_Distance_Segment_Segment");
		CGAL_2D_Location_Point_Polygon_C = (CGAL_2D_Location_Point_Polygon)GetProcAddress(hModule, "CGAL_2D_Location_Point_Polygon");
		CGAL_2D_Location_Points_Polygon_C = (CGAL_2D_Location_Points_Polygon)GetProcAddress(hModule, "CGAL_2D_Location_Points_Polygon");
		CGAL_2D_Distance_Point_Polygon_C = (CGAL_2D_Distance_Point_Polygon)GetProcAddress(hModule, "CGAL_2D_Distance_Point_Polygon");
		CGAL_2D_Intersection_Segment_Segment_C = (CGAL_2D_Intersection_Segment_Segment)GetProcAddress(hModule, "CGAL_2D_Intersection_Segment_Segment");
		CGAL_2D_Intersection_Line_Line_C = (CGAL_2D_Intersection_Line_Line)GetProcAddress(hModule, "CGAL_2D_Intersection_Line_Line");
		CGAL_2D_Intersection_Segment_Polygon_C = (CGAL_2D_Intersection_Segment_Polygon)GetProcAddress(hModule, "CGAL_2D_Intersection_Segment_Polygon");
		CGAL_2D_Polygon_Is_Clockwise_Oriented_C = (CGAL_2D_Polygon_Is_Clockwise_Oriented)GetProcAddress(hModule, "CGAL_2D_Polygon_Is_Clockwise_Oriented");
		CGAL_2D_Two_Polygons_Union_C = (CGAL_2D_Two_Polygons_Union)GetProcAddress(hModule, "CGAL_2D_Two_Polygons_Union");
		CGAL_2D_Two_Polygons_Intersection_C = (CGAL_2D_Two_Polygons_Intersection)GetProcAddress(hModule, "CGAL_2D_Two_Polygons_Intersection");
		CGAL_Decompose_Polyline_C = (CGAL_Decompose_Polyline)GetProcAddress(hModule, "CGAL_Decompose_Polyline");
		CGAL_Identify_Polycut_Extend_C = (CGAL_Identify_Polycut_Extend)GetProcAddress(hModule, "CGAL_Identify_Polycut_Extend");
		CGAL_Identify_Polycut_NotExtend_C = (CGAL_Identify_Polycut_NotExtend)GetProcAddress(hModule, "CGAL_Identify_Polycut_NotExtend");
		CGAL_Identify_Polycut_C = (CGAL_Identify_Polycut)GetProcAddress(hModule, "CGAL_Identify_Polycut");
		CGAL_2D_Polygon_One_Offsets_C = (CGAL_2D_Polygon_One_Offsets)GetProcAddress(hModule, "CGAL_2D_Polygon_One_Offsets");
		CGAL_Construct_InOutSide_Polygon_C = (CGAL_Construct_InOutSide_Polygon)GetProcAddress(hModule, "CGAL_Construct_InOutSide_Polygon");
		CGAL_2D_Intersection_Ray_Segment_C = (CGAL_2D_Intersection_Ray_Segment)GetProcAddress(hModule, "CGAL_2D_Intersection_Ray_Segment");
		CGAL_Get_Angle_Kerf_Offset_Tan_C = (CGAL_Get_Angle_Kerf_Offset_Tan)GetProcAddress(hModule, "CGAL_Get_Angle_Kerf_Offset_Tan");
		CGAL_3D_Distance_Point_Segment_C = (CGAL_3D_Distance_Point_Segment)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Segment");
		CGAL_3D_Plane_Fitting_C = (CGAL_3D_Plane_Fitting)GetProcAddress(hModule, "CGAL_3D_Plane_Fitting");
		CGAL_3D_Plane_Point_Projection_C = (CGAL_3D_Plane_Point_Projection)GetProcAddress(hModule, "CGAL_3D_Plane_Point_Projection");
		CGAL_3D_Plane_Points_Projection_C = (CGAL_3D_Plane_Points_Projection)GetProcAddress(hModule, "CGAL_3D_Plane_Points_Projection");
		CGAL_3D_Plane_3D_to_2D_Point_C = (CGAL_3D_Plane_3D_to_2D_Point)GetProcAddress(hModule, "CGAL_3D_Plane_3D_to_2D_Point");
		CGAL_3D_Plane_2D_to_3D_Point_C = (CGAL_3D_Plane_2D_to_3D_Point)GetProcAddress(hModule, "CGAL_3D_Plane_2D_to_3D_Point");
		CGAL_3D_Plane_3D_to_2D_Points_C = (CGAL_3D_Plane_3D_to_2D_Points)GetProcAddress(hModule, "CGAL_3D_Plane_3D_to_2D_Points");
		CGAL_3D_Plane_2D_to_3D_Points_C = (CGAL_3D_Plane_2D_to_3D_Points)GetProcAddress(hModule, "CGAL_3D_Plane_2D_to_3D_Points");
		CGAL_3D_Projection_Point_Segment_C = (CGAL_3D_Projection_Point_Segment)GetProcAddress(hModule, "CGAL_3D_Projection_Point_Segment");
		CGAL_3D_Distance_Point_Point_C = (CGAL_3D_Distance_Point_Point)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Point");
		CGAL_3D_Distance_Point_Polygon_C = (CGAL_3D_Distance_Point_Polygon)GetProcAddress(hModule, "CGAL_3D_Distance_Point_Polygon");
		CGAL_2D_Polygon_Triangulation_C = (CGAL_2D_Polygon_Triangulation)GetProcAddress(hModule, "CGAL_2D_Polygon_Triangulation");
		CGAL_Remesh_Surface_by_Adding_Feature_C = (CGAL_Remesh_Surface_by_Adding_Feature)GetProcAddress(hModule, "CGAL_Remesh_Surface_by_Adding_Feature");
		CGAL_3D_Read_Triangle_Mesh_C = (CGAL_3D_Read_Triangle_Mesh)GetProcAddress(hModule, "CGAL_3D_Read_Triangle_Mesh");
		CGAL_Mesh_Edges_C = (CGAL_Mesh_Edges)GetProcAddress(hModule, "CGAL_Mesh_Edges");
		CGAL_3D_Intersection_Sphere_Ray_C = (CGAL_3D_Intersection_Sphere_Ray)GetProcAddress(hModule, "CGAL_3D_Intersection_Sphere_Ray");
		CGAL_3D_Intersection_Ray_Triangle_C = (CGAL_3D_Intersection_Ray_Triangle)GetProcAddress(hModule, "CGAL_3D_Intersection_Ray_Triangle");
		CGAL_3D_Intersection_Ray_Mesh_C = (CGAL_3D_Intersection_Ray_Mesh)GetProcAddress(hModule, "CGAL_3D_Intersection_Ray_Mesh");
		CGAL_3D_Intersection_Rays_Mesh_Vector3d_C = (CGAL_3D_Intersection_Rays_Mesh_Vector3d)GetProcAddress(hModule, "CGAL_3D_Intersection_Rays_Mesh_Vector3d");
		CGAL_3D_Intersection_Rays_Mesh_C1_Bool_C = (CGAL_3D_Intersection_Rays_Mesh_C1_Bool)GetProcAddress(hModule, "CGAL_3D_Intersection_Rays_Mesh_C1_Bool");
		CGAL_3D_Intersection_Rays_Mesh_C2_Bool_C = (CGAL_3D_Intersection_Rays_Mesh_C2_Bool)GetProcAddress(hModule, "CGAL_3D_Intersection_Rays_Mesh_C2_Bool");
	};

	static PL& Instance()
	{
		static PL instance;
		return instance;
	};

	HMODULE hModule;
	CGAL_Test_PGL CGAL_Test_PGL_C;
	CGAL_Vector_Base CGAL_Vector_Base_C;
	CGAL_Export_Path_Segment CGAL_Export_Path_Segment_C;
	CGAL_Export_Path_Point CGAL_Export_Path_Point_C;
	CGAL_2D_Distance_Point_Point CGAL_2D_Distance_Point_Point_C;
	CGAL_2D_Distance_Point_Line CGAL_2D_Distance_Point_Line_C;
	CGAL_2D_Distance_Point_Segment CGAL_2D_Distance_Point_Segment_C;
	CGAL_2D_Distance_Segment_Segment CGAL_2D_Distance_Segment_Segment_C;
	CGAL_2D_Location_Point_Polygon CGAL_2D_Location_Point_Polygon_C;
	CGAL_2D_Location_Points_Polygon CGAL_2D_Location_Points_Polygon_C;
	CGAL_2D_Distance_Point_Polygon CGAL_2D_Distance_Point_Polygon_C;
	CGAL_2D_Intersection_Segment_Segment CGAL_2D_Intersection_Segment_Segment_C;
	CGAL_2D_Intersection_Line_Line CGAL_2D_Intersection_Line_Line_C;
	CGAL_2D_Intersection_Segment_Polygon CGAL_2D_Intersection_Segment_Polygon_C;
	CGAL_2D_Polygon_Is_Clockwise_Oriented CGAL_2D_Polygon_Is_Clockwise_Oriented_C;
	CGAL_2D_Two_Polygons_Union CGAL_2D_Two_Polygons_Union_C;
	CGAL_2D_Two_Polygons_Intersection CGAL_2D_Two_Polygons_Intersection_C;
	CGAL_Decompose_Polyline CGAL_Decompose_Polyline_C;
	CGAL_Identify_Polycut_Extend CGAL_Identify_Polycut_Extend_C;
	CGAL_Identify_Polycut_NotExtend CGAL_Identify_Polycut_NotExtend_C;
	CGAL_Identify_Polycut CGAL_Identify_Polycut_C;
	CGAL_2D_Polygon_One_Offsets CGAL_2D_Polygon_One_Offsets_C;
	CGAL_Construct_InOutSide_Polygon CGAL_Construct_InOutSide_Polygon_C;
	CGAL_2D_Intersection_Ray_Segment CGAL_2D_Intersection_Ray_Segment_C;
	CGAL_Get_Angle_Kerf_Offset_Tan CGAL_Get_Angle_Kerf_Offset_Tan_C;
	CGAL_3D_Distance_Point_Segment CGAL_3D_Distance_Point_Segment_C;
	CGAL_3D_Plane_Fitting CGAL_3D_Plane_Fitting_C;
	CGAL_3D_Plane_Point_Projection CGAL_3D_Plane_Point_Projection_C;
	CGAL_3D_Plane_Points_Projection CGAL_3D_Plane_Points_Projection_C;
	CGAL_3D_Plane_3D_to_2D_Point CGAL_3D_Plane_3D_to_2D_Point_C;
	CGAL_3D_Plane_2D_to_3D_Point CGAL_3D_Plane_2D_to_3D_Point_C;
	CGAL_3D_Plane_3D_to_2D_Points CGAL_3D_Plane_3D_to_2D_Points_C;
	CGAL_3D_Plane_2D_to_3D_Points CGAL_3D_Plane_2D_to_3D_Points_C;
	CGAL_3D_Projection_Point_Segment CGAL_3D_Projection_Point_Segment_C;
	CGAL_3D_Distance_Point_Point CGAL_3D_Distance_Point_Point_C;
	CGAL_3D_Distance_Point_Polygon CGAL_3D_Distance_Point_Polygon_C;
	CGAL_2D_Polygon_Triangulation CGAL_2D_Polygon_Triangulation_C;
	CGAL_Remesh_Surface_by_Adding_Feature CGAL_Remesh_Surface_by_Adding_Feature_C;
	CGAL_3D_Read_Triangle_Mesh CGAL_3D_Read_Triangle_Mesh_C;
	CGAL_Mesh_Edges CGAL_Mesh_Edges_C;
	CGAL_3D_Intersection_Sphere_Ray CGAL_3D_Intersection_Sphere_Ray_C;
	CGAL_3D_Intersection_Ray_Triangle CGAL_3D_Intersection_Ray_Triangle_C;
	CGAL_3D_Intersection_Ray_Mesh CGAL_3D_Intersection_Ray_Mesh_C;
	CGAL_3D_Intersection_Rays_Mesh_Vector3d CGAL_3D_Intersection_Rays_Mesh_Vector3d_C;
	CGAL_3D_Intersection_Rays_Mesh_C1_Bool CGAL_3D_Intersection_Rays_Mesh_C1_Bool_C;
	CGAL_3D_Intersection_Rays_Mesh_C2_Bool CGAL_3D_Intersection_Rays_Mesh_C2_Bool_C;
};
#endif
