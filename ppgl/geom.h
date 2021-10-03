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

using namespace PGL;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_2.h>


#include <CGAL/subdivision_method_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <CGAL/barycenter.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/intersections.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#include <cstdlib>
#include <iterator>
#include <CGAL/Random.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_set_2.h>



typedef CGAL::Simple_cartesian<double> KC;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<KC> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>  SDG2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Line_2 Line_2;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef K::Vector_2 Vector_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Ray_2 Ray_2;
typedef K::Direction_2 Direction_2;
typedef K::Line_3 Line_3;
typedef K::Plane_3 Plane_3;
typedef K::Vector_3 Vector_3;
typedef K::Segment_3 Segment_3;
typedef K::Direction_3 Direction_3;
typedef K::Ray_3 Ray_3;
typedef K::Triangle_3 Triangle;

typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
typedef Polyhedron_3::Facet_iterator Poly_facet_iterator;
typedef Polyhedron_3::Point_3 Poly_point_3;
typedef Polyhedron_3::HalfedgeDS Poly_halfedgeds_3;
typedef Polyhedron_3::Halfedge_handle Halfedge_handle;
typedef Polyhedron_3::Vertex_handle Vertex_handle;
typedef Polyhedron_3::Halfedge_around_vertex_const_circulator HV_circulator;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron_3> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits_poly;
typedef CGAL::AABB_tree<Traits_poly> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef boost::optional< Tree::Intersection_and_primitive_id<Ray_3>::Type> Ray_intersection;

typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef std::vector<K::Point_3> Polyline_type;
typedef std::list< Polyline_type > Polylines;
typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh> HGSP;
typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;


typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Mesh_Primitive;
typedef CGAL::AABB_traits<K, Mesh_Primitive> Mesh_Traits;
typedef CGAL::AABB_tree<Mesh_Traits> Mesh_Tree;
typedef boost::optional<Mesh_Tree::Intersection_and_primitive_id<Ray_3>::Type> Mesh_Ray_intersection;


typedef KC::Triangle_3 Triangle_3;
typedef std::list<Triangle_3>::iterator Iterator_3;
typedef CGAL::AABB_triangle_primitive<KC, Iterator_3> Primitive_3;
typedef CGAL::AABB_traits<KC, Primitive_3> AABB_triangle_traits_3;
typedef CGAL::AABB_tree<AABB_triangle_traits_3> Tree_3;



struct FaceInfo2
{
	FaceInfo2(){}
	int nesting_level;
	bool in_domain(){
		return nesting_level % 2 == 1;
	}
};


// ----------------------- A CGAL::Vertex with decoration ------------------
template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class Vertex : public  Vb {
	typedef Vb superclass;
public:
	typedef typename Vb::Vertex_handle      Vertex_handle;
	typedef typename Vb::Point              Point;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef Vertex<Gt, Vb2> Other;
	};

public:
	Vertex() : superclass() {}
	Vertex(const Point & p) : superclass(p) {}
	int index;
};

typedef CGAL::Triangulation_2<K>         Triangulation;
typedef Vertex<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;

template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:

    std::vector<double> coords;
   Vector1i1 tris;

    polyhedron_builder(std::vector<double> &c,Vector1i1 &t) {
        coords = c;
        tris = t;
    }

    void operator()(HDS &hds) override {
        typedef typename HDS::Vertex Vertex;
        typedef typename Vertex::Point Point;

        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
        B.begin_surface(coords.size() / 3, tris.size() / 3);

        for (int i = 0; i < (int) coords.size(); i += 3) {
            B.add_vertex(Point(coords[i + 0], coords[i + 1], coords[i + 2]));
        }

        for (int i = 0; i < (int) tris.size(); i += 3) {
            B.begin_facet();
            B.add_vertex_to_facet(tris[i + 0]);
            B.add_vertex_to_facet(tris[i + 1]);
            B.add_vertex_to_facet(tris[i + 2]);
            B.end_facet();
        }
    }
};

Vector3d PointVector3d(Point_3 p);
Point_3 VectorPoint3d(Vector3d p);
Vector2d PointVector2d(Point_2 p);
Point_2 VectorPoint2d(Vector2d p);

void  Construct_Polyhedron(Polyhedron_3& polyhedron, const Vector3d1& vecs, const Vector1i1& face_id_0, const Vector1i1& face_id_1, const Vector1i1& face_id_2);
void  Construct_Polyhedron(Polyhedron_3& polyhedron, std::string path);
void  Construct_Polyhedron(Polyhedron_3& polyhedron, std::string path, Vector3d1& vecs, Vector1i1& face_id_0, Vector1i1& face_id_1, Vector1i1& face_id_2);

extern "C" PPGL_EXPORT void Test_PGL(Vector3d n);


//implementation in "io.cpp"
//####################################################################################
extern "C" PPGL_EXPORT void CGAL_Vector_Base(Vector3d n, Vector3d &);
extern "C" PPGL_EXPORT void CGAL_Export_Path_Segment(std::ofstream &export_file_output, int &export_index,
                                                               std::string s_name, double r, double g, double b,
                                                               Vector3d &start,
                                                               Vector3d &end, double radius);
extern "C" PPGL_EXPORT void CGAL_Export_Path_Point(std::ofstream &export_file_output, int &export_index,
                                                             std::string s_name, double r, double g, double b,
                                                             Vector3d point,
                                                             double radius);

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
extern "C" PPGL_EXPORT double GetAngleKerfOffsetTan(const Vector2d &a, const Vector2d &b);


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

extern "C" PPGL_EXPORT void CGAL_Remesh_Surface_by_Adding_Feature(const Vector3d1 &feature,
                                                                            const Vector1i1 &face_ids,
                                                                            const Vector3d1 &vecs,
                                                                            const Vector1i1 &face_id_0,
                                                                            const Vector1i1 &face_id_1,
                                                                            const Vector1i1 &face_id_2,
                                                                           Vector1i1 &igl_cutting_0_edges,
                                                                           Vector1i1 &igl_cutting_1_edges,
                                                                            Vector3d1 &igl_cutting_points,
                                                                            Vector1i2 &cutting_faces);

extern "C" PPGL_EXPORT void CGAL_3D_Read_Triangle_Mesh(std::string path, Vector3d1 &vecs,
                                                                Vector1i1 &face_id_0,
                                                                Vector1i1 &face_id_1,
                                                                Vector1i1 &face_id_2);
extern "C" PPGL_EXPORT void CGAL_Mesh_Edges(std::string path);

extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Mesh(Vector3d p, Vector3d n, std::string path);
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh(Vector3d1 ps, Vector3d1 ns, std::string path, Vector3d1& inters);

#endif
