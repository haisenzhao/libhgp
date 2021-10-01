#ifndef mydll_hpp
#define mydll_hpp

#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <algorithm>
#include <list>
#include <carpentry_geom_export.h>
#include <math.hpp>

using namespace Math;

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
    std::vector<int> tris;

    polyhedron_builder(std::vector<double> &c, std::vector<int> &t) {
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

//implementation in "io.cpp"
//####################################################################################
extern "C" CARPENTRY_GEOM_EXPORT void CGAL_Vector_Base(Vector3d n, Vector3d &);
extern "C" CARPENTRY_GEOM_EXPORT void CGAL_Export_Path_Segment(std::ofstream &export_file_output, int &export_index,
                                                               std::string s_name, double r, double g, double b,
                                                               Vector3d &start,
                                                               Vector3d &end, double radius);
extern "C" CARPENTRY_GEOM_EXPORT void CGAL_Export_Path_Point(std::ofstream &export_file_output, int &export_index,
                                                             std::string s_name, double r, double g, double b,
                                                             Vector3d point,
                                                             double radius);

//implementation in "twoD.cpp"
//####################################################################################
extern "C" CARPENTRY_GEOM_EXPORT double CGAL_2D_Distance_Point_Point(Vector2d p_0, Vector2d p_1);
extern "C" CARPENTRY_GEOM_EXPORT double CGAL_2D_Distance_Point_Line(Vector2d v, Vector2d l_0, Vector2d l_1);
extern "C" CARPENTRY_GEOM_EXPORT double CGAL_2D_Distance_Point_Segment(Vector2d v, Vector2d s_0, Vector2d s_1);
extern "C" CARPENTRY_GEOM_EXPORT double
CGAL_2D_Distance_Segment_Segment(Vector2d s_0, Vector2d s_1, Vector2d e_0, Vector2d e_1);
extern "C" CARPENTRY_GEOM_EXPORT bool CGAL_2D_Location_Point_Polygon(Vector2d p, std::vector<Vector2d> py);
extern "C" CARPENTRY_GEOM_EXPORT bool CGAL_2D_Location_Points_Polygon(const std::vector<Vector2d> &ps,
                                                                      const std::vector<Vector2d> &py);
extern "C" CARPENTRY_GEOM_EXPORT double CGAL_2D_Distance_Point_Polygon(Vector2d p, std::vector<Vector2d> py);

extern "C" CARPENTRY_GEOM_EXPORT bool CGAL_2D_Intersection_Segment_Segment
        (Vector2d s_0_s, Vector2d s_0_e, Vector2d s_1_s, Vector2d s_1_e, Vector2d &inter);

extern "C" CARPENTRY_GEOM_EXPORT bool CGAL_2D_Intersection_Line_Line
        (const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s, const Vector2d &s_1_e, Vector2d &inter);

extern "C" CARPENTRY_GEOM_EXPORT bool
CGAL_2D_Intersection_Segment_Polygon(Vector2d s_s, Vector2d s_e, std::vector<Vector2d> &p);
extern "C" CARPENTRY_GEOM_EXPORT bool CGAL_2D_Polygon_Is_Clockwise_Oriented(std::vector<Vector2d> &ps);
extern "C" CARPENTRY_GEOM_EXPORT double
CGAL_2D_Two_Polygons_Union(std::vector<Vector2d> poly_0, std::vector<Vector2d> poly_1,
                           std::vector<std::vector<Vector2d> > &inter_polygons);

extern "C" CARPENTRY_GEOM_EXPORT double CGAL_2D_Two_Polygons_Intersection(const std::vector<Vector2d> &poly_0,
                                                                          const std::vector<Vector2d> &poly_1);

extern "C" CARPENTRY_GEOM_EXPORT std::vector<int>
CGAL_Decompose_Polyline(std::vector<Vector2d> &polyline, double threshold);
extern "C" CARPENTRY_GEOM_EXPORT bool
CGAL_Identify_Polycut_Extend(const std::vector<Vector2d> &polygon, const Vector2d &s,
                             const Vector2d &e, Vector2d &ns, Vector2d &ne);
extern "C" CARPENTRY_GEOM_EXPORT bool
CGAL_Identify_Polycut_NotExtend(const std::vector<Vector2d> &polygon, const Vector2d &s,
                                const Vector2d &e);
extern "C" CARPENTRY_GEOM_EXPORT bool CGAL_Identify_Polycut(const std::vector<Vector2d> &polygon,
                                                            const std::vector<Vector2d> &cutLine,
                                                            std::vector<std::pair<bool, bool> > &result);

extern "C" CARPENTRY_GEOM_EXPORT void CGAL_2D_Polygon_One_Offsets(std::vector<Vector2d> &poly,
                                                                  double d,
                                                                  std::vector<std::vector<Vector2d> > &offset_polys);

extern "C" CARPENTRY_GEOM_EXPORT bool
CGAL_Construct_InOutSide_Polygon(const std::vector<Vector2d> &py, const Vector2d &p, const Vector2d &q, bool &isPInside,
                                 bool &isQInside);
extern "C" CARPENTRY_GEOM_EXPORT bool
CGAL_2D_Intersection_Ray_Segment(const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s,
                                 const Vector2d &s_1_e, Vector2d &inter);
extern "C" CARPENTRY_GEOM_EXPORT double GetAngleKerfOffsetTan(const Vector2d &a, const Vector2d &b);
//implementation in "threeD.cpp"
//####################################################################################
extern "C" CARPENTRY_GEOM_EXPORT double CGAL_3D_Distance_Point_Segment(Vector3d p, Vector3d s_s, Vector3d s_e);
extern "C" CARPENTRY_GEOM_EXPORT void
CGAL_3D_Plane_Fitting(std::vector<Vector3d> &points, Vector3d &plane_p, Vector3d &plane_n);
extern "C" CARPENTRY_GEOM_EXPORT void CGAL_3D_Plane_Point_Projection(Vector3d &plane_p, Vector3d &plane_n, Vector3d &p,
                                                                     Vector3d &result);
extern "C" CARPENTRY_GEOM_EXPORT void CGAL_3D_Plane_Points_Projection(Vector3d &plane_p, Vector3d &plane_n,
                                                                      std::vector<Vector3d> &points,
                                                                      std::vector<Vector3d> &project_points);

extern "C" CARPENTRY_GEOM_EXPORT void
CGAL_3D_Plane_3D_to_2D_Point(Vector3d &plane_p, Vector3d &plane_n, Vector3d &point_3d,
                             Vector2d &result);
extern "C" CARPENTRY_GEOM_EXPORT void
CGAL_3D_Plane_2D_to_3D_Point(Vector3d &plane_p, Vector3d &plane_n, Vector2d &points_2d,
                             Vector3d &result);

extern "C" CARPENTRY_GEOM_EXPORT void CGAL_3D_Plane_3D_to_2D_Points(Vector3d &plane_p, Vector3d &plane_n,
                                                                    std::vector<Vector3d> &points_3d,
                                                                    std::vector<Vector2d> &points_2d);
extern "C" CARPENTRY_GEOM_EXPORT void CGAL_3D_Plane_2D_to_3D_Points(Vector3d &plane_p, Vector3d &plane_n,
                                                                    std::vector<Vector2d> &points_2d,
                                                                    std::vector<Vector3d> &points_3d);

extern "C" CARPENTRY_GEOM_EXPORT Vector3d CGAL_3D_Projection_Point_Segment(Vector3d p, Vector3d s_s, Vector3d s_e);
extern "C" CARPENTRY_GEOM_EXPORT double
CGAL_3D_Distance_Point_Point(double p_0_x, double p_0_y, double p_0_z, double p_1_x,
                             double p_1_y, double p_1_z);
extern "C" CARPENTRY_GEOM_EXPORT double
CGAL_3D_Distance_Point_Polygon(const std::vector<Vector3d> &py, const Vector3d &p);


extern "C" CARPENTRY_GEOM_EXPORT void CGAL_2D_Polygon_Triangulation(const std::vector<std::vector<Vector2d>> &polys, std::vector<std::vector<int>> &faces);

//implementation in "mesh.cpp"
//####################################################################################

extern "C" CARPENTRY_GEOM_EXPORT void CGAL_Remesh_Surface_by_Adding_Feature(const std::vector<Vector3d> &feature,
                                                                            const std::vector<int> &face_ids,
                                                                            const std::vector<Vector3d> &vecs,
                                                                            const std::vector<int> &face_id_0,
                                                                            const std::vector<int> &face_id_1,
                                                                            const std::vector<int> &face_id_2,
                                                                            std::vector<int> &igl_cutting_0_edges,
                                                                            std::vector<int> &igl_cutting_1_edges,
                                                                            std::vector<Vector3d> &igl_cutting_points,
                                                                            std::vector<std::vector<int> > &cutting_faces);

extern "C" CARPENTRY_GEOM_EXPORT void CGAL_3D_Read_Triangle_Mesh(std::string path, std::vector<Vector3d> &vecs,
                                                                 std::vector<int> &face_id_0,
                                                                 std::vector<int> &face_id_1,
                                                                 std::vector<int> &face_id_2);
extern "C" CARPENTRY_GEOM_EXPORT void CGAL_Mesh_Edges(std::string path);

#endif
