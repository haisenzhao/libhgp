#ifndef include_cgal
#define include_cgal
#pragma once

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

//typedef CGAL::Surface_mesh_shortest_path_traits<K, Polyhedron_3> Traits;
//typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
//typedef boost::graph_traits<Polyhedron_3> Graph_traits;
//typedef Graph_traits::vertex_iterator vertex_iterator;
//typedef Graph_traits::face_iterator face_iterator;

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
	FaceInfo2():nesting_level(0){}
	int nesting_level;
	bool in_domain() {
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
	Vertex(const Point& p) : superclass(p) {}
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

	polyhedron_builder(std::vector<double>& c, Vector1i1& t) {
		coords = c;
		tris = t;
	}

	void operator()(HDS& hds) override {
		typedef typename HDS::Vertex Vertex;
		typedef typename Vertex::Point Point;

		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
		B.begin_surface(coords.size() / 3, tris.size() / 3);

		for (int i = 0; i < (int)coords.size(); i += 3) {
			B.add_vertex(Point(coords[i + 0], coords[i + 1], coords[i + 2]));
		}

		for (int i = 0; i < (int)tris.size(); i += 3) {
			B.begin_facet();
			B.add_vertex_to_facet(tris[i + 0]);
			B.add_vertex_to_facet(tris[i + 1]);
			B.add_vertex_to_facet(tris[i + 2]);
			B.end_facet();
		}
	}
};

struct Edge
{
	int source;
	int end;
	Edge() {}
	Edge(int s, int e)
	{
		source = s;
		end = e;
	}
};

#endif
