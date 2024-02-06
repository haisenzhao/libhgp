#include "libhgp.h"




// Transformation
// Input: 3d point from cgal 
// Output: 3d vector from glm
Vector3d PointVector3d(Point_3 p) {
    return Vector3d(p[0], p[1], p[2]);
}

// Transformation
// Input: 3d vector from glm
// Output: 3d point from cgal 
Point_3 VectorPoint3d(Vector3d p) {
    return Point_3(p[0], p[1], p[2]);
}

// Transformation
// Input: 2d point from cgal 
// Output: 2d vector from glm
Vector2d PointVector2d(Point_2 p) {
    return Vector2d(p[0], p[1]);
}

// Transformation
// Input: 2d vector from glm
// Output: 2d point from cgal   
Point_2 VectorPoint2d(Vector2d p) {
    return Point_2(p[0], p[1]);
}



Point_3 point_to_3d(const Point_2& p, Plane_3& pl)
{
	Vector_3 basis[2];
	auto pop = pl.point();
	const Vector_3 vpop(pop.x(), pop.y(), pop.z());
	basis[0] = pl.base1() / CGAL::sqrt(pl.base1().squared_length());
	basis[1] = pl.base2() / CGAL::sqrt(pl.base2().squared_length());
	Vector_3 nr(pl.a(), pl.b(), pl.c());
	const Point_3 vi = pop + (p.x() * basis[0] + p.y() * basis[1]);
	return vi;
}

Point_2 point_to_2d(const Point_3& p, Plane_3& pl)
{
	Vector_3 basis[2];
	auto pop = pl.point();
	const Vector_3 vpop(pop.x(), pop.y(), pop.z());
	basis[0] = pl.base1() / CGAL::sqrt(pl.base1().squared_length());
	basis[1] = pl.base2() / CGAL::sqrt(pl.base2().squared_length());
	const Vector_3 ter(pop, p);
	return Point_2(ter * basis[0], ter * basis[1]);
}




/********************************************************
* Function name: insert_polygon
* Description: Insert a polygon into the Constrained Delaunay Triangulation (CDT) data structure.
* Parameters:
* @cdt: The Constrained Delaunay Triangulation.
* @polygon: The 2D polygon to be inserted.
* @indexInt: A vector of integers representing indices associated with the vertices of the polygon.
* Return: void
*********************************************************/
void insert_polygon(CDT& cdt, const Polygon_2& polygon, Vector1i1& indexInt)
{
	if (polygon.is_empty()) return;
	int index = 0;

	CDT::Vertex_handle v_prev = cdt.insert(*CGAL::cpp11::prev(polygon.vertices_end()));

	for (Polygon_2::Vertex_iterator vit = polygon.vertices_begin();
		vit != polygon.vertices_end(); ++vit)
	{
		CDT::Vertex_handle vh = cdt.insert(*vit);
		vh->index = indexInt[index];
		index++;
		cdt.insert_constraint(vh, v_prev);
		v_prev = vh;
	}
}
/********************************************************
* Function name: mark_domains
* Description: Mark the domains in the Constrained Delaunay Triangulation (CDT) data structure.
* Parameters:
* @ct: The Constrained Delaunay Triangulation.
* @start: The starting face for marking the domains.
* @index: The index used for marking the domains.
* @border: A list of CDT edges representing the borders of the domains.
* Return: void
*********************************************************/
void mark_domains(CDT& ct, CDT::Face_handle start, int index, std::list<CDT::Edge>& border)
{
	if (start->info().nesting_level != -1) {
		return;
	}
	std::list<CDT::Face_handle> queue;
	queue.push_back(start);
	while (!queue.empty()) {
		CDT::Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1) {
			fh->info().nesting_level = index;
			for (int i = 0; i < 3; i++) {
				CDT::Edge e(fh, i);
				CDT::Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1) {
					if (ct.is_constrained(e)) border.push_back(e);
					else queue.push_back(n);
				}
			}
		}
	}
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident 
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void mark_domains(CDT& cdt)
{
	for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) {
		it->info().nesting_level = -1;
	}
	std::list<CDT::Edge> border;
	mark_domains(cdt, cdt.infinite_face(), 0, border);
	while (!border.empty()) {
		CDT::Edge e = border.front();
		border.pop_front();
		CDT::Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1) {
			mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
		}
	}
}


/********************************************************
* Function name: Polygon2D
* Description: Create a 2D polygon from a vector of 2D points.
* Parameters:
* @py: A vector of 2D points representing the vertices of the polygon.
* Return: A 2D polygon created from the input points.
*********************************************************/
Polygon_2 Polygon2D(const Vector2d1& py)
{
	Polygon_2 polygon;
	for (auto p : py) polygon.push_back(VectorPoint2d(p));
	return polygon;
};





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

void ComputeEdgeLables(const int size_of_vertice, Halfedge_handle& start_hh, VectorPI1& edges, Vector1i1& lables);
void HGP_2D_Polygon_Simple_0(Vector2d1 points_2d);
void ComputeRemeshTriangles(const Vector3d1& vecs, const Vector1i1& face_id_0, const Vector1i1& face_id_1, const Vector1i1& face_id_2,
	const VectorPI1& edges, const Vector3d1& cutting_points, const Vector1i1& multi_cps, const Vector1i1& two_class_lables, const std::string output_path);

void  Construct_Polyhedron(Polyhedron_3& polyhedron, const Vector3d1& vecs, const Vector1i1& face_id_0, const Vector1i1& face_id_1, const Vector1i1& face_id_2)
{
	Vector1d1 coords;
	Vector1i1 tris;
	for (int i = 0; i < vecs.size(); i++)
	{
		coords.push_back(vecs[i][0]);
		coords.push_back(vecs[i][1]);
		coords.push_back(vecs[i][2]);
	}
	for (int i = 0; i < face_id_0.size(); i++)
	{
		tris.push_back(face_id_0[i]);
		tris.push_back(face_id_1[i]);
		tris.push_back(face_id_2[i]);
	}

	polyhedron_builder<Poly_halfedgeds_3> builder(coords, tris);
	polyhedron.delegate(builder);

	CGAL::set_halfedgeds_items_id(polyhedron);
	std::size_t facet_id = 0;
	for (Polyhedron_3::Facet_iterator facet_it = polyhedron.facets_begin();
		facet_it != polyhedron.facets_end(); ++facet_it, ++facet_id) {
		facet_it->id() = facet_id;
	}

	Vector1d1().swap(coords);
	Vector1i1().swap(tris);
}

void  Construct_Polyhedron(Polyhedron_3& polyhedron, const char* path_)
{
	std::string path = path_;
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		std::ifstream input(path);
		input >> polyhedron;
		input.close();

		CGAL::set_halfedgeds_items_id(polyhedron);
		std::size_t facet_id = 0;
		for (Polyhedron_3::Facet_iterator facet_it = polyhedron.facets_begin();
			facet_it != polyhedron.facets_end(); ++facet_it, ++facet_id) {
			facet_it->id() = facet_id;
		}
	}
	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		Vector3d1 vecs;
		Vector1i1 face_id_0;
		Vector1i1 face_id_1;
		Vector1i1 face_id_2;
		Functs::LoadObj3d(path_, vecs, face_id_0, face_id_1, face_id_2);
		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);
	}
}

void  Construct_Polyhedron(Polyhedron_3& polyhedron, const char* path_, Vector3d1& vecs, Vector1i1& face_id_0, Vector1i1& face_id_1, Vector1i1& face_id_2)
{
	std::string path = path_;
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		std::ifstream input(path);
		input >> polyhedron;
		input.close();

		CGAL::set_halfedgeds_items_id(polyhedron);
		std::size_t facet_id = 0;
		for (Polyhedron_3::Facet_iterator facet_it = polyhedron.facets_begin();
			facet_it != polyhedron.facets_end(); ++facet_it, ++facet_id) {
			facet_it->id() = facet_id;
		}

		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin(); iter != polyhedron.vertices_end(); iter++)
		{
			Point_3 p = iter->point();
			vecs.push_back(Vector3d(p[0], p[1], p[2]));
		}

		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			face_id_0.push_back((int)iter->halfedge()->next()->next()->vertex()->id());
			face_id_1.push_back((int)iter->halfedge()->vertex()->id());
			face_id_2.push_back((int)iter->halfedge()->next()->vertex()->id());
		}
	}
	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		Functs::LoadObj3d(path_, vecs, face_id_0, face_id_1, face_id_2);
		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);
	}
}



/********************************************************
* Function name : OutsidePointInsideTriangle
* Description  :	Determine if a point is inside a triangle by the center of gravity coordinate method.
* Parameter  :
* @face				 Triangular surface
* @p				 3d point
* Return  :true -- inside  ,  false -- outside
**********************************************************/
bool OutsidePointInsideTriangle(Poly_facet_iterator& face, Vector3d p)
{
	Point_3 p0 = face->halfedge()->next()->next()->vertex()->point();
	Point_3 p1 = face->halfedge()->vertex()->point();
	Point_3 p2 = face->halfedge()->next()->vertex()->point();
	Plane_3 plane(p1, CGAL::cross_product(p2 - p1, p0 - p1));
	Point_3 project = plane.projection(VectorPoint3d(p));

	Vector3d v0 = PointVector3d(p0);
	Vector3d v1 = PointVector3d(p1);
	Vector3d v2 = PointVector3d(p2);
	Vector3d vp = PointVector3d(project);

	double u, v, w;
	Functs::Barycentric(vp, v0, v1, v2, u, v, w);

	if ((u >= 0.0 && u <= 1.0) && (v >= 0.0 && v <= 1.0) && (w >= 0.0 && w <= 1.0)) {
		return true;
	}
	else {
		return false;
	}
}

/********************************************************
* Function name :Intersection
* Description        :Determine whether a given two-dimensional line segment intersects any half of the line.
* Parameter         :
* @hh       half edge
* @nb       integer
* @inside		3d point
* @outside		3d point
* @handle		the handle to store intersecting halves
* @intersection		3d point to store intersection
* Return          :false -- no intersection  ,  true -- intersection exists
**********************************************************/
bool Intersection(Halfedge_handle& hh, int nb, Vector3d inside, Vector3d outside, Halfedge_handle& handle,
	Vector3d& intersection) {
	Point_3 p0 = hh->next()->next()->vertex()->point();
	Point_3 p1 = hh->vertex()->point();
	Point_3 p2 = hh->next()->vertex()->point();
	Plane_3 plane(p1, CGAL::cross_product(p2 - p1, p0 - p1));

	Vector2d inside_2d = PointVector2d(plane.to_2d(VectorPoint3d(inside)));
	Vector2d outside_2d = PointVector2d(plane.to_2d(VectorPoint3d(outside)));

	for (int i = 0; i < nb; i++) {
		hh = hh->next();
		Vector2d edge_0 = PointVector2d(plane.to_2d(hh->vertex()->point()));
		Vector2d edge_1 = PointVector2d(plane.to_2d(hh->opposite()->vertex()->point()));
		Vector3d edge_3d_0 = PointVector3d(hh->vertex()->point());
		Vector3d edge_3d_1 = PointVector3d(hh->opposite()->vertex()->point());

		Vector2d iter;
		if (HGP_2D_Intersection_Segment_Segment(inside_2d, outside_2d, edge_0, edge_1, iter)) {
			intersection = PointVector3d(plane.to_3d(VectorPoint2d(iter)));
			intersection = HGP_3D_Projection_Point_Segment(intersection, edge_3d_0, edge_3d_1);
			handle = hh;
			return true;
		}
	}
	return false;
}


/********************************************************
* Function name :RelatedFaceNormal
* Description        : Calculate the normal of the triangle (or polygon) face closest to the given point.
* Parameter         :
* @polyhedron       Data structures for polyhedra or polygons
* @tree       Tree structure of spatial queries
* @normals		normal vector
* @source		the given 3d point
* Return          :a new 3d vector of the normal
**********************************************************/
Vector3d RelatedFaceNormal(Polyhedron_3& polyhedron, Tree& tree, Vector3d1& normals, Vector3d source)
{
	Point_3 query(source[0], source[1], source[2]);
	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

	Point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
	Point_3 p1 = pp.second->halfedge()->vertex()->point();
	Point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

	int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
	int point_id_1 = pp.second->halfedge()->vertex()->id();
	int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

	double u, v, w;
	Functs::Barycentric(PointVector3d(query), PointVector3d(p0), PointVector3d(p1), PointVector3d(p2), u, v, w);
	return (double)u * normals[point_id_0] + (double)v * normals[point_id_1] + (double)w * normals[point_id_2];
}






