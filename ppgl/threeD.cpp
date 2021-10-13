#include "geom.h"

Point_3 point_to_3d(const Point_2 &p, Plane_3 &pl) {
    Vector_3 basis[2];
    auto pop = pl.point();
    const Vector_3 vpop(pop.x(), pop.y(), pop.z());
    basis[0] = pl.base1() / CGAL::sqrt(pl.base1().squared_length());
    basis[1] = pl.base2() / CGAL::sqrt(pl.base2().squared_length());
    Vector_3 nr(pl.a(), pl.b(), pl.c());
    const Point_3 vi = pop + (p.x() * basis[0] + p.y() * basis[1]);
    return vi;
}

Point_2 point_to_2d(const Point_3 &p, Plane_3 &pl) {
    Vector_3 basis[2];
    auto pop = pl.point();
    const Vector_3 vpop(pop.x(), pop.y(), pop.z());
    basis[0] = pl.base1() / CGAL::sqrt(pl.base1().squared_length());
    basis[1] = pl.base2() / CGAL::sqrt(pl.base2().squared_length());
    const Vector_3 ter(pop, p);
    return Point_2(ter * basis[0], ter * basis[1]);
}

extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Segment(const Vector3d& p, const Vector3d & s_s, const Vector3d & s_e) {
    return sqrt((double) CGAL::squared_distance(VectorPoint3d(p), Segment_3(VectorPoint3d(s_s), VectorPoint3d(s_e))));
}

extern "C" PPGL_EXPORT void CGAL_3D_Plane_Fitting(const Vector3d1 &points, Vector3d &plane_p, Vector3d &plane_n) {
    // centroid of 3D points
    std::vector<Point_3> points_3;
    for (int i = 0; i < points.size(); i++)
        points_3.push_back(VectorPoint3d(points[i]));
    Point_3 center = CGAL::centroid(points_3.begin(), points_3.end(), CGAL::Dimension_tag<0>());

    std::vector<Triangle> triangles;
    for (int i = 0; i < points_3.size() - 1; i++)
        triangles.emplace_back(points_3[i], points_3[i + 1], center);

    Plane_3 plane;
    // fit plane to whole triangles
    linear_least_squares_fitting_3(triangles.begin(), triangles.end(), plane, CGAL::Dimension_tag<2>());

    //set plane point
    plane_p = PointVector3d(plane.projection(center));

    //set plane vector
    plane_n[0] = plane.orthogonal_vector().x();
    plane_n[1] = plane.orthogonal_vector().y();
    plane_n[2] = plane.orthogonal_vector().z();
}

extern "C" PPGL_EXPORT void CGAL_3D_Plane_Point_Projection(const Vector3d &plane_p, const Vector3d &plane_n, const Vector3d &p, Vector3d &result) 
{
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    result = PointVector3d(plane.projection(VectorPoint3d(p)));
}

extern "C" PPGL_EXPORT void CGAL_3D_Plane_Points_Projection
(const Vector3d &plane_p, const Vector3d &plane_n, const Vector3d1 &points,Vector3d1 &project_points) 
{
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    for (int i = 0; i < points.size(); i++)
        project_points.push_back(PointVector3d(plane.projection(VectorPoint3d(points[i]))));
}

extern "C" PPGL_EXPORT void
CGAL_3D_Plane_3D_to_2D_Point(const Vector3d &plane_p, const Vector3d &plane_n, const Vector3d &point_3d, Vector2d &result) 
{
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    result = PointVector2d(point_to_2d(VectorPoint3d(point_3d), plane));
}

extern "C" PPGL_EXPORT void
CGAL_3D_Plane_2D_to_3D_Point(const Vector3d &plane_p, const Vector3d &plane_n, const Vector2d &points_2d, Vector3d &result) 
{
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    result = PointVector3d(point_to_3d(VectorPoint2d(points_2d), plane));
}

extern "C" PPGL_EXPORT void CGAL_3D_Plane_3D_to_2D_Points
(const Vector3d &plane_p, const Vector3d &plane_n, const Vector3d1 &points_3d, Vector2d1 &points_2d) 
{
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    for (int i = 0; i < points_3d.size(); i++)
        points_2d.push_back(PointVector2d(plane.to_2d(VectorPoint3d(points_3d[i]))));
}

extern "C" PPGL_EXPORT void CGAL_3D_Plane_2D_to_3D_Points
(const Vector3d &plane_p, const Vector3d &plane_n, const Vector2d1 &points_2d, Vector3d1 &points_3d) 
{
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    for (int i = 0; i < points_2d.size(); i++)
        points_3d.push_back(PointVector3d(plane.to_3d(VectorPoint2d(points_2d[i]))));
}

extern "C" PPGL_EXPORT
double CGAL_3D_Distance_Point_Point(const Vector3d& v0, const Vector3d& v1)
{
	return sqrt((double)CGAL::squared_distance(Point_3(v0[0], v0[1], v0[2]), Point_3(v1[0], v1[1], v1[2])));
}


extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Segment(const Vector3d& p, const Vector3d & s_s, const Vector3d & s_e) {
    Line_3 l(VectorPoint3d(s_s), VectorPoint3d(s_e));
    auto m_p = PointVector3d(l.projection(VectorPoint3d(p)));
	double d_m_s = CGAL_3D_Distance_Point_Point(m_p, s_s);
	double d_m_e = CGAL_3D_Distance_Point_Point(m_p, s_e);
	double d_s_e = CGAL_3D_Distance_Point_Point(s_s, s_e);

    if (d_m_s >= d_s_e)
        return s_e;
    if (d_m_e >= d_s_e)
        return s_s;
    return m_p;
}

extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Segment_Ref(const Vector3d &v, const Vector3d &s_0,
                                                                           const Vector3d &s_1) {
    return sqrt((double) CGAL::squared_distance(Point_3(v[0], v[1], v[2]),
                                                Segment_3(Point_3(s_0[0], s_0[1], s_0[2]),
                                                          Point_3(s_1[0], s_1[1], s_1[2]))));
}

extern "C" PPGL_EXPORT double
CGAL_3D_Distance_Point_Polygon(const Vector3d1 &py, const Vector3d &p) {
    double distance = 1000000000000.0;
    for (int i = 0; i < py.size(); i++)
        distance = std::min(distance, CGAL_3D_Distance_Point_Segment_Ref(p, py[i], py[(i + 1) % py.size()]));
    return distance;
}


void insert_polygon(CDT& cdt, const Polygon_2& polygon,Vector1i1 &indexInt){
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

void
mark_domains(CDT& ct,
CDT::Face_handle start,
int index,
std::list<CDT::Edge>& border)
{
	if (start->info().nesting_level != -1){
		return;
	}
	std::list<CDT::Face_handle> queue;
	queue.push_back(start);
	while (!queue.empty()){
		CDT::Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1){
			fh->info().nesting_level = index;
			for (int i = 0; i < 3; i++){
				CDT::Edge e(fh, i);
				CDT::Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1){
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
void
mark_domains(CDT& cdt)
{
	for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
		it->info().nesting_level = -1;
	}
	std::list<CDT::Edge> border;
	mark_domains(cdt, cdt.infinite_face(), 0, border);
	while (!border.empty()){
		CDT::Edge e = border.front();
		border.pop_front();
		CDT::Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1){
			mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
		}
	}
}

extern "C" PPGL_EXPORT
void CGAL_2D_Polygon_Triangulation(const Vector2d2 &polys, Vector1i2 &faces)
{
	int nb = 0;
	CDT cdt;
	for (int i = 0; i < polys.size(); i++)
	{
		Polygon_2 polygon;
		std::vector<int> indexInt;
		for (int j = 0; j < polys[i].size(); j++)
		{
			polygon.push_back(Point_2(polys[i][j][0], polys[i][j][1]));
			indexInt.emplace_back(j + nb);
		}
		nb += (int)polys[i].size();
		insert_polygon(cdt, polygon, indexInt);
	}

	//Mark facets that are inside the domain bounded by the polygon
	mark_domains(cdt);

	for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
		fit != cdt.finite_faces_end(); ++fit)
	if (fit->info().in_domain())
		faces.emplace_back(std::vector<int>{fit->vertex(2)->index, fit->vertex(1)->index, fit->vertex(0)->index});
}


//IO mesh
/***************************************************************************************************/


//void CGAL_Load_Obj(std::string path, std::vector<double> &coords,Vector1i1 &tris) {
//    double x, y, z;
//    char line[1024], v0[1024], v1[1024], v2[1024];
//
//    // open the file, return if open fails
//    if (!Functs::DetectExisting(path))
//        Functs::MAssert("Path does not exist: "+ path);
//    FILE *fp = fopen(path.c_str(), "r");
//    if (!fp) return;
//
//    while (fgets(line, 1024, fp)) {
//        if (line[0] == 'v') {
//            sscanf(line, "%*s%lf%lf%lf", &x, &y, &z);
//            coords.push_back(x);
//            coords.push_back(y);
//            coords.push_back(z);
//        } else if (line[0] == 'f') {
//            sscanf(line, "%*s%s%s%s", v0, v1, v2);
//            tris.push_back(get_first_integer(v0) - 1);
//            tris.push_back(get_first_integer(v1) - 1);
//            tris.push_back(get_first_integer(v2) - 1);
//        }
//    }
//    fclose(fp);
//}

extern "C" PPGL_EXPORT void CGAL_3D_Output_Triangle_Mesh
(const std::string& path, const Vector3d1 & vecs, const Vector1i1 & face_id_0, const Vector1i1 & face_id_1, const Vector1i1 & face_id_2)
{
	std::ofstream file(path);

	for (auto& vec : vecs)
		file << "v " << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;

    for (int i = 0; i < face_id_0.size(); i++)
        file << "f " << face_id_0[i]+1 << " " << face_id_1[i]+1 << " " << face_id_2[i]+1 << std::endl;

	file.clear();
	file.close();
}

extern "C" PPGL_EXPORT void CGAL_3D_Read_Triangle_Mesh(const std::string& path, Vector3d1 &vecs,
                        Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2)
{
    //if (path.substr(path.size() - 3, path.size()) == "off")
    //{
    //	Polyhedron_3 polyhedron;
    //	Construct_Polyhedron(polyhedron, path);

    //	for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
    //		iter != polyhedron.vertices_end(); iter++)
    //	{
    //		Poly_point_3 p = iter->point();
    //		vecs.push_back(Vector3d(p[0], p[1], p[2]));
    //	}

    //	for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
    //	{
    //		//Poly_point_3 p0 = iter->halfedge()->next()->next()->vertex()->point();
    //		//Poly_point_3 p1 = iter->halfedge()->vertex()->point();
    //		//Poly_point_3 p2 = iter->halfedge()->next()->vertex()->point();
    //		face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
    //		face_id_1.push_back(iter->halfedge()->vertex()->id());
    //		face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
    //	}
    //}

    if (path.substr(path.size() - 3, path.size()) == "obj") {
        std::vector<double> coords;
       Vector1i1 tris;
        CGAL_Load_Obj(path.c_str(), coords, tris);
        if (coords.size() == 0)
            return;

        //std::cout << "Size: " << coords.size() / 3 << " " << tris.size() / 3 << std::endl;

        for (int i = 0; i < (int) coords.size(); i += 3) {
            vecs.push_back(Vector3d(coords[i + 0], coords[i + 1], coords[i + 2]));
        }

        for (int i = 0; i < (int) tris.size(); i += 3) {
            face_id_0.push_back(tris[i + 0]);
            face_id_1.push_back(tris[i + 1]);
            face_id_2.push_back(tris[i + 2]);
        }
        /*********************************************************************************/
    }
}

extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Line(const Vector3d & p, const Vector3d & l_s, const Vector3d & l_e)
{
	return sqrt((double)CGAL::squared_distance(VectorPoint3d(p), Line_3(VectorPoint3d(l_s), VectorPoint3d(l_e))));
}

extern "C" PPGL_EXPORT Vector3d CGAL_3D_Projection_Point_Line(const Vector3d & p, const Vector3d & l_s, const Vector3d & l_e)
{
	Line_3 l(VectorPoint3d(l_s), VectorPoint3d(l_e));
	Point_3 o_p = l.projection(VectorPoint3d(p));
	return PointVector3d(o_p);
}
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Segment_Segment(const Vector3d & s_0_s, const Vector3d & s_0_e, const Vector3d & s_1_s, const Vector3d & s_1_e)
{
	return sqrt((double)CGAL::squared_distance(Segment_3(VectorPoint3d(s_0_s), VectorPoint3d(s_0_e)), Segment_3(VectorPoint3d(s_1_s), VectorPoint3d(s_1_e))));
}

extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Plane(const Vector3d & v, const Vector3d & plane_p, const Vector3d & plane_n)
{
	return sqrt((double)CGAL::squared_distance(VectorPoint3d(v), Plane_3(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]))));
}



extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Line(const Vector3d & s_s, const Vector3d & s_e, const Vector3d & l_s, const Vector3d & l_e, Vector3d & inter)
{
	double d = CGAL_3D_Distance_Point_Point(s_s, s_e);

	if (!Functs::IsAlmostZero(d))
	{
		CGAL::Object result = CGAL::intersection(Segment_3(VectorPoint3d(s_s), VectorPoint3d(s_e)),
			Line_3(VectorPoint3d(l_s), VectorPoint3d(l_e)));

		if (const Point_3* ipoint = CGAL::object_cast<Point_3>(&result))
		{
			inter[0] = ipoint->x();
			inter[1] = ipoint->y();
			inter[2] = ipoint->z();
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Segment(const Vector3d & s_0_s, const Vector3d & s_0_e, const Vector3d & s_1_s, const Vector3d & s_1_e, Vector3d & iter)
{
	double d0 = CGAL_3D_Distance_Point_Point(s_0_s, s_0_e);
	double d1 = CGAL_3D_Distance_Point_Point(s_1_s, s_1_e);

	if (!Functs::IsAlmostZero(d0) && !Functs::IsAlmostZero(d1))
	{
		CGAL::Object result = CGAL::intersection(Segment_3(VectorPoint3d(s_0_s), VectorPoint3d(s_0_e)),
			Segment_3(VectorPoint3d(s_1_s), VectorPoint3d(s_1_e)));

		if (const Point_3* ipoint = CGAL::object_cast<Point_3>(&result))
		{
			iter = PointVector3d(*ipoint);
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Plane(const Vector3d & s_s, const Vector3d & s_e, const Vector3d & plane_p, const Vector3d & plane_n, Vector3d & inter)
{
	Segment_3 s(VectorPoint3d(s_s), VectorPoint3d(s_e));
	Plane_3 p(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	CGAL::Object result = CGAL::intersection(s, p);
	if (const Point_3* ipoint = CGAL::object_cast<Point_3>(&result))
	{
		inter[0] = ipoint->x();
		inter[1] = ipoint->y();
		inter[2] = ipoint->z();
		return true;
	}
	else
	{
		return false;
	}
}

extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Line_Plane(const Vector3d & l_s, const Vector3d & l_e, const Vector3d & plane_p, const Vector3d & plane_n, Vector3d & inter)
{
	Line_3 s(VectorPoint3d(l_s), VectorPoint3d(l_e));
	Plane_3 p(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));

	CGAL::Object result = CGAL::intersection(s, p);

	if (const Point_3* ipoint = CGAL::object_cast<Point_3>(&result))
	{
		inter[0] = ipoint->x();
		inter[1] = ipoint->y();
		inter[2] = ipoint->z();
		return true;
	}
	else
	{
		return false;
	}
}