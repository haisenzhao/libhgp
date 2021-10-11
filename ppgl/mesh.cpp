#include "geom.h"

//extern "C" PPGL_EXPORT


void  Construct_Polyhedron(Polyhedron_3& polyhedron, const Vector3d1& vecs, const Vector1i1& face_id_0, const Vector1i1& face_id_1, const Vector1i1& face_id_2)
{
    Vector1d1 coords;
    Vector1i1    tris;
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

	std::vector<double>().swap(coords);
	std::vector<int>().swap(tris);
}

void  Construct_Polyhedron(Polyhedron_3& polyhedron, const std::string& path)
{
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
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);
	}
}

void  Construct_Polyhedron(Polyhedron_3& polyhedron, const std::string& path, Vector3d1& vecs, Vector1i1& face_id_0, Vector1i1& face_id_1, Vector1i1& face_id_2)
{
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
			Poly_point_3 p = iter->point();
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
		CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);
	}
}


Poly_point_3 Minus(Poly_point_3 a, Poly_point_3 b) {
    return Poly_point_3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

double Dot(Poly_point_3 a, Poly_point_3 b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void Barycentric(Poly_point_3 p, Poly_point_3 a, Poly_point_3 b, Poly_point_3 c, double &u, double &v, double &w) {
    Poly_point_3 v0 = Minus(b, a), v1 = Minus(c, a), v2 = Minus(p, a);

    double d00 = Dot(v0, v0);
    double d01 = Dot(v0, v1);
    double d11 = Dot(v1, v1);
    double d20 = Dot(v2, v0);
    double d21 = Dot(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
}

void CGAL_Barycentric(Vector3d p, Vector3d a, Vector3d b, Vector3d c, double &u, double &v, double &w) {
    Barycentric(Poly_point_3(p[0], p[1], p[2]), Poly_point_3(a[0], a[1], a[2]), Poly_point_3(b[0], b[1], b[2]),
                Poly_point_3(c[0], c[1], c[2]), u, v, w);
}

//Project p onto the planar surface of 3d triangle
//Checking the position relationship between the p and 3d triangle
//face: 3d triangle
//p: 3d point
//return true: inside
//return false: outside
bool OutsidePointInsideTriangle(Poly_facet_iterator &face, Vector3d p) {
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
    CGAL_Barycentric(vp, v0, v1, v2, u, v, w);

    if ((u >= 0.0 && u <= 1.0) && (v >= 0.0 && v <= 1.0) && (w >= 0.0 && w <= 1.0)) {
        return true;
    } else {
        return false;
    }
}


bool Intersection(Halfedge_handle &hh, int nb, Vector3d inside, Vector3d outside, Halfedge_handle &handle,
                  Vector3d &intersection) {
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
        if (CGAL_2D_Intersection_Segment_Segment(inside_2d, outside_2d, edge_0, edge_1, iter)) {
            intersection = PointVector3d(plane.to_3d(VectorPoint2d(iter)));
            intersection = CGAL_3D_Projection_Point_Segment(intersection, edge_3d_0, edge_3d_1);
            handle = hh;
            return true;
        }
    }
    return false;
}

extern "C" PPGL_EXPORT  void
CGAL_Remesh_Surface_by_Adding_Feature(const Vector3d1 &feature, const Vector1i1 &face_ids,
                                      const Vector3d1 &vecs,
                                      const Vector1i1 &face_id_0, const Vector1i1 &face_id_1,
                                      const Vector1i1 &face_id_2,
                                     Vector1i1 &igl_cutting_0_edges,Vector1i1 &igl_cutting_1_edges,
                                      Vector3d1 &igl_cutting_points,
                                      Vector1i2 &cutting_faces) {
    //polyhedron
    Polyhedron_3 polyhedron;
    Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

    //iterations
    std::vector<Poly_facet_iterator> all_faces_iters;//related faces of projecting points
    std::vector<Poly_facet_iterator> project_faces;//related faces of projecting points
    for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++) {
        all_faces_iters.push_back(iter);

    }
    for (int i = 0; i < face_ids.size(); i++) {
        project_faces.push_back(all_faces_iters[face_ids[i]]);

    }
    std::vector<Poly_facet_iterator>().swap(all_faces_iters);

    /*******************************************/
    //searching for all of the cutting points on edges
    Vector3d1 cutting_points;
    /*******************************************/
    int cur_face_id = (int)project_faces[0]->id();
    Poly_facet_iterator cur_face = project_faces[0];
    Halfedge_handle cur_handle = cur_face->halfedge();
    Vector3d inside = feature[0];

    int next_index = 1;
    std::vector<Halfedge_handle> handles;
    int iteration = 0;
    while (true) {
        std::cout << iteration << " / " << feature.size() << std::endl;
        //searching for the outside point of the current triangle
        bool goon = false;
        while (true) {
            if (cur_face_id == project_faces[next_index]->id()) {
                inside = feature[next_index];
				next_index = (next_index + 1);
                next_index = next_index % feature.size();
                if (next_index == 0) break;
            } else {
                if (OutsidePointInsideTriangle(cur_face, feature[next_index])) {
					next_index = next_index + 1;
                    next_index = next_index% feature.size();
                    if (next_index == 0) break;
                } else {
                    goon = true;
                    break;
                }
            }
        }

        if (!goon) break;

        Halfedge_handle handle;
        Vector3d intersection;
        bool b;
        if (iteration == 0)
            b = Intersection(cur_handle, 3, inside, feature[next_index], handle, intersection);
        else
            b = Intersection(cur_handle, 2, inside, feature[next_index], handle, intersection);

        if (b) {
            cutting_points.push_back(intersection);
            handles.push_back(cur_handle);

            //move next step
            inside = intersection;
            cur_handle = handle->opposite();
            cur_face_id = (int)cur_handle->face()->id();
            cur_face = cur_handle->face();

            if (cur_face_id == project_faces[0]->id()) {
                break;
            }
        } else {
            cutting_points.erase(cutting_points.begin() + cutting_points.size() - 1);
            handles.erase(handles.begin() + handles.size() - 1);
			next_index = next_index + 1;
            next_index = next_index % feature.size();
            if (next_index == 0) break;

            inside = cutting_points[cutting_points.size() - 1];
            cur_handle = handles[handles.size() - 1]->opposite();
            cur_face_id = (int)cur_handle->face()->id();
            cur_face = cur_handle->face();
        }
        iteration++;
    }

    for (int i = 0; i < handles.size(); i++) {
        cutting_faces[(int)handles[i]->face()->id()].push_back((int)igl_cutting_0_edges.size());
        cutting_faces[(int)handles[i]->opposite()->face()->id()].push_back((int)igl_cutting_0_edges.size());
        igl_cutting_0_edges.push_back((int)handles[i]->vertex()->id());
        igl_cutting_1_edges.push_back((int)handles[i]->opposite()->vertex()->id());
    }
    igl_cutting_points = cutting_points;
}


extern "C" PPGL_EXPORT void CGAL_Mesh_Edges(const std::string& path) {
}



extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Sphere_Ray(const double& center_x, const double& center_y, const double& center_z, const double& radius,
	const double& ray_origin_x, const double& ray_origin_y, const double& ray_origin_z, const double& ray_direction_x, const double& ray_direction_y, const double& ray_direction_z,
	std::vector<double>& i_x, std::vector<double>& i_y, std::vector<double>& i_z)
{
	//Wm5::Sphere3<double> sphere(Wm5::Vector3d(center_x, center_y, center_z), radius);
	//Wm5::Ray3<double> ray(Wm5::Vector3d(ray_origin_x, ray_origin_y, ray_origin_z), Wm5::Vector3d(ray_direction_x, ray_direction_y, ray_direction_z));

	//Wm5::IntrRay3Sphere3d intr(ray,sphere);

	//intr.Test();
	//intr.Find();

	//int nb = intr.GetQuantity();

	//for (int i = 0; i < nb; i++)
	//{
	//	Wm5::Vector3d p = intr.GetPoint(i);
	//	i_x.push_back(p[0]);
	//	i_y.push_back(p[1]);
	//	i_z.push_back(p[2]);
	//}

	//return nb>0;

	return false;
}

extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Triangle(const Vector3d& p, const Vector3d & n, const Vector3d & p0, const Vector3d & p1, const Vector3d & p2)
{
	Ray_3 ray(VectorPoint3d(p), Vector_3(n[0], n[1], n[2]));

	K::Triangle_3 tri(VectorPoint3d(p0), VectorPoint3d(p1), VectorPoint3d(p2));
	CGAL::Object result = CGAL::intersection(ray, tri);
	if (const Point_2* ipoint = CGAL::object_cast<Point_2>(&result))
	{
		return true;
	}
	else
	{
		return false;
	}
}

extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Mesh(const Vector3d& p, const Vector3d & n, const std::string& path)
{
	std::cerr << "CGAL_3D_Intersection_Ray_Mesh..." << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Ray_3 ray(Point_3(p[0], p[1], p[2]), Vector_3(n[0], n[1], n[2]));

	if (tree.do_intersect(ray))
		return true;
	else
		return false;
}

//test each group directions (nes[i]) for each point in ps
extern "C" PPGL_EXPORT void  CGAL_3D_Intersection_Rays_Mesh_C1_Bool(const Vector3d1& ps, const Vector3d2& nes, const std::string& path, Vector1b2& inters)
{
	//input validation
	if (ps.size() != nes.size()  || !Functs::DetectExisting(path))
	{
		Functs::MAssert("ps.size() != nes.size() || inters.size() != 0 || !Functs::DetectExisting(path)");
	}

	//construct polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	//intersection
	inters = Vector1b2(ps.size(), Vector1b1());
	for (int i = 0; i < ps.size(); i++)
	{
		Functs::OutputIterInfo("CGAL_3D_Intersection_Rays_Mesh_Bool", (int)ps.size(), (int)i,(int)100);
		Point_3 p3=VectorPoint3d(ps[i]);
		for (int j = 0; j < nes[i].size(); j++)
		{
			Ray_3 ray(p3, Vector_3(nes[i][j][0], nes[i][j][1], nes[i][j][2]));
			inters[i].push_back(tree.do_intersect(ray));
		}
	}
}

//test all directions (ns) for each point in ps
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_C2_Bool(const Vector3d1& ps, const Vector3d1& ns, const std::string & path, Vector1b2& inters)
{
	//input validation
	if (!Functs::DetectExisting(path))
		Functs::MAssert("!Functs::DetectExisting(path)");

	//construct polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	//intersection
	inters = Vector1b2(ps.size(), Vector1b1());
	for (int i = 0; i < ps.size(); i++)
	{
		Functs::OutputIterInfo("CGAL_3D_Intersection_Rays_Mesh_Bool", (int)ps.size(), (int)i, (int)10);
		Point_3 p3 = VectorPoint3d(ps[i]);

		for (int j = 0; j < ns.size(); j++)
		{
			Ray_3 ray(p3, Vector_3(ns[j][0], ns[j][1], ns[j][2]));
			inters[i].push_back(tree.do_intersect(ray));
		}
	}
}

extern "C" PPGL_EXPORT void  CGAL_3D_Points_Inside_Triangles_C1_Bool(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2, const Vector3d1& points, std::vector<bool>& insides)
{
	//build polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> inside(polyhedron);

	for (int i = 0; i < points.size(); i++)
	{
		CGAL::Bounded_side res = inside(Point_3(points[i][0], points[i][1], points[i][2]));
		if (res == CGAL::ON_BOUNDED_SIDE)
			insides.push_back(true);
		else
			insides.push_back(false);
	}
}

extern "C" PPGL_EXPORT void CGAL_3D_Points_Inside_Triangles_C2_Bool(const std::string& path, const Vector3d1& points, std::vector<bool>& insides)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> inside(polyhedron);

	for (int i = 0; i < points.size(); i++)
	{
		CGAL::Bounded_side res = inside(Point_3(points[i][0], points[i][1], points[i][2]));
		if (res == CGAL::ON_BOUNDED_SIDE)
			insides.push_back(true);
		else
			insides.push_back(false);
	}
}

//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Dart_Sampling_C1(const std::string & outside_path, const double& d, Vector3d1 & sampling_points, const int& total_iter)
{
	if (!(d > 0 && d < 1.0))
		Functs::MAssert("CGAL_3D_Mesh_Dart_Sampling_C1 if (!(d > 0 && d < 1.0))");

	Polyhedron_3 out_polyhedron;
	Vector3d1 out_vecs;
	Vector1i1 out_face_id_0, out_face_id_1, out_face_id_2;
	Construct_Polyhedron(out_polyhedron, outside_path, out_vecs, out_face_id_0, out_face_id_1, out_face_id_2);
	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> out_checker(out_polyhedron);
	Vector3d out_minC, out_maxC;
	Functs::GetBoundingBox(out_vecs, out_minC, out_maxC);

	double diagonal_length = CGAL_3D_Distance_Point_Point(out_minC, out_maxC);
	double minimal_d = d * diagonal_length;

	int run = 0;
	while (run < total_iter)
	{
		run++;
		double x = rand() / double(RAND_MAX);
		double y = rand() / double(RAND_MAX);
		double z = rand() / double(RAND_MAX);
		x = (out_maxC[0] - out_minC[0]) * x + out_minC[0];
		y = (out_maxC[1] - out_minC[1]) * y + out_minC[1];
		z = (out_maxC[2] - out_minC[2]) * z + out_minC[2];

		CGAL::Bounded_side res = out_checker(Point_3(x, y, z));
		if (res == CGAL::ON_BOUNDED_SIDE)
		{
			double distance = CGAL_IA_MAX_DOUBLE;
			for (int i = 0; i < sampling_points.size(); i++)
				distance = std::min(distance, CGAL_3D_Distance_Point_Point(sampling_points[i], Vector3d(x, y, z)));

			if (distance > minimal_d)
			{
				sampling_points.push_back(Vector3d(x, y, z));
				if (sampling_points.size() % 100 == 0) std::cerr << sampling_points.size() << " ";
				run = 0;
			}
		}
	}
	std::cerr << std::endl;
}

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Dart_Sampling_C2(const std::string & outside_path, const std::string & inside_path, const double& d, Vector3d1 & sampling_points, const int& total_iter)
{
	if (!(d > 0 && d < 1.0))
		Functs::MAssert("CGAL_3D_Mesh_Dart_Sampling_C1 if (!(d > 0 && d < 1.0))");

	//outside
	Polyhedron_3 out_polyhedron;
	Vector3d1 out_vecs;
	Vector1i1 out_face_id_0, out_face_id_1, out_face_id_2;
	Construct_Polyhedron(out_polyhedron, outside_path, out_vecs, out_face_id_0, out_face_id_1, out_face_id_2);
	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> out_checker(out_polyhedron);
	Vector3d out_minC, out_maxC;
	Functs::GetBoundingBox(out_vecs, out_minC, out_maxC);

	//inside
	Polyhedron_3 in_polyhedron;
	Construct_Polyhedron(in_polyhedron, inside_path);
	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> in_checker(in_polyhedron);

	double diagonal_length = CGAL_3D_Distance_Point_Point(out_minC, out_maxC);
	double minimal_d = d * diagonal_length;

	int run = 0;
	while (run < total_iter)
	{
		run++;
		double x = rand() / double(RAND_MAX);
		double y = rand() / double(RAND_MAX);
		double z = rand() / double(RAND_MAX);
		x = (out_maxC[0] - out_minC[0]) * x + out_minC[0];
		y = (out_maxC[1] - out_minC[1]) * y + out_minC[1];
		z = (out_maxC[2] - out_minC[2]) * z + out_minC[2];

		CGAL::Bounded_side out_res = out_checker(Point_3(x, y, z));
		CGAL::Bounded_side in_res = in_checker(Point_3(x, y, z));
		if (out_res == CGAL::ON_BOUNDED_SIDE && in_res == CGAL::ON_UNBOUNDED_SIDE)
		{
			double distance = CGAL_IA_MAX_DOUBLE;
			for (int i = 0; i < sampling_points.size(); i++)
				distance = std::min(distance, CGAL_3D_Distance_Point_Point(sampling_points[i], Vector3d(x, y, z)));

			if (distance > minimal_d)
			{
				sampling_points.push_back(Vector3d(x, y, z));

				if (sampling_points.size() % 100 == 0) std::cerr<< sampling_points.size()<<" ";

				run = 0;
			}
		}
	}
	std::cerr << std::endl;
}


extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_Vector3d(const Vector3d1& ps, const Vector3d1& ns, const std::string& path, Vector3d1& inters)
{
	std::ifstream input(path);
	Mesh mesh;
	input >> mesh;
	Mesh_Tree tree(faces(mesh).first, faces(mesh).second, mesh);

	for (int i = 0; i < ps.size(); i++)
	{
		Ray_3 ray(Point_3(ps[i][0], ps[i][1], ps[i][2]), Vector_3(ns[i][0], ns[i][1], ns[i][2]));

		std::list<Mesh_Ray_intersection> intersections;
		tree.all_intersections(ray, std::back_inserter(intersections));

		bool b = false;

		double min_d = 100000000000.0;

		Vector3d near_p;
		for (auto iter = intersections.begin(); iter != intersections.end(); iter++)
		{
			const Point_3* p = boost::get<Point_3>(&(iter->value().first));

			if (p)
			{
				Vector3d v(p->x(), p->y(), p->z());
				double d = CGAL_3D_Distance_Point_Point(v, ps[i]);

				if (min_d > d)
				{
					min_d = d;
					near_p = v;
				}
				b = true;
			}
		}
		if (!b)
		{
			inters.push_back(ps[i]);
		}
		else
		{
			inters.push_back(near_p);
		}
	}
}
