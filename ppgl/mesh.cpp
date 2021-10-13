#include "geom.h"

#include <Mathematics/MeshCurvature.h>

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


Poly_point_3 Minus(Poly_point_3 a, Poly_point_3 b)
{
	return Poly_point_3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
float Dot(Poly_point_3 a, Poly_point_3 b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void Barycentric(Poly_point_3 p, Poly_point_3 a, Poly_point_3 b, Poly_point_3 c, double& u, double& v, double& w)
{
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

void CGAL_Barycentric(Vector3d p, Vector3d a, Vector3d b, Vector3d c, double &u, double &v, double &w) 
{
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
	//gte::Sphere3<double> sphere(gte::Vector3d(center_x, center_y, center_z), radius);
	//gte::Ray3<double> ray(gte::Vector3d(ray_origin_x, ray_origin_y, ray_origin_z), gte::Vector3d(ray_direction_x, ray_direction_y, ray_direction_z));

	//gte::IntrRay3Sphere3d intr(ray,sphere);

	//intr.Test();
	//intr.Find();

	//int nb = intr.GetQuantity();

	//for (int i = 0; i < nb; i++)
	//{
	//	gte::Vector3d p = intr.GetPoint(i);
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


//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C1(const std::string & outside_path, const double& d, Vector3d1 & sampling_points)
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

	double x(out_minC[0]);
	while (x < out_maxC[0])
	{
		double y(out_minC[1]);
		while (y < out_maxC[1])
		{
			double z(out_minC[2]);
			while (z < out_maxC[2])
			{
				CGAL::Bounded_side res = out_checker(Point_3(x, y, z));
				if (res == CGAL::ON_BOUNDED_SIDE)
				{
					sampling_points.push_back(Vector3d(x, y, z));
					if (sampling_points.size() % 100 == 0) std::cerr << sampling_points.size() << " ";
				}

				z += minimal_d;
			}
			y += minimal_d;
		}

		x += minimal_d;
	}

	std::cerr << std::endl;
}

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C2(const std::string & outside_path, const std::string & inside_path, const double& d, Vector3d1 & sampling_points)
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

	double x(out_minC[0]);
	while (x < out_maxC[0])
	{
		double y(out_minC[1]);
		while (y < out_maxC[1])
		{
			double z(out_minC[2]);
			while (z < out_maxC[2])
			{
				CGAL::Bounded_side out_res = out_checker(Point_3(x, y, z));
				CGAL::Bounded_side in_res = in_checker(Point_3(x, y, z));
				if (out_res == CGAL::ON_BOUNDED_SIDE && in_res == CGAL::ON_UNBOUNDED_SIDE)
				{
					sampling_points.push_back(Vector3d(x, y, z));
					if (sampling_points.size() % 100 == 0) std::cerr << sampling_points.size() << " ";
				}

				z += minimal_d;
			}
			y += minimal_d;
		}

		x += minimal_d;
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




extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Triangle(const Vector3d & p, const Vector3d & t_0, const Vector3d & t_1, const Vector3d & t_2)
{
	KC::Point_3 a(t_0[0], t_0[1], t_0[2]);
	KC::Point_3 b(t_1[0], t_1[1], t_1[2]);
	KC::Point_3 c(t_2[0], t_2[1], t_2[2]);

	std::list<Triangle_3> triangles;
	triangles.push_back(Triangle_3(a, b, c));

	Tree_3 tree(triangles.begin(), triangles.end());
	KC::Point_3 point_query(p[0], p[1], p[2]);

	return sqrt(tree.squared_distance(point_query));
}
extern "C" PPGL_EXPORT double CGAL_3D_Distance_Point_Triangles(const Vector3d & p, const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2)
{
	//build polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	//build tree
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	return sqrt(tree.squared_distance(Point_3(p[0], p[1], p[2])));
}
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Nearest_Point_Triangles(const Vector3d & p, const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2)
{
	//build polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	//build tree
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Point_3 c_p = tree.closest_point(Point_3(p[0], p[1], p[2]));
	return Vector3d(c_p[0], c_p[1], c_p[2]);
}

extern "C" PPGL_EXPORT void CGAL_3D_Distance_Point_Mesh(const std::string & path, const Vector3d1 & query_points, std::vector<double>&distances)
{
	std::cout << "CGAL_3D_Distance_Point_Mesh" << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < query_points.size(); i++)
	{
		if (i % (query_points.size() / 10) == 0)
		{
			std::cout << (double)i / (double)query_points.size() << std::endl;
		}
		distances.push_back(sqrt(tree.squared_distance(VectorPoint3d(query_points[i]))));
	}
}

extern "C" PPGL_EXPORT void CGAL_3D_Neareast_Point_Mesh(const std::string & path, const Vector3d1 & ves, Vector3d1 & ners)
{
	std::cout << "CGAL_3D_Distance_Point_Mesh" << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < ves.size(); i++)
	{
		if (i % (ves.size() / 10) == 0)
		{
			std::cout << (double)i / (double)ves.size() << std::endl;
		}
		Point_3 p = tree.closest_point(Point_3(ves[i][0], ves[i][1], ves[i][2]));

		ners.push_back(Vector3d(p[0], p[1], p[2]));
	}
}


extern "C" PPGL_EXPORT void  CGAL_3D_Mesh_Near_Triangles(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const Vector3d1 & points, const double& d, std::vector<std::vector<int>>&triangles)
{
	//build polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	//build tree
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < points.size(); i++)
	{
		std::vector<int> triangle;

		Poly_point_3 query(points[i][0], points[i][1], points[i][2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

		std::priority_queue<Polyhedron_3::Facet_handle> facets;
		std::vector<int> save_index;
		facets.push(pp.second);
		save_index.push_back(pp.second->id());

		while (facets.size() != 0)
		{
			Polyhedron_3::Facet_handle fh = facets.top();
			triangle.push_back(fh->id());
			facets.pop();

			std::vector<Polyhedron_3::Facet_handle> neighbors;

			if (!fh->halfedge()->is_border_edge())
				neighbors.push_back(fh->halfedge()->opposite()->face());
			if (!fh->halfedge()->next()->is_border_edge())
				neighbors.push_back(fh->halfedge()->next()->opposite()->face());
			if (!fh->halfedge()->next()->next()->is_border_edge())
				neighbors.push_back(fh->halfedge()->next()->next()->opposite()->face());

			for (int j = 0; j < neighbors.size(); j++)
			{
				Polyhedron_3::Facet_handle n_fh = neighbors[j];
				std::vector<Poly_point_3> n_fh_vecs;
				n_fh_vecs.push_back(n_fh->halfedge()->vertex()->point());
				n_fh_vecs.push_back(n_fh->halfedge()->next()->vertex()->point());
				n_fh_vecs.push_back(n_fh->halfedge()->next()->next()->vertex()->point());

				bool add_bool = false;
				for (int k = 0; k < 3; k++)
				{
					double distance = sqrt((double)CGAL::squared_distance(pp.first, n_fh_vecs[k]));
					if (distance < d) {
						add_bool = true;
						break;
					}
				}
				if (add_bool && !(std::find(save_index.begin(), save_index.end(), n_fh->id()) != save_index.end()))
				{
					facets.push(n_fh);
					save_index.push_back(n_fh->id());
				}
			}
		}

		triangles.push_back(triangle);
	}
}



extern "C" PPGL_EXPORT void CGAL_3D_Points_inside_Triangles_C1(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const Vector3d1 & points, std::vector<bool>&insides)
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

extern "C" PPGL_EXPORT void CGAL_3D_Points_inside_Triangles_C2(const std::string & path, const Vector3d1 & points, std::vector<bool>&insides)
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



extern "C" PPGL_EXPORT void CGAL_Mesh_Loop_Subdivision_One_Step(Vector3d1 & vecs, std::vector<int>&face_id_0, std::vector<int>&face_id_1, std::vector<int>&face_id_2)
{
	Vector3d1 loop_vecs = vecs;
	std::vector<int> loop_face_id_0;
	std::vector<int> loop_face_id_1;
	std::vector<int> loop_face_id_2;

	//edges
	std::vector<Edge> edges;

	std::vector<std::vector<int>> vecs_neighbors(vecs.size(), std::vector<int>());
	std::vector<std::vector<int>> vecs_neighbors_labels(vecs.size(), std::vector<int>());

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		edges.push_back(Edge(index_0, index_1));
		edges.push_back(Edge(index_1, index_0));
		edges.push_back(Edge(index_0, index_2));
		edges.push_back(Edge(index_2, index_0));
		edges.push_back(Edge(index_1, index_2));
		edges.push_back(Edge(index_2, index_1));
	}

	for (int i = 0; i < edges.size(); i++)
	{
		int source = edges[i].source;
		int end = edges[i].end;

		bool b = false;
		for (int j = 0; j < vecs_neighbors[source].size(); j++)
		{
			if (vecs_neighbors[source][j] == end)
			{
				b = true;
				break;
			}
		}
		if (!b)
		{
			vecs_neighbors[source].push_back(end);
			vecs_neighbors_labels[source].push_back(-1);
		}
	}
	std::vector<Edge>().swap(edges);

	for (int i = 0; i < vecs.size(); i++)
	{
		int  source = i;
		for (int j = 0; j < vecs_neighbors[i].size(); j++)
		{
			int end = vecs_neighbors[i][j];
			if (vecs_neighbors_labels[i][j] < 0)
			{
				edges.push_back(Edge(source, end));
				vecs_neighbors_labels[i][j] = edges.size() - 1;

				for (int k = 0; k < vecs_neighbors[end].size(); k++)
				{
					if (vecs_neighbors[end][k] == source)
					{
						vecs_neighbors_labels[end][k] = edges.size() - 1;
					}
				}
			}
		}
	}

	//loop_vecs
	Vector3d1 edge_middle_points;
	for (int i = 0; i < edges.size(); i++)
		edge_middle_points.push_back((vecs[edges[i].source] + vecs[edges[i].end]) / (double)2.0);
	for (int i = 0; i < edge_middle_points.size(); i++) loop_vecs.push_back(edge_middle_points[i]);

	//loop faces
	std::vector<Edge> face_edges;
	std::vector<int> face_edges_id;

	for (int i = 0; i < face_id_0.size(); i++) {
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		face_edges.push_back(Edge(index_0, index_1));
		face_edges.push_back(Edge(index_1, index_2));
		face_edges.push_back(Edge(index_2, index_0));
	}
	for (int i = 0; i < face_edges.size(); i++)
	{
		int source = face_edges[i].source;
		int end = face_edges[i].end;

		for (int j = 0; j < vecs_neighbors[source].size(); j++)
		{
			if (vecs_neighbors[source][j] == end)
			{
				face_edges_id.push_back(vecs_neighbors_labels[source][j]);
				break;
			}
		}
	}

	//     0
	//    2 0
	//  2  1  1
	for (int i = 0; i < face_edges.size(); i = i + 3)
	{
		int index_0 = face_edges[i].source;
		int index_1 = face_edges[i + 1].source;
		int index_2 = face_edges[i + 2].source;

		int edge_id_0 = face_edges_id[i] + vecs.size();
		int edge_id_1 = face_edges_id[i + 1] + vecs.size();
		int edge_id_2 = face_edges_id[i + 2] + vecs.size();

		loop_face_id_0.push_back(index_0);
		loop_face_id_1.push_back(edge_id_0);
		loop_face_id_2.push_back(edge_id_2);

		loop_face_id_0.push_back(edge_id_0);
		loop_face_id_1.push_back(index_1);
		loop_face_id_2.push_back(edge_id_1);

		loop_face_id_0.push_back(edge_id_2);
		loop_face_id_1.push_back(edge_id_1);
		loop_face_id_2.push_back(index_2);

		loop_face_id_0.push_back(edge_id_2);
		loop_face_id_1.push_back(edge_id_0);
		loop_face_id_2.push_back(edge_id_1);
	}

	//release
	Vector3d1().swap(vecs);
	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);

	vecs = loop_vecs;
	face_id_0 = loop_face_id_0;
	face_id_1 = loop_face_id_1;
	face_id_2 = loop_face_id_2;

	Vector3d1().swap(loop_vecs);
	std::vector<int>().swap(loop_face_id_0);
	std::vector<int>().swap(loop_face_id_1);
	std::vector<int>().swap(loop_face_id_2);

	Vector3d1().swap(edge_middle_points);
	std::vector<Edge>().swap(face_edges);
	std::vector<int>().swap(face_edges_id);
	std::vector<Edge>().swap(edges);
}


extern "C" PPGL_EXPORT void CGAL_Mesh_Subdivision(const std::string & in_path, const std::string & sub, const int& step, const std::string & out_path)
{
	//load the input obj
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, in_path);
	//subdivision

	if (sub == "Loop" || sub == "loop" || sub == "l" || sub == "L")
		CGAL::Subdivision_method_3::Loop_subdivision(polyhedron, step);

	if (sub == "Sqrt" || sub == "sqrt" || sub == "s" || sub == "S")
		CGAL::Subdivision_method_3::Sqrt3_subdivision(polyhedron, step);

	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
		iter != polyhedron.vertices_end(); iter++)
	{
		Poly_point_3 p = iter->point();
		vecs.push_back(Vector3d(p[0], p[1], p[2]));
	}

	for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
	{
		face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
		face_id_1.push_back(iter->halfedge()->vertex()->id());
		face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
	}

	CGAL_Output_Obj_C2(out_path, vecs, face_id_0, face_id_1, face_id_2);
}



extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C1(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, std::vector<double>&max_curs, std::vector<double>&min_curs)
{
	int verticeSize = vecs.size();
	int faceindiceSize = face_id_0.size() * 3;

	gte::Vector3<double>* points = new gte::Vector3<double>[verticeSize];
	unsigned int* indices = new unsigned int[faceindiceSize];

	for (int i = 0; i < verticeSize; i++)
	{
		points[i][0] = vecs[i][0];
		points[i][1] = vecs[i][1];
		points[i][2] = vecs[i][2];
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		indices[3 * i] = face_id_0[i];
		indices[3 * i + 1] = face_id_1[i];
		indices[3 * i + 2] = face_id_2[i];
	}
	gte::MeshCurvature<double> meshCurvature;

	meshCurvature((size_t)verticeSize, points, (size_t)(faceindiceSize / 3), indices, 1e-06f);

	for (int i = 0; i < verticeSize; i++)
	{
		max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
		min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);
	}
}
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C2(const Vector3d1 & vecs, const std::vector<std::vector<int>>&face_ids, std::vector<double>&max_curs, std::vector<double>&min_curs)
{
	std::vector<int> face_id_0, face_id_1, face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}

	CGAL_3D_Mesh_Curvature_C1(vecs, face_id_0, face_id_1, face_id_2, max_curs, min_curs);

	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);

}
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C3(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, std::vector<double>&max_curs, std::vector<double>&min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions)
{
	int verticeSize = vecs.size();
	int faceindiceSize = face_id_0.size() * 3;

	gte::Vector3<double>* points = new gte::Vector3<double>[verticeSize];
	unsigned int* indices = new unsigned int[faceindiceSize];

	for (int i = 0; i < verticeSize; i++)
	{
		points[i][0] = vecs[i][0];
		points[i][1] = vecs[i][1];
		points[i][2] = vecs[i][2];
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		indices[3 * i] = face_id_0[i];
		indices[3 * i + 1] = face_id_1[i];
		indices[3 * i + 2] = face_id_2[i];
	}

	gte::MeshCurvature<double> meshCurvature;
	meshCurvature((size_t)verticeSize, points, (size_t)(faceindiceSize / 3), indices, 1e-06f);

	for (int i = 0; i < verticeSize; i++)
	{
		max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
		min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);

		gte::Vector3<double> max_curs_direction = meshCurvature.GetMaxDirections()[i];
		gte::Vector3<double> min_curs_direction = meshCurvature.GetMinDirections()[i];

		max_curs_directions.push_back(Vector3d(max_curs_direction[0], max_curs_direction[1], max_curs_direction[2]));
		min_curs_directions.push_back(Vector3d(min_curs_direction[0], min_curs_direction[1], min_curs_direction[2]));
	}
}
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C4(const Vector3d1 & vecs, const std::vector<std::vector<int>>&face_ids, std::vector<double>&max_curs, std::vector<double>&min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions)
{
	std::vector<int> face_id_0, face_id_1, face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}

	CGAL_3D_Mesh_Curvature_C3(vecs, face_id_0, face_id_1, face_id_2, max_curs, min_curs, max_curs_directions, min_curs_directions);

	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);
}
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C5(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, std::vector<double>&max_curs, std::vector<double>&min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions, Vector3d1 & normals)
{
	int verticeSize = vecs.size();
	int faceindiceSize = face_id_0.size() * 3;

	gte::Vector3<double>* points = new gte::Vector3<double>[verticeSize];
	unsigned int* indices = new unsigned int[faceindiceSize];

	for (int i = 0; i < verticeSize; i++)
	{
		points[i][0] = vecs[i][0];
		points[i][1] = vecs[i][1];
		points[i][2] = vecs[i][2];
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		indices[3 * i] = face_id_0[i];
		indices[3 * i + 1] = face_id_1[i];
		indices[3 * i + 2] = face_id_2[i];
	}
	gte::MeshCurvature<double> meshCurvature;
	meshCurvature((size_t)verticeSize, points, (size_t)(faceindiceSize / 3), indices, 1e-06f);

	for (int i = 0; i < verticeSize; i++)
	{
		max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
		min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);

		gte::Vector3<double> max_curs_direction = meshCurvature.GetMaxDirections()[i];
		gte::Vector3<double> min_curs_direction = meshCurvature.GetMinDirections()[i];

		max_curs_directions.push_back(Vector3d(max_curs_direction[0], max_curs_direction[1], max_curs_direction[2]));
		min_curs_directions.push_back(Vector3d(min_curs_direction[0], min_curs_direction[1], min_curs_direction[2]));

		gte::Vector3<double> normal = meshCurvature.GetNormals()[i];
		normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
	}
}

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Curvature_C6(const Vector3d1 & vecs, const std::vector<std::vector<int>>&face_ids, std::vector<double>&max_curs, std::vector<double>&min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions, Vector3d1 & normals)
{
	std::vector<int> face_id_0, face_id_1, face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}

	CGAL_3D_Mesh_Curvature_C5(vecs, face_id_0, face_id_1, face_id_2, max_curs, min_curs, max_curs_directions, min_curs_directions, normals);

	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);
}

extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C1(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, std::vector<bool>&bools)
{
	std::vector<bool>().swap(bools);

	std::vector<std::vector<int>> vecs_neigbor(vecs.size(), std::vector<int>());
	std::vector<std::vector<int>> vecs_neigbor_lable(vecs.size(), std::vector<int>());
	std::vector<int> edges;
	for (int i = 0; i < face_id_0.size(); i++) {
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		edges.push_back(index_0);
		edges.push_back(index_1);
		edges.push_back(index_1);
		edges.push_back(index_2);
		edges.push_back(index_2);
		edges.push_back(index_0);

		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
	}

	for (int i = 0; i < edges.size(); i = i + 2) {
		int index_0 = edges[i];
		int index_1 = edges[i + 1];

		int search_0 = -1;
		for (int j = 0; j < vecs_neigbor[index_0].size() && search_0 < 0; j++) {
			if (vecs_neigbor[index_0][j] == index_1) {
				search_0 = j;
				vecs_neigbor_lable[index_0][j]++;
			}
		}

		if (search_0 < 0) {
			vecs_neigbor[index_0].push_back(index_1);
			vecs_neigbor_lable[index_0].push_back(1);
		}
	}

	for (int i = 0; i < vecs.size(); i++) {
		bool b = false;
		for (int j = 0; j < vecs_neigbor_lable[i].size() & !b; j++) {
			if (vecs_neigbor_lable[i][j] == 1) {
				b = true;
			}
		}
		bools.push_back(b);
	}

	std::vector<std::vector<int>>().swap(vecs_neigbor);
	std::vector<std::vector<int>>().swap(vecs_neigbor_lable);
	std::vector<int>().swap(edges);
}
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C2(const std::string & path, std::vector<bool>&bools)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_3D_Triangle_Mesh_Boundary_C1(vecs, face_id_0, face_id_1, face_id_2, bools);
}

extern "C" PPGL_EXPORT void CGAL_3D_Connecting_Segments_C1(Vector2d2 & segments, Vector2d2 & lines)
{
	//save connecting relations
	std::vector<bool> used(segments.size(), false);

	std::vector<int> relations;
#pragma region get_relations
	for (int i = 0; i < segments.size(); i++)
	{
		for (int j = 0; j < segments.size(); j++)
		{
			if (i != j && !used[i] && !used[j])
			{
				bool b_0_0 = Functs::IsAlmostZero_Double(CGAL_2D_Distance_Point_Point(segments[i][0], segments[j][0]), 1.0E-09);
				bool b_0_1 = Functs::IsAlmostZero_Double(CGAL_2D_Distance_Point_Point(segments[i][0], segments[j][1]), 1.0E-09);
				bool b_1_0 = Functs::IsAlmostZero_Double(CGAL_2D_Distance_Point_Point(segments[i][1], segments[j][0]), 1.0E-09);
				bool b_1_1 = Functs::IsAlmostZero_Double(CGAL_2D_Distance_Point_Point(segments[i][1], segments[j][1]), 1.0E-09);

				if ((b_0_0 && b_1_1) || (b_0_1 && b_1_0))
				{
					used[j] = true;
					continue;
				}

				if (b_0_0)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_0_1)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
				if (b_1_0)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_1_1)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
			}
		}
	}
#pragma endregion

	std::vector<std::vector<int>> ones;


	while (true)
	{
		int index = -1;
		int end = -1;

		for (int i = 0; i < segments.size(); i++)
		{
			if (!used[i]) {
				index = i;
				end = 0;
				used[i] = true;
				break;
			}
		}

		if (index < 0)break;

		Vector2d1 line(1, segments[index][end]);

		std::vector<int> one(1, index);

		while (true)
		{
			end = 1 - end;
			bool search = false;
			for (int i = 0; i < relations.size(); i = i + 4)
			{
				if (relations[i] == index && relations[i + 1] == end && !used[relations[i + 2]])
				{
					line.push_back(segments[relations[i + 2]][relations[i + 3]]);
					one.push_back(relations[i + 2]);
					index = relations[i + 2];
					end = relations[i + 3];
					used[index] = true;
					search = true;
					break;
				}
				if (relations[i + 2] == index && relations[i + 3] == end && !used[relations[i]])
				{
					line.push_back(segments[relations[i]][relations[i + 1]]);
					one.push_back(relations[i]);
					index = relations[i];
					end = relations[i + 1];
					used[index] = true;
					search = true;
					break;
				}
			}
			if (!search) { break; }
		}

		ones.push_back(one);
		lines.push_back(line);
	}
}
extern "C" PPGL_EXPORT void CGAL_3D_Connecting_Segments_C2(Vector3d2 & segments, Vector3d2 & lines)
{
	//save connecting relations
	std::vector<bool> used(segments.size(), false);

	std::vector<int> relations;
#pragma region get_relations
	for (int i = 0; i < segments.size(); i++)
	{
		for (int j = 0; j < segments.size(); j++)
		{
			if (i != j && !used[i] && !used[j])
			{
				bool b_0_0 = Functs::IsAlmostZero_Double(CGAL_3D_Distance_Point_Point(segments[i][0], segments[j][0]), 1.0E-09);
				bool b_0_1 = Functs::IsAlmostZero_Double(CGAL_3D_Distance_Point_Point(segments[i][0], segments[j][1]), 1.0E-09);
				bool b_1_0 = Functs::IsAlmostZero_Double(CGAL_3D_Distance_Point_Point(segments[i][1], segments[j][0]), 1.0E-09);
				bool b_1_1 = Functs::IsAlmostZero_Double(CGAL_3D_Distance_Point_Point(segments[i][1], segments[j][1]), 1.0E-09);

				if ((b_0_0 && b_1_1) || (b_0_1 && b_1_0))
				{
					used[j] = true;
					continue;
				}

				if (b_0_0)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_0_1)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
				if (b_1_0)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_1_1)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
			}
		}
	}
#pragma endregion

	std::vector<std::vector<int>> ones;


	while (true)
	{
		int index = -1;
		int end = -1;

		for (int i = 0; i < segments.size(); i++)
		{
			if (!used[i]) {
				index = i;
				end = 0;
				used[i] = true;
				break;
			}
		}

		if (index < 0)break;

		Vector3d1 line(1, segments[index][end]);

		std::vector<int> one(1, index);

		while (true)
		{
			end = 1 - end;
			bool search = false;
			for (int i = 0; i < relations.size(); i = i + 4)
			{
				if (relations[i] == index && relations[i + 1] == end && !used[relations[i + 2]])
				{
					line.push_back(segments[relations[i + 2]][relations[i + 3]]);
					one.push_back(relations[i + 2]);
					index = relations[i + 2];
					end = relations[i + 3];
					used[index] = true;
					search = true;
					break;
				}
				if (relations[i + 2] == index && relations[i + 3] == end && !used[relations[i]])
				{
					line.push_back(segments[relations[i]][relations[i + 1]]);
					one.push_back(relations[i]);
					index = relations[i];
					end = relations[i + 1];
					used[index] = true;
					search = true;
					break;
				}
			}
			if (!search) { break; }
		}

		ones.push_back(one);
		lines.push_back(line);
	}
}

extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C3(Vector3d1 & vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, Vector3d2 & boundaries)
{
	Vector3d2 segments;

	std::vector<std::vector<int>> vecs_neigbor(vecs.size(), std::vector<int>());
	std::vector<std::vector<int>> vecs_neigbor_lable(vecs.size(), std::vector<int>());
	std::vector<int> edges;
	for (int i = 0; i < face_id_0.size(); i++) {
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		edges.push_back(index_0);
		edges.push_back(index_1);
		edges.push_back(index_1);
		edges.push_back(index_2);
		edges.push_back(index_2);
		edges.push_back(index_0);

		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
	}

	for (int i = 0; i < edges.size(); i = i + 2) {
		int index_0 = edges[i];
		int index_1 = edges[i + 1];

		int search_0 = -1;
		for (int j = 0; j < vecs_neigbor[index_0].size() && search_0 < 0; j++) {
			if (vecs_neigbor[index_0][j] == index_1) {
				search_0 = j;
				vecs_neigbor_lable[index_0][j]++;
			}
		}

		if (search_0 < 0) {
			vecs_neigbor[index_0].push_back(index_1);
			vecs_neigbor_lable[index_0].push_back(1);
		}
	}

	for (int i = 0; i < vecs.size(); i++)
	{
		int index_0 = i;

		for (int j = 0; j < vecs_neigbor_lable[i].size(); j++)
		{
			if (vecs_neigbor_lable[i][j] == 1)
			{
				int index_1 = vecs_neigbor[i][j];

				Vector3d1 segment;
				segment.push_back(vecs[index_0]);
				segment.push_back(vecs[index_1]);

				segments.push_back(segment);

				//delete
				for (int k = 0; k < vecs_neigbor[index_1].size(); k++)
				{
					if (vecs_neigbor[index_1][k] == index_0)
					{
						vecs_neigbor_lable[index_1][k] = 0;
					}
				}

			}
		}
	}

	CGAL_3D_Connecting_Segments_C2(segments, boundaries);

	for (int i = 0; i < segments.size(); i++)
		Vector3d1().swap(segments[i]);
	Vector3d2().swap(segments);
	std::vector<std::vector<int>>().swap(vecs_neigbor);
	std::vector<std::vector<int>>().swap(vecs_neigbor_lable);
	std::vector<int>().swap(edges);
}
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C4(Vector3d1 & vecs, std::vector<std::vector<int>>&face_ids, Vector3d2 & boundaries)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (auto& face : face_ids)
	{
		face_id_0.emplace_back(face[0]);
		face_id_1.emplace_back(face[1]);
		face_id_2.emplace_back(face[2]);
	}

	CGAL_3D_Triangle_Mesh_Boundary_C3(vecs, face_id_0, face_id_1, face_id_2, boundaries);
}
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C5(std::string path, Vector3d2 & boundaries)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_3D_Triangle_Mesh_Boundary_C3(vecs, face_id_0, face_id_1, face_id_2, boundaries);
}


extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth_C1(const std::string & in_path, const std::string & out_path, const int laplace_nb)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(in_path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_Mesh_Laplace_Smooth_C2(vecs, face_id_0, face_id_1, face_id_2, laplace_nb);
	CGAL_Output_Obj_C2(out_path, vecs, face_id_0, face_id_1, face_id_2);
}

extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Vecs_Neighbors(Vector3d1 & vecs, std::vector<int>&face_id_0, std::vector<int>&face_id_1, std::vector<int>&face_id_2, std::vector<std::vector<int>>&neighs)
{
	for (int i = 0; i < vecs.size(); i++)
		neighs.push_back(std::vector<int>());

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int id_0 = face_id_0[i];
		int id_1 = face_id_1[i];
		int id_2 = face_id_2[i];

		neighs[id_0].push_back(id_1);
		neighs[id_0].push_back(id_2);

		neighs[id_1].push_back(id_0);
		neighs[id_1].push_back(id_2);

		neighs[id_2].push_back(id_0);
		neighs[id_2].push_back(id_1);
	}

	for (int i = 0; i < neighs.size(); i++)
	{
		sort(neighs[i].begin(), neighs[i].end());
		neighs[i].erase(unique(neighs[i].begin(), neighs[i].end()), neighs[i].end());
	}
}

extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth_C2(Vector3d1 & vecs, std::vector<int>&face_id_0, std::vector<int>&face_id_1, std::vector<int>&face_id_2, const int laplace_nb)
{
	std::vector<bool> vertices_boundary;
	CGAL_3D_Triangle_Mesh_Boundary_C1(vecs, face_id_0, face_id_1, face_id_2, vertices_boundary);
	std::vector<std::vector<int>> neighs;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(vecs, face_id_0, face_id_1, face_id_2, neighs);

	for (int iter = 0; iter < laplace_nb; iter++)
	{
		Vector3d1 new_vecs;
		for (int i = 0; i < vecs.size(); i++)
		{
			if (!vertices_boundary[i])
			{
				Vector3d v(0.0, 0.0, 0.0);
				double w = 0.0;
				for (int j = 0; j < neighs[i].size(); j++)
				{
					double d = CGAL_3D_Distance_Point_Point(vecs[neighs[i][j]], vecs[i]);
					v += vecs[neighs[i][j]] * (double)(1.0 / d);
					w += (1.0 / d);
				}
				v = vecs[i] * (double)0.5 + (double)0.5 * v / (double)w;

				new_vecs.push_back(v);
			}
			else
			{
				new_vecs.push_back(vecs[i]);
			}

		}
		vecs = new_vecs;
	}
}


extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Vecs_Faces(Vector3d1 & vecs, std::vector<int>&face_id_0, std::vector<int>&face_id_1, std::vector<int>&face_id_2,
	std::vector<std::vector<int>>&surface_vectices_to_face)
{
	surface_vectices_to_face = std::vector<std::vector<int>>(vecs.size(), std::vector<int>());

	std::vector<std::vector<int>> sets(vecs.size(), std::vector<int>());
	for (int i = 0; i < face_id_0.size(); i++)
	{
		//surface_vectices_to_face
		sets[face_id_0[i]].emplace_back(i);
		sets[face_id_1[i]].emplace_back(i);
		sets[face_id_2[i]].emplace_back(i);
	}

	for (int i = 0; i < vecs.size(); i++)
	{
		set<int>s(sets[i].begin(), sets[i].end());
		vector<int> vec;
		vec.assign(s.begin(), s.end());
		surface_vectices_to_face[i] = vec;
	}

}

extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges(Vector3d1 & vecs, std::vector<int>&face_id_0, std::vector<int>&face_id_1, std::vector<int>&face_id_2,
	std::vector<std::vector<std::vector<int>>>&surface_vectices_to_neighbor_edges)
{
	std::vector<std::vector<int>> surface_vectices_to_face;
	CGAL_3D_Triangle_Mesh_Vecs_Faces(vecs, face_id_0, face_id_1, face_id_2, surface_vectices_to_face);

	for (int i = 0; i < vecs.size(); i++)
	{
		int vertice_id = i;

		std::vector<std::vector<int>> edges;

		for (int j = 0; j < surface_vectices_to_face[i].size(); j++)
		{
			int surface_id = surface_vectices_to_face[i][j];

			std::vector<int> face;
			face.push_back(face_id_0[surface_id]);
			face.push_back(face_id_1[surface_id]);
			face.push_back(face_id_2[surface_id]);

			for (int k = 0; k < face.size(); k++)
			{
				if (face[k] == vertice_id)
				{
					int vertice_id_0 = face[(k + 1) % 3];
					int vertice_id_1 = face[(k + 2) % 3];

					std::vector<int> edge;
					edge.push_back(vertice_id_0);
					edge.push_back(vertice_id_1);
					edges.push_back(edge);
					break;
				}
			}
		}
		surface_vectices_to_neighbor_edges.push_back(edges);
	}
}


extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth_by_Curvature(Vector3d1 & vecs, std::vector<int>&face_id_0, std::vector<int>&face_id_1, std::vector<int>&face_id_2, double& low_curvature)
{
	std::vector<double> max_curvature;
	std::vector<double> min_curvature;
	Vector3d1 max_curvature_direction;
	Vector3d1 min_curvature_direction;

	std::vector<bool> boundary;
	CGAL_3D_Triangle_Mesh_Boundary_C1(vecs, face_id_0, face_id_1, face_id_2, boundary);

	Vector3d1 vecs_normals;

	int last_nb = 0;

	std::vector<std::vector<int>> vecs_neighbors;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(vecs, face_id_0, face_id_1, face_id_2, vecs_neighbors);

	std::vector<std::vector<std::vector<int>>> surface_vectices_to_neighbor_edges;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges(vecs, face_id_0, face_id_1, face_id_2, surface_vectices_to_neighbor_edges);

	int stop = 0;

	double par_0 = 0.1;
	double par_1 = 0.6;
	double par_2 = 0.3;

	for (int iter = 0; iter < 500; iter++)
		//while (true)
	{
		//compute vertices curvature
		std::vector<double>().swap(max_curvature);
		std::vector<double>().swap(min_curvature);
		Vector3d1().swap(max_curvature_direction);
		Vector3d1().swap(min_curvature_direction);
		Vector3d1().swap(vecs_normals);

		CGAL_3D_Mesh_Curvature_C5(vecs, face_id_0, face_id_1, face_id_2, max_curvature, min_curvature, max_curvature_direction, min_curvature_direction, vecs_normals);

		int nb = 0;
		double minimal_cur = 100000.0;
		for (int i = 0; i < vecs.size(); i++)
		{
			if (min_curvature[i] < low_curvature && !boundary[i])
			{
				nb++;
				minimal_cur = std::min(minimal_cur, min_curvature[i]);
			}
		}


		//terminal condition

		if (nb < 5)break;

		if (abs(last_nb - nb) < 2)
			stop++;
		else
			stop = 0;

		if (stop == 60)
		{
			//break;
			par_0 = 0.10;
			par_1 = 0.65;
			par_2 = 0.25;
		}

		if (stop == 100)
		{
			par_0 = 0.10;
			par_1 = 0.70;
			par_2 = 0.20;
		}

		if (stop == 300)
		{
			break;
		}

		std::cout << "Current low curvature points number: " << stop << " " << nb << " " << minimal_cur << std::endl;

		last_nb = nb;

		//one iteration
		Vector3d1 iteration_vecs;

		for (int i = 0; i < vecs.size(); i++)
		{
			if (boundary[i])
			{
				iteration_vecs.push_back(vecs[i]);
			}
			else
			{
				bool run = min_curvature[i] < low_curvature + 0.1;

				for (int j = 0; j < vecs_neighbors[i].size() && !run; j++)
				{
					if (min_curvature[vecs_neighbors[i][j]] < low_curvature)
					{
						run = true;
					}

					if (true)
					{
						int vertice_id = vecs_neighbors[i][j];

						for (int k = 0; k < vecs_neighbors[vertice_id].size() && !run; k++)
						{
							if (min_curvature[vecs_neighbors[vertice_id][k]] < low_curvature)
							{
								run = true;
							}
						}
					}
				}

				if (run)
				{
					Vector3d cur_v = vecs[i] + Functs::SetVectorLength(vecs_normals[i], 0.001 * min_curvature[i] / low_curvature);

					Vector3d smooth_v(0.0, 0.0, 0.0);

					double weight = 0.0;
					for (int j = 0; j < vecs_neighbors[i].size(); j++)
					{
						double l = CGAL_3D_Distance_Point_Point(vecs[vecs_neighbors[i][j]], vecs[i]);
						smooth_v += vecs[vecs_neighbors[i][j]] * (double)l;
						weight += l;
					}
					smooth_v = smooth_v / (double)weight;

					if (min_curvature[i] < 0.0 && max_curvature[i]>0.0)
					{
						//smooth_v = vecs[i];
						iteration_vecs.push_back((double)par_0 * vecs[i] + (double)par_1 * cur_v + (double)par_2 * smooth_v);
					}
					else
						iteration_vecs.push_back((double)par_0 * vecs[i] + (double)par_1 * cur_v + (double)par_2 * smooth_v);

					//iteration_vecs.push_back(LaplaceMeshSmoothForOnePoint(i, vecs, min_curvature_direction, vecs_normals, surface_vectices_to_neighbor_edges));
				}
				else
				{
					iteration_vecs.push_back(vecs[i]);
				}
			}
		}

		Vector3d1().swap(vecs);
		vecs = iteration_vecs;
		Vector3d1().swap(iteration_vecs);
	}
}

extern "C" PPGL_EXPORT void CGAL_Mesh_Loop_Subdivision_Own_Version(const std::string & in_path, const int& step, const std::string & out_path, const int& laplace_nb)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	CGAL_3D_Read_Triangle_Mesh(in_path, vecs, face_id_0, face_id_1, face_id_2);

	for (int i = 0; i < step; i++)
	{
		CGAL_Mesh_Loop_Subdivision_One_Step(vecs, face_id_0, face_id_1, face_id_2);
		CGAL_Mesh_Laplace_Smooth_C2(vecs, face_id_0, face_id_1, face_id_2, laplace_nb);
	}

	CGAL_Output_Obj_C2(out_path, vecs, face_id_0, face_id_1, face_id_2);

	Vector3d1().swap(vecs);
	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);
}

extern "C" PPGL_EXPORT void CGAL_Rotation_Obj(const std::string & path, const double& angle, const Vector3d & axis)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	for (int i = 0; i < vecs.size(); i++)
	{
		Vector3d v = Functs::RotationAxis(vecs[i], angle, axis);
		vecs[i] = v;
	}
	CGAL_Output_Obj_C2(path, vecs, face_id_0, face_id_1, face_id_2);
}

extern "C" PPGL_EXPORT void CGAL_Slicer_Mesh(const std::string & path, const Vector3d & plane_normal, const std::vector<double> & plane_d, Vector3d3 & offsetses, Vector3d2 & offsets)
{
	std::ifstream input(path.c_str());
	Mesh mesh;
	if (!input || !(input >> mesh) || mesh.is_empty()) {
		std::cerr << "Not a valid off file." << std::endl;
		return;
	}
	// Slicer constructor from the mesh
	CGAL::Polygon_mesh_slicer<Mesh, K> slicer(mesh);
	Polylines polylines;

	// Use the Slicer constructor from a pre-built AABB_treen
	AABB_tree tree(edges(mesh).first, edges(mesh).second, mesh);

	CGAL::Polygon_mesh_slicer<Mesh, K> slicer_aabb(mesh, tree);

	for (int i = 0; i < plane_d.size(); i++)
	{
		std::cout << i << "/" << plane_d.size() << std::endl;

		slicer_aabb(K::Plane_3(plane_normal[0], plane_normal[1], plane_normal[2], -plane_d[i]), std::back_inserter(polylines));

		std::vector<std::vector<double>> xs;
		std::vector<std::vector<double>> ys;
		std::vector<std::vector<double>> zs;

		Vector3d2 circles;

		Polylines::iterator iter;
		for (iter = polylines.begin(); iter != polylines.end(); iter++)
		{
			Vector3d1 one_circle;

			Polyline_type::iterator p_iter;
			for (p_iter = iter->begin(); p_iter != iter->end(); p_iter++)
			{
				one_circle.push_back(Vector3d(p_iter->x(), p_iter->y(), p_iter->z()));
			}
			circles.push_back(one_circle);
			offsets.push_back(one_circle);
		}
		polylines.clear();
		offsetses.push_back(circles);
	}
}

//
//
///***************************************************************************************************/
//
////shortest geodesic path
///***************************************************************************************************/
//extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C1(const std::string & path, Vector3d1 & xyzs)
//{
//	// read input polyhedron
//	Polyhedron_3 polyhedron;
//	Construct_Polyhedron(polyhedron, path);
//	// pick up a random face
//	const size_t randSeed = 7915421;
//	CGAL::Random rand(randSeed);
//	const int target_face_index = rand.get_int(0, num_faces(polyhedron));
//	face_iterator face_it = faces(polyhedron).first;
//	std::advance(face_it, target_face_index);
//	// ... and define a barycentric coordinate inside the face
//	Traits::Barycentric_coordinate face_location = { { 0.25, 0.5, 0.25 } };
//	// construct a shortest path query object and add a source point
//	Surface_mesh_shortest_path shortest_paths(polyhedron);
//	shortest_paths.add_source_point(*face_it, face_location);
//
//	vertex_iterator vit = polyhedron.vertices_begin();
//	std::vector<Traits::Point_3> points;
//	shortest_paths.shortest_path_points_to_source_points(*vit, std::back_inserter(points));
//
//	for (int i = 0; i < points.size(); i++)
//	{
//		xyzs.push_back(Vector3d(points[i][0], points[i][1], points[i][2]));
//	}
//
//}
//
//

//

//
//
//void RelatedFaceAndBarycentric(const Polyhedron_3& polyhedron, const Tree& tree,
//	const Vector3d& source, double& u, double& v, double& w, Poly_point_3& nearest_point, face_iterator& face_it)
//{
//	Poly_point_3 query(source[0], source[1], source[2]);
//	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
//	nearest_point = pp.first;
//	face_it = pp.second;
//
//	Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
//	Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
//	Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();
//
//	Barycentric(pp.first, p0, p1, p2, u, v, w);
//}
//
//
//extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C2(Polyhedron_3 & polyhedron, const Tree & tree,
//	Vector3d & source, Vector3d & target, Vector3d1 & xyzs)
//{
//	////////////////////////////////////////
//	face_iterator source_face, target_face;
//	double source_x_w, source_y_w, source_z_w;
//	double target_x_w, target_y_w, target_z_w;
//	Poly_point_3 source_nearest_point, target_nearest_point;
//
//	RelatedFaceAndBarycentric(polyhedron, tree, source, source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
//	RelatedFaceAndBarycentric(polyhedron, tree, target, target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);
//
//	Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
//	Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
//	//////////////////////////////////////////////////////////////
//
//	Surface_mesh_shortest_path shortest_paths(polyhedron);
//	shortest_paths.add_source_point(*source_face, source_face_location);
//
//	std::vector<Traits::Point_3> points;
//	shortest_paths.shortest_path_points_to_source_points(*target_face, target_face_location, std::back_inserter(points));
//
//	for (int i = points.size() - 1; i >= 0; i--)
//	{
//		xyzs.push_back(Vector3d(points[i][0], points[i][1], points[i][2]));
//	}
//}
//
//extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C3(std::string path, Vector3d source, Vector3d target, Vector3d1 & output)
//{
//	Polyhedron_3 polyhedron;
//	Construct_Polyhedron(polyhedron, path);
//
//	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
//	tree.accelerate_distance_queries();
//
//	//////////////////////////////////////////////////////////////
//	face_iterator source_face, target_face;
//	double source_x_w, source_y_w, source_z_w;
//	double target_x_w, target_y_w, target_z_w;
//	Poly_point_3 source_nearest_point, target_nearest_point;
//
//	RelatedFaceAndBarycentric(polyhedron, tree, source, source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
//	RelatedFaceAndBarycentric(polyhedron, tree, target, target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);
//
//	Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
//	Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
//	//////////////////////////////////////////////////////////////
//
//	Surface_mesh_shortest_path shortest_paths(polyhedron);
//	shortest_paths.add_source_point(*source_face, source_face_location);
//
//	std::vector<Traits::Point_3> points;
//	shortest_paths.shortest_path_points_to_source_points(*target_face, target_face_location, std::back_inserter(points));
//
//	for (int i = points.size() - 1; i >= 0; i--)
//	{
//		output.push_back(Vector3d(points[i][0], points[i][1], points[i][2]));
//	}
//}
//
//extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C4(std::string path, Vector3d1 sources,
//	Vector3d1 targets, Vector3d2 & xyzes)
//{
//	Polyhedron_3 polyhedron;
//
//	Construct_Polyhedron(polyhedron, path);
//
//	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
//	tree.accelerate_distance_queries();
//
//	std::cout << "Start to compute the geodesic path..." << std::endl;
//
//	for (int i = 0; i < sources.size(); i++)
//	{
//		std::cout << "Path: " << i << std::endl;
//		Vector3d1 xyzs;
//		CGAL_Shortest_Geodesic_Path_C2(polyhedron, tree, sources[i], targets[i], xyzs);
//		xyzes.push_back(xyzs);
//	}
//}
//
//
//extern "C" PPGL_EXPORT double CGAL_Geodesic_Distance(const std::string & path, const Vector3d & source, const Vector3d & target)
//{
//	std::cout << "one time geodesic computing.." << std::endl;
//
//	Polyhedron_3 polyhedron;
//	Construct_Polyhedron(polyhedron, path);
//
//	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
//	tree.accelerate_distance_queries();
//
//	//////////////////////////////////////////////////////////////
//	face_iterator source_face, target_face;
//	double source_x_w, source_y_w, source_z_w;
//	double target_x_w, target_y_w, target_z_w;
//	Poly_point_3 source_nearest_point, target_nearest_point;
//
//	RelatedFaceAndBarycentric(polyhedron, tree, source, source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
//	RelatedFaceAndBarycentric(polyhedron, tree, target, target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);
//
//	Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
//	Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
//	//////////////////////////////////////////////////////////////
//
//	Surface_mesh_shortest_path shortest_paths(polyhedron);
//	shortest_paths.add_source_point(*source_face, source_face_location);
//
//	return shortest_paths.shortest_distance_to_source_points(*target_face, target_face_location).first;
//}
//
//Vector3d NearestPoint(const Polyhedron_3& polyhedron, const Tree& tree, const Vector3d& source)
//{
//	Poly_point_3 query(source[0], source[1], source[2]);
//	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
//	return Vector3d(pp.first.x(), pp.first.y(), pp.first.z());
//}
//
//extern "C" PPGL_EXPORT Vector3d1 CGAL_Project_Points_Onto_Surface_C1(const Vector3d1 & vecs, const std::vector<int> & face_id_0, const std::vector<int> & face_id_1, const std::vector<int> & face_id_2, const Vector3d1 & points)
//{
//	//construct polyhedron
//	Polyhedron_3 polyhedron;
//	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);
//
//	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
//	tree.accelerate_distance_queries();
//
//	Vector3d1 temp;
//	for (int i = 0; i < points.size(); i++)
//	{
//		temp.push_back(NearestPoint(polyhedron, tree, points[i]));
//	}
//	return temp;
//}
//
//extern "C" PPGL_EXPORT Vector3d1 CGAL_Project_Points_Onto_Surface_C2(const std::string & path, const Vector3d1 & points)
//{
//	Polyhedron_3 polyhedron;
//	Construct_Polyhedron(polyhedron, path);
//
//	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
//	tree.accelerate_distance_queries();
//
//	Vector3d1 temp;
//	for (int i = 0; i < points.size(); i++)
//	{
//		temp.push_back(NearestPoint(polyhedron, tree, points[i]));
//	}
//	return temp;
//}
