#include "geom.h"

#include <Mathematics/MeshCurvature.h>
#include <Mathematics/BSplineCurveFit.h>
#include "NewtonApple_hull3D.h"
#include "kdtree.h"

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
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
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



//Project p onto the planar surface of 3d triangle
//Checking the position relationship between the p and 3d triangle
//face: 3d triangle
//p: 3d point
//return true: inside
//return false: outside
bool OutsidePointInsideTriangle(Poly_facet_iterator &face, Vector3d p) 
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


extern "C" PPGL_EXPORT void CGAL_Mesh_Edges(const char* path) {
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

extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Segment_Mesh(const Vector3d& s, const Vector3d& e, const char* path)
{
	std::cerr << "CGAL_3D_Intersection_Segment_Mesh..." << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Segment_3 seg(VectorPoint3d(s),VectorPoint3d(e));

	if (tree.do_intersect(seg))
		return true;
	else
		return false;
}

extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Segments_Mesh(const Vector3d1& ss, const Vector3d1& ee, const char* path, Vector1b1& inters)
{
	Functs::MAssert(ss.size()==ee.size(),"ss.size()!=ee.size()");

	std::cerr << "CGAL_3D_Intersection_Segment_Mesh..." << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < ss.size(); i++)
	{
		Segment_3 seg(VectorPoint3d(ss[i]), VectorPoint3d(ee[i]));
		inters.push_back(tree.do_intersect(seg));
	}

}

extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Polygons_Mesh(const Vector3d2& polygons, const char* path, Vector1b1& inters)
{
	for (auto& polygon : polygons)
		Functs::MAssert(polygon.size()==4,"Input polygons are not valid: polygon.size()!=4");

	std::cerr << "CGAL_3D_Intersection_Polygons_Mesh..." << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for(int i=0;i<polygons.size();i++)
	{
		Functs::OutputIterInfo("Polygons:", polygons.size(), i, 10);
		auto& polygon = polygons[i];

		int num_inters = inters.size();
		for (int j = 0; j < polygon.size(); j++)
		{
			auto s = polygon[j];
			auto e = polygon[(j+1)%polygon.size()];
			Segment_3 seg(VectorPoint3d(s), VectorPoint3d(e));
			if (tree.do_intersect(seg))
			{
				inters.push_back(true);
				break;
			}
		}

		if (inters.size() == num_inters)
		{
			Triangle triangle_0(VectorPoint3d(polygon[0]), VectorPoint3d(polygon[1]), VectorPoint3d(polygon[2]));
			Triangle triangle_1(VectorPoint3d(polygon[2]), VectorPoint3d(polygon[3]), VectorPoint3d(polygon[0]));
			int num_0 = tree.number_of_intersected_primitives(triangle_0);
			int num_1 = tree.number_of_intersected_primitives(triangle_1);
			inters.push_back(num_0 != 0 || num_1 != 0);
		}
	}
}

//check whether there is a polygon intersected with the input mesh
extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Polygons_Mesh_Bool(const Vector3d2& polygons, const char* path)
{
	for (auto& polygon : polygons)
		Functs::MAssert(polygon.size() == 4, "Input polygons are not valid: polygon.size()!=4");

	std::cerr << "CGAL_3D_Intersection_Polygons_Mesh..." << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < polygons.size(); i++)
	{
		Functs::OutputIterInfo("Polygons:", polygons.size(), i, 10);
		auto& polygon = polygons[i];

		for (int j = 0; j < polygon.size(); j++)
		{
			auto s = polygon[j];
			auto e = polygon[(j + 1) % polygon.size()];
			Segment_3 seg(VectorPoint3d(s), VectorPoint3d(e));
			if (tree.do_intersect(seg))
			{
				return true;
			}
		}

		Triangle triangle_0(VectorPoint3d(polygon[0]), VectorPoint3d(polygon[1]), VectorPoint3d(polygon[2]));
		Triangle triangle_1(VectorPoint3d(polygon[2]), VectorPoint3d(polygon[3]), VectorPoint3d(polygon[0]));
		if (tree.number_of_intersected_primitives(triangle_0) != 0) return true;
		if (tree.number_of_intersected_primitives(triangle_1) != 0) return true;
	}
	return false;
}

extern "C" PPGL_EXPORT bool CGAL_3D_Intersection_Ray_Mesh(const Vector3d& p, const Vector3d & n, const char* path)
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
extern "C" PPGL_EXPORT void  CGAL_3D_Intersection_Rays_Mesh_C1_Bool(const Vector3d1& ps, const Vector3d2& nes, const char* path, Vector1b2& inters)
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
extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_C2_Bool(const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector1b2& inters)
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
			//Segment_3 ray(p3, Point_3(ns[j][0], ns[j][1], ns[j][2]));
			inters[i].push_back(tree.do_intersect(ray));

			//std::list<Ray_intersection> intersections;
			//tree.all_intersections(ray, std::back_inserter(intersections));
		}
	}
}

extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d(
	const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector1d2& inters)
{
	auto IntersectionVector3d = [](const Vector3d& origin, std::list<Ray_intersection>& intersections)
	{
		bool b = false;

		double min_d = MAXDOUBLE;

		Vector3d near_p;
		for (auto iter = intersections.begin(); iter != intersections.end(); iter++)
		{
			const Point_3* p = boost::get<Point_3>(&(iter->value().first));

			if (p)
			{
				Vector3d v(p->x(), p->y(), p->z());
				double d = CGAL_3D_Distance_Point_Point(v, origin);

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
			return -1.0;
		}
		else
		{
			return min_d;
		}
	};

	//input validation
	if (!Functs::DetectExisting(path))
		Functs::MAssert("!Functs::DetectExisting(path)");

	//construct polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	//intersection
	inters = Vector1d2(ps.size(), Vector1d1());
	for (int i = 0; i < ps.size(); i++)
	{
		Functs::OutputIterInfo("CGAL_3D_Intersection_Rays_Mesh_C2_Vector3d", (int)ps.size(), (int)i, (int)10);
		Point_3 p3 = VectorPoint3d(ps[i]);



		for (int j = 0; j < ns.size(); j++)
		{
			Ray_3 ray(p3, Vector_3(ns[j][0], ns[j][1], ns[j][2]));
			//inters[i].push_back(tree.do_intersect(ray));

			auto inter = tree.first_intersection(ray);

			if (inter)
			{
				if (boost::get<Point_3>(&(inter->first))) {
					const Point_3* p = boost::get<Point_3>(&(inter->first));
					//std::cout << *p << std::endl;
					inters[i].push_back(Functs::GetDistance(PointVector3d(*p), ps[i]));
				}
				else
				{
					inters[i].push_back(-1.0);
				}
			}
			else
			{
				inters[i].push_back(-1.0);
			}

			if (tree.do_intersect(ray))
			{
				if (inters[i].back() < 0.0)
				{
					Functs::MAssert("inters[i].back() < 0.0: "+std::to_string(i)+":"+std::to_string(j));
				}
			}

			//std::list<Ray_intersection> intersections;
			//tree.all_intersections(ray, std::back_inserter(intersections));
			//inters[i].push_back(IntersectionVector3d(ps[i], intersections));
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

extern "C" PPGL_EXPORT void CGAL_3D_Points_Inside_Triangles_C2_Bool(const char* path, const Vector3d1& points, std::vector<bool>& insides)
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
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Dart_Sampling_C1(const char* outside_path, const double& d, Vector3d1 & sampling_points, const int& total_iter)
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

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Dart_Sampling_C2(const char* outside_path, const char* inside_path, const double& d, Vector3d1 & sampling_points, const int& total_iter)
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
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C1(const char* outside_path, const double& d, Vector3d1 & sampling_points)
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

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C2(const char* outside_path, const char* inside_path, const double& d, Vector3d1 & sampling_points)
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

//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Cube_Surface_Sampling_C1(const double& cube_size, const double& d, Vector3d2& sampling_points, VectorPI2& neighbors, const bool& compute_neighbors)
{
	auto half_cube = cube_size / 2.0;

	VectorPI1 sneighbors;
	Vector3d1 z3d = Functs::Vector2d3d(CGAL_2D_Square_Regular_Sampling_C3(d, sneighbors, compute_neighbors));
	z3d = Functs::PosApplyM(z3d, Functs::ScaleMatrix(Vector3d(cube_size, cube_size, 0.0)));//-z +z
	auto y3d = Functs::PosApplyM(z3d, Functs::RotationMatrix(Vector3d(1.0, 0.0, 0.0), Math_PI / 2.0));//-y +y
	auto x3d = Functs::PosApplyM(z3d, Functs::RotationMatrix(Vector3d(0.0, 1.0, 0.0), Math_PI / 2.0));//-x +x

	sampling_points.push_back(Functs::PosApplyM(x3d, Functs::TranslationMatrix(Vector3d(half_cube, 0.0, 0.0))));//x
	sampling_points.push_back(Functs::PosApplyM(x3d, Functs::TranslationMatrix(Vector3d(-half_cube, 0.0, 0.0))));//x
	sampling_points.push_back(Functs::PosApplyM(y3d, Functs::TranslationMatrix(Vector3d(0.0, half_cube, 0.0))));//y
	sampling_points.push_back(Functs::PosApplyM(y3d, Functs::TranslationMatrix(Vector3d(0.0, -half_cube, 0.0))));//y
	sampling_points.push_back(Functs::PosApplyM(z3d, Functs::TranslationMatrix(Vector3d(0.0, 0.0, half_cube))));//z
	sampling_points.push_back(Functs::PosApplyM(z3d, Functs::TranslationMatrix(Vector3d(0.0, 0.0, -half_cube))));//z

	if (compute_neighbors)
	{
		neighbors.push_back(sneighbors);
		neighbors.push_back(sneighbors);
		neighbors.push_back(sneighbors);
		neighbors.push_back(sneighbors);
		neighbors.push_back(sneighbors);
		neighbors.push_back(sneighbors);
	}



}
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Cube_Surface_Sampling_C2(const double& cube_size, const double& d, Vector3d2& sampling_points)
{
	VectorPI2 neighbors;
	CGAL_3D_Cube_Surface_Sampling_C1(cube_size, d, sampling_points, neighbors, false);
}
//d: percentage value of the length of the diagonal of the bounding box.
extern "C" PPGL_EXPORT void CGAL_3D_Cube_Surface_Sampling_C3(const double& cube_size, const double& d, Vector3d2& sampling_points, VectorPI2& neighbors)
{
	CGAL_3D_Cube_Surface_Sampling_C1(cube_size, d, sampling_points, neighbors, true);

}

//extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Segments_Mesh(const Vector3d1& ss, const Vector3d1& ee, const char* path, Vector1b1& inters)
//{
//	Functs::MAssert(ss.size() != ee.size(), "ss.size()!=ee.size()");
//
//	std::cerr << "CGAL_3D_Intersection_Segment_Mesh..." << std::endl;
//
//	Polyhedron_3 polyhedron;
//	Construct_Polyhedron(polyhedron, path);
//	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
//	tree.accelerate_distance_queries();
//
//	for (int i = 0; i < ss.size(); i++)
//	{
//		Segment_3 ray(VectorPoint3d(ss[i]), VectorPoint3d(ee[i]));
//		inters.push_back(tree.do_intersect(ray));
//	}
//}

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Regular_Sampling_C3(const char* outside_path, const char* inside_path, const double& d, Vector3d1 & sampling_points,VectorPI1& neighbors)
{
	Functs::MAssert((d > 0 && d < 1.0), "CGAL_3D_Mesh_Dart_Sampling_C1 if (!(d > 0 && d < 1.0))");

	//outside
	Polyhedron_3 out_polyhedron;
	Vector3d1 out_vecs;
	Vector1i1 out_face_id_0, out_face_id_1, out_face_id_2;
	Vector3d out_minC, out_maxC;
	Construct_Polyhedron(out_polyhedron, outside_path, out_vecs, out_face_id_0, out_face_id_1, out_face_id_2);
	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> out_checker(out_polyhedron);
	Functs::GetBoundingBox(out_vecs, out_minC, out_maxC);

	//inside
	Polyhedron_3 in_polyhedron;
	Construct_Polyhedron(in_polyhedron, inside_path);
	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> in_checker(in_polyhedron);

	//tree
	Tree in_tree(faces(in_polyhedron).first, faces(in_polyhedron).second, in_polyhedron);
	in_tree.accelerate_distance_queries();

	double diagonal_length = CGAL_3D_Distance_Point_Point(out_minC, out_maxC);
	double minimal_d = d * diagonal_length;

	double x(out_minC[0]);
	Vector1i3 xyzb;
	while (x < out_maxC[0])
	{
		double y(out_minC[1]);
		Vector1i2 yzb;
		while (y < out_maxC[1])
		{
			double z(out_minC[2]);
			Vector1i1 zb;
			while (z < out_maxC[2])
			{
				CGAL::Bounded_side out_res = out_checker(Point_3(x, y, z));
				CGAL::Bounded_side in_res = in_checker(Point_3(x, y, z));
				if (out_res == CGAL::ON_BOUNDED_SIDE && in_res == CGAL::ON_UNBOUNDED_SIDE)
				{
					zb.push_back(sampling_points.size());
					sampling_points.push_back(Vector3d(x, y, z));
					if (sampling_points.size() % 100 == 0) std::cerr << sampling_points.size() << " ";
				}
				else
				{
					zb.push_back(-1);
				}

				z += minimal_d;
			}
			yzb.push_back(zb);
			y += minimal_d;
		}
		xyzb.push_back(yzb);
		x += minimal_d;
	}

	// get the inside or outside relation of a
	auto GetPB = [](const Vector1i3& xyzb, const Vector3i& a)
	{
		const int xe = xyzb.size();
		const int ye = xyzb.front().size();
		const int ze = xyzb.front().front().size();
		if (!(a[0] >= 0 && a[1] >= 0 && a[2] >= 0 && a[0] < xe && a[1] < ye && a[2] < ze)) return -1;
		return xyzb[a[0]][a[1]][a[2]];
	};

	auto CheckIntersection = [](const Tree& in_tree, const Vector3d1& sampling_points, const int& i, const int& j)
	{
		Segment_3 seg(VectorPoint3d(sampling_points[i]), VectorPoint3d(sampling_points[j]));
		return in_tree.do_intersect(seg);
	};

	auto PushEdge = [&](const Vector3i& xyz3i, const Vector3i& v3i)
	{
		if (GetPB(xyzb, v3i) >= 0)
		{
			if (GetPB(xyzb, xyz3i) >= 0)
			{
				if (!CheckIntersection(in_tree, sampling_points, xyzb[xyz3i[0]][xyz3i[1]][xyz3i[2]], xyzb[v3i[0]][v3i[1]][v3i[2]]))
					neighbors.push_back(std::pair<int, int>(xyzb[xyz3i[0]][xyz3i[1]][xyz3i[2]], xyzb[v3i[0]][v3i[1]][v3i[2]]));
			}
		}
	};

	for (int xi = 0; xi < xyzb.size(); xi++)
	{
		for (int yi = 0; yi < xyzb[xi].size(); yi++)
		{
			for (int zi = 0; zi < xyzb[xi][yi].size(); zi++)
			{
				if (xyzb[xi][yi][zi] >= 0)
				{
					PushEdge(Vector3i(xi, yi, zi), Vector3i(xi + 1, yi, zi));
					PushEdge(Vector3i(xi, yi, zi), Vector3i(xi, yi + 1, zi));
					PushEdge(Vector3i(xi, yi, zi), Vector3i(xi, yi, zi + 1));

					PushEdge(Vector3i(xi, yi, zi), Vector3i(xi + 1, yi + 1, zi));
					PushEdge(Vector3i(xi, yi, zi), Vector3i(xi, yi + 1, zi + 1));
					PushEdge(Vector3i(xi, yi, zi), Vector3i(xi + 1, yi, zi + 1));

					PushEdge(Vector3i(xi + 1, yi, zi), Vector3i(xi, yi + 1, zi));//xy
					PushEdge(Vector3i(xi, yi + 1, zi), Vector3i(xi, yi, zi + 1));//yz
					PushEdge(Vector3i(xi+1, yi, zi), Vector3i(xi, yi, zi + 1));//xz

					PushEdge(Vector3i(xi, yi, zi), Vector3i(xi + 1, yi + 1, zi + 1));
					PushEdge(Vector3i(xi + 1, yi, zi), Vector3i(xi, yi + 1, zi + 1));
					PushEdge(Vector3i(xi, yi, zi + 1), Vector3i(xi + 1, yi + 1, zi));
					PushEdge(Vector3i(xi+1, yi, zi + 1), Vector3i(xi, yi+1, zi));
				}
			}
		}
	}

	std::cerr << std::endl;
}

extern "C" PPGL_EXPORT void CGAL_3D_Intersection_Rays_Mesh_Vector3d(const Vector3d1& ps, const Vector3d1& ns, const char* path, Vector3d1& inters)
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

extern "C" PPGL_EXPORT void CGAL_3D_Distance_Point_Mesh(const char* path, const Vector3d1 & query_points, std::vector<double>&distances)
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

extern "C" PPGL_EXPORT void CGAL_3D_Neareast_Point_Mesh(const char* path, const Vector3d1 & ves, Vector3d1 & ners)
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

		Point_3 query(points[i][0], points[i][1], points[i][2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

		std::priority_queue<Polyhedron_3::Facet_handle> facets;
		std::vector<int> save_index;
		facets.push(pp.second);
		save_index.push_back((int)pp.second->id());

		while (facets.size() != 0)
		{
			Polyhedron_3::Facet_handle fh = facets.top();
			triangle.push_back((int)fh->id());
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
				std::vector<Point_3> n_fh_vecs;
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
					save_index.push_back((int)n_fh->id());
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

extern "C" PPGL_EXPORT void CGAL_3D_Points_inside_Triangles_C2(const char* path, const Vector3d1 & points, std::vector<bool>&insides)
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
				vecs_neighbors_labels[i][j] =(int)edges.size() - 1;

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


extern "C" PPGL_EXPORT void CGAL_Mesh_Subdivision(const char* in_path, const char* sub, const int& step, const char* out_path)
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
		Point_3 p = iter->point();
		vecs.push_back(Vector3d(p[0], p[1], p[2]));
	}

	for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
	{
		face_id_0.push_back((int)iter->halfedge()->next()->next()->vertex()->id());
		face_id_1.push_back((int)iter->halfedge()->vertex()->id());
		face_id_2.push_back((int)iter->halfedge()->next()->vertex()->id());
	}

	Functs::OutputObj3d(out_path, vecs, face_id_0, face_id_1, face_id_2);
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
		for (int j = 0; j < vecs_neigbor_lable[i].size() && !b; j++) {
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
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C2(const char* path, std::vector<bool>&bools)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	Functs::LoadObj3d(path, vecs, face_id_0, face_id_1, face_id_2);
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
extern "C" PPGL_EXPORT void CGAL_3D_Triangle_Mesh_Boundary_C5(const char* path, Vector3d2 & boundaries)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	Functs::LoadObj3d(path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_3D_Triangle_Mesh_Boundary_C3(vecs, face_id_0, face_id_1, face_id_2, boundaries);
}


extern "C" PPGL_EXPORT void CGAL_Mesh_Laplace_Smooth_C1(const char* in_path, const char* out_path, const int laplace_nb)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	Functs::LoadObj3d(in_path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_Mesh_Laplace_Smooth_C2(vecs, face_id_0, face_id_1, face_id_2, laplace_nb);
	Functs::OutputObj3d(out_path, vecs, face_id_0, face_id_1, face_id_2);
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
					Vector3d n = vecs_normals[i];
					Functs::SetVectorLength(n, 0.001 * min_curvature[i] / low_curvature);
					Vector3d cur_v = vecs[i] + n;

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

extern "C" PPGL_EXPORT void CGAL_Mesh_Loop_Subdivision_Own_Version(const char* in_path, const int& step, const char* out_path, const int& laplace_nb)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	Functs::LoadObj3d(in_path, vecs, face_id_0, face_id_1, face_id_2);

	for (int i = 0; i < step; i++)
	{
		CGAL_Mesh_Loop_Subdivision_One_Step(vecs, face_id_0, face_id_1, face_id_2);
		CGAL_Mesh_Laplace_Smooth_C2(vecs, face_id_0, face_id_1, face_id_2, laplace_nb);
	}

	Functs::OutputObj3d(out_path, vecs, face_id_0, face_id_1, face_id_2);

	Vector3d1().swap(vecs);
	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);
}

extern "C" PPGL_EXPORT void CGAL_Rotation_Obj(const char* path, const double& angle, const Vector3d & axis)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	Functs::LoadObj3d(path, vecs, face_id_0, face_id_1, face_id_2);
	for (int i = 0; i < vecs.size(); i++)
	{
		Vector3d v = Functs::RotationAxis(vecs[i], angle, axis);
		vecs[i] = v;
	}
	Functs::OutputObj3d(path, vecs, face_id_0, face_id_1, face_id_2);
}

extern "C" PPGL_EXPORT void CGAL_Slicer_Mesh(const char* path, const Vector3d & plane_normal, const std::vector<double> & plane_d, Vector3d3 & offsetses, Vector3d2 & offsets)
{
	std::ifstream input(path);
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



/***************************************************************************************************/
//shortest geodesic path
/***************************************************************************************************/
extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C1(const char* path, Vector3d1 & xyzs)
{
	// read input polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	// pick up a random face
	const size_t randSeed = 7915421;
	CGAL::Random rand(randSeed);
	const int target_face_index = rand.get_int(0, num_faces(polyhedron));
	face_iterator face_it = faces(polyhedron).first;
	std::advance(face_it, target_face_index);
	// ... and define a barycentric coordinate inside the face
	Traits::Barycentric_coordinate face_location = { { 0.25, 0.5, 0.25 } };
	// construct a shortest path query object and add a source point
	Surface_mesh_shortest_path shortest_paths(polyhedron);
	shortest_paths.add_source_point(*face_it, face_location);

	vertex_iterator vit = polyhedron.vertices_begin();
	std::vector<Traits::Point_3> points;
	shortest_paths.shortest_path_points_to_source_points(*vit, std::back_inserter(points));

	for (int i = 0; i < points.size(); i++)
	{
		xyzs.push_back(Vector3d(points[i][0], points[i][1], points[i][2]));
	}
}


void RelatedFaceAndBarycentric(const Polyhedron_3& polyhedron, const Tree& tree, const Vector3d& source, double& u, double& v, double& w, Point_3& nearest_point, face_iterator& face_it)
{
	Point_3 query(source[0], source[1], source[2]);
	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
	nearest_point = pp.first;
	face_it = pp.second;

	Point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
	Point_3 p1 = pp.second->halfedge()->vertex()->point();
	Point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

	Functs::Barycentric(PointVector3d(pp.first), PointVector3d(p0), PointVector3d(p1), PointVector3d(p2), u,v , w);
}


extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C3(const char* path, Vector3d source, Vector3d target, Vector3d1 & output)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	//////////////////////////////////////////////////////////////
	face_iterator source_face, target_face;
	double source_x_w, source_y_w, source_z_w;
	double target_x_w, target_y_w, target_z_w;
	Point_3 source_nearest_point, target_nearest_point;

	RelatedFaceAndBarycentric(polyhedron, tree, source, source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
	RelatedFaceAndBarycentric(polyhedron, tree, target, target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);

	Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
	Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
	//////////////////////////////////////////////////////////////

	Surface_mesh_shortest_path shortest_paths(polyhedron);
	shortest_paths.add_source_point(*source_face, source_face_location);

	std::vector<Traits::Point_3> points;
	shortest_paths.shortest_path_points_to_source_points(*target_face, target_face_location, std::back_inserter(points));

	for (int i = points.size() - 1; i >= 0; i--)
	{
		output.push_back(Vector3d(points[i][0], points[i][1], points[i][2]));
	}
}

extern "C" PPGL_EXPORT void CGAL_Shortest_Geodesic_Path_C4(const char* path_, Vector3d1 sources,
	Vector3d1 targets, Vector3d2 & xyzes)
{
	auto CGAL_Shortest_Geodesic_Path_C2=[](Polyhedron_3 & polyhedron, const Tree & tree, Vector3d & source, Vector3d & target, Vector3d1 & xyzs)
	{
		////////////////////////////////////////
		face_iterator source_face, target_face;
		double source_x_w, source_y_w, source_z_w;
		double target_x_w, target_y_w, target_z_w;
		Point_3 source_nearest_point, target_nearest_point;

		RelatedFaceAndBarycentric(polyhedron, tree, source, source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
		RelatedFaceAndBarycentric(polyhedron, tree, target, target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);

		Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
		Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
		//////////////////////////////////////////////////////////////

		Surface_mesh_shortest_path shortest_paths(polyhedron);
		shortest_paths.add_source_point(*source_face, source_face_location);

		std::vector<Traits::Point_3> points;
		shortest_paths.shortest_path_points_to_source_points(*target_face, target_face_location, std::back_inserter(points));

		for (int i = points.size() - 1; i >= 0; i--)
		{
			xyzs.push_back(Vector3d(points[i][0], points[i][1], points[i][2]));
		}
	};

	std::string path = path_;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path_);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();
	std::cout << "Start to compute the geodesic path..." << std::endl;

	for (int i = 0; i < sources.size(); i++)
	{
		std::cout << "Path: " << i << std::endl;
		Vector3d1 xyzs;
		CGAL_Shortest_Geodesic_Path_C2(polyhedron, tree, sources[i], targets[i], xyzs);
		xyzes.push_back(xyzs);
	}
}


extern "C" PPGL_EXPORT double CGAL_Geodesic_Distance(const char* path, const Vector3d & source, const Vector3d & target)
{
	std::cout << "one time geodesic computing.." << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	//////////////////////////////////////////////////////////////
	face_iterator source_face, target_face;
	double source_x_w, source_y_w, source_z_w;
	double target_x_w, target_y_w, target_z_w;
	Point_3 source_nearest_point, target_nearest_point;

	RelatedFaceAndBarycentric(polyhedron, tree, source, source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
	RelatedFaceAndBarycentric(polyhedron, tree, target, target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);

	Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
	Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
	//////////////////////////////////////////////////////////////

	Surface_mesh_shortest_path shortest_paths(polyhedron);
	shortest_paths.add_source_point(*source_face, source_face_location);

	return shortest_paths.shortest_distance_to_source_points(*target_face, target_face_location).first;
}

extern "C" PPGL_EXPORT Vector3d1 CGAL_Project_Points_Onto_Surface_C1(const Vector3d1 & vecs, const std::vector<int> & face_id_0, const std::vector<int> & face_id_1, const std::vector<int> & face_id_2, const Vector3d1 & points)
{
	auto NearestPoint=[](const Polyhedron_3 & polyhedron, const Tree & tree, const Vector3d & source)
	{
		Point_3 query(source[0], source[1], source[2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
		return Vector3d(pp.first.x(), pp.first.y(), pp.first.z());
	};

	//construct polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Vector3d1 temp;
	for (int i = 0; i < points.size(); i++)
	{
		temp.push_back(NearestPoint(polyhedron, tree, points[i]));
	}
	return temp;
}

extern "C" PPGL_EXPORT Vector3d1 CGAL_Project_Points_Onto_Surface_C2(const char* path, const Vector3d1 & points)
{
	auto NearestPoint = [](const Polyhedron_3& polyhedron, const Tree& tree, const Vector3d& source)
	{
		Point_3 query(source[0], source[1], source[2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
		return Vector3d(pp.first.x(), pp.first.y(), pp.first.z());
	};

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Vector3d1 temp;
	for (int i = 0; i < points.size(); i++)
	{
		temp.push_back(NearestPoint(polyhedron, tree, points[i]));
	}
	return temp;
}



extern "C" PPGL_EXPORT void CGAL_3D_Triangel_Mesh_Most_Inside_Point(const Vector3d1 & vecs_, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, Vector3d & inside)
{
	//inside
	std::vector<bool> vec_boundary;
	CGAL_3D_Triangle_Mesh_Boundary_C1(vecs_, face_id_0, face_id_1, face_id_2, vec_boundary);

	kdtree* tree = kd_create(3);

	auto vecs = vecs_;
	for (int i = 0; i < vecs.size(); i++)
	{
		if (vec_boundary[i])
		{
			void* val = &vecs[i];
			kd_insert3(tree, vecs[i][0], vecs[i][1], vecs[i][2], val);
		}
	}

	double max_d = 0.0;
	for (int i = 0; i < vecs.size(); i++)
	{
		if (!vec_boundary[i])
		{
			double* pos = new double[3];
			pos[0] = vecs[i][0];
			pos[1] = vecs[i][1];
			pos[2] = vecs[i][2];
			struct kdres* r = kd_nearest(tree, pos);
			double position[3];
			*(int*)kd_res_item(r, position);
			double d = CGAL_3D_Distance_Point_Point(Vector3d(pos[0], pos[1], pos[2]),Vector3d(position[0], position[1], position[2]));
			if (d > max_d)
			{
				max_d = d;
				inside = vecs[i];
			}
		}
	}

	std::vector<bool>().swap(vec_boundary);
}

extern "C" PPGL_EXPORT double CGAL_3D_One_Triangle_Area(const Vector3d & v0, const Vector3d & v1, const Vector3d & v2)
{
	return Functs::GetTriangleArea(v0, v1, v2);
}

extern "C" PPGL_EXPORT double CGAL_3D_Triangle_Mesh_Area(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2)
{
	double area = 0.0;

	for (int i = 0; i < face_id_0.size(); i++) {
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		Vector3d v0 = vecs[index_0];
		Vector3d v1 = vecs[index_1];
		Vector3d v2 = vecs[index_2];
		area += CGAL_3D_One_Triangle_Area(v0, v1, v2);
	}

	return area;
}



extern "C" PPGL_EXPORT void CGAL_3D_Convex_Hulls_C1(const Vector3d1 & vec, Vector3d1 & hull_points)
{
	std::vector<Point_3> points;
	for (int i = 0; i < vec.size(); i++)
		points.push_back(VectorPoint3d(vec[i]));

	if (points.size() <= 3)
	{
		hull_points = vec;
		return;
	}
	// define polyhedron to hold convex hull
	Polyhedron_3 poly;
	// compute convex hull of non-collinear points
	CGAL::convex_hull_3(points.begin(), points.end(), poly);
	std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;
	for (Polyhedron_3::Vertex_iterator iter = poly.vertices_begin(); iter != poly.vertices_end(); iter++)
	{
		Point_3 p = iter->point();
		hull_points.push_back(PointVector3d(p));
	}
}

extern "C" PPGL_EXPORT void CGAL_3D_Convex_Hulls_C2(const Vector3d1 & vec, Vector3d1 & hull_points, std::vector<int>&hulls_surface_0, std::vector<int>&hulls_surface_1, std::vector<int>&hulls_surface_2)
{
	std::vector<R3> pts, pts2;
	R3 pt;
	pts.clear();

	for (int i = 0; i < vec.size(); i++)
	{
		//points.push_back(Point_3(xs[i], ys[i], zs[i]));
		pt.id = i;
		pt.r = (float)vec[i][0];
		pt.c = (float)vec[i][1];
		pt.z = (float)vec[i][2];
		pts.push_back(pt);
	}

	std::vector<Tri> tris;

	std::vector<int> outx;
	int nx = de_duplicateR3(pts, outx, pts2);
	pts.clear();
	int ts = NewtonApple_hull_3D(pts2, tris);

	for (int i = 0; i < (int)tris.size(); i++)
	{
		pts.push_back(pts2[tris[i].a]);
		pts.push_back(pts2[tris[i].b]);
		pts.push_back(pts2[tris[i].c]);
	}
	pts2.clear();
	outx.clear();
	tris.clear();
	nx = de_duplicateR3(pts, outx, pts2);
	ts = NewtonApple_hull_3D(pts2, tris);


	for (int i = 0; i < pts2.size(); i++)
	{
		hull_points.push_back(Vector3d(pts2[i].r, pts2[i].c, pts2[i].z));
	}

	for (int i = 0; i < (int)tris.size(); i++)
	{
		hulls_surface_0.push_back(tris[i].a);
		hulls_surface_1.push_back(tris[i].b);
		hulls_surface_2.push_back(tris[i].c);
	}
	pts2.clear();

}
extern "C" PPGL_EXPORT void CGAL_3D_Convex_Hulls_C3(const Vector3d1 & vec, Vector3d1 & hull_points, Vector3d1 & plane_p, Vector3d1 & plane_n)
{
	std::vector<R3> pts, pts2;
	R3 pt;
	pts.clear();

	for (int i = 0; i < vec.size(); i++)
	{
		//points.push_back(Point_3(xs[i], ys[i], zs[i]));
		pt.id = i;
		pt.r = (float)vec[i][0];
		pt.c = (float)vec[i][1];
		pt.z = (float)vec[i][2];
		pts.push_back(pt);
	}

	std::vector<Tri> tris;

	std::vector<int> outx;
	int nx = de_duplicateR3(pts, outx, pts2);
	pts.clear();

	//int ts = NewtonApple_Delaunay( pts2, tris);
	int ts = NewtonApple_hull_3D(pts2, tris);

	int nr = (int)tris.size();

	for (int i = 0; i < nr; i++)
	{
		pts.push_back(pts2[tris[i].a]);
		pts.push_back(pts2[tris[i].b]);
		pts.push_back(pts2[tris[i].c]);

		Point_3 p_0 = Point_3(pts2[tris[i].a].r, pts2[tris[i].a].c, pts2[tris[i].a].z);
		Point_3 p_1 = Point_3(pts2[tris[i].b].r, pts2[tris[i].b].c, pts2[tris[i].b].z);
		Point_3 p_2 = Point_3(pts2[tris[i].c].r, pts2[tris[i].c].c, pts2[tris[i].c].z);

		Plane_3 plane(p_0, p_1, p_2);

		Point_3 p((p_0[0] + p_1[0] + p_2[0]) / 3.0, (p_0[1] + p_1[1] + p_2[1]) / 3.0, (p_0[2] + p_1[2] + p_2[2]) / 3.0);
		plane_p.push_back(PointVector3d(p));

		Vector_3 v = plane.orthogonal_vector();
		plane_n.push_back(Vector3d(v[0], v[1], v[2]));
	}
	pts2.clear();
	nx = de_duplicateR3(pts, outx, pts2);
}

extern "C" PPGL_EXPORT void CGAL_3D_Convex_Hulls_C4(const Vector3d1 & vec, Vector3d1 & hull_points, std::vector<int>&hulls_surface_0, std::vector<int>&hulls_surface_1, std::vector<int>&hulls_surface_2, Vector3d1 & plane_p, Vector3d1 & plane_n)
{
	std::vector<R3> pts, pts2;
	R3 pt;
	pts.clear();

	for (int i = 0; i < vec.size(); i++)
	{
		//points.push_back(Point_3(xs[i], ys[i], zs[i]));
		pt.id = i;
		pt.r = (float)vec[i][0];
		pt.c = (float)vec[i][1];
		pt.z = (float)vec[i][2];
		pts.push_back(pt);
	}

	std::vector<Tri> tris;

	std::vector<int> outx;
	int nx = de_duplicateR3(pts, outx, pts2);
	pts.clear();
	int ts = NewtonApple_hull_3D(pts2, tris);

	for (int i = 0; i < (int)tris.size(); i++)
	{
		pts.push_back(pts2[tris[i].a]);
		pts.push_back(pts2[tris[i].b]);
		pts.push_back(pts2[tris[i].c]);
	}
	pts2.clear();
	outx.clear();
	tris.clear();
	nx = de_duplicateR3(pts, outx, pts2);
	ts = NewtonApple_hull_3D(pts2, tris);


	for (int i = 0; i < pts2.size(); i++)
	{
		hull_points.push_back(Vector3d(pts2[i].r, pts2[i].c, pts2[i].z));
	}

	for (int i = 0; i < (int)tris.size(); i++)
	{
		hulls_surface_0.push_back(tris[i].a);
		hulls_surface_1.push_back(tris[i].b);
		hulls_surface_2.push_back(tris[i].c);

		Point_3 p_0 = Point_3(pts2[tris[i].a].r, pts2[tris[i].a].c, pts2[tris[i].a].z);
		Point_3 p_1 = Point_3(pts2[tris[i].b].r, pts2[tris[i].b].c, pts2[tris[i].b].z);
		Point_3 p_2 = Point_3(pts2[tris[i].c].r, pts2[tris[i].c].c, pts2[tris[i].c].z);

		Plane_3 plane(p_0, p_1, p_2);

		Point_3 p((p_0[0] + p_1[0] + p_2[0]) / 3.0, (p_0[1] + p_1[1] + p_2[1]) / 3.0, (p_0[2] + p_1[2] + p_2[2]) / 3.0);
		//Point_3 p = plane.point();

		plane_p.push_back(PointVector3d(p));

		Vector_3 v = plane.orthogonal_vector();
		plane_n.push_back(Vector3d(v[0], v[1], v[2]));
	}
	pts2.clear();
}



extern "C" PPGL_EXPORT void CGAL_Mesh_Field_Query_C1(const char* path, const Vector3d1 & gradients, const Vector3d1 & input_points, Vector3d1 & points_gradients)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < input_points.size(); i++)
	{
		Point_3 query(input_points[i][0], input_points[i][1], input_points[i][2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

		Point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
		Point_3 p1 = pp.second->halfedge()->vertex()->point();
		Point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

		int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
		int point_id_1 = pp.second->halfedge()->vertex()->id();
		int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

		double u, v, w;
		Functs::Barycentric(PointVector3d(query), PointVector3d(p0), PointVector3d(p1), PointVector3d(p2), u, v, w);

		points_gradients.push_back((double)u * gradients[point_id_0] + (double)v * gradients[point_id_1] + (double)w * gradients[point_id_2]);
	}
}
extern "C" PPGL_EXPORT void CGAL_Mesh_Field_Query_C2(const char* path, const std::vector<double>&gradient_values, const Vector3d1 & input_points, std::vector<double>&points_gradient_values)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < input_points.size(); i++)
	{
		Point_3 query(input_points[i][0], input_points[i][1], input_points[i][2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

		Point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
		Point_3 p1 = pp.second->halfedge()->vertex()->point();
		Point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

		int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
		int point_id_1 = pp.second->halfedge()->vertex()->id();
		int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

		double u, v, w;
		Functs::Barycentric(PointVector3d(query), PointVector3d(p0), PointVector3d(p1), PointVector3d(p2), u, v, w);
		points_gradient_values.push_back((double)u * gradient_values[point_id_0] + (double)v * gradient_values[point_id_1] + (double)w * gradient_values[point_id_2]);
	}
}
extern "C" PPGL_EXPORT void CGAL_Mesh_Field_Query_C3(const char* path, const std::vector<double>&gradient_values, const Vector3d2 & input_point_es, std::vector<std::vector<double>>&points_gradient_value_es)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < input_point_es.size(); i++)
	{
		std::vector<double> values;
		for (int j = 0; j < input_point_es[i].size(); j++)
		{
			Point_3 query(input_point_es[i][j][0], input_point_es[i][j][1], input_point_es[i][j][2]);
			Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

			Point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
			Point_3 p1 = pp.second->halfedge()->vertex()->point();
			Point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

			int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
			int point_id_1 = pp.second->halfedge()->vertex()->id();
			int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

			double u, v, w;
			Functs::Barycentric(PointVector3d(query), PointVector3d(p0), PointVector3d(p1), PointVector3d(p2), u, v, w);

			values.push_back((double)u * gradient_values[point_id_0] + (double)v * gradient_values[point_id_1] + (double)w * gradient_values[point_id_2]);
		}
		points_gradient_value_es.push_back(values);
	}
}

extern "C" PPGL_EXPORT void CGAL_Curvature_Mesh(const char* path_, const Vector3d1 & input_points, std::vector<double>&max_curs, std::vector<double>&min_curs, Vector3d1 & max_curs_directions, Vector3d1 & min_curs_directions)
{
	auto Curvature_Direction_Adjusting=[](Vector3d & cur_0, Vector3d & cur_1, Vector3d & cur_2)
	{
		int f_i = 0;
		int f_j = 0;
		int f_k = 0;

		double angle = 1000000000.0;
		for (int i = -1; i <= 1; i = i + 2)
		{
			for (int j = -1; j <= 1; j = j + 2)
			{
				for (int k = -1; k <= 1; k = k + 2)
				{
					double angle_0 = Functs::GetAngleBetween(cur_0 * (double)i, cur_1 * (double)j);
					double angle_1 = Functs::GetAngleBetween(cur_0 * (double)i, cur_2 * (double)k);
					double angle_2 = Functs::GetAngleBetween(cur_2 * (double)k, cur_1 * (double)j);

					if (angle > angle_0 + angle_1 + angle_2)
					{
						angle = angle_0 + angle_1 + angle_2;
						f_i = i;
						f_j = j;
						f_k = k;
					}
				}
			}
		}

		cur_0 = cur_0 * (double)f_i;
		cur_1 = cur_1 * (double)f_j;
		cur_2 = cur_2 * (double)f_k;
	};

	std::cout << "CGAL_Curvature_Mesh.." << std::endl;

	std::string path = path_;
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		Polyhedron_3 polyhedron;
		Construct_Polyhedron(polyhedron, path_);
		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//get mesh vertices and surface
		//std::vector<double> mesh_xs, mesh_ys, mesh_zs;
		Vector3d1 vecs;
		std::vector<int> mesh_face_id_0, mesh_face_id_1, mesh_face_id_2;

		//Functs::LoadObj3d(path, vecs, mesh_face_id_0, mesh_face_id_1, mesh_face_id_2);

		//Polyhedron_3 polyhedron;
		//Construct_Polyhedron(polyhedron, vecs, mesh_face_id_0, mesh_face_id_1, mesh_face_id_2);
		//Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		//tree.accelerate_distance_queries();

		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
			iter != polyhedron.vertices_end(); iter++)
		{
			Point_3 p = iter->point();
			vecs.push_back(Vector3d(p[0], p[1], p[2]));
		}
		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			mesh_face_id_0.push_back((int)iter->halfedge()->next()->next()->vertex()->id());
			mesh_face_id_1.push_back((int)iter->halfedge()->vertex()->id());
			mesh_face_id_2.push_back((int)iter->halfedge()->next()->vertex()->id());
		}

		//compute surface normals
		int verticeSize = vecs.size();
		int faceindiceSize = mesh_face_id_0.size() * 3;

		gte::Vector3<double>* points = new gte::Vector3<double>[verticeSize];
		unsigned int* indices = new unsigned int[faceindiceSize];

		for (int i = 0; i < verticeSize; i++)
		{
			points[i][0] = vecs[i][0];
			points[i][1] = vecs[i][1];
			points[i][2] = vecs[i][2];
		}

		for (int i = 0; i < mesh_face_id_0.size(); i++)
		{
			indices[3 * i] = mesh_face_id_0[i];
			indices[3 * i + 1] = mesh_face_id_1[i];
			indices[3 * i + 2] = mesh_face_id_2[i];
		}
		gte::MeshCurvature<double> meshCurvature;

		meshCurvature((size_t)verticeSize, points, (size_t)(faceindiceSize / 3), indices, 1e-06f);

		//meshCurvature((size_t)verticeSize, points, (size_t)(faceindiceSize / 3), indices, 1e-06f);

		std::vector<double> mesh_max_curs;
		std::vector<double> mesh_min_curs;

		Vector3d1 mesh_max_curs_directions;
		Vector3d1 mesh_min_curs_directions;

		for (int i = 0; i < verticeSize; i++)
		{
			mesh_max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
			mesh_min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);

			gte::Vector3<double> max_curs_direction = meshCurvature.GetMaxDirections()[i];
			gte::Vector3<double> min_curs_direction = meshCurvature.GetMinDirections()[i];

			mesh_max_curs_directions.push_back(Vector3d(max_curs_direction[0], max_curs_direction[1], max_curs_direction[2]));
			mesh_min_curs_directions.push_back(Vector3d(min_curs_direction[0], min_curs_direction[1], min_curs_direction[2]));
		}

		for (int i = 0; i < input_points.size(); i++)
		{
			Point_3 query(input_points[i][0], input_points[i][1], input_points[i][2]);
			Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

			Point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
			Point_3 p1 = pp.second->halfedge()->vertex()->point();
			Point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

			int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
			int point_id_1 = pp.second->halfedge()->vertex()->id();
			int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

			double u, v, w;
			Functs::Barycentric(PointVector3d(query), PointVector3d(p0), PointVector3d(p1), PointVector3d(p2), u, v, w);

			max_curs.push_back(u * mesh_max_curs[point_id_0] + v * mesh_max_curs[point_id_1] + w * mesh_max_curs[point_id_2]);
			min_curs.push_back(u * mesh_min_curs[point_id_0] + v * mesh_min_curs[point_id_1] + w * mesh_min_curs[point_id_2]);

			Curvature_Direction_Adjusting(mesh_max_curs_directions[point_id_0], mesh_max_curs_directions[point_id_1], mesh_max_curs_directions[point_id_2]);
			Curvature_Direction_Adjusting(mesh_min_curs_directions[point_id_0], mesh_min_curs_directions[point_id_1], mesh_min_curs_directions[point_id_2]);

			Vector3d max_curs_direction = (double)u * mesh_max_curs_directions[point_id_0] + (double)v * mesh_max_curs_directions[point_id_1] + (double)w * mesh_max_curs_directions[point_id_2];
			Vector3d min_curs_direction = (double)u * mesh_min_curs_directions[point_id_0] + (double)v * mesh_min_curs_directions[point_id_1] + (double)w * mesh_min_curs_directions[point_id_2];

			max_curs_directions.push_back(max_curs_direction);
			min_curs_directions.push_back(min_curs_direction);
		}
	}
}



extern "C" PPGL_EXPORT void CGAL_Normal_Mesh_C1(const char* path_, const Vector3d1 & mesh_points, Vector3d1 & mesh_normals)
{
	std::string path = path_;
	std::cout << "CGAL_Normal_Mesh.." << std::endl;
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		Polyhedron_3 polyhedron;
		Construct_Polyhedron(polyhedron, path_);

		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//get mesh vertices and surface
		std::vector<double> mesh_xs, mesh_ys, mesh_zs;
		std::vector<int> mesh_face_id_0, mesh_face_id_1, mesh_face_id_2;
		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
			iter != polyhedron.vertices_end(); iter++)
		{
			Point_3 p = iter->point();
			mesh_xs.push_back(p[0]);
			mesh_ys.push_back(p[1]);
			mesh_zs.push_back(p[2]);
		}
		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			mesh_face_id_0.push_back((int)iter->halfedge()->next()->next()->vertex()->id());
			mesh_face_id_1.push_back((int)iter->halfedge()->vertex()->id());
			mesh_face_id_2.push_back((int)iter->halfedge()->next()->vertex()->id());
		}

		//compute surface normals
		int verticeSize = mesh_xs.size();
		int faceindiceSize = mesh_face_id_0.size() * 3;

		gte::Vector3<double>* points = new gte::Vector3<double>[verticeSize];
		unsigned int* indices = new unsigned int[faceindiceSize];

		for (int i = 0; i < verticeSize; i++)
		{
			points[i][0] = mesh_xs[i];
			points[i][1] = mesh_ys[i];
			points[i][2] = mesh_zs[i];
		}

		for (int i = 0; i < mesh_face_id_0.size(); i++)
		{
			indices[3 * i] = mesh_face_id_0[i];
			indices[3 * i + 1] = mesh_face_id_1[i];
			indices[3 * i + 2] = mesh_face_id_2[i];
		}
		gte::MeshCurvature<double> meshCurvature;
		meshCurvature((size_t)verticeSize, points, (size_t)(faceindiceSize / 3), indices, 1e-06f);

		Vector3d1 normals;
		for (int i = 0; i < verticeSize; i++)
		{
			gte::Vector3<double> normal = meshCurvature.GetNormals()[i];
			normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
		}

		for (int i = 0; i < mesh_points.size(); i++)
		{
			Vector3d n = RelatedFaceNormal(polyhedron, tree, normals, mesh_points[i]);
			mesh_normals.push_back(n);
		}
	}

	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		Vector3d1 vecs;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		Functs::LoadObj3d(path_, vecs, face_id_0, face_id_1, face_id_2);

		Polyhedron_3 polyhedron;

		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//compute normals
		/**********************************************************/
		Vector3d1 normals(vecs.size(), Vector3d(0.0, 0.0, 0.0));
		std::vector<double> areas(vecs.size(), 0.0);
		for (int i = 0; i < face_id_0.size(); i++)
		{
			double area = CGAL_3D_One_Triangle_Area(vecs[face_id_0[i]], vecs[face_id_1[i]], vecs[face_id_2[i]]);
			Point_3 p0(vecs[face_id_0[i]][0], vecs[face_id_0[i]][1], vecs[face_id_0[i]][2]);
			Point_3 p1(vecs[face_id_1[i]][0], vecs[face_id_1[i]][1], vecs[face_id_1[i]][2]);
			Point_3 p2(vecs[face_id_2[i]][0], vecs[face_id_2[i]][1], vecs[face_id_2[i]][2]);
			Vector_3   n = CGAL::cross_product(p2 - p1, p0 - p1);

			Vector3d n0(n[0], n[1], n[2]);

			normals[face_id_0[i]] += (double)area * n0;
			normals[face_id_1[i]] += (double)area * n0;
			normals[face_id_2[i]] += (double)area * n0;
			areas[face_id_0[i]] += area;
			areas[face_id_1[i]] += area;
			areas[face_id_2[i]] += area;
		}

		for (int i = 0; i < vecs.size(); i++)
		{
			normals[i] = normals[i] / (double)areas[i];
		}
		/**********************************************************/

		for (int i = 0; i < mesh_points.size(); i++)
		{
			Vector3d n = RelatedFaceNormal(polyhedron, tree, normals, mesh_points[i]);
			mesh_normals.push_back(n);
		}
	}
}
extern "C" PPGL_EXPORT void CGAL_Normal_Mesh_C2(const char* path_, const Vector3d2 & mesh_pointses, Vector3d2 & mesh_normalses)
{
	std::string path = path_;
	std::cout << "CGAL_Normal_Mesh.." << std::endl;
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		Polyhedron_3 polyhedron;
		Construct_Polyhedron(polyhedron, path_);

		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//get mesh vertices and surface
		std::vector<double> mesh_xs, mesh_ys, mesh_zs;
		std::vector<int> mesh_face_id_0, mesh_face_id_1, mesh_face_id_2;
		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
			iter != polyhedron.vertices_end(); iter++)
		{
			Point_3 p = iter->point();
			mesh_xs.push_back(p[0]);
			mesh_ys.push_back(p[1]);
			mesh_zs.push_back(p[2]);
		}
		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			mesh_face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
			mesh_face_id_1.push_back(iter->halfedge()->vertex()->id());
			mesh_face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
		}

		//compute surface normals
		int verticeSize = mesh_xs.size();
		int faceindiceSize = mesh_face_id_0.size() * 3;

		gte::Vector3<double>* points = new gte::Vector3<double>[verticeSize];
		unsigned int* indices = new unsigned int[faceindiceSize];

		for (int i = 0; i < verticeSize; i++)
		{
			points[i][0] = mesh_xs[i];
			points[i][1] = mesh_ys[i];
			points[i][2] = mesh_zs[i];
		}

		for (int i = 0; i < mesh_face_id_0.size(); i++)
		{
			indices[3 * i] = mesh_face_id_0[i];
			indices[3 * i + 1] = mesh_face_id_1[i];
			indices[3 * i + 2] = mesh_face_id_2[i];
		}
		gte::MeshCurvature<double> meshCurvature;
		meshCurvature((size_t)verticeSize, points, (size_t)(faceindiceSize / 3), indices, 1e-06f);

		Vector3d1 normals;
		for (int i = 0; i < verticeSize; i++)
		{
			gte::Vector3<double> normal = meshCurvature.GetNormals()[i];
			normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
		}

		for (int i = 0; i < mesh_pointses.size(); i++)
		{
			Vector3d1 mesh_normals;
			for (int j = 0; j < mesh_pointses[i].size(); j++)
			{
				Vector3d n = RelatedFaceNormal(polyhedron, tree, normals, mesh_pointses[i][j]);
				mesh_normals.push_back(n);
			}
			mesh_normalses.push_back(mesh_normals);
			Vector3d1().swap(mesh_normals);
		}
	}

	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		Vector3d1 vecs;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		Functs::LoadObj3d(path_, vecs, face_id_0, face_id_1, face_id_2);

		Polyhedron_3 polyhedron;

		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//compute normals
		/**********************************************************/
		Vector3d1 normals(vecs.size(), Vector3d(0.0, 0.0, 0.0));
		std::vector<double> areas(vecs.size(), 0.0);
		for (int i = 0; i < face_id_0.size(); i++)
		{
			double area = CGAL_3D_One_Triangle_Area(vecs[face_id_0[i]], vecs[face_id_1[i]], vecs[face_id_2[i]]);
			Point_3 p0(vecs[face_id_0[i]][0], vecs[face_id_0[i]][1], vecs[face_id_0[i]][2]);
			Point_3 p1(vecs[face_id_1[i]][0], vecs[face_id_1[i]][1], vecs[face_id_1[i]][2]);
			Point_3 p2(vecs[face_id_2[i]][0], vecs[face_id_2[i]][1], vecs[face_id_2[i]][2]);
			Vector_3   n = CGAL::cross_product(p2 - p1, p0 - p1);

			Vector3d n0(n[0], n[1], n[2]);

			normals[face_id_0[i]] += (double)area * n0;
			normals[face_id_1[i]] += (double)area * n0;
			normals[face_id_2[i]] += (double)area * n0;
			areas[face_id_0[i]] += area;
			areas[face_id_1[i]] += area;
			areas[face_id_2[i]] += area;
		}

		for (int i = 0; i < vecs.size(); i++)
		{
			normals[i] = normals[i] / (double)areas[i];
		}
		/**********************************************************/

		for (int i = 0; i < mesh_pointses.size(); i++)
		{
			Vector3d1 mesh_normals;
			for (int j = 0; j < mesh_pointses[i].size(); j++)
			{
				Vector3d n = RelatedFaceNormal(polyhedron, tree, normals, mesh_pointses[i][j]);
				mesh_normals.push_back(n);
			}
			mesh_normalses.push_back(mesh_normals);
			Vector3d1().swap(mesh_normals);
		}

	}
}
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Normal_C1(const Vector3d1 & ps, const std::vector<std::vector<int>>&face_ids, Vector3d1 & normals)
{
	int verticeSize = ps.size();
	int faceindiceSize = face_ids.size() * 3;

	gte::Vector3<double>* points = new gte::Vector3<double>[verticeSize];
	unsigned int* indices = new unsigned int[faceindiceSize];

	for (int i = 0; i < verticeSize; i++)
	{
		points[i][0] = ps[i][0];
		points[i][1] = ps[i][1];
		points[i][2] = ps[i][2];
	}
	for (int i = 0; i < face_ids.size(); i++)
	{
		indices[3 * i] = face_ids[i][0];
		indices[3 * i + 1] = face_ids[i][1];
		indices[3 * i + 2] = face_ids[i][2];
	}
	gte::MeshCurvature<double> meshCurvature;
	meshCurvature((size_t)verticeSize, points, (size_t)(faceindiceSize / 3), indices, 1e-06f);

	for (int i = 0; i < verticeSize; i++)
	{
		gte::Vector3<double> normal = meshCurvature.GetNormals()[i];
		normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
	}
}
extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Normal_C2(const Vector3d1 & ps, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, Vector3d1 & normals)
{
	int verticeSize = ps.size();
	int faceindiceSize = face_id_0.size() * 3;

	gte::Vector3<double>* points = new gte::Vector3<double>[verticeSize];
	unsigned int* indices = new unsigned  int[faceindiceSize];

	for (int i = 0; i < verticeSize; i++)
	{
		points[i][0] = ps[i][0];
		points[i][1] = ps[i][1];
		points[i][2] = ps[i][2];
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
		gte::Vector3<double> normal = meshCurvature.GetNormals()[i];
		normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
	}
}


extern "C" PPGL_EXPORT Vector3d CGAL_3D_Mesh_Center_C1(const Vector3d2 & ps)
{
	Vector3d1 ps1;
	for (auto ps_ : ps)for (auto p : ps_)ps1.emplace_back(p);
	return CGAL_3D_Mesh_Center_C2(ps1);
}
extern "C" PPGL_EXPORT Vector3d CGAL_3D_Mesh_Center_C2(const Vector3d1 & ps)
{
	Vector3d center(0.0, 0.0, 0.0);

	for (int i = 0; i < ps.size(); i++)
	{
		center += ps[i];
	}

	center = center / (double)ps.size();
	return center;
}

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Boundingbox_C1(const Vector3d2 & ps, Vector3d & min_corner, Vector3d & max_corner)
{
	Vector3d1 ps1;
	for (auto ps_ : ps)for (auto p : ps_)ps1.emplace_back(p);
	CGAL_3D_Mesh_Boundingbox_C2(ps1, min_corner, max_corner);
}

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Boundingbox_C2(const Vector3d1 & ps, Vector3d & min_corner, Vector3d & max_corner)
{
	min_corner = ps[0];
	max_corner = ps[0];
	for (int i = 0; i < ps.size(); i++)
	{
		min_corner[0] = std::min(min_corner[0], ps[i][0]);
		min_corner[1] = std::min(min_corner[1], ps[i][1]);
		min_corner[2] = std::min(min_corner[2], ps[i][2]);
		max_corner[0] = std::max(max_corner[0], ps[i][0]);
		max_corner[1] = std::max(max_corner[1], ps[i][1]);
		max_corner[2] = std::max(max_corner[2], ps[i][2]);
	}
}

extern "C" PPGL_EXPORT void CGAL_Surface_Decomposition(const char* path, std::vector<double>&face_sdf, int& regions_nb, std::vector<int>&face_segments)
{
	Polyhedron_3 polyhedron;

	Construct_Polyhedron(polyhedron, path);

	// create a property-map for SDF values
	typedef std::map<Polyhedron_3::Facet_const_handle, double> Facet_double_map;
	Facet_double_map internal_sdf_map;
	boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

	// compute SDF values using default parameters for number of rays, and cone angle
	//CGAL::sdf_values(polyhedron, sdf_property_map);

	int index = 0;
	for (Polyhedron_3::Facet_const_iterator facet_it = polyhedron.facets_begin();
		facet_it != polyhedron.facets_end(); ++facet_it) {
		sdf_property_map[facet_it] = face_sdf[index];
		index++;
	}

	// create a property-map for segment-ids
	typedef std::map<Polyhedron_3::Facet_const_handle, std::size_t> Facet_int_map;
	Facet_int_map internal_segment_map;
	boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

	// segment the mesh using default parameters for number of levels, and smoothing lambda
	// Any other scalar values can be used instead of using SDF values computed using the CGAL function
	//std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(polyhedron, sdf_property_map, segment_property_map);

	//regions_nb = number_of_segments;

	const std::size_t number_of_clusters = 100;       // use 4 clusters in soft clustering
	const double smoothing_lambda = 0.08;  // importance of surface features, suggested to be in-between [0,1]
	regions_nb = CGAL::segmentation_from_sdf_values(
		polyhedron, sdf_property_map, segment_property_map, (int)number_of_clusters, smoothing_lambda);

	for (Polyhedron_3::Facet_const_iterator facet_it = polyhedron.facets_begin();
		facet_it != polyhedron.facets_end(); ++facet_it) {
		face_segments.push_back(segment_property_map[facet_it]);
	}
}

extern "C" PPGL_EXPORT void CGAL_Intergral_Curvature(const Vector2d1 & input_points, const int& sampling_points_nb, const double& radius, const double& thresholder, Vector2d1 & output_points, std::vector<double>&output_rates)
{
	auto Strip_Get_Total_length=[](const Vector2d1 & input_points)
	{
		double total_length = 0.0;
		if (input_points.size() >= 2) {
			for (int i = 0; i < input_points.size() - 1; i++) {
				total_length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[i + 1]);
			}
		}
		return total_length;
	};

	auto Strip_Get_One_Point_From_Strip=[](double d, const Vector2d1 & input_points, int& index)
	{
		Vector2d v;

		double length = 0.0;
		for (int i = 0; i < input_points.size() - 1; i++)
		{
			length += sqrt((double)CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]));
		}
		double total_length = length;

		length = 0.0;

		for (int i = 0; i < input_points.size() - 1; i++)
		{
			double l = sqrt((double)CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]));

			if (d * total_length >= length && d * total_length <= length + l)
			{
				double ll = (d - length / total_length) * total_length / l;
				v[0] = input_points[i].x + (input_points[(i + 1) % input_points.size()].x - input_points[i].x) * ll;
				v[1] = input_points[i].y + (input_points[(i + 1) % input_points.size()].y - input_points[i].y) * ll;
				index = i;
				break;
			}
			length += l;
		}
		return v;
	};

	auto GenerateACircle=[](int divided_nb, Vector2d v, double distance, Vector2d1 & circle_points)
	{
		for (int i = 0; i < divided_nb; i++)
		{
			Vector2d p0(v[0] + abs(distance) * sin(i * 2 *Math_PI / (double)divided_nb), v[1] + abs(distance) * cos(i * 2 * Math_PI / (double)divided_nb));
			circle_points.push_back(p0);
		}
	};

	auto CircleIntersectWithSegment=[](Vector2d1 & circle_points, Vector2d v0, Vector2d v1)
	{
		for (int i = 0; i < circle_points.size(); i++)
		{
			Vector2d inter;
			if (CGAL_2D_Intersection_Segment_Segment(circle_points[i], circle_points[(i + 1) % circle_points.size()], v0, v1, inter))
			{
				return inter;
			}
		}
		Vector2d v(0.0, 0.0);
		return v;
	};

	auto FindLastIndexPoint=[&](const Vector2d1 & input_points, Vector2d center, double radius, Vector2d1 & circle_points, int start_index, Vector2d1 & output_points)
	{
		for (int i = start_index; i >= 0; i--)
		{
			if (i >= 0 && i < input_points.size())
			{
				double distance = CGAL_2D_Distance_Point_Point(input_points[i], center);

				if (CGAL_2D_Location_Point_Polygon(input_points[i], circle_points))
				{
					output_points.push_back(input_points[i]);
				}
				else
				{
					Vector2d v;
					if (output_points.size() > 0)
					{
						v = CircleIntersectWithSegment(circle_points, input_points[i], output_points[output_points.size() - 1]);
					}
					else
					{
						v = CircleIntersectWithSegment(circle_points, input_points[i], center);
					}
					output_points.push_back(v);
					break;
				}
			}
		}
	};


	auto FindNextIndexPoint=[&](const Vector2d1& input_points, Vector2d center, double radius, Vector2d1& circle_points, int start_index, Vector2d1& output_points)
	{
		for (int i = start_index; i < input_points.size(); i++)
		{
			if (i >= 0 && i < input_points.size())
			{
				if (CGAL_2D_Location_Point_Polygon(input_points[i], circle_points))
				{
					output_points.push_back(input_points[i]);
				}
				else
				{
					Vector2d v;

					if (output_points.size() > 0)
					{
						v = CircleIntersectWithSegment(circle_points, input_points[i], output_points[output_points.size() - 1]);
					}
					else
					{
						v = CircleIntersectWithSegment(circle_points, input_points[i], center);
					}
					output_points.push_back(v);
					break;
				}
			}
		}
	};

	auto Circuit_Get_Total_length=[](Vector2d1& input_points)
	{
		double total_length = 0.0;
		for (int i = 0; i < input_points.size(); i++) {
			total_length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
		}
		return total_length;
	};

	auto Circuit_Find_Nearest_Point_Par=[&](Vector2d v, Vector2d1& input_points)
	{
		Vector2d n_p;

		double total_length = Circuit_Get_Total_length(input_points);

		double min_d = DBL_MAX;
		int min_i = -1;
		for (int i = 0; i < input_points.size(); i++)
		{
			Vector2d p0(input_points[i].x, input_points[i].y);
			Vector2d p1(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y);

			double l = CGAL_2D_Distance_Point_Segment(v, p0, p1);

			if (l < min_d)
			{
				min_d = l;
				min_i = i;
			}
		}

		if (min_i >= 0)
		{
			double length = 0.0;
			for (int i = 0; i < min_i; i++)
			{
				length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
			}

			Vector2d p0(input_points[min_i].x, input_points[min_i].y);
			Vector2d p1(input_points[(min_i + 1) % input_points.size()].x, input_points[(min_i + 1) % input_points.size()].y);

			double l0 = CGAL_2D_Distance_Point_Point(v, p0);
			double l1 = CGAL_2D_Distance_Point_Point(v, p1);

			if (min_d < l0 && min_d < l1 && abs(min_d - l0)>0.0001 && abs(min_d - l1)>0.0001)
			{
				double l = CGAL_2D_Distance_Point_Point(p0, p1);
				if (l < 0.00001)
				{
					v[0] = p0[0];
					v[1] = p0[1];
				}
				else
				{
					Vector2d vec(p1[0] - p0[0], p1[1] - p0[1]);
					Vector2d r_vec(vec[1], -vec[0]);
					if (vec[0] < 0.00001)
					{
						r_vec[0] = -vec[1];
						r_vec[1] = vec[0];
					}

					if (vec[1] < 0.00001)
					{
						r_vec[0] = vec[1];
						r_vec[1] = -vec[0];
					}

					if (CGAL_2D_Intersection_Segment_Line(p0, p1, v, v + r_vec, n_p))
					{

					}
					else
					{
						assert(false);
					}
				}
			}
			else
			{
				if (l0 < l1)
				{
					n_p[0] = p0[0];
					n_p[1] = p0[1];
				}
				else
				{
					n_p[0] = p1[0];
					n_p[1] = p1[1];
				}

			}

			length += CGAL_2D_Distance_Point_Point(input_points[min_i], n_p);

			return length / total_length;
		}
		else
		{
			assert(false);
		}

		return -1.0;
	};


	auto Circuit_Get_One_Point_From_Offset=[](double d, Vector2d1& input_points)
	{
		Vector2d v;
		double length = 0.0;
		for (int i = 0; i < input_points.size(); i++)
		{
			length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
		}

		double total_length = length;
		length = 0.0;

		for (int i = 0; i < input_points.size(); i++)
		{
			double l = CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);

			if (d * total_length >= length && d * total_length <= length + l)
			{
				double ll = (d - length / total_length) * total_length / l;
				v[0] = input_points[i].x + (input_points[(i + 1) % input_points.size()].x - input_points[i].x) * ll;
				v[1] = input_points[i].y + (input_points[(i + 1) % input_points.size()].y - input_points[i].y) * ll;
				break;
			}
			length += l;
		}

		return v;
	};



	auto Circuit_Select_One_Part_Offset=[&](Vector2d1& input_points, double d0, double d1, Vector2d1& vecs)
	{
		if (abs(d0 - d1) < 0.0000001)
		{
			Vector2d v = Circuit_Get_One_Point_From_Offset(d0, input_points);
			vecs.push_back(v);
		}
		else
		{
			Vector2d v = Circuit_Get_One_Point_From_Offset(d0, input_points);
			vecs.push_back(v);

			if (d1 >= 0 && d1 <= 1.0)
			{
				double total_length = Circuit_Get_Total_length(input_points);

				std::vector<double> vec_ds;

				vec_ds.push_back(0.0);
				double length = 0.0;
				for (int i = 0; i < input_points.size(); i++)
				{
					length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
					vec_ds.push_back(length / total_length);
				}


				if (d0 > d1)
				{
					for (int i = vec_ds.size() - 1; i >= 0; i--)
					{
						if (vec_ds[i]<d0 && vec_ds[i]>d1)
						{
							v = Circuit_Get_One_Point_From_Offset(vec_ds[i], input_points);

							if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001 && abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
							{
								vecs.push_back(v);
							}
						}

						if (vec_ds[i] < d1)
						{
							break;
						}
					}
				}
				else
				{
					for (int i = vec_ds.size() - 1; i > 0; i--)
					{
						if (vec_ds[i] < d0)
						{
							v = Circuit_Get_One_Point_From_Offset(vec_ds[i], input_points);

							if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001 && abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
							{
								vecs.push_back(v);
							}
						}
					}

					for (int i = vec_ds.size() - 1; i > 0; i--)
					{
						if (vec_ds[i] > d1)
						{
							v = Circuit_Get_One_Point_From_Offset(vec_ds[i], input_points);
							if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001 && abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
							{
								vecs.push_back(v);
							}
						}
					}
				}

				v = Circuit_Get_One_Point_From_Offset(d1, input_points);
				if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001 && abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
				{
					vecs.push_back(v);
				}

				if (abs(vecs[1][0] - vecs[0][0]) < 0.000001 && abs(vecs[1][1] - vecs[0][1]) < 0.000001)
				{
					vecs.erase(vecs.begin());
				}

				std::vector<double>().swap(vec_ds);
			}
		}
	};

	auto GetHalfCircleArea=[&](Vector2d1& last_points, Vector2d1& next_points, Vector2d center, double radius, Vector2d1& circle_points)
	{
		Vector2d1 polygon_points;
		polygon_points.push_back(center);

		for (int i = 0; i < last_points.size() - 1; i++)
			polygon_points.push_back(last_points[i]);

		double par_0 = Circuit_Find_Nearest_Point_Par(last_points[last_points.size() - 1], circle_points);
		double par_1 = Circuit_Find_Nearest_Point_Par(next_points[next_points.size() - 1], circle_points);
		Circuit_Select_One_Part_Offset(circle_points, par_0, par_1, polygon_points);

		for (int i = next_points.size() - 2; i >= 0; i--)
			polygon_points.push_back(next_points[i]);

		return CGAL_2D_Polygon_Area(polygon_points);
	};

	double total_length = Strip_Get_Total_length(input_points);

	for (int i = 0; i < sampling_points_nb; i++)
	{
		if (i % (int)(sampling_points_nb / 100.0) == 0)
		{
			std::cout << "Intergral: " << i << " / " << sampling_points_nb << std::endl;
		}

		double par = i * 1.0 / (double)sampling_points_nb;

		if (par > radius / total_length && par < 1.0 - radius / total_length)
		{

			int last_index;
			Vector2d center = Strip_Get_One_Point_From_Strip(par, input_points, last_index);

			Vector2d1 circle_points;
			GenerateACircle(30, center, radius, circle_points);

			Vector2d1 last_points;
			FindLastIndexPoint(input_points, center, radius, circle_points, last_index, last_points);
			Vector2d1 next_points;
			FindNextIndexPoint(input_points, center, radius, circle_points, last_index + 1, next_points);

			if (last_points.size() > 0 && next_points.size() > 0)
			{
				double area = GetHalfCircleArea(last_points, next_points, center, radius, circle_points);

				if (area > Math_PI * radius * radius * 0.5)
				{
					area = Math_PI * radius * radius - area;
				}

				area = area / (Math_PI * radius * radius);

				if (area < thresholder)
				{
					output_points.push_back(center);
					output_rates.push_back(area);
				}
			}
		}

	}
}

extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Gradient(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const std::vector<double>&psd, Vector3d1 & vecs_gradients, Vector3d1 & faces_gradients)
{
	Vector3d1().swap(vecs_gradients);
	Vector3d1().swap(faces_gradients);

	for (int i = 0; i < vecs.size(); i++)
		vecs_gradients.push_back(Vector3d(0.0, 0.0, 0.0));

	std::vector<double> areas(vecs.size(), 0.0);

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		Vector3d v_0 = vecs[index_0];
		Vector3d v_1 = vecs[index_1];
		Vector3d v_2 = vecs[index_2];

		Vector3d   n = Functs::GetCrossproduct(v_0 - v_2, v_1 - v_2);
		Functs::SetVectorLength(n, 1.0);

		double area = CGAL_3D_One_Triangle_Area(v_0, v_1, v_2);

		Vector3d gradient(0.0, 0.0, 0.0);

		double d_0 = (double)psd[index_0];
		double d_1 = (double)psd[index_1];
		double d_2 = (double)psd[index_2];

		gradient += (double)psd[index_0] * Functs::GetCrossproduct(n, v_2 - v_1);
		gradient += (double)psd[index_1] * Functs::GetCrossproduct(n, v_0 - v_2);
		gradient += (double)psd[index_2] * Functs::GetCrossproduct(n, v_1 - v_0);

		Vector3d face_gradient = gradient / (double)(2.0 * area);

		faces_gradients.push_back(face_gradient);

		vecs_gradients[index_0] += (double)area * face_gradient;
		areas[index_0] += area;

		vecs_gradients[index_1] += (double)area * face_gradient;
		areas[index_1] += area;

		vecs_gradients[index_2] += (double)area * face_gradient;
		areas[index_2] += area;

	}

	for (int i = 0; i < vecs.size(); i++)
	{
		vecs_gradients[i] = vecs_gradients[i] / (double)areas[i];
	}
}


extern "C" PPGL_EXPORT bool CGAL_3D_Mesh_Extract_Isoline(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const std::vector<double>&psd, const double& d, Vector3d2 & isolines)
{
	auto InterpolationPoint=[](Vector3d v0, Vector3d v1, double psd0, double  psd1, double d, Vector3d & inter)
	{
		if (Functs::IsAlmostZero(d - psd0)) psd0 = d;
		if (Functs::IsAlmostZero(psd1 - d)) psd1 = d;

		if (!Functs::IsAlmostZero(psd1 - psd0) && ((d >= psd0 && d <= psd1) || (d >= psd1 && d <= psd0)))
		{
			if (d >= psd1 && d <= psd0)
			{
				Vector3d v = v0;
				v0 = v1;
				v1 = v;

				double d = psd0;
				psd0 = psd1;
				psd1 = d;
			}
			inter = (double)((d - psd0) / (psd1 - psd0)) * v1 + (double)((psd1 - d) / (psd1 - psd0)) * v0;
			return true;
		}
		else
		{
			return false;
		}
	};

	//compute gredient
	Vector3d1 vecs_gradients;
	Vector3d1 faces_gradients;
	CGAL_3D_Mesh_Gradient(vecs, face_id_0, face_id_1, face_id_2, psd, vecs_gradients, faces_gradients);

	Vector3d2 segments;
	VectorPI2 edges;

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		Vector3d v_0 = vecs[index_0];
		Vector3d v_1 = vecs[index_1];
		Vector3d v_2 = vecs[index_2];

		Vector3d  n = Functs::GetCrossproduct(v_0 - v_2, v_1 - v_2);
		Vector3d gradient = faces_gradients[i];

		Vector3d inter_0_1;
		Vector3d inter_1_2;
		Vector3d inter_2_0;

		bool b_0 = InterpolationPoint(v_0, v_1, psd[index_0], psd[index_1], d, inter_0_1);
		bool b_1 = InterpolationPoint(v_1, v_2, psd[index_1], psd[index_2], d, inter_1_2);
		bool b_2 = InterpolationPoint(v_2, v_0, psd[index_2], psd[index_0], d, inter_2_0);

		if (b_0 && b_1)
		{
			Vector3d1 segment;
			segment.push_back(inter_0_1);
			segment.push_back(inter_1_2);

			VectorPI1 endedge;
			endedge.push_back(std::pair<int,int>(index_0, index_1));
			endedge.push_back(std::pair<int, int>(index_1, index_2));

			if (!Functs::IsAlmostZero(CGAL_3D_Distance_Point_Point(segment[0], segment[1])))
			{
				Vector3d  gradient_n = Functs::GetCrossproduct(segment[0] - segment[1], gradient);
				double angle = Functs::GetAngleBetween(gradient_n, n);
				if (angle > Math_PI / 2.0)
					std::reverse(segment.begin(), segment.end());

				segments.push_back(segment);
				edges.push_back(endedge);
			}
		}
		if (b_1 && b_2)
		{
			Vector3d1 segment;
			segment.push_back(inter_1_2);
			segment.push_back(inter_2_0);

			VectorPI1 endedge;
			endedge.push_back(std::pair<int,int>(index_1, index_2));
			endedge.push_back(std::pair<int, int>(index_2, index_0));

			if (!Functs::IsAlmostZero(CGAL_3D_Distance_Point_Point(segment[0], segment[1])))
			{
				Vector3d  gradient_n = Functs::GetCrossproduct(segment[0] - segment[1], gradient);
				double angle = Functs::GetAngleBetween(gradient_n, n);
				if (angle > Math_PI / 2.0)
					std::reverse(segment.begin(), segment.end());

				segments.push_back(segment);
				edges.push_back(endedge);
			}
		}
		if (b_2 && b_0)
		{
			Vector3d1 segment;
			segment.push_back(inter_2_0);
			segment.push_back(inter_0_1);

			VectorPI1 endedge;
			endedge.push_back(std::pair<int,int>(index_2, index_0));
			endedge.push_back(std::pair<int,int>(index_0, index_1));

			if (!Functs::IsAlmostZero(CGAL_3D_Distance_Point_Point(segment[0], segment[1])))
			{
				Vector3d  gradient_n = Functs::GetCrossproduct(segment[0] - segment[1], gradient);
				double angle = Functs::GetAngleBetween(gradient_n, n);
				if (angle > Math_PI / 2.0)
					std::reverse(segment.begin(), segment.end());

				segments.push_back(segment);
				edges.push_back(endedge);
			}
		}
	}

	std::cout << "CGAL_3D_Mesh_Extract_Isoline: " << segments.size() << std::endl;


	//std::ofstream export_fie("G:\\segments.obj");
	//int export_int = 1;

	//for (int i = 0; i < segments.size(); i++)
	//{
	//	Functs::Export_Segment(export_fie, export_int, "segment_" + Functs::IntString(i), 1.0, 0.0, 0.0, segments[i][0], segments[i][1], 0.002);
	//}

	//export_fie.clear();
	//export_fie.close();


	if (segments.size() == 0) return false;

	CGAL_3D_Connecting_Segments_C2(segments, isolines);

	//std::ofstream export_fie_0("G:\\isolines.obj");
	//export_int = 1;

	//for (int i = 0; i < isolines.size(); i++)
	//{
	//	for (int j = 0; j < isolines.size()-1; j++)
	//	{
	//		Functs::Export_Segment(export_fie_0, export_int, "isolines_" + Functs::IntString(i) + "_" + Functs::IntString(j), 1.0, 0.0, 0.0, isolines[i][j], isolines[i][(j + 1)], 0.002);
	//	}

	//}

	//export_fie_0.clear();
	//export_fie_0.close();


	return true;
}
extern "C" PPGL_EXPORT void CGAL_BSplineCurveFit(const Vector3d1 & samples, Vector3d1 & output)
{
	/********************************/
	//float mult = 2.0f / (1000 - 1), t;
	//int i;
	//for (i = 0; i < 1000; ++i)
	//{
	//	t = -1.0f + mult*i;
	//	float angle = 2.0f*gte::Mathf::TWO_PI*t;
	//	float amplitude = 1.0f - t*t;
	//	Vector3d v(amplitude*gte::Mathf::Cos(angle), amplitude*gte::Mathf::Sin(angle),t);
	//	samples.push_back(v);
	//}
	/********************************/

	//int dimension = 3;
	//int numSamples = samples.size();

	//gte::Vector3<double>;

	//gte::Vector3d* mSamples;
	//mSamples = new1<gte::Vector3d>(numSamples);
	//for (int i = 0; i < numSamples; ++i)
	//{
	//	mSamples[i].X() = samples[i][0];
	//	mSamples[i].Y() = samples[i][1];
	//	mSamples[i].Z() = samples[i][2];
	//}

	//int  mNumCtrlPoints = numSamples / 2;
	//int mDegree = 6;

	//gte::BSplineCurveFitd* mSpline;

	//mSpline = new0 gte::BSplineCurveFitd(dimension, numSamples, (const double*)mSamples, mDegree, mNumCtrlPoints);

	//// Sample it the same number of times as the original data.

	//float mult = 1.0 / (numSamples * 2 - 1);
	//for (int i = 0; i < numSamples * 2; ++i)
	//{
	//	double* pos = new double[3];
	//	mSpline->GetPosition(mult * i, pos);
	//	output.push_back(Vector3d(pos[0], pos[1], pos[2]));
	//}
}

//assumption: one projecting line can intersect with the three edges at most twice times
//it's impossible to meet one edge for more than one times
extern "C" PPGL_EXPORT void CGAL_Cut_Surface(const Vector3d1 & boundary, const Vector3d & inside_point, const char* full_path, char* output_path)
{
	Vector3d2 multi_boundary(1, boundary);
	CGAL_Cut_Surface_by_Multi_Boundaries(multi_boundary, inside_point, full_path, output_path);
	return;
}

bool DebugInformation()
{
	return true;
}


void ComputeEdgeLables(const int size_of_vertice, Halfedge_handle& start_hh, std::vector<std::pair<int, int>>& edges, std::vector<int>& lables)
{
	std::vector<int> two_class_lables(size_of_vertice, 0);
	std::vector<bool> edge_lables(size_of_vertice, false);
	std::vector<bool> over(size_of_vertice, false);
	for (int i = 0; i < edges.size(); i++)
	{
		edge_lables[edges[i].first] = true;
		edge_lables[edges[i].second] = true;
	}

	std::queue<Halfedge_handle> queue;
	std::queue<int> queue_1;

	queue.push(start_hh);
	queue_1.push(-1);
	over[start_hh->vertex()->id()] = true;

	std::vector<std::pair<Vector3d, int>> dsd;

	int iter = 0;
	while (queue.size() != 0)
	{
		Halfedge_handle hh = queue.front();
		int id = hh->vertex()->id();

		two_class_lables[hh->vertex()->id()] = 1;
		Vector3d v = PointVector3d(hh->vertex()->point());

		if (iter % 100 == 0)
			std::cerr << iter << ": " << hh->vertex()->id() << " / " << size_of_vertice << std::endl;
		//////////

		dsd.emplace_back(std::pair<Vector3d, int>(v, queue_1.front()));

		if (false)
		{
			std::ofstream  export_file_output_0("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting\\" + std::to_string(iter) + "_" + std::to_string(hh->vertex()->id()) + "_" + std::to_string(queue_1.front()) + ".obj");
			int export_index_0 = 1;
			std::string str = "inside" + std::to_string(iter) + "_" + std::to_string(hh->vertex()->id()) + "_" + std::to_string(queue_1.front());
			Functs::Export_Point(export_file_output_0, export_index_0,
				str.c_str(), v, 0.1);
			export_file_output_0.clear();
			export_file_output_0.close();
		}


		queue.pop();
		queue_1.pop();

		Halfedge_handle iter_hh = hh;
		do
		{
			int opposite_id = iter_hh->opposite()->vertex()->id();

			if (two_class_lables[opposite_id] == 0 && !over[opposite_id])
			{
				if (!(edge_lables[id] && edge_lables[opposite_id]))
				{
					queue.push(iter_hh->opposite());
					queue_1.push(iter);
					over[opposite_id] = true;
				}
				else
				{
					bool b = false;
					for (int i = 0; i < edges.size(); i++)
					{
						if ((edges[i].first == id && edges[i].second == opposite_id) || (edges[i].second == id && edges[i].first == opposite_id))
						{
							b = true;
							break;
						}
					}
					if (!b)
					{
						queue.push(iter_hh->opposite());
						queue_1.push(iter);
						over[opposite_id] = true;
					}
				}
			}

			if (!(edge_lables[id] && edge_lables[opposite_id]) && two_class_lables[opposite_id] == 0 && !over[opposite_id])
			{
				queue.push(iter_hh->opposite());
				queue_1.push(iter);
				over[opposite_id] = true;
			}

			iter_hh = iter_hh->opposite()->prev();
		} while (iter_hh != hh);

		iter++;
	}


	auto AAAA = [](std::vector<std::pair<Vector3d, int>>& dsd)
	{
		std::ofstream  export_file_output_0("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\debug.obj");
		int export_index_0 = 1;
		int start_search = 22756;
		while (true)
		{
			if (start_search >= 0)
				Functs::Export_Point(export_file_output_0, export_index_0, "point", dsd[start_search].first, 0.1);
			else
				break;

			start_search = dsd[start_search].second;
		}

		export_file_output_0.clear();
		export_file_output_0.close();
	};

	//AAAA(dsd);

	lables = two_class_lables;
}


void CGAL_2D_Polygon_Simple_0(Vector2d1 points_2d)
{
	for (int i = 0; i < points_2d.size(); i++)
	{
		//i//i+1
		auto s0 = points_2d[i];
		auto e0 = points_2d[(i + 1) % points_2d.size()];

		for (int j = i + 2; j < points_2d.size() + i - 1; j++)
		{
			auto s1 = points_2d[(j) % points_2d.size()];
			auto e1 = points_2d[(j + 1) % points_2d.size()];
			Vector2d inter;
			if (CGAL_2D_Intersection_Segment_Segment(s0, e0, s1, e1, inter))
			{
				std::cerr << i << " " << (j) % points_2d.size() << std::endl;
			}
		}
	}
}

void ComputeRemeshTriangles(const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2,
	const std::vector<std::pair<int, int>>& edges, const Vector3d1& cutting_points, const std::vector<int>& multi_cps, const std::vector<int>& two_class_lables, const std::string output_path)
{
	auto HandleBoundaryTriangle = [](const Vector3d1& vecs, const std::vector<int>& remesh_lables, const std::vector<std::pair<int, int>>& edges, const Vector3d1& cutting_points,
		std::vector<std::vector<int>>& remesh_triangles, std::vector<int>& remesh_triangles_lables, Vector3d1& new_points, std::vector<int>& cutting_point_ids,
		const int index_0, const int index_1, const int index_2)
	{
		auto GetEdgeMiddlePoint = [](const int index_0, const int index_1, const std::vector<std::pair<int, int>>& edges)
		{
			for (int i = 0; i < edges.size(); i++)
				if ((edges[i].first == index_0 && edges[i].second == index_1) || (edges[i].second == index_0 && edges[i].first == index_1))
					return i;
			std::cerr << "GetEdgeMiddlePoint Error..." << std::endl;
			system("pause");
			return -1;
		};

		Vector3d v_0 = vecs[index_0];
		Vector3d v_1 = vecs[index_1];
		Vector3d v_2 = vecs[index_2];

		int lable_0 = remesh_lables[index_0];
		int lable_1 = remesh_lables[index_1];
		int lable_2 = remesh_lables[index_2];

		int cut_0_1 = GetEdgeMiddlePoint(index_0, index_1, edges);
		int cut_0_2 = GetEdgeMiddlePoint(index_0, index_2, edges);

		Vector3d v_0_1 = cutting_points[cut_0_1];
		Vector3d v_0_2 = cutting_points[cut_0_2];

		bool b_0 = CGAL_3D_Distance_Point_Point(v_0, v_0_1) <= 0.001;
		bool b_1 = CGAL_3D_Distance_Point_Point(v_1, v_0_1) <= 0.001;
		bool b_2 = CGAL_3D_Distance_Point_Point(v_0, v_0_2) <= 0.001;
		bool b_3 = CGAL_3D_Distance_Point_Point(v_2, v_0_2) <= 0.001;

		std::vector<int> triangle;
		triangle.push_back(index_0);
		triangle.push_back(index_1);
		triangle.push_back(index_2);

		//A: 0.0-0.0
		if (b_1 && b_3)
		{
			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_0);

			cutting_point_ids[cut_0_1] = index_1;
			cutting_point_ids[cut_0_2] = index_2;

			return;
		}

		//B: 0.0-0.5
		if (b_1 && !b_2 && !b_3)
		{
			new_points.push_back(v_0_2);
			int index_0_2 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_1, index_0_2 };
			std::vector<int> triangle_1 = { index_2, index_0_2, index_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_0);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);


			cutting_point_ids[cut_0_1] = index_1;
			cutting_point_ids[cut_0_2] = index_0_2;

			return;
		}

		//C: 0.5 0.0
		if (!b_0 && !b_1 && b_3)
		{
			new_points.push_back(v_0_1);
			int index_0_1 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_0_1, index_2 };
			std::vector<int> triangle_1 = { index_1, index_2, index_0_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_0);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);


			cutting_point_ids[cut_0_1] = index_0_1;
			cutting_point_ids[cut_0_2] = index_2;

			return;
		}

		//D: 0.5 0.5
		if (!b_0 && !b_1 && !b_2 && !b_3)
		{
			new_points.push_back(v_0_1);
			new_points.push_back(v_0_2);

			int index_0_1 = new_points.size() - 2 + vecs.size();
			int index_0_2 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_0_1, index_0_2 };
			std::vector<int> triangle_1 = { index_1, index_0_2, index_0_1 };
			std::vector<int> triangle_2 = { index_2, index_0_2, index_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_0);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);

			remesh_triangles.push_back(triangle_2);
			remesh_triangles_lables.push_back(lable_2);

			cutting_point_ids[cut_0_1] = index_0_1;
			cutting_point_ids[cut_0_2] = index_0_2;
			return;
		}

		//E: 0.5 1.0
		if (!b_0 && !b_1 && b_2 && !b_3)
		{
			new_points.push_back(v_0_1);
			int index_0_1 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_0_1, index_2 };
			std::vector<int> triangle_1 = { index_1, index_2, index_0_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_1);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);

			cutting_point_ids[cut_0_1] = index_0_1;
			cutting_point_ids[cut_0_2] = index_0;

			return;
		}

		//F: 1.0 0.5 
		if (b_0 && !b_1 && !b_2 && !b_3)
		{

			new_points.push_back(v_0_2);
			int index_0_2 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_1, index_0_2 };
			std::vector<int> triangle_1 = { index_2, index_0_2, index_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_1);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);

			cutting_point_ids[cut_0_1] = index_0;
			cutting_point_ids[cut_0_2] = index_0_2;

			return;
		}

		//G: 1.0 1.0
		if (b_0 && b_2)
		{
			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_1);


			cutting_point_ids[cut_0_1] = index_0;
			cutting_point_ids[cut_0_2] = index_0;

		}

		//H: 1.0 0.0
		if (b_0 && b_3)
		{
			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_1);

			cutting_point_ids[cut_0_1] = index_0;
			cutting_point_ids[cut_0_2] = index_2;
		}

		//I: 0.0 1.0
		if (b_1 && b_2)
		{
			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_1);

			cutting_point_ids[cut_0_1] = index_1;
			cutting_point_ids[cut_0_2] = index_0;
		}

		remesh_triangles.push_back(triangle);
		remesh_triangles_lables.push_back(lable_1);
		return;
	};

	Vector3d1 	remesh_points = vecs;
	Vector3d1 add_points;
	std::vector<std::vector<int>> remesh_triangles;
	std::vector<int> remesh_triangles_lables;

	std::vector<int> cutting_point_ids(cutting_points.size(), -1);

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		int lable_0 = two_class_lables[index_0];
		int lable_1 = two_class_lables[index_1];
		int lable_2 = two_class_lables[index_2];

		if (lable_0 == lable_1 && lable_0 == lable_2)
		{
			std::vector<int> triangle;
			triangle.push_back(index_0);
			triangle.push_back(index_1);
			triangle.push_back(index_2);

			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_0);
			continue;
		}
		if (lable_1 == lable_2) HandleBoundaryTriangle(vecs, two_class_lables, edges, cutting_points, remesh_triangles, remesh_triangles_lables, add_points, cutting_point_ids, index_0, index_1, index_2);
		if (lable_0 == lable_1) HandleBoundaryTriangle(vecs, two_class_lables, edges, cutting_points, remesh_triangles, remesh_triangles_lables, add_points, cutting_point_ids, index_2, index_0, index_1);
		if (lable_0 == lable_2) HandleBoundaryTriangle(vecs, two_class_lables, edges, cutting_points, remesh_triangles, remesh_triangles_lables, add_points, cutting_point_ids, index_1, index_2, index_0);
	}


	auto GetAddPointIndex = [](Vector3d1& add_points, Vector3d v)
	{
		for (int i = 0; i < add_points.size(); i++)
		{
			double distance = CGAL_3D_Distance_Point_Point(add_points[i], v);
			if (Functs::IsAlmostZero(distance))
			{
				return i;
			}
		}
		add_points.push_back(v);
		return (int)(add_points.size() - 1);
	};


	std::vector<int> new_index;
	Vector3d1 new_points;
	for (int i = 0; i < add_points.size(); i++)
		new_index.push_back(GetAddPointIndex(new_points, add_points[i]));

	for (int i = 0; i < new_points.size(); i++)
		remesh_points.push_back(new_points[i]);

	for (int i = 0; i < remesh_triangles.size(); i++)
	{
		for (int j = 0; j < remesh_triangles[i].size(); j++)
		{
			if (remesh_triangles[i][j] >= vecs.size())
				remesh_triangles[i][j] = vecs.size() + new_index[remesh_triangles[i][j] - vecs.size()];
		}
	}

	///////////////////////////////////////////////////////
	int nb = multi_cps[0];
	std::vector<std::vector<int>>  cutting_id_points;
	std::vector<int> points;
	for (int i = 0; i < cutting_point_ids.size(); i++)
	{
		//if (cutting_point_ids[i] < 0)
		//{
		//	std::cerr << "cutting_point_ids[i] < 0" << std::endl;
		//	system("pause");
		//}

		if (cutting_point_ids[i] >= 0)
		{
			if (cutting_point_ids[i] >= vecs.size())
				cutting_point_ids[i] = vecs.size() + new_index[cutting_point_ids[i] - vecs.size()];

			if (i >= nb)
			{
				//remove duplicate points
				auto RmoveDuplicatePoints = [](const std::vector<int>& points)
				{
					std::vector<int> ndps;
					for (int i = 0; i < points.size(); i++)
					{
						if (i == 0)
							ndps.emplace_back(points[i]);
						else
						{
							if (points[i] != ndps.back())
								ndps.emplace_back(points[i]);
						}
					}

					if (ndps.front() == ndps.back()) ndps.erase(ndps.begin());
					return ndps;
				};

				auto ndps = RmoveDuplicatePoints(points);
				if (ndps.size() < 3)
				{
					std::cerr << "if (ndps.size() < 3)" << std::endl;
					system("pause");
				}

				cutting_id_points.emplace_back(points);
				nb = multi_cps[cutting_id_points.size()];
				points.clear();
			}

			points.emplace_back(cutting_point_ids[i]);
		}
	}
	cutting_id_points.emplace_back(points);
	///////////////////////////////////////////////////////


	if (true)
	{
		for (int i = 0; i < cutting_id_points.size(); i++)
		{
			std::ofstream  export_file_output("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting_point_ids_" + std::to_string(i) + ".obj");
			int export_index = 1;
			for (int j = 0; j < cutting_id_points[i].size(); j++)
			{
				Functs::Export_Segment(export_file_output, export_index, "intersection",Vector3d(1, 0, 0),
					remesh_points[cutting_id_points[i][j]], remesh_points[cutting_id_points[i][(j + 1) % cutting_id_points[i].size()]], 0.05);
			}
			export_file_output.clear();
			export_file_output.close();
		}
	}

	auto TriangleMesh = [](const std::vector<int>& ids, const Vector3d1& remesh_points,
		std::vector<std::vector<int>>& remesh_triangles, std::vector<int>& remesh_triangles_lables)
	{
		std::vector<int> nids;
		for (int i = 0; i < ids.size(); i++)
		{
			if (i == 0)
				nids.emplace_back(ids[i]);
			else
				if (ids[i] != nids.back())
					nids.emplace_back(ids[i]);
		}

		Vector3d1 points_3d;
		for (auto& id : nids) points_3d.emplace_back(remesh_points[id]);
		Vector2d1 points_2d = Functs::Vector3d2d(points_3d);

		if (!CGAL_2D_Polygon_Simple(points_2d))
		{
			CGAL_2D_Polygon_Simple_0(points_2d);

			std::cerr << "if (!CGAL_2D_Polygon_Simple(points))" << std::endl;
			system("pause");
		}

		auto triangles = CGAL_2D_Polygon_Triangulation_C3(points_2d);

		for (auto& triangle : triangles)
		{
			remesh_triangles.emplace_back(std::vector<int>{nids[triangle[2]], nids[triangle[1]], nids[triangle[0]]});
			remesh_triangles_lables.emplace_back(1);
		}
	};

	for (auto ids : cutting_id_points)
	{
		TriangleMesh(ids, remesh_points, remesh_triangles, remesh_triangles_lables);
	}


	//CGAL_Output_Boundary("G:\\abcd0.obj", remesh_points, remesh_triangles, remesh_triangles_lables);
	//CGAL_Output_Obj("G:\\remesh_full.obj", remesh_points, remesh_triangles);

	//CGAL_Output_Obj(path_0, remesh_points, remesh_triangles, remesh_triangles_lables, 0);
	Functs::OutputObj3d(output_path.c_str(), remesh_points, remesh_triangles, remesh_triangles_lables, 1);
}
//std::vector<std::vector<int>> surface_vectices_to_vectices;


//Cut a closed manifold mesh with boundaries
//multi_boundary: input boundaries
//inside_point: assign a point as the desired cutting surface
//full_path: file path of the input closed manifold mesh
//output_path: file path of output mesh
//assumption: one projecting line can intersect with the three edges of one triangle at most twice times
//            it's impossible to meet one edge for more than one times
extern "C" PPGL_EXPORT void CGAL_Cut_Surface_by_Multi_Boundaries(const Vector3d2 & multi_boundary, const Vector3d & inside_point, const char* full_path, char* output_path)
{
	auto Intersection = [](const Halfedge_handle& hh, const int nb, const Vector3d inside, const Vector3d outside, Halfedge_handle& handle, Vector3d& intersection)
	{
		std::vector<Vector3d> result_vecs;
		std::vector<Halfedge_handle> result_handles;

		Halfedge_handle hh_3d = hh;
		for (int i = 0; i < nb; i++)
		{
			hh_3d = hh_3d->next();
			Vector3d edge_0 = PointVector3d(hh_3d->vertex()->point());
			Vector3d edge_1 = PointVector3d(hh_3d->opposite()->vertex()->point());

			double inside_d = CGAL_3D_Distance_Point_Segment(inside, edge_0, edge_1);
			if (Functs::IsAlmostZero(inside_d))
			{
				result_vecs.emplace_back(inside);
				result_handles.emplace_back(hh_3d);
			}

			double outside_d = CGAL_3D_Distance_Point_Segment(outside, edge_0, edge_1);
			if (Functs::IsAlmostZero(outside_d))
			{
				result_vecs.emplace_back(outside);
				result_handles.emplace_back(hh_3d);
			}
		}


		Point_3 p0 = hh->next()->next()->vertex()->point();
		Point_3 p1 = hh->vertex()->point();
		Point_3 p2 = hh->next()->vertex()->point();
		Vector_3 n = CGAL::cross_product(p2 - p1, p0 - p1);
		Vector3d nd(n.x(), n.y(), n.z());
		Functs::SetVectorLength(nd, 1.0);
		Plane_3 plane(p1, Vector_3(nd[0], nd[1], nd[2]));

		//Point_2 point_to_2d(VectorPoint3d(inside), plane);

		Vector2d inside_2d = PointVector2d(point_to_2d(VectorPoint3d(inside), plane));
		Vector2d outside_2d = PointVector2d(point_to_2d(VectorPoint3d(outside), plane));

		hh_3d = hh;
		for (int i = 0; i < nb; i++)
		{
			hh_3d = hh_3d->next();
			Vector2d edge_0 = PointVector2d(point_to_2d(hh_3d->vertex()->point(), plane));
			Vector2d edge_1 = PointVector2d(point_to_2d(hh_3d->opposite()->vertex()->point(), plane));

			Vector3d edge_3d_0 = PointVector3d(hh_3d->vertex()->point());
			Vector3d edge_3d_1 = PointVector3d(hh_3d->opposite()->vertex()->point());

			Vector2d iter;

			if (CGAL_2D_Intersection_Segment_Segment(edge_0, edge_1, inside_2d, outside_2d, iter))
			{
				intersection = PointVector3d(point_to_3d(VectorPoint2d(iter), plane));
				intersection = CGAL_3D_Projection_Point_Segment(intersection, edge_3d_0, edge_3d_1);
				result_vecs.emplace_back(intersection);
				result_handles.emplace_back(hh_3d);
			}
		}

		if (result_vecs.empty())
		{
			return false;
		}
		else
		{
			double min_d = 1000000000000.0;
			for (int i = 0; i < result_vecs.size(); i++)
			{
				double dis = CGAL_3D_Distance_Point_Point(result_vecs[i], outside);
				if (min_d > dis)
				{
					min_d = dis;
					intersection = result_vecs[i];
					handle = result_handles[i];
				}
			}

			return true;
		}


	};

	//kd close query
	auto KD_Close_Query = [](kdtree* kd_tree, Vector3d query, const std::vector<Vector3d>& full_vecs,
		const std::vector<int>& full_face_id_0, const std::vector<int>& full_face_id_1, const std::vector<int>& full_face_id_2,
		const std::vector<std::vector<int>>& surface_vectices_to_face)
	{
		double* pos = new double[3];
		pos[0] = query[0];
		pos[1] = query[1];
		pos[2] = query[2];
		struct kdres* r = kd_nearest(kd_tree, pos);
		double position[3];
		int index = *(int*)kd_res_item(r, position);

		std::vector<int> faces = surface_vectices_to_face[index];
		std::vector<int> all_faces;
		int iter = 0;
		while (true)
		{
			for (int i = 0; i < faces.size(); i++) all_faces.emplace_back(faces[i]);

			std::vector<int> next_faces;
			for (int i = 0; i < faces.size(); i++)
			{
				auto vec_0 = full_face_id_0[faces[i]];
				auto vec_1 = full_face_id_1[faces[i]];
				auto vec_2 = full_face_id_2[faces[i]];
				auto d = CGAL_3D_Distance_Point_Triangle(query, full_vecs[vec_0], full_vecs[vec_1], full_vecs[vec_2]);
				if (Functs::IsAlmostZero(d))return faces[i];

				for (int j = 0; j < surface_vectices_to_face[vec_0].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_0][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_0][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_0][j]);
					}
				}

				for (int j = 0; j < surface_vectices_to_face[vec_1].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_1][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_1][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_1][j]);
					}
				}

				for (int j = 0; j < surface_vectices_to_face[vec_2].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_2][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_2][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_2][j]);
					}
				}
			}

			if (next_faces.empty()) break;

			if (iter > 5)break;


			faces = next_faces;
			iter++;
		}

		std::cerr << "KD_Close_Query" << std::endl;
		system("pause");

		return -1;
	};


	auto KD_Close_Query_0 = [](kdtree* kd_tree, Vector3d query, const std::vector<Vector3d>& full_vecs,
		const std::vector<int>& full_face_id_0, const std::vector<int>& full_face_id_1, const std::vector<int>& full_face_id_2,
		const std::vector<std::vector<int>>& surface_vectices_to_face)
	{
		double* pos = new double[3];
		pos[0] = query[0];
		pos[1] = query[1];
		pos[2] = query[2];
		struct kdres* r = kd_nearest(kd_tree, pos);
		double position[3];
		int index = *(int*)kd_res_item(r, position);

		std::vector<int> faces = surface_vectices_to_face[index];
		std::vector<int> all_faces;
		int iter = 0;
		int return_face_id = -1;
		double return_face_d = 1000000000000000.0;
		while (true)
		{
			for (int i = 0; i < faces.size(); i++) all_faces.emplace_back(faces[i]);

			std::vector<int> next_faces;
			for (int i = 0; i < faces.size(); i++)
			{
				auto vec_0 = full_face_id_0[faces[i]];
				auto vec_1 = full_face_id_1[faces[i]];
				auto vec_2 = full_face_id_2[faces[i]];
				auto d = CGAL_3D_Distance_Point_Triangle(query, full_vecs[vec_0], full_vecs[vec_1], full_vecs[vec_2]);
				//if (Functs::IsAlmostZero(d))return faces[i];

				if (d < return_face_d)
				{
					return_face_d = d;
					return_face_id = faces[i];
				}

				for (int j = 0; j < surface_vectices_to_face[vec_0].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_0][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_0][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_0][j]);
					}
				}

				for (int j = 0; j < surface_vectices_to_face[vec_1].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_1][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_1][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_1][j]);
					}
				}

				for (int j = 0; j < surface_vectices_to_face[vec_2].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_2][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_2][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_2][j]);
					}
				}
			}

			if (next_faces.empty()) break;

			if (iter > 5)break;


			faces = next_faces;
			iter++;
		}

		if (return_face_id >= 0)return return_face_id;

		std::cerr << "KD_Close_Query" << std::endl;
		system("pause");

		return -1;
	};


	auto AAA = [](const char* path, Vector3d2& multi_projects, std::vector<std::vector<Poly_facet_iterator>>& multi_project_faces,
		const std::vector<Vector3d>& full_vecs, const std::vector<int>& full_face_id_0, const std::vector<int>& full_face_id_1, const std::vector<int>& full_face_id_2)
	{
		std::ofstream  export_file_output("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\projects.obj");
		int export_index_0 = 1;
		for (int i = 0; i < multi_projects.size(); i++)
		{
			for (int j = 0; j < multi_projects[i].size(); j++)
			{
				Vector3d end_0 = multi_projects[i][j];
				Vector3d end_1 = multi_projects[i][(j + 1) % multi_projects[i].size()];
				std::string str = "projection_" + Functs::IntString(i);
				Functs::Export_Segment(export_file_output, export_index_0, str.c_str(),  end_0, end_1, 0.002);
			}
		}
		export_file_output.clear();
		export_file_output.close();

		std::vector<Vector3d> vecs;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		for (int i = 0; i < multi_project_faces.size(); i++)
		{
			for (int j = 0; j < multi_project_faces[i].size(); j++)
			{
				Poly_facet_iterator cur_face = multi_project_faces[i][j];
				int face_id = cur_face->id();
				face_id_0.emplace_back(vecs.size());
				vecs.emplace_back(full_vecs[full_face_id_0[face_id]]);
				face_id_1.emplace_back(vecs.size());
				vecs.emplace_back(full_vecs[full_face_id_1[face_id]]);
				face_id_2.emplace_back(vecs.size());
				vecs.emplace_back(full_vecs[full_face_id_2[face_id]]);

			}
		}

		Functs::OutputObj3d("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\triangles.obj", vecs, face_id_0, face_id_1, face_id_2);
	};

	auto BBBB = [](Halfedge_handle& cur_handle, Halfedge_handle& handle, Vector3d intersection,
		bool b, int iteration, int next_index, Vector3d1& cutting_points, Vector3d2& multi_projects, std::vector<std::vector<Poly_facet_iterator>>& multi_project_faces, int i, Vector3d inside)
	{

		Point_3 p0 = cur_handle->next()->next()->vertex()->point();
		Point_3 p1 = cur_handle->vertex()->point();
		Point_3 p2 = cur_handle->next()->vertex()->point();
		Point_3 p3 = cur_handle->opposite()->vertex()->point();
		Vector3d v0 = PointVector3d(p0);
		Vector3d v1 = PointVector3d(p1);
		Vector3d v2 = PointVector3d(p2);
		Vector3d v3 = PointVector3d(p3);
		Vector3d center = (v0 + v1 + v2) / (double)3.0;

		//edge.push_back(handles[i]->vertex()->id());
		//edge.push_back(handles[i]->opposite()->vertex()->id());
		//handles
		std::string file_path;

		if (b)
			file_path = "Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting\\cutting_points_" + Functs::IntString(iteration) + "_" + Functs::IntString(next_index) + "_1.obj";
		else
			file_path = "Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting\\cutting_points_" + Functs::IntString(iteration) + "_" + Functs::IntString(next_index) + "_0.obj";


		std::ofstream  export_file_output(file_path);
		int export_index = 1;
		if (false)
			for (int j = 0; j < cutting_points.size() - 1; j++)
			{
				Vector3d end_0 = cutting_points[j];
				Vector3d end_1 = cutting_points[j + 1];
				std::string str = "cutting_points_" + Functs::IntString(j);
				Functs::Export_Segment(export_file_output, export_index, str.c_str(),  end_0, end_1, 0.05);
			}
		//Functs::Export_Segment(export_file_output, export_index, "inside_outside_segment", 1.0, 0.0, 0.0, inside, multi_projects[i][next_index], 0.01);
		Functs::Export_Segment(export_file_output, export_index, "cur_tri_edge",  v1, v3, 0.001);
		std::string str = "cur_tri_center_" + std::to_string(cur_handle->face()->id());
		Functs::Export_Point(export_file_output, export_index, str.c_str(),  center, 0.002);
		Functs::Export_Point(export_file_output, export_index, "cur_tri_inside",  inside, 0.002);
		str = "cur_tri_outside_" + std::to_string(multi_project_faces[i][next_index]->id());
		Functs::Export_Point(export_file_output, export_index, str.c_str(), multi_projects[i][next_index], 0.002);

		if (b)
		{
			Vector3d vv0 = PointVector3d(handle->vertex()->point());
			Vector3d vv1 = PointVector3d(handle->opposite()->vertex()->point());
			Functs::Export_Segment(export_file_output, export_index, "next_edge",  vv0, vv1, 0.001);

			Functs::Export_Point(export_file_output, export_index, "inters_point",  intersection, 0.002);
		}

		export_file_output.clear();
		export_file_output.close();
	};

	auto CCCC = [](Vector3d1& cutting_points)
	{
		std::ofstream  export_file_output("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting_points.obj");
		int export_index = 1;
		for (int i = 0; i < cutting_points.size() - 1; i++)
		{
			std::string str = "intersection_" + Functs::IntString(i);
			Functs::Export_Segment(export_file_output, export_index, str.c_str(),  cutting_points[i], cutting_points[i + 1], 0.05);
		}
		export_file_output.clear();
		export_file_output.close();
	};

	auto DDDD = [](const char* path, std::vector<std::pair<int, int>>& edges,
		std::vector<Halfedge_handle>& handles, Vector3d1& cutting_points, Vector3d1& full_vecs)
	{
		std::ofstream  export_file_output_0(path);
		int export_index_0 = 1;

		for (int i = 0; i < handles.size(); i++)
		{
			auto name = "edge_" + std::to_string(i) + "_" + std::to_string(handles[i]->face()->id());
			Functs::Export_Segment(export_file_output_0, export_index_0, name.c_str(),  full_vecs[edges[i].first], full_vecs[edges[i].second], 0.03);
			Functs::Export_Point(export_file_output_0, export_index_0, name.c_str(),  cutting_points[i], 0.05);

			Point_3 p0 = handles[i]->next()->next()->vertex()->point();
			Point_3 p1 = handles[i]->vertex()->point();
			Point_3 p2 = handles[i]->next()->vertex()->point();
			Vector3d v0 = PointVector3d(p0);
			Vector3d v1 = PointVector3d(p1);
			Vector3d v2 = PointVector3d(p2);
			Vector3d center = (v0 + v1 + v2) / (double)3.0;
			Functs::Export_Point(export_file_output_0, export_index_0, name.c_str(),  center, 0.05);

		}

		export_file_output_0.clear();
		export_file_output_0.close();
	};

	auto EEEE = [](Vector3d inside)
	{
		std::ofstream  export_file_output_0("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\inside.obj");
		int export_index_0 = 1;
		Functs::Export_Point(export_file_output_0, export_index_0, "start_point",  inside, 0.05);
		export_file_output_0.clear();
		export_file_output_0.close();
	};

	auto FFFF = [](const char* path, Vector3d1 insides)
	{
		std::ofstream  export_file_output_0(path);
		int export_index_0 = 1;
		for (auto inside : insides)
			Functs::Export_Point(export_file_output_0, export_index_0, "start_point",  inside, 0.05);
		export_file_output_0.clear();
		export_file_output_0.close();
	};


	//build polyhedron from the full_path obj
	Polyhedron_3 polyhedron;
	Vector3d1 full_vecs;
	std::vector<int> full_face_id_0;
	std::vector<int> full_face_id_1;
	std::vector<int> full_face_id_2;
	Construct_Polyhedron(polyhedron, full_path, full_vecs, full_face_id_0, full_face_id_1, full_face_id_2);

	if (full_vecs.empty() || full_face_id_0.empty() || full_face_id_1.empty() || full_face_id_2.empty())
	{
		std::cerr << "if (full_vecs.empty() || full_face_id_0.empty() || full_face_id_1.empty() || full_face_id_2.empty())" << std::endl;
		system("pause");
	}

	//surface normals
	//full_face_iters
	Vector3d1 surface_normals;
	std::vector<Poly_facet_iterator> full_face_iters;
	for (Poly_facet_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
	{
		full_face_iters.emplace_back(iter);

		Point_3 p0 = iter->halfedge()->next()->next()->vertex()->point();
		Point_3 p1 = iter->halfedge()->vertex()->point();
		Point_3 p2 = iter->halfedge()->next()->vertex()->point();
		Vector_3 n = CGAL::cross_product(p2 - p1, p0 - p1);
		Vector3d nd(n.x(), n.y(), n.z());
		Functs::SetVectorLength(nd, 1.0);
		surface_normals.emplace_back(nd);
	}

	//surface_vectices_to_face and surface_vectices_to_vectices
	std::vector<std::vector<int>> surface_vectices_to_face;
	CGAL_3D_Triangle_Mesh_Vecs_Faces(full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, surface_vectices_to_face);
	std::vector<std::vector<int>> surface_vectices_to_vectices;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, surface_vectices_to_vectices);

	//kd tree
	std::vector<int> full_vec_ids(full_vecs.size(), 0);
	for (int i = 0; i < full_vecs.size(); i++) full_vec_ids[i] = i;
	kdtree* kd_tree = kd_create(3);
	for (int i = 0; i < full_vecs.size(); i++)
	{
		void* val = &full_vec_ids[i];
		kd_insert3(kd_tree, full_vecs[i][0], full_vecs[i][1], full_vecs[i][2], val);
	}

	//project points along input boundaries to the full_path obj
	Vector3d2 multi_projects;// projecting points
	std::vector<std::vector<Poly_facet_iterator>> multi_project_faces;//related faces of projecting points
	for (int i = 0; i < multi_boundary.size(); i++)
	{
		Vector3d1 one_projects;
		std::vector<Poly_facet_iterator> one_project_faces;
		for (int j = 0; j < multi_boundary[i].size(); j++)
		{
			std::cerr << i << " " << j << std::endl;
			Point_3 query(multi_boundary[i][j][0], multi_boundary[i][j][1], multi_boundary[i][j][2]);
			auto fid = KD_Close_Query_0(kd_tree, multi_boundary[i][j], full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, surface_vectices_to_face);
			one_projects.push_back(multi_boundary[i][j]);
			one_project_faces.push_back(full_face_iters[fid]);
		}
		multi_projects.push_back(one_projects);
		multi_project_faces.push_back(one_project_faces);
	}
	//output the projecting points
	if (DebugInformation())
		AAA("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\projects.obj", multi_projects, multi_project_faces, full_vecs, full_face_id_0, full_face_id_1, full_face_id_2);

	Functs::ClearFolder("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting\\");

	//searching for all of the cutting points on edges
	Vector3d1 cutting_points;
	std::vector<int> multi_cps;
	std::vector<Halfedge_handle> handles;
	std::vector<bool> face_used(full_face_id_0.size(), false);
	for (int i = 0; i < multi_projects.size(); i++)
	{
		std::cerr << "multi_projects: " << i << " / " << multi_projects.size() << std::endl;

		int cur_face_id = multi_project_faces[i][0]->id();
		face_used[cur_face_id] = true;
		Poly_facet_iterator cur_face = multi_project_faces[i][0];
		Halfedge_handle cur_handle = cur_face->halfedge();
		Vector3d inside = multi_projects[i][0];

		int next_index = 1;
		int iteration = 0;
		while (true)
		{
			if (multi_projects.size() == 1) std::cout << iteration << " / " << multi_projects[i].size() << " " << next_index << std::endl;

			//searching for the outside point of the current triangle
			bool goon = false;
			while (true)
			{
				if (cur_face_id == multi_project_faces[i][next_index]->id())
				{
					inside = multi_projects[i][next_index];
					if (next_index == 0) break;
					next_index = (next_index + 1) % multi_projects[i].size();
				}
				else
				{
					auto cur_face_normal = surface_normals[cur_face_id];
					auto next_face_normal = surface_normals[multi_project_faces[i][next_index]->id()];
					double angle = Functs::GetAngleBetween(cur_face_normal, next_face_normal);

					if (angle < Math_PI / 2.0 && OutsidePointInsideTriangle(cur_face, multi_projects[i][next_index]))
					{
						inside = multi_projects[i][next_index];
						if (next_index == 0) break;
						next_index = (next_index + 1) % multi_projects[i].size();
					}
					else
					{
						//if (face_used[multi_project_faces[i][next_index]->id()])
						//{
						//	inside = multi_projects[i][next_index];
						//	if (next_index == 0) break;
						//	next_index = (next_index + 1) % multi_projects[i].size();
						//}
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
				b = Intersection(cur_handle, 3, inside, multi_projects[i][next_index], handle, intersection);
			else
				b = Intersection(cur_handle, 2, inside, multi_projects[i][next_index], handle, intersection);

			//output the cutting points of each iteration
			if (true)
				BBBB(cur_handle, handle, intersection, b, iteration, next_index, cutting_points, multi_projects, multi_project_faces, i, inside);

			if (b)
			{
				cutting_points.push_back(intersection);
				handles.push_back(handle);

				//move next step
				inside = intersection;
				cur_handle = handle->opposite();
				cur_face_id = (int)cur_handle->face()->id();
				face_used[cur_face_id] = true;
				cur_face = cur_handle->face();

				if (cur_face_id == multi_project_faces[i][0]->id())
				{
					break;
				}
			}
			else
			{
				std::cerr << "if (b)" << std::endl;
				system("pause");

				cutting_points.erase(cutting_points.begin() + cutting_points.size() - 1);
				handles.erase(handles.begin() + handles.size() - 1);

				next_index = (next_index + 1) % multi_projects[i].size();
				if (next_index == 0)
					break;

				inside = cutting_points[cutting_points.size() - 1];
				cur_handle = handles[handles.size() - 1]->opposite();
				cur_face_id = (int)cur_handle->face()->id();
				cur_face = cur_handle->face();
			}
			iteration++;
		}

		multi_cps.emplace_back(cutting_points.size());
	}

	//output the cutting points
	if (false)
	{
		CCCC(cutting_points);
	}

	/*******************************************/
	std::vector<std::pair<int, int>> edges;
	/*******************************************/
	for (int i = 0; i < handles.size(); i++)
		edges.emplace_back(handles[i]->vertex()->id(), handles[i]->opposite()->vertex()->id());

	//smoothing
	//for (int i = 0; i < 3; i++) OneIterationSmoothBoundary(full_vecs, edges, cutting_points);

	auto fid = KD_Close_Query(kd_tree, inside_point, full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, surface_vectices_to_face);

	if (true)
	{
		DDDD("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\edges.obj", edges, handles, cutting_points, full_vecs);
		EEEE(PointVector3d(full_face_iters[fid]->halfedge()->vertex()->point()));
	}

	//lables
	std::vector<int> two_class_lables;
	ComputeEdgeLables((int)full_vecs.size(), full_face_iters[fid]->halfedge(), edges, two_class_lables);

	Vector3d1 lable_0_vecs;
	Vector3d1 lable_1_vecs;
	for (int i = 0; i < two_class_lables.size(); i++)
	{
		if (two_class_lables[i] == 0) lable_0_vecs.emplace_back(full_vecs[i]);
		if (two_class_lables[i] == 1) lable_1_vecs.emplace_back(full_vecs[i]);
	}
	std::cerr << lable_0_vecs.size() << " " << lable_1_vecs.size() << std::endl;

	FFFF("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\lable_0_vecs.obj", lable_0_vecs);
	FFFF("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\lable_1_vecs.obj", lable_1_vecs);

	//remesh triangles
	ComputeRemeshTriangles(full_vecs, full_face_id_0, full_face_id_1, full_face_id_2,
		edges, cutting_points, multi_cps, two_class_lables, output_path);

	std::cout << "over" << std::endl;
	return;
}


//extern "C" PPGL_EXPORT void CGAL_3D_Mesh_Gradient(const Vector3d1 & vecs, const std::vector<int>&face_id_0, const std::vector<int>&face_id_1, const std::vector<int>&face_id_2, const std::vector<double>&psd, Vector3d1 & vecs_gradients, Vector3d1 & faces_gradients);
