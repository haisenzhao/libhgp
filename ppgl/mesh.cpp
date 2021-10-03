#include "geom.h"

//extern "C" PPGL_EXPORT
void Construct_Polyhedron(Polyhedron_3 &polyhedron,
                          const Vector3d1 &vecs,
                          const Vector1i1 &face_id_0, const Vector1i1 &face_id_1,
                          const Vector1i1 &face_id_2) {
    std::vector<double> coords;
   Vector1i1 tris;
    for (int i = 0; i < vecs.size(); i++) {
        coords.push_back(vecs[i][0]);
        coords.push_back(vecs[i][1]);
        coords.push_back(vecs[i][2]);
    }
    for (int i = 0; i < face_id_0.size(); i++) {
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
   Vector1i1().swap(tris);
}


Poly_point_3 Minus(Poly_point_3 a, Poly_point_3 b) {
    return Poly_point_3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

float Dot(Poly_point_3 a, Poly_point_3 b) {
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
    int cur_face_id = project_faces[0]->id();
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
                next_index = (next_index + 1) % feature.size();
                if (next_index == 0) break;
            } else {
                if (OutsidePointInsideTriangle(cur_face, feature[next_index])) {
                    next_index = (next_index + 1) % feature.size();
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
            cur_face_id = cur_handle->face()->id();
            cur_face = cur_handle->face();

            if (cur_face_id == project_faces[0]->id()) {
                break;
            }
        } else {
            cutting_points.erase(cutting_points.begin() + cutting_points.size() - 1);
            handles.erase(handles.begin() + handles.size() - 1);

            next_index = (next_index + 1) % feature.size();
            if (next_index == 0) break;

            inside = cutting_points[cutting_points.size() - 1];
            cur_handle = handles[handles.size() - 1]->opposite();
            cur_face_id = cur_handle->face()->id();
            cur_face = cur_handle->face();
        }
        iteration++;
    }

    for (int i = 0; i < handles.size(); i++) {
        cutting_faces[handles[i]->face()->id()].push_back(igl_cutting_0_edges.size());
        cutting_faces[handles[i]->opposite()->face()->id()].push_back(igl_cutting_0_edges.size());
        igl_cutting_0_edges.push_back(handles[i]->vertex()->id());
        igl_cutting_1_edges.push_back(handles[i]->opposite()->vertex()->id());
    }
    igl_cutting_points = cutting_points;
}

