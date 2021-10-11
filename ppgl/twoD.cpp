#include "geom.h"
#include "clipper/clipper.hpp"

extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Point(const Vector2d& p_0, const Vector2d& p_1)
{
    return sqrt(pow((p_0[0] - p_1[0]), 2.0) + pow((p_0[1] - p_1[1]), 2.0));
}

extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Segment(const Vector2d& v, const Vector2d & s_0, const Vector2d & s_1) {
    return sqrt((double) CGAL::squared_distance(Point_2(v[0], v[1]),
                                                Segment_2(Point_2(s_0[0], s_0[1]), Point_2(s_1[0], s_1[1]))));
}

extern "C" PPGL_EXPORT double
CGAL_2D_Distance_Segment_Segment(const Vector2d& s_0, const Vector2d & s_1, const Vector2d & e_0, const Vector2d & e_1) {
    double d0 = CGAL_2D_Distance_Point_Segment(s_0, e_0, e_1);
    double d1 = CGAL_2D_Distance_Point_Segment(s_1, e_0, e_1);
    double d2 = CGAL_2D_Distance_Point_Segment(e_0, s_0, s_1);
    double d3 = CGAL_2D_Distance_Point_Segment(e_1, s_0, s_1);
    double min_d = d0;
    min_d = std::min(min_d, d1);
    min_d = std::min(min_d, d2);
    min_d = std::min(min_d, d3);
    return min_d;
}

extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Line(const Vector2d& v, const Vector2d & l_0, const Vector2d & l_1) {
    return sqrt((double) CGAL::squared_distance(Point_2(v[0], v[1]),
                                                Line_2(Point_2(l_0[0], l_0[1]), Point_2(l_1[0], l_1[1]))));
}

extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Polygon(const Vector2d& p, const Vector2d1 & py) {
    double distance = 1000000000000.0;
    for (int i = 0; i < py.size(); i++)
    {
        int ii = (static_cast<int>(i) + 1) % py.size();
		distance = std::min(distance, CGAL_2D_Distance_Point_Segment(p, py[i], py[ii]));
    }
    return distance;
}

extern "C" PPGL_EXPORT bool CGAL_2D_Is_Point_OutCGALPolygon(const Vector2d &p, const Polygon_2 &py) {
    return py.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE;
}

extern "C" PPGL_EXPORT bool CGAL_Construct_Polygon(const Vector2d1 &py, Polygon_2 &poly) {
    poly.clear();
    for (auto i : py)
        poly.push_back(Point_2(i[0], i[1]));
    return poly.is_simple();
}

extern "C" PPGL_EXPORT bool
CGAL_Construct_InOutSide_Polygon(const Vector2d1 &py, const Vector2d &p, const Vector2d &q, bool &isPInside,
                                 bool &isQInside) {
    Polygon_2 poly;
    poly.clear();
    for (auto i : py)
        poly.push_back(Point_2(i[0], i[1]));
    if (!poly.is_simple()) return false;

    isPInside = CGAL_2D_Is_Point_OutCGALPolygon(p, poly);
    isQInside = CGAL_2D_Is_Point_OutCGALPolygon(q, poly);
    return true;
}

extern "C" PPGL_EXPORT bool CGAL_2D_Location_Point_Polygon(const Vector2d & p, const Vector2d1 & py) {
    Polygon_2 poly;
    for (int i = 0; i < py.size(); i++)
        poly.push_back(Point_2(py[i][0], py[i][1]));

    return poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_BOUNDED_SIDE;
}

extern "C" PPGL_EXPORT bool CGAL_2D_Is_Point_OutPolygon(Vector2d p, Vector2d1 py) {
    Polygon_2 poly;
    for (int i = 0; i < py.size(); i++)
        poly.push_back(Point_2(py[i][0], py[i][1]));

    return poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE;
}


extern "C" PPGL_EXPORT bool CGAL_2D_Location_Points_Polygon(const Vector2d1 &ps,
                                                                      const Vector2d1 &py) {
    Polygon_2 poly;
    for (int i = 0; i < py.size(); i++)
        poly.push_back(Point_2(py[i][0], py[i][1]));

    for (auto p : ps) {
        if (poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE)
            return false;
    }
    return true;
}

extern "C" PPGL_EXPORT void CGAL_2D_Polygon_Dart_Sampling(const Vector2d1 & py, const double& d, Vector2d1 & sampling_points, const int& total_iter)
{
	if (!(d > 0 && d < 1.0))
		Functs::MAssert("CGAL_3D_Mesh_Dart_Sampling_C1 if (!(d > 0 && d < 1.0))");

	Polygon_2 poly;
	for (int i = 0; i < py.size(); i++)
		poly.push_back(Point_2(py[i][0], py[i][1]));

	double xmin = poly.bbox().xmin();
	double ymin = poly.bbox().ymin();
	double xmax = poly.bbox().xmax();
	double ymax = poly.bbox().ymax();

	double diagonal_length = CGAL_2D_Distance_Point_Point(Vector2d(xmin,ymin), Vector2d(xmax,ymax));
	double minimal_d = d * diagonal_length;

	int run = 0;
	Vector2d1 insert_points;

	while (run < total_iter)
	{
		run++;
		double x = rand() / double(RAND_MAX);
		double y = rand() / double(RAND_MAX);
		x = (xmax - xmin) * x + xmin;
		y = (ymax - ymin) * y + ymin;

		if (poly.bounded_side(Point_2(x, y)) == CGAL::ON_BOUNDED_SIDE)
		{
			double distance = CGAL_IA_MAX_DOUBLE;
			for (int i = 0; i < insert_points.size(); i++)
				distance = std::min(distance, CGAL_2D_Distance_Point_Point(insert_points[i], Vector2d(x, y)));

			if (distance > minimal_d)
			{
				insert_points.push_back(Vector2d(x, y));
				run = 0;
			}
		}
	}

	for (int i = 0; i < insert_points.size(); i++)
	{
		double distance = CGAL_IA_MAX_DOUBLE;
		for (int j = 0; j < py.size(); j++)
			distance = std::min(distance, CGAL_2D_Distance_Point_Point(py[j], insert_points[i]));
		if (distance > d / 2.0)
			sampling_points.push_back(insert_points[i]);
	}
}

extern "C" bool CGAL_2D_Intersection_Segment_Segment
        (const Vector2d & s_0_s, const Vector2d & s_0_e, const Vector2d & s_1_s, const Vector2d & s_1_e, Vector2d &inter) {
    CGAL::Object result = intersection(Segment_2(Point_2(s_0_s[0], s_0_s[1]), Point_2(s_0_e[0], s_0_e[1])),
                                       Segment_2(Point_2(s_1_s[0], s_1_s[1]), Point_2(s_1_e[0], s_1_e[1])));

    if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        inter[0] = ipoint->x();
        inter[1] = ipoint->y();
        return true;
    }
    return false;
}

extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Ray_Segment
        (const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s, const Vector2d &s_1_e, Vector2d &inter) {
    Point_2 st(s_0_s[0], s_0_s[1]);

    /*std::cerr << "CGAL_2D_Intersection_Ray_Segment " << st << std::endl;
    std::cerr << "CGAL_2D_Intersection_Ray_Segment " << s_0_e[0] << " " << s_0_e[1] << std::endl;
    std::cerr << "Ray_Segment " << s_1_s[0] << " " << s_1_s[1] << std::endl;
    std::cerr << "Ray_Segment " << s_1_e[0] << " " << s_1_e[1] << std::endl;*/

    Vector_2 dir(s_0_e[0], s_0_e[1]);
    dir = dir / CGAL::sqrt(dir * dir);
    Point_2 nd = st + 100000.0 * dir;
    CGAL::Object result = intersection(Segment_2(st, nd),
                                       Segment_2(Point_2(s_1_s[0], s_1_s[1]), Point_2(s_1_e[0], s_1_e[1])));

    if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        //if (!Math::IsAlmostZero((st - *ipoint).squared_length()))
        {
            inter[0] = ipoint->x();
            inter[1] = ipoint->y();
            //			std::cerr << "has intersection" << std::endl;
            return true;
        }
        //else return false;
    }
    return false;

    //Point_2 st(s_0_s[0]+ s_0_e[0]*10000.0, s_0_s[1] + s_0_e[1] * 10000.0);
    /*
    Ray_2 ray(Point_2(s_0_s[0], s_0_s[1]), Vector_2(s_0_e[0], s_0_e[1]));
    Segment_2 seg(Point_2(s_1_s[0], s_1_s[1]), Point_2(s_1_e[0], s_1_e[1]));
    auto result = CGAL::intersection(ray, seg);

    if (result) std::cerr << "not empty" << std::endl;
    if (const Point_2 * ipoint = boost::get<Point_2>(&*result))
    {
        std::cerr << *ipoint << std::endl;
        inter[0] = ipoint->x();
        inter[1] = ipoint->y();
        std::cerr << "has intersection" << std::endl;
        return true;
    }
    else
    {
        std::cerr << "no intersection" << std::endl;

        return false;
    }*/
}

extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Ray_Segment_Ignore_Endpoints
        (const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s, const Vector2d &s_1_e, Vector2d &inter) {
    Point_2 st(s_0_s[0], s_0_s[1]);
    Vector_2 dir(s_0_e[0], s_0_e[1]);
    dir = dir / CGAL::sqrt(dir * dir);
    Point_2 nd = st + 100000.0 * dir;
    CGAL::Object result = intersection(Segment_2(st, nd),
                                       Segment_2(Point_2(s_1_s[0], s_1_s[1]), Point_2(s_1_e[0], s_1_e[1])));

    if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        // TODO: Major change: if intersected point is one of the end-points of segment, then ignore it
        // TODO: Finished on 08-12-2019
        if (Functs::IsAlmostZero((Point_2(s_1_e[0], s_1_e[1]) - *ipoint).squared_length()) ||
            Functs::IsAlmostZero((Point_2(s_1_e[0], s_1_e[1]) - *ipoint).squared_length())) {
            return false;
        }

        inter[0] = ipoint->x();
        inter[1] = ipoint->y();
        return true;
    }
    return false;
}


extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Ray_Polygon(
        const Vector2d &r_s,
        const Vector2d &r_d,
        const Vector2d1 &poly,
        Vector2d &pnt) {
    int nIntersect = 0;
    for (int i = 0; i < poly.size(); i++) {
        Vector2d inter;
        int ii = (static_cast<int>(i) + 1) % poly.size();
        if (CGAL_2D_Intersection_Ray_Segment(r_s, r_d, poly[i], poly[ii], inter)) {
            ++nIntersect;
            return true;
        }
    }
    return false;
}


extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Line_Line
        (const Vector2d &s_0_s, const Vector2d &s_0_e, const Vector2d &s_1_s, const Vector2d &s_1_e, Vector2d &inter) {
    CGAL::Object result = intersection(Line_2(Point_2(s_0_s[0], s_0_s[1]), Point_2(s_0_e[0], s_0_e[1])),
                                       Line_2(Point_2(s_1_s[0], s_1_s[1]), Point_2(s_1_e[0], s_1_e[1])));

    if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        inter[0] = ipoint->x();
        inter[1] = ipoint->y();
        return true;
    }
    return false;
}

extern "C" PPGL_EXPORT bool
CGAL_2D_Intersection_Segment_Polygon(const Vector2d & s_s, const Vector2d & s_e, Vector2d1 &p) {
    for (int i = 0; i < p.size(); i++) {
        Vector2d inter;
        int ii = (static_cast<int>(i) + 1) % p.size();
        if (CGAL_2D_Intersection_Segment_Segment(s_s, s_e, p[i], p[ii], inter)) {
            return true;
        }
    }
    return false;
}

extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Is_Clockwise_Oriented(const Vector2d1 &ps) {
    Polygon_2 poly;
    for (int i = 0; i < ps.size(); i++)
        poly.push_back(Point_2(ps[i][0], ps[i][1]));

    return poly.is_clockwise_oriented();
}

extern "C" PPGL_EXPORT double CGAL_2D_Two_Polygons_Intersection(const Vector2d1 &poly_0,
                                                                          const Vector2d1 &poly_1) {
    double scale = 1000000.0;

    ClipperLib::Paths subj(1);
    for (int i = 0; i < poly_0.size(); i++)
        subj[0] << ClipperLib::IntPoint((ClipperLib::cInt)(poly_0[i][0] * scale), (ClipperLib::cInt)(poly_0[i][1] * scale));

    ClipperLib::Paths cliper(1);
    for (int i = 0; i < poly_1.size(); i++)
        cliper[0] << ClipperLib::IntPoint((ClipperLib::cInt)(poly_1[i][0] * scale), (ClipperLib::cInt)(poly_1[i][1] * scale));

    ClipperLib::Paths solution;
    ClipperLib::Clipper c;
    c.AddPaths(subj, ClipperLib::ptSubject, true);
    c.AddPaths(cliper, ClipperLib::ptClip, true);
    c.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

    double area = 0.0;

    for (int i = 0; i < solution.size(); i++) {
        Polygon_2 poly_2;
        for (int j = 0; j < solution[i].size(); j++) {
            poly_2.push_back(Point_2(((double) solution[i][j].X) / scale, ((double) solution[i][j].Y) / scale));
        }
        area += poly_2.area();
    }

    return area;
}

extern "C" PPGL_EXPORT double
CGAL_2D_Two_Polygons_Union(const Vector2d1& poly_0, const Vector2d1 & poly_1, Vector2d2  &inter_polygons) {
    double scale = 1000000.0;

    ClipperLib::Paths subj(1);
    for (int i = 0; i < poly_0.size(); i++)
        subj[0] << ClipperLib::IntPoint((ClipperLib::cInt)(poly_0[i][0] * scale), (ClipperLib::cInt)(poly_0[i][1] * scale));

    ClipperLib::Paths cliper(1);
    for (int i = 0; i < poly_1.size(); i++)
        cliper[0] << ClipperLib::IntPoint((ClipperLib::cInt)(poly_1[i][0] * scale), (ClipperLib::cInt)(poly_1[i][1] * scale));

    ClipperLib::Paths solution;
    ClipperLib::Clipper c;
    c.AddPaths(subj, ClipperLib::ptSubject, true);
    c.AddPaths(cliper, ClipperLib::ptClip, true);
    c.Execute(ClipperLib::ctUnion, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

    double area = 0.0;

    for (int i = 0; i < solution.size(); i++) {
        Polygon_2 poly_2;

        std::vector<double> xs;
        std::vector<double> ys;

        Vector2d1 polygon;

        for (int j = 0; j < solution[i].size(); j++) {
            poly_2.push_back(Point_2(((double) solution[i][j].X) / scale, ((double) solution[i][j].Y) / scale));
            polygon.push_back(Vector2d(((double) solution[i][j].X) / scale, ((double) solution[i][j].Y) / scale));
        }

        if (poly_2.area() > 0.0) {
            inter_polygons.push_back(polygon);
        }

        area += poly_2.area();
    }
    return area;
}


static void RemoveClosePoints(Vector2d1 &poly) {
   Vector1i1 remove_int;

    if (poly.size() > 2) {
        for (int i = 0; i < poly.size() - 1; i++) {
            int ii = (static_cast<int>(i) + 1) % poly.size();
            double d = CGAL_2D_Distance_Point_Point(poly[i], poly[ii]);
            if (d < 0.00001) remove_int.push_back(ii);
        }

        for (int i = (int)remove_int.size() - 1; i >= 0; i--) {
            poly.erase(poly.begin() + remove_int[i]);
        }
    }
}

extern "C" PPGL_EXPORT void CGAL_2D_Polygon_One_Offsets(const Vector2d1 &poly, const double& d, Vector2d2  &offset_polys) {
    if (!(poly.size() > 0)) return;

    double scale = 1000000.0;

    ClipperLib::ClipperOffset co;

    co.MiterLimit = 1000.0;


    co.ArcTolerance = co.ArcTolerance * scale / 1000.0;

    ClipperLib::Path subj;
    ClipperLib::Paths solution;

    //build the most outside path
    for (int i = 0; i < poly.size(); i++)
        subj << ClipperLib::IntPoint((ClipperLib::cInt)(poly[i].x * scale), (ClipperLib::cInt)(poly[i].y * scale));
    co.AddPath(subj, ClipperLib::jtMiter, ClipperLib::etClosedPolygon);

    // execute
    co.Execute(solution, -d * scale);

    //output
    for (int i = 0; i < solution.size(); i++) {
        Vector2d1 one_offset;

        Polygon_2 poly_2;
        for (int j = 0; j < solution[i].size(); j++) {
            double x = ((double) solution[i][j].X) / scale;
            double y = ((double) solution[i][j].Y) / scale;
            one_offset.push_back(Vector2d(x, y));
            poly_2.push_back(Point_2(x, y));
        }

        if (poly_2.is_clockwise_oriented()) {
            std::reverse(one_offset.begin(), one_offset.end());
        }

        //remove closed points
        RemoveClosePoints(one_offset);
        offset_polys.push_back(one_offset);
        Vector2d1().swap(one_offset);
    }
}

extern "C" PPGL_EXPORT void
CGAL_Decompose_Polyline(const Vector2d1 &polyline, const double& threshold, Vector1i1& result) {
    for (auto &p : polyline) 
    {

    }
}

// This one is used to intersect a polygon with a segment
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_NotExtend(
        const Vector2d1 &polygon,
        const Vector2d &s, const Vector2d &e) {
    const Vector2d dir = normalize(e - s);
    const int polySize = (int)polygon.size();

    for (int i = 0; i < polySize; i++) {
        Vector2d inter1;
        const auto &pop1 = polygon[i];
        const auto &pop2 = polygon[(i + 1) % polySize];
        const Vector2d segDir = normalize(pop1 - pop2);
        const Vector2d conDir = normalize(pop1 - s);
        //std::cerr << "Va = " << glm::abs(dot(dir, segDir)) << " Vb = " << glm::abs(dot(dir, conDir)) << std::endl;

        // 1. If a segment is parallel to an edge, then igonre
        if (glm::abs(dot(dir, segDir)) > 0.9999 &&
            glm::abs(dot(dir, conDir)) > 0.9999) {
            continue;
        }
        // 2. if segment point is on the end-points of segment, then ignore it as well
        if (CGAL_2D_Intersection_Ray_Segment_Ignore_Endpoints(s, dir, pop1, pop2, inter1)) {
            //const Vector2d intDir(inter1 - s);
            //if (glm::abs(glm::dot(intDir, dir)) < 0.99)
            if (!Functs::IsAlmostZero(length(inter1 - s)) &&
                !Functs::IsAlmostZero(length(inter1 - pop1)) &&
                !Functs::IsAlmostZero(length(inter1 - pop2))) {
                //std::cerr << "intersected at " << inter1[0] << " " << inter1[1] << std::endl;
                return false;
            }
        }
    }

    return true;
}

void OutputRectangle(std::string path, const Vector2d2  &points) {
    std::ofstream file(path);

    for (int i = 0; i < points.size(); i++) {
        for (int j = 0; j < points[i].size(); j++) {
            file << "v " << points[i][j][0] << " " << points[i][j][1] << " " << 0 << std::endl;
        }
    }

    int nb = 1;

    for (int i = 0; i < points.size(); i++) {
        file << "f ";
        for (int j = 0; j < points[i].size(); j++) {
            file << Functs::IntString(nb) << " ";
            nb++;
        }
        file << "" << std::endl;
    }


    file.clear();
    file.close();
}

extern "C" PPGL_EXPORT double CGAL_Get_Angle_Kerf_Offset_Tan(const Vector2d &a, const Vector2d &b) {
    auto na = normalize(a);
    auto nb = normalize(b);
    return glm::tan(glm::acos(glm::abs(dot(na, nb))));
}

// This one is used to intersect a polygon with a line
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_Extend(
        const Vector2d1 &polygon,
        const Vector2d &s, const Vector2d &e,
        Vector2d &ns, Vector2d &ne) {
    const double eps = 0.1;
    //std::cerr << "input s = " << s[0] << " " << s[1] << std::endl;
    //std::cerr << "input e = " << e[0] << " " << e[1] << std::endl;
    const Vector2d cutDir = normalize(e - s);

    //std::cerr << "input cut = " << cutDir[0] << " " << cutDir[1] << std::endl;
    ns = s, ne = e;
    const int polySize = (int)polygon.size();
    Vector2d1 rayD1Int, rayD2Int;

    const auto pd1 = CGAL_2D_Distance_Point_Polygon(s, polygon);
    const auto pd2 = CGAL_2D_Distance_Point_Polygon(e, polygon);

    // if the points are in the outside of polygon, then we intersect them with polygon
    Polygon_2 cgalPoly;

    if (!CGAL_Construct_Polygon(polygon, cgalPoly)) {
        std::cerr << "Fatal error: CGAL_Construct_Polygon" << std::endl;
        return false;
    }


    bool isoutside1 = false, isoutside2 = false;
    isoutside1 = CGAL_2D_Is_Point_OutCGALPolygon(s, cgalPoly);
    isoutside2 = CGAL_2D_Is_Point_OutCGALPolygon(e, cgalPoly);
    //std::cerr << "inside1 = " << inside1 << " inside2 = " << inside2 << std::endl;

    if (Functs::IsAlmostZero(pd1)) {
    } else if (isoutside1) {
        Vector2d1 raySnap;
        for (int i = 0; i < polySize; i++) {
            Vector2d inter;
            const auto &pop1 = polygon[i];
            const auto &pop2 = polygon[(i + 1) % polySize];
            const Vector2d segDir = normalize(pop1 - pop2);
            const Vector2d conDir = normalize(pop1 - s);
            if (glm::abs(dot(cutDir, segDir)) > 0.9999 &&
                glm::abs(dot(cutDir, conDir)) > 0.9999) {
                raySnap.push_back(pop1);
                raySnap.push_back(pop2);
            }
            if (CGAL_2D_Intersection_Ray_Segment(s - eps * cutDir, cutDir, pop1, pop2, inter)) {
                auto tanAngle = CGAL_Get_Angle_Kerf_Offset_Tan(segDir, cutDir);
                inter -= (1.5875 / tanAngle * cutDir);
                //std::cerr << "ray snap 1 angle" << tanAngle << std::endl;
                //std::cerr << "inter = " << inter[0] << " " << inter[1] << std::endl;
                raySnap.push_back(inter);
            }
        }
        //std::cerr << "raysnap1 = " << raySnap.size() << std::endl;
        if (raySnap.empty())
            return false;

        if (raySnap.size() != 1) {
            double minDist = DBL_MAX;
            for (const auto &v : raySnap) {
                const double tmpDist2 = length2(v - s);
                if (tmpDist2 < minDist) {
                    ns = v;
                    minDist = tmpDist2;
                }
            }
        } else {
            std::cerr << "raysnap1 % 2 != 0" << " size = " << raySnap.size() << std::endl;
        }
    } else {
        for (int i = 0; i < polySize; i++) {
            Vector2d inter;
            const auto &pop1 = polygon[i];
            const auto &pop2 = polygon[(i + 1) % polySize];
            Vector2d segDir = normalize(pop1 - pop2);

            if (CGAL_2D_Intersection_Ray_Segment(s + eps * cutDir, -cutDir, pop1, pop2, inter)) {
                auto tanAngle = CGAL_Get_Angle_Kerf_Offset_Tan(segDir, cutDir);
                //std::cerr << "rayd1int" << tanAngle << std::endl;
                inter -= (1.5875 / tanAngle * cutDir);
                rayD1Int.push_back(inter);
            }
        }
        //std::cerr << "rayd1int = " << rayD1Int.size() << std::endl;
        if (rayD1Int.empty()) {
            return false;
        }

        double minDist = DBL_MAX;
        for (const auto &v : rayD1Int) {
            const double tmpDist2 = length2(v - s);
            if (tmpDist2 < minDist) {
                ns = v;
                minDist = tmpDist2;
            }
        }
    }

    if (Functs::IsAlmostZero(pd2)) {
    } else if (isoutside2) {
        Vector2d1 raySnap;
        Vector2d inter;
        for (int i = 0; i < polySize; i++) {
            const auto &pop1 = polygon[i];
            const auto &pop2 = polygon[(i + 1) % polySize];
            const Vector2d segDir = normalize(pop1 - pop2);
            const Vector2d conDir = normalize(pop1 - s);
            if (glm::abs(dot(cutDir, segDir)) > 0.9999 &&
                glm::abs(dot(cutDir, conDir)) > 0.9999) {
                raySnap.push_back(pop1);
                raySnap.push_back(pop2);
            } else if (CGAL_2D_Intersection_Ray_Segment(e + eps * cutDir, -cutDir, pop1, pop2, inter)) {
                auto tanAngle = CGAL_Get_Angle_Kerf_Offset_Tan(segDir, cutDir);
                inter += (1.5875 / tanAngle * cutDir);
                raySnap.push_back(inter);
            }
        }
        //std::cerr << "raysnap2 = " << raySnap.size() << std::endl;
        if (raySnap.empty())
            return false;
        if (raySnap.size() % 2 == 0) {
            double minDist = DBL_MAX;
            for (const auto &v : raySnap) {
                //std::cerr << "update = " << v[0] << " " << v[1] << std::endl;
                const double tmpDist2 = length2(v - e);
                if (tmpDist2 < minDist) {
                    ne = v;
                    minDist = tmpDist2;
                }
            }
        } else {
            std::cerr << "raysnap2 % 2 != 0" << " size = " << raySnap.size() << std::endl;
            //std::string filename = "D:\\cgaldebug\\" + std::to_string(rand()).append(".obj");
            //std::cerr << "raysnap2 % 2 != 0 and see " << filename << std::endl;
            // 			Vector2d2  poly = { polygon };
            // 			OutputRectangle(filename, poly);
            // 			auto fromS = e + eps * cutDir;
            // 			auto toE = -cutDir;
            // 			std::cerr << fromS[0] << " " << fromS[1] << " -> " << toE[0] << " " << toE[1]  << std::endl;
            // 			system("pause");
        }
    } else {
        for (int i = 0; i < polySize; i++) {
            Vector2d inter;
            const auto &pop1 = polygon[i];
            const auto &pop2 = polygon[(i + 1) % polySize];
            if (CGAL_2D_Intersection_Ray_Segment(e - eps * cutDir, cutDir, pop1, pop2, inter)) {
                Vector2d segDir = normalize(pop1 - pop2);
                auto tanAngle = CGAL_Get_Angle_Kerf_Offset_Tan(segDir, cutDir);
                //std::cerr << "angle = " << tanAngle << std::endl;
                inter += (1.5875 / tanAngle * cutDir);
                rayD2Int.push_back(inter);
            }
        }
        //std::cerr << "rayd2int = " << rayD2Int.size() << std::endl;
        if (rayD2Int.empty()) return false;

        double minDist = DBL_MAX;
        for (const auto &v : rayD2Int) {
            //std::cerr << "update = " << v[0] << " " << v[1] << std::endl;
            const double tmpDist2 = length2(v - e);
            if (tmpDist2 < minDist) {
                ne = v;
                minDist = tmpDist2;
            }
        }
    }

    //std::cerr << "output s = " << ns[0] << " " << ns[1] << std::endl;
    //std::cerr << "output e = " << ne[0] << " " << ne[1] << std::endl;

    //std::cerr << "pd1 = " << pd1 << " pd2 = " << pd2 << std::endl;
    return true;
}

// This one is used to intersect a polygon with a line
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut_ExtendOld(
        const Vector2d1 &polygon,
        const Vector2d &s, const Vector2d &e,
        Vector2d &ns, Vector2d &ne) {
    const double eps = 1e-4;
    Vector2d cutDir = normalize(e - s);
    ns = s, ne = e;
    Vector2d ts = s, te = e;
    const int polySize = (int)polygon.size();
    Vector2d1 rayD1Int, rayD2Int;

    const auto pd1 = CGAL_2D_Distance_Point_Polygon(s, polygon);
    const auto pd2 = CGAL_2D_Distance_Point_Polygon(e, polygon);
    bool isoutside1 = false, isoutside2 = false;

    if (!Functs::IsAlmostZero(pd1)) {
        isoutside1 = CGAL_2D_Is_Point_OutPolygon(s, polygon);
        //std::cerr << "inside1 = " << inside1 << " inside2 = " << inside2 << std::endl;
        if (isoutside1) {
            Vector2d1 raySnap;
            Vector2d inter;
            for (int i = 0; i < polySize; i++) {
                if (CGAL_2D_Intersection_Ray_Segment(s - eps * cutDir, cutDir, polygon[i],
                                                     polygon[(i + 1) % polySize], inter))
                    raySnap.push_back(inter);
            }
            //std::cerr << "raysnap1 = " << raySnap.size() << std::endl;
            if (raySnap.empty())
                return false;
            double minDist = DBL_MAX;
            for (const auto &v : raySnap) {
                const double tmpDist2 = length2(v - s);
                if (tmpDist2 < minDist) {
                    ts = v;
                    minDist = tmpDist2;
                }
            }
        }
    }

    if (!Functs::IsAlmostZero(pd2)) {
        isoutside2 = CGAL_2D_Is_Point_OutPolygon(e, polygon);
        if (isoutside2) {
            Vector2d1 raySnap;
            Vector2d inter;
            for (int i = 0; i < polySize; i++) {
                if (CGAL_2D_Intersection_Ray_Segment(e + eps * cutDir, -cutDir, polygon[i],
                                                     polygon[(i + 1) % polySize], inter))
                    raySnap.push_back(inter);
            }
            //std::cerr << "raysnap2 = " << raySnap.size() << std::endl;
            if (raySnap.empty())
                return false;
            double minDist = DBL_MAX;
            for (const auto &v : raySnap) {
                const double tmpDist2 = length2(v - e);
                if (tmpDist2 < minDist) {
                    te = v;
                    minDist = tmpDist2;
                }
            }
        }
    }

    //std::cerr << "ts = " << ts[0] << " " << ts[1] << std::endl;
    //std::cerr << "te = " << te[0] << " " << te[1] << std::endl;
    if (isoutside1 && isoutside2) {
        ns = ts;
        ne = te;
        return true;
    }

    //std::cerr << "pd1 = " << pd1 << " pd2 = " << pd2 << std::endl;
    if (!Functs::IsAlmostZero(pd1)) {
        for (int i = 0; i < polySize; i++) {
            Vector2d inter1;
            if (CGAL_2D_Intersection_Ray_Segment(ts + eps * cutDir, -cutDir, polygon[i],
                                                 polygon[(i + 1) % polySize], inter1))
                rayD1Int.push_back(inter1);
        }
        //std::cerr << "rayd1int = " << rayD1Int.size() << std::endl;
        if (rayD1Int.empty()) {
            // 			Vector2d2  debug = { polygon };
            // 			OutputRectangle("D:\\final.obj", debug);
            return false;
        }

        double minDist = DBL_MAX;
        for (const auto &v : rayD1Int) {
            const double tmpDist2 = length2(v - s);
            if (tmpDist2 < minDist) {
                ns = v;
                minDist = tmpDist2;
            }
        }
    }

    if (!Functs::IsAlmostZero(pd2)) {
        for (int i = 0; i < polySize; i++) {
            Vector2d inter2;
            if (CGAL_2D_Intersection_Ray_Segment(te - eps * cutDir, cutDir, polygon[i],
                                                 polygon[(i + 1) % polySize], inter2))
                rayD2Int.push_back(inter2);
        }
        //std::cerr << "rayd2int = " << rayD2Int.size() << std::endl;
        if (rayD2Int.empty()) return false;

        double minDist = DBL_MAX;
        for (const auto &v : rayD2Int) {
            const double tmpDist2 = length2(v - e);
            if (tmpDist2 < minDist) {
                ne = v;
                minDist = tmpDist2;
            }
        }
    }

    return true;
}


extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut(const Vector2d1 &polygon,const Vector2d1 &cutLine, VectorPB1 &result) {
    // N-1 edges, default 0 -> cant be fabricated
    result = VectorPB1(cutLine.size() - 1, std::make_pair<bool, bool>(false, false));

    Polygon_2 poly;
    for (auto &p : polygon)
        poly.push_back(Point_2(p[0], p[1]));

    if (!poly.is_simple()) {
        std::cerr << "Dynamic CGAL: Polygon is not simple" << std::endl;
        return false;
    }

    result.front().first = true;
    result.back().second = true;

    const double polyArea = poly.area();
    const double offLine = polyArea / 1e4;

    for (auto i = 1; i < cutLine.size(); ++i) {
        int ii = i - 1;
        Point_2 curA(cutLine[ii][0], cutLine[ii][1]);
        Point_2 curB(cutLine[i][0], cutLine[i][1]);
        Vector_2 cutDir = curB - curA;
        //cutDir /= CGAL::sqrt(cutDir.squared_length());
        cutDir = Vector_2(cutDir.x() / CGAL::sqrt(cutDir.squared_length()),
                          cutDir.y() / CGAL::sqrt(cutDir.squared_length()));

        //cutDir *= offLine;
        cutDir = Vector_2(cutDir.x() * offLine, cutDir.y() * offLine);

        Vector_2 invCutDir = -cutDir;

        //curA += invCutDir;
        curA = Point_2(curA.x() + invCutDir.x(), curA.y() + invCutDir.y());
        //curB += cutDir;
        curB = Point_2(curB.x() + cutDir.x(), curB.y() + cutDir.y());

        if (poly.bounded_side(curA) != CGAL::ON_BOUNDED_SIDE) result[ii].first = true;
        if (poly.bounded_side(curB) != CGAL::ON_BOUNDED_SIDE) result[ii].second = true;
    }
    //
    // 	for (auto& p : cutLine)
    // 	{
    // 		switch (poly.bounded_side(Point_2(p[0], p[1]))) {
    // 		case CGAL::ON_BOUNDED_SIDE:
    // 			std::cerr << " is inside the polygon.\n";
    // 			break;
    // 		case CGAL::ON_BOUNDARY:
    // 			std::cerr << " is on the polygon boundary.\n";
    // 			break;
    // 		case CGAL::ON_UNBOUNDED_SIDE:
    // 			std::cerr << " is outside the polygon.\n";
    // 			break;
    // 		}
    // 	}
    return true;
}
