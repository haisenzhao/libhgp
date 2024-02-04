#include "geom.h"
#include "clipper/clipper.hpp"
#include <Mathematics/Vector2.h>

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

/********************************************************
* Function name: CGAL_2D_Distance_Point_Point
* Description: Calculate the Euclidean distance between two 2D points.
* Parameters:
* @p_0: The first 2D point.
* @p_1: The second 2D point.
* Return: The Euclidean distance between the two input points.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Point(const Vector2d& p_0, const Vector2d& p_1)
{
    return sqrt(pow((p_0[0] - p_1[0]), 2.0) + pow((p_0[1] - p_1[1]), 2.0));
}

/********************************************************
* Function name: CGAL_2D_Distance_Point_Segment
* Description: Calculate the distance between a 2D point and a line segment.
* Parameters:
* @v: The 2D point for which the distance is calculated.
* @s_0: The starting point of the line segment.
* @s_1: The ending point of the line segment.
* Return: The distance between the point and the line segment.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Segment(const Vector2d& v, const Vector2d & s_0, const Vector2d & s_1) {
    return sqrt((double) CGAL::squared_distance(Point_2(v[0], v[1]), Segment_2(Point_2(s_0[0], s_0[1]), Point_2(s_1[0], s_1[1]))));
}

/********************************************************
* Function name: CGAL_2D_Distance_Segment_Segment
* Description: Calculate the minimum distance between two 2D line segments.
* Parameters:
* @s_0: The starting point of the first line segment.
* @s_1: The ending point of the first line segment.
* @e_0: The starting point of the second line segment.
* @e_1: The ending point of the second line segment.
* Return: The minimum distance between the two line segments.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Segment_Segment(const Vector2d& s_0, const Vector2d & s_1, const Vector2d & e_0, const Vector2d & e_1) {
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

/********************************************************
* Function name: CGAL_2D_Distance_Point_Line
* Description: Calculate the distance between a 2D point and a 2D line.
* Parameters:
* @v: The 2D point for which you want to calculate the distance.
* @l_0: The first point on the line.
* @l_1: The second point on the line.
* Return: The distance between the point and the line.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Line(const Vector2d& v, const Vector2d & l_0, const Vector2d & l_1) {
    return sqrt((double) CGAL::squared_distance(Point_2(v[0], v[1]),
                                                Line_2(Point_2(l_0[0], l_0[1]), Point_2(l_1[0], l_1[1]))));
}

/********************************************************
* Function name: CGAL_2D_Distance_Point_Polygon
* Description: Calculate the minimum distance between a 2D point and a 2D polygon.
* Parameters:
* @p: The 2D point for which you want to calculate the distance.
* @py: A vector of 2D points representing the vertices of the polygon.
* Return: The minimum distance between the point and the polygon.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Polygon(const Vector2d& p, const Vector2d1 & py) {
    double distance = 1000000000000.0;
    for (int i = 0; i < py.size(); i++)
    {
        int ii = (static_cast<int>(i) + 1) % py.size();
		distance = std::min(distance, CGAL_2D_Distance_Point_Segment(p, py[i], py[ii]));
    }
    return distance;
}

/********************************************************
* Function name: CGAL_2D_Distance_Point_Polygons
* Description: Calculate the minimum distance between a 2D point and a collection of 2D polygons.
* Parameters:
* @p: The 2D point for which you want to calculate the distance.
* @pys: A vector of vectors, where each inner vector represents a 2D polygon defined by its vertices.
* Return: The minimum distance between the point and the collection of polygons.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Point_Polygons(const Vector2d & p, const Vector2d2 & pys)
{
	double distance = 1000000000000.0;
	for (int i = 0; i < pys.size(); i++)
		distance = std::min(distance, CGAL_2D_Distance_Point_Polygon(p, pys[i]));
	return distance;
}

/********************************************************
* Function name: CGAL_2D_Is_Point_OutCGALPolygon
* Description: Check if a 2D point is outside a CGAL polygon.
* Parameters:
* @p: The 2D point you want to check.
* @py: The CGAL polygon (Polygon_2) against which you want to check if the point is outside.
* Return:
* - true if the point is outside the polygon (on the unbounded side).
* - false if the point is inside or on the boundary of the polygon.
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Is_Point_OutCGALPolygon(const Vector2d &p, const Polygon_2 &py) {
    return py.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE;
}

/********************************************************
* Function name: CGAL_Construct_Polygon
* Description: Construct a CGAL polygon from a list of 2D points and check if it is simple.
* Parameters:
* @py: A list of 2D points (Vector2d1) that represent the vertices of the polygon.
* @poly: The CGAL polygon (Polygon_2) to be constructed.
* Return:
* - true if the constructed polygon is simple (has no self-intersections).
* - false if the constructed polygon is not simple (has self-intersections).
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_Construct_Polygon(const Vector2d1 &py, Polygon_2 &poly) {
    poly.clear();
    for (auto i : py)
        poly.push_back(Point_2(i[0], i[1]));
    return poly.is_simple();
}

/********************************************************
* Function name: CGAL_Construct_InOutSide_Polygon
* Description: Construct a CGAL polygon from a list of 2D points, check if it is simple,
* and determine if two points are inside or outside the polygon.
* Parameters:
* @py: A list of 2D points (Vector2d1) that represent the vertices of the polygon.
* @p: The first 2D point to be checked.
* @q: The second 2D point to be checked.
* @isPInside: A reference to a boolean variable that will be set to true if point 'p' is inside the polygon.
* @isQInside: A reference to a boolean variable that will be set to true if point 'q' is inside the polygon.
* Return:- true if the polygon is successfully constructed, is simple, and points 'p' and 'q' are checked; - false if the polygon construction fails or if it is not simple.
*********************************************************/
extern "C" PPGL_EXPORT bool
CGAL_Construct_InOutSide_Polygon(const Vector2d1 &py, const Vector2d &p, const Vector2d &q, bool &isPInside,
                                 bool &isQInside) {
    Polygon_2 poly = Polygon2D(py);
    if (!poly.is_simple()) return false;

    isPInside = CGAL_2D_Is_Point_OutCGALPolygon(p, poly);
    isQInside = CGAL_2D_Is_Point_OutCGALPolygon(q, poly);
    return true;
}

/********************************************************
* Function name: CGAL_2D_Location_Point_Polygon
* Description: Check the location of a 2D point relative to a CGAL polygon.
* Parameters:
* @p: The 2D point (Vector2d) to be checked for its location.
* @py: A list of 2D points (Vector2d1) that represent the vertices of the CGAL polygon.
* Return: - true if the point 'p' is located inside the polygon (bounded side). - false if the point 'p' is located outside the polygon (unbounded side).
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Location_Point_Polygon(const Vector2d & p, const Vector2d1 & py) {
	Polygon_2 poly = Polygon2D(py);

    return poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_BOUNDED_SIDE;
}

/********************************************************
* Function name: CGAL_2D_Is_Point_OutPolygon
* Description: Check if a 2D point is outside a CGAL polygon.
* Parameters:
* @p: The 2D point (Vector2d) to be checked for its location.
* @py: A list of 2D points (Vector2d1) that represent the vertices of the CGAL polygon.
* Return: - true if the point 'p' is located outside the polygon (unbounded side). - false if the point 'p' is located inside the polygon (bounded side).
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Is_Point_OutPolygon(Vector2d p, Vector2d1 py) {
    Polygon_2 poly = Polygon2D(py);
    return poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE;
}

/********************************************************
* Function name: CGAL_2D_Location_Points_Polygon
* Description: Check if a list of 2D points are all located inside a CGAL polygon.
* Parameters:
* @ps: A list of 2D points (Vector2d1) to be checked for their locations.
* @py: A list of 2D points (Vector2d1) that represent the vertices of the CGAL polygon.
* Return: - true if all points in the list 'ps' are located inside the polygon (bounded side). - false if at least one point in the list 'ps' is located outside the polygon (unbounded side).
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Location_Points_Polygon(const Vector2d1 &ps,
                                                                      const Vector2d1 &py) {
	Polygon_2 poly = Polygon2D(py);

    for (auto p : ps) {
        if (poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE)
            return false;
    }
    return true;
}

/********************************************************
* Function name: CGAL_2D_Polygon_Dart_Sampling
* Description: Sample points inside a CGAL polygon using the dart throwing algorithm.
* Parameters:
* @py: A list of 2D points (Vector2d1) representing the vertices of the CGAL polygon.
* @d: The minimum distance between sampled points, as a fraction of the diagonal length of the polygon's bounding box.
* @sampling_points: An output list (Vector2d1) that will contain the sampled points.
* @total_iter: The maximum number of iterations to attempt for point sampling.
* Return: This function doesn't return a value, but it populates the 'sampling_points' list with the sampled points.
*********************************************************/
extern "C" PPGL_EXPORT void CGAL_2D_Polygon_Dart_Sampling(const Vector2d1 & py, const double& d, Vector2d1 & sampling_points, const int& total_iter)
{
	Functs::MAssert(d>0&&d<1.0,"CGAL_2D_Polygon_Dart_Sampling if (!(d > 0 && d < 1.0))");

	Polygon_2 poly = Polygon2D(py);

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

extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Polygon_Regular_Sampling_C1(const Vector2d1& py, const double& d)
{
	VectorPI1 neighbors;
	return CGAL_2D_Polygon_Regular_Sampling_C3(py, d, neighbors, false);
}


extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Polygon_Regular_Sampling_C2(const Vector2d1& py, const double& d, VectorPI1& neighbors)
{
	return CGAL_2D_Polygon_Regular_Sampling_C3(py, d, neighbors, true);
}

/********************************************************
* Function name: CGAL_2D_Polygon_Regular_Sampling_C3
* Description: Computes a regular sampling of a 2D polygon.
* Parameters:
* @py: A vector of 2D points representing a polygon.
* @d: percentage value of the length of the diagonal of the bounding box.
* @neighbors: A vector to store neighbor relationships (optional).
* @compute_neighbors: A flag to indicate whether to compute neighbor relationships (optional).
* Return: A vector of regularly sampled 2D points within the polygon.
* *********************************************************/
extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Polygon_Regular_Sampling_C3(const Vector2d1& py, const double& d, VectorPI1& neighbors, const bool& compute_neighbors)
{
	//check input
	Functs::MAssert(d > 0 && d < 1.0, "CGAL_2D_Polygon_Regular_Sampling_C1 if (!(d > 0 && d < 1.0))");
	Functs::MAssert(py.size() >= 3, "CGAL_2D_Polygon_Regular_Sampling_C1: py.size()>=3");
	Functs::MAssert(CGAL_2D_Polygon_Area(py) > 0.0, "CGAL_2D_Polygon_Regular_Sampling_C1: CGAL_2D_Polygon_Area(py)>0.0");

	Polygon_2 poly = Polygon2D(py);

	Vector2d minC, maxC;
	Functs::GetBoundingBox(py, minC, maxC);

	double minimal_d = d * Functs::GetDistance(minC, maxC);

	Vector2d1 sampling_points;
	double x(minC[0]);
	Vector1i2 xyb;
	while (x < maxC[0])
	{
		double y(minC[1]);
		Vector1i1 yb;
		while (y < maxC[1])
		{
			if (poly.bounded_side(Point_2(x, y)) == CGAL::ON_BOUNDED_SIDE)
			{
				yb.push_back(sampling_points.size());
				sampling_points.push_back(Vector2d(x, y));
			}
			else
			{
				yb.push_back(-1);
			}
			y += minimal_d;
		}
		x += minimal_d;
		xyb.push_back(yb);
	}

	neighbors.clear();
	if(compute_neighbors)
	{
		// get the inside or outside relation of a
		auto GetPB = [](const Vector1i2& xyb, const Vector2i& a)
		{
			const int xe = xyb.size();
			const int ye = xyb.front().size();
			if (!(a[0] >= 0 && a[1] >= 0 && a[0] < xe && a[1] < ye)) return -1;
			return xyb[a[0]][a[1]];
		};

		auto PushEdge = [&](const Vector2i& xy2i, const Vector2i& v2i)
		{
			if (GetPB(xyb, v2i) >= 0)
			{
				if (GetPB(xyb, xy2i) >= 0)
				{
					int index_0 = xyb[xy2i[0]][xy2i[1]];
					int index_1 = xyb[v2i[0]][v2i[1]];
					if (!CGAL_2D_Intersection_Segment_Polygon(sampling_points[index_0], sampling_points[index_1], py))
					{
						neighbors.push_back(std::pair<int, int>(index_0, index_1));
					}
				}
			}
		};

		for (int xi = 0; xi < xyb.size(); xi++)
		{
			for (int yi = 0; yi < xyb[xi].size(); yi++)
			{
				if (xyb[xi][yi] >= 0)
				{
					PushEdge(Vector2i(xi, yi), Vector2i(xi + 1, yi));
					PushEdge(Vector2i(xi, yi), Vector2i(xi, yi + 1));
					PushEdge(Vector2i(xi, yi), Vector2i(xi + 1, yi + 1));
					PushEdge(Vector2i(xi + 1, yi), Vector2i(xi, yi + 1));
				}
			}
		}
	}

	return sampling_points;
}


extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Square_Regular_Sampling_C1(const double& d)
{
	return CGAL_2D_Polygon_Regular_Sampling_C1(Functs::GetUnitSquare(), d);
}
extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Square_Regular_Sampling_C2(const double& d, VectorPI1& neighbors)
{
	return CGAL_2D_Polygon_Regular_Sampling_C2(Functs::GetUnitSquare(), d, neighbors);
}

extern "C" PPGL_EXPORT Vector2d1 CGAL_2D_Square_Regular_Sampling_C3(const double& d, VectorPI1& neighbors, const bool& compute_neighbors)
{
	return CGAL_2D_Polygon_Regular_Sampling_C3(Functs::GetUnitSquare(), d, neighbors, compute_neighbors);
}

/********************************************************
* Function name: CGAL_2D_Intersection_Segment_Segment
* Description: Computes the intersection point between two line segments in 2D, if it exists.
* Parameters:
* @s_0_s: Starting point of the first line segment.
* @s_0_e: Ending point of the first line segment.
* @s_1_s: Starting point of the second line segment.
* @s_1_e: Ending point of the second line segment.
* @inter: If an intersection point exists, it is stored in this vector (x, y).
* Return: True if an intersection point exists, false otherwise.
* *********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Segment
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

/********************************************************
* Function name: CGAL_2D_Intersection_Ray_Segment
* Description: Computes the intersection point between a ray and a line segment in 2D, if it exists.
* Parameters:
* @s_0_s: Starting point of the ray.
* @s_0_e: Endpoint of the ray (direction vector).
* @s_1_s: Starting point of the line segment.
* @s_1_e: Ending point of the line segment.
* @inter: If an intersection point exists, it is stored in this vector (x, y).
* Return: True if an intersection point exists, false otherwise.
* *********************************************************/
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

/********************************************************
* Function name: CGAL_2D_Intersection_Ray_Segment_Ignore_Endpoints
* Description: Computes the intersection point between a ray and a line segment in 2D, excluding endpoints of the segment.
* Parameters:
* @s_0_s: Starting point of the ray.
* @s_0_e: Endpoint of the ray (direction vector).
* @s_1_s: Starting point of the line segment.
* @s_1_e: Ending point of the line segment.
* @inter: If an intersection point exists and is not an endpoint of the segment, it is stored in this vector (x, y).
* Return: True if a valid intersection point exists, false otherwise.
* *********************************************************/
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


/********************************************************
* Function name: CGAL_2D_Intersection_Ray_Polygon
* Description: Determines if a ray intersects a 2D polygon and, if so, returns the intersection point.
* Parameters:
* @r_s: Starting point of the ray.
* @r_d: Direction vector of the ray.
* @poly: Points defining the vertices of the polygon as a list of 2D vectors.
* @pnt: If an intersection point exists, it is stored in this vector (x, y).
* Return: True if the ray intersects the polygon and a valid intersection point is found, false otherwise.
* *********************************************************/
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

/********************************************************
* Function name: CGAL_2D_Intersection_Line_Line
* Description: Determines if two 2D lines intersect and, if so, returns the intersection point.
* Parameters:
* @s_0_s: Starting point of the first line.
* @s_0_e: Ending point of the first line.
* @s_1_s: Starting point of the second line.
* @s_1_e: Ending point of the second line.
* @inter: If lines intersect, the intersection point (x, y) is stored in this vector.
* Return: True if the lines intersect and a valid intersection point is found, false otherwise.
* *********************************************************/
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

/********************************************************
* Function name: CGAL_2D_Intersection_Segment_Line
* Description: Determines if a 2D line segment and a 2D line intersect, and if so, returns the intersection point.
* Parameters:
* @s_s: Starting point of the line segment.
* @s_e: Ending point of the line segment.
* @l_s: Starting point of the line.
* @l_e: Ending point of the line.
* @inter: If the line segment and line intersect, the intersection point (x, y) is stored in this vector.
* Return: True if they intersect and a valid intersection point is found, false otherwise.
* *********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Segment_Line(const Vector2d & s_s, const Vector2d & s_e, const Vector2d & l_s, const Vector2d & l_e, Vector2d & inter)
{
	CGAL::Object result = CGAL::intersection(Segment_2(Point_2(s_s[0], s_s[1]), Point_2(s_e[0], s_e[1])),
		Line_2(Point_2(l_s[0], l_s[1]), Point_2(l_e[0], l_e[1])));

	if (const Point_2* ipoint = CGAL::object_cast<Point_2>(&result))
	{
		inter[0] = ipoint->x();
		inter[1] = ipoint->y();
		return true;
	}
	else
	{
		return false;
	}
}

/********************************************************
* Function name: CGAL_2D_Intersection_Segment_Polygon
* Description: Determines if a 2D line segment intersects with a 2D polygon by checking for segment-polygon intersection.
* Parameters:
* @s_s: Starting point of the line segment.
* @s_e: Ending point of the line segment.
* @p: A vector of points defining the vertices of the polygon.
* Return: True if the line segment intersects with the polygon, false otherwise.
* *********************************************************/
extern "C" PPGL_EXPORT bool
CGAL_2D_Intersection_Segment_Polygon(const Vector2d & s_s, const Vector2d & s_e, const Vector2d1 &p) {
    for (int i = 0; i < p.size(); i++) {
        Vector2d inter;
        int ii = (static_cast<int>(i) + 1) % p.size();
        if (CGAL_2D_Intersection_Segment_Segment(s_s, s_e, p[i], p[ii], inter)) {
            return true;
        }
    }
    return false;
}

//extern "C" PPGL_EXPORT bool CGAL_2D_Intersection_Polygon_Polygon(const Vector2d1 & p1, const Vector2d1 & p2);

extern "C" PPGL_EXPORT bool
CGAL_2D_Intersection_Polygon_Polygon(const Vector2d1 & p1, const Vector2d1 & p2) {

	Polygon_2 poly1 = Polygon2D(p1);
	Polygon_2 poly2 = Polygon2D(p2);

	return CGAL::do_intersect(poly1, poly2);
}


extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Is_Clockwise_Oriented(const Vector2d1 &ps) {
	Polygon_2 poly = Polygon2D(ps);
    return poly.is_clockwise_oriented();
}

/********************************************************
* Function name: CGAL_2D_Two_Polygons_Intersection
* Description: Calculates the intersection area of two 2D polygons.
* Parameters:
* @poly_0: First polygon defined by a vector of 2D points.
* @poly_1: Second polygon defined by a vector of 2D points.
* Return: The intersection area of the two polygons.
* *********************************************************/
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

/********************************************************
* Function name: CGAL_2D_Two_Polygons_Union
* Description: Calculates the union area of two 2D polygons and provides the resulting polygons.
* Parameters:
* @poly_0: First polygon defined by a vector of 2D points.
* @poly_1: Second polygon defined by a vector of 2D points.
* @inter_polygons: Output parameter for resulting polygons after union.
* Return: The union area of the two polygons.
* *********************************************************/
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

        Vector1d1 xs;
        Vector1d1 ys;

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

/********************************************************
* Function name: CGAL_2D_Polygon_One_Offsets
* Description: Computes one-side offsets of a 2D polygon.
* Parameters:
* @poly: Input polygon defined by a vector of 2D points.
* @d: Offset distance.
* @offset_polys: Output parameter for resulting offset polygons.
* *********************************************************/
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
        one_offset= Functs::RemoveClosePoints(one_offset);
        offset_polys.push_back(one_offset);
        Vector2d1().swap(one_offset);
    }
}


/********************************************************
* Function name: CGAL_2D_Polygons_One_Offsets
* Description: Computes one-side offsets of a collection of 2D polygons.
* Parameters:
* @polys: Input polygons defined by a vector of vector of 2D points.
* @d: Offset distance.
* @offset_polys: Output parameter for resulting offset polygons.
* *********************************************************/
extern "C" PPGL_EXPORT void CGAL_2D_Polygons_One_Offsets(const Vector2d2 & polys, const double& d, Vector2d2 & offset_polys)
{
	if (!(polys.size() > 0 && polys[0].size() > 0)) return;

	double scale = 1000000.0;

	ClipperLib::ClipperOffset co;
	co.ArcTolerance = co.ArcTolerance * scale / 1000.0;

	ClipperLib::Path subj;
	ClipperLib::Paths solution;

	//build the most outside path
	for (int i = 0; i < polys[0].size(); i++)
		subj << ClipperLib::IntPoint(polys[0][i][0] * scale, polys[0][i][1] * scale);
	co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);

	//build the following paths
	for (int i = 1; i < polys.size(); i++)
	{
		subj.clear();
		for (int j = 0; j < polys[i].size(); j++)
			subj << ClipperLib::IntPoint(polys[i][j][0] * scale, polys[i][j][1] * scale);
		co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
	}

	// execute
	co.Execute(solution, -d * scale);

	//output
	for (int i = 0; i < solution.size(); i++)
	{
        Vector2d1 one_offset_xys;

		Polygon_2 poly_2;
		for (int j = 0; j < solution[i].size(); j++)
		{
            one_offset_xys.push_back(Vector2d(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
			poly_2.push_back(Point_2(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
		}

		if (poly_2.is_clockwise_oriented())
		{
			std::reverse(one_offset_xys.begin(), one_offset_xys.end());
		}

		//remove closed points
        one_offset_xys = Functs::RemoveClosePoints(one_offset_xys);

        offset_polys.push_back(one_offset_xys);
        Vector2d1().swap(one_offset_xys);
	}
}

extern "C" PPGL_EXPORT bool CGAL_2D_Polygons_Simple(const Vector2d2 & poly)
{
	for (auto& py : poly)
		if (!CGAL_2D_Polygon_Simple(py))
			return false;
	return true;
}

extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Simple(const Vector2d1 & py)
{
	Polygon_2 poly = Polygon2D(py);
	return poly.is_simple();
}

/********************************************************
* Function name: CGAL_2D_Polygon_Simple_Inter
* Description: Checks if a 2D polygon is simple (non-self-intersecting) by testing all non-adjacent edge pairs for intersections.
* Parameters:
* @poly: Input polygon defined by a vector of 2D points.
* Return: True if the polygon is simple; false otherwise.
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Simple_Inter(const Vector2d1 & poly)
{
	for (int i = 0; i < poly.size(); i++)
	{
		auto s_0 = poly[i];
		auto e_0 = poly[(i + 1) % poly.size()];

		for (int j = 0; j < poly.size(); j++)
		{
			auto s_1 = poly[j];
			auto e_1 = poly[(j + 1) % poly.size()];

			if (i != j && i != (j + 1) % poly.size() && (i + 1) % poly.size() != j && (i + 1) % poly.size() != (j + 1) % poly.size())
			{
				Vector2d inter;
				bool b = CGAL_2D_Intersection_Segment_Segment(s_0, e_0, s_1, e_1, inter);

				if (b)
				{
					return false;
				}
			}
		}
	}
	return true;
}

/********************************************************
* Function name: CGAL_2D_Convex_Hulls
* Description: Computes the convex hull of a set of 2D points using CGAL's convex_hull_2 algorithm.
* Parameters:
* @vec: Input vector of 2D points.
* @hull_points: Output vector to store the points on the convex hull.
*********************************************************/
extern "C" PPGL_EXPORT void CGAL_2D_Convex_Hulls(const Vector2d1 & vec, Vector2d1 & hull_points)
{
	std::vector<Point_2> points;
	std::vector<Point_2> results;

	for (int i = 0; i < vec.size(); i++)
		points.push_back(VectorPoint2d(vec[i]));

	CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(results));
	std::cout << results.size() << " points on the convex hull" << std::endl;

	for (int i = 0; i < results.size(); i++)
		hull_points.push_back(PointVector2d(results[i]));
}

/********************************************************
* Function name: CGAL_2D_OBB_Box
* Description: Computes the oriented bounding box (OBB) of a set of 2D points.
* Parameters:
* @vec: Input vector of 2D points.
* @center: Output vector to store the OBB center.
* @axis_0: Output vector to store the first OBB axis.
* @axis_1: Output vector to store the second OBB axis.
* @extent_0: Output variable to store the extent along the first axis.
* @extent_1: Output variable to store the extent along the second axis.
*********************************************************/
extern "C" PPGL_EXPORT void CGAL_2D_OBB_Box(const Vector2d1 & vec, Vector2d & center, Vector2d & axis_0, Vector2d & axis_1, double& entent_0, double& entent_1)
{
	gte::Vector2<double>* points = new gte::Vector2<double>[vec.size()];

	for (int i = 0; i < vec.size(); i++)
	{
		points[i][0] = vec[i][0];
		points[i][1] = vec[i][1];
	}

    gte::IntrinsicsVector2<double> iv(vec.size(), points, static_cast<double>(0));

    center[0] = iv.origin[0];
	center[1] = iv.origin[1];
	axis_0[0] = iv.direction[0][0];
	axis_0[1] = iv.direction[0][1];
	axis_1[0] = iv.direction[1][0];
	axis_1[1] = iv.direction[1][1];

	entent_0 = iv.max[0 ]- iv.min[0];
	entent_1 = iv.max[1] - iv.min[1];
}

extern "C" PPGL_EXPORT void
CGAL_Decompose_Polyline(const Vector2d1 &polyline, const double& threshold, Vector1i1& result) {
    for (auto &p : polyline) 
    {

    }
}

/********************************************************
* Function name: CGAL_Identify_Polycut_NotExtend
* Description: Identifies if a line segment does not extend into a polygon.
* Parameters:
* @polygon: Input vector representing the polygon vertices.
* @s: Start point of the segment.
* @e: End point of the segment.
* Returns: True if the segment does not extend into the polygon, false otherwise.
*********************************************************/
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

/********************************************************
* Function name: OutputRectangle
* Description: Writes 2D points to an .obj file format.
* Parameters:
* @path: Path to the output .obj file.
* @points: 2D points to be written to the file.
*********************************************************/
void OutputRectangle(const char* path, const Vector2d2  &points) {
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
            file << Functs::Int2String(nb) << " ";
            nb++;
        }
        file << "" << std::endl;
    }


    file.clear();
    file.close();
}

/********************************************************
* Function name: CGAL_Get_Angle_Kerf_Offset_Tan
* Description: Calculate the tangent of the angle between two 2D vectors.
* Parameters:
* @a: First 2D vector.
* @b: Second 2D vector.
* Return: The tangent of the angle between the vectors.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_Get_Angle_Kerf_Offset_Tan(const Vector2d &a, const Vector2d &b) {
    auto na = normalize(a);
    auto nb = normalize(b);
    return glm::tan(glm::acos(glm::abs(dot(na, nb))));
}

/********************************************************
* Function name: CGAL_2D_Projection_Point_Segment
* Description: Project a 2D point onto a line segment defined by two endpoints.
* Parameters:
* @p: The 2D point to be projected.
* @s: The starting endpoint of the line segment.
* @e: The ending endpoint of the line segment.
* Return: Vector2d representing the projected point on the line segment.
*********************************************************/
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Projection_Point_Segment(const Vector2d & p, const Vector2d & s, const Vector2d & e)
{
	Line_2 l(VectorPoint2d(s), VectorPoint2d(e));
	Point_2 m_p = l.projection(VectorPoint2d(p));

	double d_m_s = CGAL_2D_Distance_Point_Point(Vector2d(m_p[0], m_p[1]), s);
	double d_m_e = CGAL_2D_Distance_Point_Point(Vector2d(m_p[0], m_p[1]), e);
	double d_s_e = CGAL_2D_Distance_Point_Point(s,e);

	if (d_m_s >= d_s_e) return e;
	if (d_m_e >= d_s_e) return s;
    return PointVector2d(m_p);
}

/********************************************************
* Function name: CGAL_2D_Detect_Polygon_Inside_C1
* Description: Detect if a 2D point is inside or near the boundary of a polygon.
* Parameters:
* @outside_py: A vector of 2D points defining the outer boundary of the polygon.
* @p: The 2D point to be tested.
* Return: A boolean indicating if the point is inside or near the boundary of the polygon.
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C1(const Vector2d1 & outside_py, const Vector2d & p)
{
	Vector2d outside_max, outside_min;
    Functs::GetBoundingBox(outside_py, outside_min, outside_max);

	if (outside_max[0] < p[0])return false;
	if (outside_max[1] < p[1])return false;
	if (p[0] < outside_min[0])return false;
	if (p[1] < outside_min[1])return false;

	Polygon_2 poly = Polygon2D(outside_py);
	if (poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE)
	{
		double dis = CGAL_2D_Distance_Point_Polygon(p, outside_py);
		if (dis > 0.1)
			return false;
	}

	return true;
}

/********************************************************
* Function name: CGAL_2D_Detect_Polygon_Inside_C2
* Description: Detect if a set of 2D points are inside or near the boundary of a polygon.
* Parameters:
* @outside_py: A vector of 2D points defining the outer boundary of the polygon.
* @inside_py: A vector of 2D points to be tested for being inside or near the polygon.
* Return: A boolean indicating if all points are inside or near the boundary of the polygon.
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C2(const Vector2d1 & outside_py, const Vector2d1 & inside_py)
{
	Vector2d outside_max, outside_min, inside_max, inside_min;
    Functs::GetBoundingBox(outside_py, outside_min, outside_max);
	Functs::GetBoundingBox(inside_py, inside_min, inside_max);


	if (outside_max[0] < inside_min[0])return false;
	if (outside_max[1] < inside_min[1])return false;
	if (inside_max[0] < outside_min[0])return false;
	if (inside_max[1] < outside_min[1])return false;


	Polygon_2 poly = Polygon2D(outside_py);

	for (auto p : inside_py)
	{
		if (poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE)
		{
			double dis = CGAL_2D_Distance_Point_Polygon(p, outside_py);
			if (dis > 0.1)
				return false;
		}
	}
	return true;
}

/********************************************************
* Function name: CGAL_2D_Detect_Polygon_Inside_C3
* Description: Detect if a 2D point is inside or near the boundary of multiple polygons.
* Parameters:
* @outside_pys: A vector of vectors of 2D points, each defining the outer boundary of a polygon.
* @p: A 2D point to be tested for being inside or near any of the polygons.
* Return: A boolean indicating if the point is inside or near the boundary of any of the polygons.
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C3(const Vector2d2 & outside_pys, const Vector2d & p)
{
	if (!CGAL_2D_Detect_Polygon_Inside_C1(outside_pys[0], p)) return false;

	for (int i = 1; i < outside_pys.size(); i++)
	{
		Polygon_2 poly = Polygon2D(outside_pys[i]);
		if (poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_BOUNDED_SIDE)
		{
			double dis = CGAL_2D_Distance_Point_Polygon(p, outside_pys[i]);
			if (dis > 0.1)
				return false;
		}
	}

	return true;
}

/********************************************************
* Function name: CGAL_2D_Detect_Polygon_Inside_C4
* Description: Detect if a set of 2D points are inside or near the boundary of multiple polygons.
* Parameters:
* @outside_pys: A vector of vectors of 2D points, each defining the outer boundary of a polygon.
* @inside_py: A vector of 2D points defining a polygon that should be inside or near the boundary of the outer polygons.
* Return: A boolean indicating if the inner polygon is inside or near the boundary of any of the outer polygons.
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C4(const Vector2d2 & outside_pys, const Vector2d1 & inside_py)
{
	if (!CGAL_2D_Detect_Polygon_Inside_C2(outside_pys[0], inside_py)) return false;

	for (int i = 1; i < outside_pys.size(); i++)
	{
		Polygon_2 poly = Polygon2D(outside_pys[i]);

		for (auto p : inside_py)
		{
			if (poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_BOUNDED_SIDE)
			{
				double dis = CGAL_2D_Distance_Point_Polygon(p, outside_pys[i]);
				if (dis > 0.1)
					return false;
			}
		}
	}

	return true;
}

extern "C" PPGL_EXPORT bool CGAL_2D_Detect_Polygon_Inside_C5(const Vector2d2 & outside_pys, const Vector2d2 & inside_pys)
{
	for (auto inside_py : inside_pys)
	{
		if (!CGAL_2D_Detect_Polygon_Inside_C4(outside_pys, inside_py))
		{
			return false;
		}
	}
	return true;
}

/********************************************************
* Function name: CGAL_2D_Distance_Polygon_Polygon
* Description: Calculate the minimum distance between two 2D polygons.
* Parameters:
* @poly_0: A vector of 2D points defining the first polygon.
* @poly_1: A vector of 2D points defining the second polygon.
* Return: The minimum distance between the two polygons.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Polygon_Polygon(const Vector2d1 & poly_0, const Vector2d1 & poly_1)
{
	double min_d = 1000000000.0;
	for (int i = 0; i < poly_0.size(); i++)
	{
		Vector2d v = CGAL_2D_Nearest_Point_Polygon_C1(poly_0[i], poly_1);
		double l = CGAL_2D_Distance_Point_Point(poly_0[i], v);
		if (l < min_d)
		{
			min_d = l;
		}
	}
	return min_d;
}

/********************************************************
* Function name: CGAL_2D_Distance_Polygons_Polygons
* Description: Calculate the minimum distance between sets of 2D polygons.
* Parameters:
* @poly_0: A vector of vectors of 2D points defining the first set of polygons.
* @poly_1: A vector of vectors of 2D points defining the second set of polygons.
* Return: The minimum distance between the two sets of polygons.
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Distance_Polygons_Polygons(const Vector2d2 & poly_0, const Vector2d2 & poly_1)
{
	double min_d = 1000000000.0;
	for (auto poly_ : poly_0)
	{
		for (auto poly__ : poly_1)
		{
			double l = CGAL_2D_Distance_Polygon_Polygon(poly_, poly__);
			if (l < min_d)
			{
				min_d = l;
			}
		}
	}
	return min_d;
}

/********************************************************
* Function name: CGAL_2D_Nearest_Point_Polygon_C1
* Description: Find the nearest point on a 2D polygon to a given point.
* Parameters:
* @v: The 2D point for which the nearest point on the polygon is to be found.
* @poly: A vector of 2D points defining the polygon.
* Return: The nearest point on the polygon to the input point.
*********************************************************/
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Nearest_Point_Polygon_C1(const Vector2d & v, const Vector2d1 & poly)
{
	double min_d = 1000000000.0;
	int min_i = -1;

	for (int i = 0; i < poly.size(); i++)
	{
		double l = CGAL_2D_Distance_Point_Segment(v, poly[i], poly[(i + 1) % poly.size()]);

		if (l < min_d)
		{
			min_d = l;
			min_i = i;
		}
	}

	return CGAL_2D_Projection_Point_Segment(v, poly[min_i], poly[(min_i + 1) % poly.size()]);
}

/********************************************************
* Function name: CGAL_2D_Nearest_Point_Polygon_C2
* Description: Find the nearest point on a 2D polygon to a given point and return both the nearest point and the minimum distance.
* Parameters:
* @v: The 2D point for which the nearest point on the polygon is to be found.
* @poly: A vector of 2D points defining the polygon.
* @p: The nearest point on the polygon to the input point.
* @min_d: The minimum distance from the input point to the polygon.
*********************************************************/
extern "C" PPGL_EXPORT void CGAL_2D_Nearest_Point_Polygon_C2(const Vector2d & v, const Vector2d1 & poly, Vector2d & p, double& min_d)
{
	min_d = 1000000000.0;
	int min_i = -1;

	for (int i = 0; i < poly.size(); i++)
	{
		double l = CGAL_2D_Distance_Point_Segment(v, poly[i], poly[(i + 1) % poly.size()]);

		if (l < min_d)
		{
			min_d = l;
			min_i = i;
		}
	}

    p= CGAL_2D_Projection_Point_Segment(v, poly[min_i], poly[(min_i + 1) % poly.size()]);
}

/********************************************************
* Function name: CGAL_2D_Nearest_Point_Polygons
* Description: Find the nearest point on a collection of 2D polygons to a given point.
* Parameters:
* @v: The 2D point for which the nearest point on the polygons is to be found.
* @polys: A vector of vectors of 2D points representing a collection of polygons.
* Return: The nearest 2D point on any of the polygons to the input point.
*********************************************************/
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Nearest_Point_Polygons(const Vector2d & v, const Vector2d2 & polys)
{
	Vector2d result;
	double min_d = 1000000000.0;
	for (int i = 0; i < polys.size(); i++)
	{
		Vector2d p;
		double p_d;
		CGAL_2D_Nearest_Point_Polygon_C2(v, polys[i], p, p_d);

		if (p_d < min_d)
		{
			min_d = p_d;
			result = p;
		}
	}
	return result;
}

/********************************************************
* Function name: CGAL_2d_Polygon_Boundingbox
* Description: Compute the bounding box (axis-aligned) of a collection of 2D points.
* Parameters:
* @ps: A vector of 2D points for which the bounding box is to be computed.
* @min_corner: A vector to store the minimum corner (bottom-left) of the bounding box.
* @max_corner: A vector to store the maximum corner (top-right) of the bounding box.
*********************************************************/
extern "C" PPGL_EXPORT void CGAL_2d_Polygon_Boundingbox(const Vector2d1 & ps, Vector2d & min_corner, Vector2d & max_corner)
{
    Functs::GetBoundingBox(ps, min_corner, max_corner);
}

/********************************************************
* Function name: CGAL_2D_Polygon_Area
* Description: Compute the signed area of a 2D polygon defined by a vector of 2D points.
* Parameters:
* @py: A vector of 2D points representing the vertices of the polygon.
* Return: The signed area of the polygon (positive if vertices are in counterclockwise order, negative if clockwise).
*********************************************************/
extern "C" PPGL_EXPORT double CGAL_2D_Polygon_Area(const Vector2d1 & py)
{
	Polygon_2 poly = Polygon2D(py);
	return abs(poly.area());
}

/********************************************************
* Function name: CGAL_2D_Polygon_Inside_Point_C1
* Description: Find the approximate center of a convex or concave 2D polygon defined by a vector of 2D points.
* Parameters:
* @py: A vector of 2D points representing the vertices of the polygon.
* Return: A 2D point representing the approximate center of the polygon.
*********************************************************/
extern "C" PPGL_EXPORT Vector2d CGAL_2D_Polygon_Inside_Point_C1(const Vector2d1 & py)
{
	double inside_x, inside_y;
	Polygon_2 poly = Polygon2D(py);

	if (poly.is_clockwise_oriented()) poly.reverse_orientation();

	double max_dis = -10000.0;
	for (int i = 0; i < py.size(); i++) {
		for (int j = i + 2; j < py.size(); j++) {
			double p_x = (py[i][0] + py[j][0]) / 2.0;
			double p_y = (py[i][0] + py[j][0]) / 2.0;
			if (poly.bounded_side(Point_2(p_x, p_y)) == CGAL::ON_BOUNDED_SIDE)
			{
				double dis = CGAL_2D_Distance_Point_Point(py[i], py[j]);
				if (dis > max_dis)
				{
					max_dis = dis;
					inside_x = p_x;
					inside_y = p_y;
				}
			}
		}
	}
	return Vector2d(inside_x, inside_y);
}

/********************************************************
* Function name: CGAL_2D_Polygon_Inside_Point_C2
* Description: Find a point inside multiple 2D polygons. The function generates random points and checks if they are inside all given polygons.
* Parameters:
* @polys: A vector of vectors, where each inner vector represents the vertices of a polygon.
* @inner_vec: A 2D point representing the found point inside the polygons (if successful).
* Return: A boolean indicating if a point inside all polygons was found.
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_2D_Polygon_Inside_Point_C2(const Vector2d2 & polys, Vector2d & inner_vec)
{
	auto CheckValid = [](std::vector<Polygon_2>& cgal_polys, double p_x, double p_y)
	{
		if (cgal_polys.front().bounded_side(Point_2(p_x, p_y)) == CGAL::ON_BOUNDED_SIDE)
		{
			for (int i = 1; i < cgal_polys.size(); i++)
			{
				if (cgal_polys[i].bounded_side(Point_2(p_x, p_y)) == CGAL::ON_BOUNDED_SIDE)
				{
					return false;
				}
			}
			return true;
		}
		return false;
	};

	std::vector<Polygon_2> cgal_polys;
	for (int i = 0; i < polys.size(); i++)
	{
		Polygon_2 poly = Polygon2D(polys[i]);

		if (poly.is_clockwise_oriented()) poly.reverse_orientation();
		cgal_polys.emplace_back(poly);
	}

	auto outside_poly = cgal_polys.front();

	double xmin = outside_poly.bbox().xmin();
	double ymin = outside_poly.bbox().ymin();
	double xmax = outside_poly.bbox().xmax();
	double ymax = outside_poly.bbox().ymax();
	int inter = 0;
	int success_iter = 0;
	double dis_success;
	while (true)
	{
		double p_x = rand() / double(RAND_MAX) * (xmax - xmin) + xmin;
		double p_y = rand() / double(RAND_MAX) * (ymax - ymin) + ymin;

		if (CheckValid(cgal_polys, p_x, p_y))
		{
			if (success_iter == 0)
			{
				dis_success = CGAL_2D_Distance_Point_Polygons(Vector2d(p_x, p_y), polys);
				inner_vec[0] = p_x;
				inner_vec[1] = p_y;
			}
			else
			{
				//double distance = CGAL_2D_Distance_Point_Polygon(Vector2d(p_x, p_y), polys);
				double distance = p_y;
				if (distance > dis_success)
				{
					dis_success = distance;
					inner_vec[0] = p_x;
					inner_vec[1] = p_y;
				}
			}
			success_iter++;
		}

		if (inter > 10000 || success_iter > 50)
		{
			break;
		}
		inter++;
	}

	return success_iter != 0;
}

/********************************************************
* Function name: CGAL_Identify_Polycut_Extend
* Description: Identifies and extends the boundary of a polygon cut. Given a polygon, a starting point (s), and an ending point (e), this function calculates the extended points (ns and ne) for the cut.
* Parameters:
* @polygon: A vector of 2D points representing the vertices of the polygon.
* @s: The starting point of the cut.
* @e: The ending point of the cut.
* @ns: Output parameter for the starting point of the extended cut.
* @ne: Output parameter for the ending point of the extended cut.
* Return: A boolean indicating if the operation was successful.
*********************************************************/
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

/********************************************************
* Function name: CGAL_Identify_Polycut_ExtendOld
* Description: An alternative implementation to identify and extend the boundary of a polygon cut. Given a polygon, a starting point (s), and an ending point (e), this function calculates the extended points (ns and ne) for the cut.
* Parameters:
* @polygon: A vector of 2D points representing the vertices of the polygon.
* @s: The starting point of the cut.
* @e: The ending point of the cut.
* @ns: Output parameter for the starting point of the extended cut.
* @ne: Output parameter for the ending point of the extended cut.
* Return: A boolean indicating if the operation was successful.
*********************************************************/
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

/********************************************************
* Function name: CGAL_Identify_Polycut
* Description: Identifies the intersection points of a cutting line with a polygon's boundary.
* Parameters:
* @polygon: A vector of 2D points representing the vertices of the polygon.
* @cutLine: A vector of 2D points representing the vertices of the cutting line.
* @result: Output parameter, a vector of pairs of boolean values indicating whether the cutting line enters and exits the polygon at each segment. The vector's size is one less than the size of cutLine.
* Return: A boolean indicating if the operation was successful.
*********************************************************/
extern "C" PPGL_EXPORT bool CGAL_Identify_Polycut(const Vector2d1 &polygon,const Vector2d1 &cutLine, VectorPB1 &result) {
    // N-1 edges, default 0 -> cant be fabricated
    result = VectorPB1(cutLine.size() - 1, std::make_pair<bool, bool>(false, false));

    Polygon_2 poly =Polygon2D(polygon);

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

/********************************************************
* Function name: GetIndex
* Description: Searches for an edge in the vector of GridEdge structures that matches the specified start and end points.
* Parameters:
* @grid_edges: A vector of GridEdge structures representing edges in a grid.
* @i_0: A pair of integers representing the start point coordinates.
* @i_1: A pair of integers representing the end point coordinates.
* Return: The index of the edge that matches the specified start and end points, or -1 if no matching edge is found.
*********************************************************/
int GetIndex(std::vector<GridEdge>& grid_edges, std::pair<int,int> i_0, std::pair<int,int> i_1)
{
	int s_i = i_0.first;
	int s_j = i_0.second;
	int e_i = i_1.first;
	int e_j = i_1.second;
	for (int i = 0; i < grid_edges.size(); i++) {
		if (grid_edges[i].s_i == s_i && grid_edges[i].s_j == s_j && grid_edges[i].e_i == e_i && grid_edges[i].e_j == e_j) {
			return i;
		}
		if (grid_edges[i].s_i == e_i && grid_edges[i].s_j == e_j && grid_edges[i].e_i == s_i && grid_edges[i].e_j == s_j) {
			return i;
		}
	}
	return -1;
}

/********************************************************
* Function name: Decomposition_Mapping
* Description: Decomposes a boundary represented by a vector of coordinate pairs into multiple parts and stores them in separate vectors.
* Parameters:
* @one_boundary: A vector of coordinate pairs representing a boundary.
* @boundary_xs: A vector of vectors that will store the x-coordinates of the decomposed boundaries.
* @boundary_ys: A vector of vectors that will store the y-coordinates of the decomposed boundaries.
* @remove_b: A boolean flag indicating whether to remove duplicate entries from one_boundary.
*********************************************************/
void Decomposition_Mapping(VectorPI1& one_boundary,
	Vector1d2& boundary_xs, Vector1d2& boundary_ys, bool remove_b = false)
{
	//remove multi grid index
	if (remove_b) {
		VectorPI1 new_one_boundary(1, one_boundary[0]);
		for (int i = 0; i < one_boundary.size(); i++) {
			std::pair<int,int>& last = new_one_boundary[new_one_boundary.size() - 1];
			if (!(one_boundary[i].first == last.first && one_boundary[i].second == last.second)) {
				new_one_boundary.push_back(one_boundary[i]);
			}
		}
		VectorPI1().swap(one_boundary);
		one_boundary = new_one_boundary;
		VectorPI1().swap(new_one_boundary);
	}

	if (one_boundary.size() <= 2) return;

	//find repeating index
	bool run = false;
	int index_0 = -1;
	int index_1 = -1;
	for (int i = 0; i < one_boundary.size() - 1 && !run; i++) {
		std::pair<int,int>& current = one_boundary[i];
		for (int j = i + 1; j < one_boundary.size() && !run; j++) {
			std::pair<int,int>& next = one_boundary[j];
			if (current.first == next.first && current.second == next.second) {
				index_0 = i;
				index_1 = j;
				run = true;
			}
		}
	}

	if (!run) {
		Vector1d1 xs;
		Vector1d1 ys;
		for (int i = 0; i < one_boundary.size(); i++) {
			xs.push_back(one_boundary[i].first - 1.0);
			ys.push_back(one_boundary[i].second - 1.0);
		}
		boundary_xs.push_back(xs);
		boundary_ys.push_back(ys);
		Vector1d1().swap(xs);
		Vector1d1().swap(ys);
	}
	else {
		VectorPI1 one_boundary_0;
		VectorPI1 one_boundary_1;
		for (int i = index_0; i < index_1; i++) {
			one_boundary_0.push_back(one_boundary[i]);
		}
		for (int i = index_1; i < one_boundary.size(); i++) {
			one_boundary_1.push_back(one_boundary[i]);
		}
		for (int i = 0; i < index_0; i++) {
			one_boundary_1.push_back(one_boundary[i]);
		}
		Decomposition_Mapping(one_boundary_0, boundary_xs, boundary_ys);
		Decomposition_Mapping(one_boundary_1, boundary_xs, boundary_ys);
	}
}

/********************************************************
* Function name: CGAL_Image_Grid_Decomposition_C1
* Description: Decomposes an image represented as a grid into boundaries and stores the x and y coordinates of these boundaries.
* Parameters:
* @image: A 2D vector representing the image grid.
* @boundary_xs: A vector of vectors that will store the x-coordinates of the detected boundaries.
* @boundary_ys: A vector of vectors that will store the y-coordinates of the detected boundaries.
*********************************************************/
extern "C" PPGL_EXPORT void CGAL_Image_Grid_Decomposition_C1(Vector1i2&image, Vector1d2&boundary_xs, Vector1d2&boundary_ys)
{
	//refine  lables
	Vector1i2 grid;
	for (int i = 0; i < image.size() + 2; i++) {
		Vector1i1 one_grid_raw;
		for (int j = 0; j < image[0].size() + 2; j++) {
			if (i >= 1 && i < image.size() + 1 && j >= 1 && j < image.size() + 1) {
				one_grid_raw.push_back(image[i - 1][j - 1]);
			}
			else {
				one_grid_raw.push_back(0);
			}
		}
		grid.push_back(one_grid_raw);
	}


	//std::ofstream lable_file("D:\\123.lable");
	//for (int j = 0; j < grid.size(); j++)
	//{
	//	for (int k = 0; k < grid[j].size(); k++)
	//	{
	//		lable_file << grid[j][k] << " ";
	//	}
	//	lable_file << "" << std::endl;
	//}
	//lable_file.clear();
	//lable_file.close();


	//get grid edges
	std::vector<GridEdge> grid_edges;
	for (int i = 0; i < grid.size(); i++) {
		for (int j = 0; j < grid[i].size(); j++) {
			if (i + 1 < grid.size()) {
				if (grid[i][j] != grid[i + 1][j]) {
					grid_edges.push_back(GridEdge(i, j, i + 1, j, 0));
				}
			}
			if (j + 1 < grid[i].size()) {
				if (grid[i][j] != grid[i][j + 1]) {
					grid_edges.push_back(GridEdge(i, j, i, j + 1, 1));
				}
			}
		}
	}

	//connect cut edges
	//std::pair<int,int> 1 2
	//std::pair<int,int> s e
	//std::pair<int,int> 3 4
	Vector1i2 grid_relations;
	Vector1b1 grid_edges_used;
	for (int i = 0; i < grid_edges.size(); i++) {
		grid_edges_used.push_back(false);
		std::pair<int,int> index_1;
		std::pair<int,int> index_2;
		std::pair<int,int> index_3;
		std::pair<int,int> index_4;
		std::pair<int,int> index_s(grid_edges[i].s_i, grid_edges[i].s_j);
		std::pair<int,int> index_e(grid_edges[i].e_i, grid_edges[i].e_j);

		if (grid_edges[i].type == 0) {
			index_1.first = grid_edges[i].s_i;
			index_1.second = grid_edges[i].s_j - 1;
			index_2.first = grid_edges[i].e_i;
			index_2.second = grid_edges[i].e_j - 1;
			index_3.first = grid_edges[i].s_i;
			index_3.second = grid_edges[i].s_j + 1;
			index_4.first = grid_edges[i].e_i;
			index_4.second = grid_edges[i].e_j + 1;
		}

		if (grid_edges[i].type == 1) {
			index_1.first = grid_edges[i].s_i - 1;
			index_1.second = grid_edges[i].s_j;
			index_2.first = grid_edges[i].e_i - 1;
			index_2.second = grid_edges[i].e_j;
			index_3.first = grid_edges[i].s_i + 1;
			index_3.second = grid_edges[i].s_j;
			index_4.first = grid_edges[i].e_i + 1;
			index_4.second = grid_edges[i].e_j;
		}

		int lable_0 = GetIndex(grid_edges, index_1, index_2);
		int lable_1 = GetIndex(grid_edges, index_1, index_s);
		int lable_2 = GetIndex(grid_edges, index_2, index_e);
		int lable_3 = GetIndex(grid_edges, index_3, index_4);
		int lable_4 = GetIndex(grid_edges, index_3, index_s);
		int lable_5 = GetIndex(grid_edges, index_4, index_e);

		Vector1i1 gr;
		if (lable_1 >= 0 && (grid[index_s.first][index_s.second] == 0 || (grid[index_s.first][index_s.second] == 1 && grid[index_2.first][index_2.second] == 0))) {
			gr.push_back(lable_1);
		}
		if (lable_2 >= 0 && (grid[index_e.first][index_e.second] == 0 || (grid[index_e.first][index_e.second] == 1 && grid[index_1.first][index_1.second] == 0))) {
			gr.push_back(lable_2);
		}
		if (lable_4 >= 0 && (grid[index_s.first][index_s.second] == 0 || (grid[index_s.first][index_s.second] == 1 && grid[index_4.first][index_4.second] == 0))) {
			gr.push_back(lable_4);
		}
		if (lable_5 >= 0 && (grid[index_e.first][index_e.second] == 0 || (grid[index_e.first][index_e.second] == 1 && grid[index_3.first][index_3.second] == 0))) {
			gr.push_back(lable_5);
		}
		if (lable_0 >= 0 && lable_1 < 0 && lable_2 < 0) {
			gr.push_back(lable_0);
		}
		if (lable_3 >= 0 && lable_4 < 0 && lable_5 < 0) {
			gr.push_back(lable_3);
		}
		grid_relations.push_back(gr);
		//if ((grid_edges[i].s_i == 11+1 && grid_edges[i].s_j == 15+1) ||
		//	(grid_edges[i].e_i == 11+1 && grid_edges[i].e_j == 15+1))
		//{
		//	std::cout << "Edge: (" << grid_edges[i].s_i - 1 << ", " << grid_edges[i].s_j - 1 << ")_("
		//		<< grid_edges[i].e_i - 1 << ", " << grid_edges[i].e_j - 1 << ")" << std::endl;
		//	for (int j = 0; j < gr.size(); j++)
		//	{
		//		int index = gr[j];
		//		std::cout << "(" << grid_edges[index].s_i - 1 << ", " << grid_edges[index].s_j - 1 << ")_("
		//			<< grid_edges[index].e_i - 1 << ", " << grid_edges[index].e_j - 1 << ")" << std::endl;
		//	}
		//}
	}

	//std::vector<GridEdgeRelation> grid_relations;
	//Vector1b1 grid_edges_used;
	//for (int i = 0; i < grid_edges.size(); i++)
	Vector1i2 grid_boundaries;
	while (true) {
		int start_edge_index = -1;

		for (int i = 0; i < grid_edges_used.size(); i++) {
			if (!grid_edges_used[i]) {
				start_edge_index = i;
				break;
			}
		}
		if (start_edge_index < 0) break;
		grid_edges_used[start_edge_index] = true;
		Vector1i1 one_boundary(1, start_edge_index);
		while (true) {
			int next_edge_index = -1;
			for (int i = 0; i < grid_relations[start_edge_index].size(); i++) {
				if (!grid_edges_used[grid_relations[start_edge_index][i]]) {
					next_edge_index = grid_relations[start_edge_index][i];
					break;
				}
			}
			if (next_edge_index < 0) break;
			grid_edges_used[next_edge_index] = true;
			one_boundary.push_back(next_edge_index);
			start_edge_index = next_edge_index;
		}
		grid_boundaries.push_back(one_boundary);
	}

	//transform back

#if 0
	{
		for (int i = 0; i < grid_boundaries.size(); i++) {
			Vector1d1 xs;
			Vector1d1 ys;
			for (int j = 0; j < grid_boundaries[i].size(); j++) {
				std::pair<int,int> index_s(grid_edges[grid_boundaries[i][j]].s_i, grid_edges[grid_boundaries[i][j]].s_j);
				std::pair<int,int> index_e(grid_edges[grid_boundaries[i][j]].e_i, grid_edges[grid_boundaries[i][j]].e_j);
				xs.push_back((index_s.first + index_e.first) / 2.0 - 1.0);
				ys.push_back((index_s.second + index_e.second) / 2.0 - 1.0);
			}
			boundary_xs.push_back(xs);
			boundary_ys.push_back(ys);
		}
	}
#else
	{
		for (int i = 0; i < grid_boundaries.size(); i++) {
			VectorPI1 one_boundary;
			for (int j = 0; j < grid_boundaries[i].size(); j++) {
				std::pair<int,int> index_s(grid_edges[grid_boundaries[i][j]].s_i, grid_edges[grid_boundaries[i][j]].s_j);
				std::pair<int,int> index_e(grid_edges[grid_boundaries[i][j]].e_i, grid_edges[grid_boundaries[i][j]].e_j);
				if (grid[index_s.first][index_s.second] == 1) {
					one_boundary.push_back(index_s);
				}
				else {
					one_boundary.push_back(index_e);
				}
			}
			Decomposition_Mapping(one_boundary, boundary_xs, boundary_ys, true);
			VectorPI1().swap(one_boundary);
		}
	}
#endif

	Vector1i2().swap(grid);
	std::vector<GridEdge>().swap(grid_edges);
	Vector1i2().swap(grid_relations);
	Vector1b1().swap(grid_edges_used);
	Vector1i2().swap(grid_boundaries);
}


/********************************************************
* Function name: CGAL_Image_Grid_Decomposition_Conservative_C1
* Description: Decomposes an image represented as a grid into conservative boundaries and stores the x and y coordinates of these boundaries.
* Parameters:
* @image: A 2D vector representing the image grid.
* @boundary_xs: A vector of vectors that will store the x-coordinates of the detected boundaries.
* @boundary_ys: A vector of vectors that will store the y-coordinates of the detected boundaries.
*********************************************************/

extern "C" PPGL_EXPORT void CGAL_Image_Grid_Decomposition_Conservative_C1(Vector1i2&image, Vector1d2&boundary_xs, Vector1d2&boundary_ys)
{
	VectorPI2 boundaries;

	//#pragma omp parallel 
	{
		//2
		for (int i = 0; i < image.size(); i++)
		{
			for (int j = 0; j < image[i].size(); j++)
			{
				if (image[i][j] == 1)
				{
					bool run = true;
					for (int k = -1; k < 2 && run; k++)
					{
						for (int m = -1; m < 2 && run; m++)
						{
							int index_0 = i + k;
							int index_1 = j + m;
							if (index_0 < 0 || index_1 < 0 || index_0 >= image.size() || index_1 >= image[i].size() || image[index_0][index_1] == 0)
							{
								image[i][j] = 2;
								run = false;
							}
						}
					}
				}
			}
		}

		//3
		VectorPI1 lables_3;
		Vector1b1 lables_3_used;
		for (int i = 0; i < image.size(); i++)
		{
			for (int j = 0; j < image[i].size(); j++)
			{
				if (image[i][j] == 2)
				{
					bool run = true;
					for (int k = -1; k < 2 && run; k++)
					{
						for (int m = -1; m < 2 && run; m++)
						{
							int index_0 = i + k;
							int index_1 = j + m;
							if (!(index_0 < 0 || index_1 < 0 || index_0 >= image.size() || index_1 >= image[i].size()) && image[index_0][index_1] == 1)
							{
								image[i][j] = 3 + lables_3.size();
								lables_3.push_back(std::pair<int,int>(i, j));
								lables_3_used.push_back(true);
								run = false;
							}
						}
					}

				}
			}
		}

		//order
		while (true)
		{
			//start_index
			std::pair<int, int> start_index(-1, -1);
			for (int i = 0; i < lables_3_used.size(); i++)
			{
				if (lables_3_used[i])
				{
					start_index.first = lables_3[i].first;
					start_index.second = lables_3[i].second;
					break;
				}
			}

			if (start_index.first < 0 || start_index.second < 0) break;

			VectorPI1 one_boundary;
			one_boundary.push_back(start_index);
			lables_3_used[image[start_index.first][start_index.second] - 3] = false;
			image[start_index.first][start_index.second] = 2;


			while (true)
			{
				std::pair<int,int> next_index(-1, -1);

				//next_index
				bool run = true;
				for (int k = -1; k < 2 && run; k++)
				{
					for (int m = -1; m < 2 && run; m++)
					{
						if (abs(k) + abs(m) == 1)
						{
							int index_0 = start_index.first + k;
							int index_1 = start_index.second + m;
							if (!(index_0 < 0 || index_1 < 0 || index_0 >= image.size() || index_1 >= image[0].size()))
							{
								if (image[index_0][index_1] >= 3)
								{
									next_index.first = index_0;
									next_index.second = index_1;
									run = false;
								}
							}
						}

					}
				}

				if (next_index.first < 0 || next_index.second < 0) break;

				one_boundary.push_back(next_index);
				lables_3_used[image[next_index.first][next_index.second] - 3] = false;
				image[next_index.first][next_index.second] = 2;

				start_index.first = next_index.first;
				start_index.second = next_index.second;
			}

			if (abs(one_boundary[0].first - one_boundary[one_boundary.size() - 1].first) <= 1 && abs(one_boundary[0].second - one_boundary[one_boundary.size() - 1].second) <= 1 && one_boundary.size() >= 3)
				boundaries.push_back(one_boundary);
		}
	}

	for (int i = 0; i < boundaries.size(); i++)
	{
		Vector1d1 xs;
		Vector1d1 ys;
		for (int j = 0; j < boundaries[i].size(); j++)
		{
			xs.push_back(boundaries[i][j].first);
			ys.push_back(boundaries[i][j].second);
		}

		boundary_xs.push_back(xs);
		boundary_ys.push_back(ys);
	}

	VectorPI2().swap(boundaries);
}

/********************************************************
* Function name: CGAL_Image_Grid_Decomposition_C2
* Description: Decomposes an image represented as a grid into boundaries and stores the boundary point coordinates.
* Parameters:
* @image: A 2D vector representing the image grid.
* @boundaries: A vector of vectors that will store the boundary point coordinates (as pairs of x and y coordinates).
*********************************************************/
extern "C" PPGL_EXPORT void CGAL_Image_Grid_Decomposition_C2(Vector1i2&image, Vector2d2 & boundaries)
{
	Vector1d2 xs;
	Vector1d2 ys;

	CGAL_Image_Grid_Decomposition_C1(image, xs, ys);

	for (int i = 0; i < xs.size(); i++)
	{
		Vector2d1 one_boundary;
		for (int j = 0; j < xs[i].size(); j++)
		{
			one_boundary.push_back(Vector2d(xs[i][j], ys[i][j]));
		}
		boundaries.push_back(one_boundary);
	}
}

/********************************************************
* Function name: CGAL_Image_Grid_Decomposition_Conservative_C2
* Description: Decomposes an image represented as a grid using conservative method into boundaries and stores the boundary point coordinates.
* Parameters:
* @image: A 2D vector representing the image grid.
* @boundaries: A vector of vectors that will store the boundary point coordinates (as pairs of x and y coordinates).
*********************************************************/
extern "C" PPGL_EXPORT void CGAL_Image_Grid_Decomposition_Conservative_C2(Vector1i2&image, Vector2d2 & boundaries)
{
	Vector1d2 xs;
	Vector1d2 ys;

	CGAL_Image_Grid_Decomposition_Conservative_C1(image, xs, ys);

	for (int i = 0; i < xs.size(); i++)
	{
		Vector2d1 one_boundary;
		for (int j = 0; j < xs[i].size(); j++)
		{
			one_boundary.push_back(Vector2d(xs[i][j], ys[i][j]));
		}
		boundaries.push_back(one_boundary);
	}
}
