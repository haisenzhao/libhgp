#include "libhgp.h"






/********************************************************
* Function name: HGP_3D_Distance_Point_Segment
* Description: Calculate the Euclidean distance between a 3D point and a line segment.
* Parameters:
* @p: 3D point for which the distance is to be calculated.
* @s_s: Start point of the line segment.
* @s_e: End point of the line segment.
* Return: Euclidean distance between the 3D point @p and the line segment defined by @s_s and @s_e.
*********************************************************/
extern "C" LIBHGP_EXPORT double HGP_3D_Distance_Point_Segment(const Vector3d & p, const Vector3d & s_s, const Vector3d & s_e) {
	return sqrt((double)CGAL::squared_distance(VectorPoint3d(p), Segment_3(VectorPoint3d(s_s), VectorPoint3d(s_e))));
}
/********************************************************
* Function name: HGP_3D_Plane_Fitting
* Description: Fit a plane to a set of 3D points using linear least squares fitting.
* Parameters:
* @points: A vector of 3D points to fit the plane to.
* @plane_p: The calculated point on the fitted plane.
* @plane_n: The normal vector of the fitted plane.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void HGP_3D_Plane_Fitting(const Vector3d1 & points, Vector3d & plane_p, Vector3d & plane_n) {
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
/********************************************************
* Function name: HGP_3D_Plane_Point_Projection
* Description: Project a 3D point onto a plane defined by a point on the plane and its normal vector.
* Parameters:
* @plane_p: A point on the plane.
* @plane_n: The normal vector of the plane.
* @p: The 3D point to be projected onto the plane.
* @result: The projected point on the plane.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void HGP_3D_Plane_Point_Projection(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d & p, Vector3d & result)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	result = PointVector3d(plane.projection(VectorPoint3d(p)));
}
/********************************************************
* Function name: HGP_3D_Plane_Points_Projection
* Description: Project a list of 3D points onto a plane defined by a point on the plane and its normal vector.
* Parameters:
* @plane_p: A point on the plane.
* @plane_n: The normal vector of the plane.
* @points: A list of 3D points to be projected onto the plane.
* @project_points: The list of projected points on the plane.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void HGP_3D_Plane_Points_Projection
(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d1 & points, Vector3d1 & project_points)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	for (int i = 0; i < points.size(); i++)
		project_points.push_back(PointVector3d(plane.projection(VectorPoint3d(points[i]))));
}
/********************************************************
* Function name: HGP_3D_Plane_3D_to_2D_Point
* Description: Project a 3D point onto a 2D plane defined by a point on the plane and its normal vector.
* Parameters:
* @plane_p: A point on the 2D plane.
* @plane_n: The normal vector of the 2D plane.
* @point_3d: The 3D point to be projected onto the 2D plane.
* @result: The 2D point result of projecting 'point_3d' onto the 2D plane.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void
HGP_3D_Plane_3D_to_2D_Point(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d & point_3d, Vector2d & result)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	result = PointVector2d(point_to_2d(VectorPoint3d(point_3d), plane));
}

/********************************************************
* Function name: HGP_3D_Plane_2D_to_3D_Point
* Description: Project a 2D point onto a 3D plane defined by a point on the plane and its normal vector.
* Parameters:
* @plane_p: A point on the 3D plane.
* @plane_n: The normal vector of the 3D plane.
* @points_2d: The 2D point to be projected onto the 3D plane.
* @result: The 3D point result of projecting 'points_2d' onto the 3D plane.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void
HGP_3D_Plane_2D_to_3D_Point(const Vector3d & plane_p, const Vector3d & plane_n, const Vector2d & points_2d, Vector3d & result)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	result = PointVector3d(point_to_3d(VectorPoint2d(points_2d), plane));
}
/********************************************************
* Function name: HGP_3D_Plane_3D_to_2D_Points
* Description: Project multiple 3D points onto a 2D plane defined by a point on the plane and its normal vector.
* Parameters:
* @plane_p: A point on the 3D plane.
* @plane_n: The normal vector of the 3D plane.
* @points_3d: A vector of 3D points to be projected onto the 2D plane.
* @points_2d: A vector to store the resulting 2D points after projection.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void HGP_3D_Plane_3D_to_2D_Points
(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d1 & points_3d, Vector2d1 & points_2d)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	for (int i = 0; i < points_3d.size(); i++)
		points_2d.push_back(PointVector2d(plane.to_2d(VectorPoint3d(points_3d[i]))));
}
/********************************************************
* Function name: HGP_3D_Plane_3Ds_to_2Ds_Points
* Description: Project multiple sets of 3D points onto a 2D plane defined by a point on the plane and its normal vector.
* Parameters:
* @plane_p: A point on the 3D plane.
* @plane_n: The normal vector of the 3D plane.
* @points_3d: A vector of vectors of 3D points to be projected onto the 2D plane.
* @points_2d: A vector of vectors to store the resulting 2D points after projection for each set of 3D points.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void HGP_3D_Plane_3Ds_to_2Ds_Points(const Vector3d & plane_p, const Vector3d & plane_n, const Vector3d2 & points_3d, Vector2d2 & points_2d)
{
	for (auto& p3ds : points_3d)
	{
		Vector2d1 points_2d_;
		HGP_3D_Plane_3D_to_2D_Points(plane_p, plane_n, p3ds, points_2d_);
		points_2d.push_back(points_2d_);
	}
}
/********************************************************
* Function name: HGP_3D_Plane_2D_to_3D_Points
* Description: Project a set of 2D points onto a 3D plane.
* Parameters:
* @plane_p: A point on the 3D plane where the 2D points will be projected.
* @plane_n: The normal vector of the 3D plane.
* @points_2d: A vector of 2D points to be projected onto the 3D plane.
* @points_3d: A vector to store the resulting 3D points after projection.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void HGP_3D_Plane_2D_to_3D_Points
(const Vector3d & plane_p, const Vector3d & plane_n, const Vector2d1 & points_2d, Vector3d1 & points_3d)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	for (int i = 0; i < points_2d.size(); i++)
		points_3d.push_back(PointVector3d(plane.to_3d(VectorPoint2d(points_2d[i]))));
}
/********************************************************
* Function name: HGP_3D_Distance_Point_Point
* Description: Calculate the 3D Euclidean distance between two points.
* Parameters:
* @v0: The first 3D point.
* @v1: The second 3D point.
* Return: The Euclidean distance between the two points.
*********************************************************/
extern "C" LIBHGP_EXPORT
double HGP_3D_Distance_Point_Point(const Vector3d & v0, const Vector3d & v1)
{
	return sqrt((double)CGAL::squared_distance(Point_3(v0[0], v0[1], v0[2]), Point_3(v1[0], v1[1], v1[2])));
}

/********************************************************
* Function name: HGP_3D_Projection_Point_Segment
* Description: Project a 3D point onto a line segment defined by two 3D points.
* Parameters:
* @p: The 3D point to be projected.
* @s_s: The start point of the line segment.
* @s_e: The end point of the line segment.
* Return: The projected point on the line segment.
*********************************************************/
extern "C" LIBHGP_EXPORT Vector3d HGP_3D_Projection_Point_Segment(const Vector3d & p, const Vector3d & s_s, const Vector3d & s_e) {
	Line_3 l(VectorPoint3d(s_s), VectorPoint3d(s_e));
	auto m_p = PointVector3d(l.projection(VectorPoint3d(p)));
	double d_m_s = HGP_3D_Distance_Point_Point(m_p, s_s);
	double d_m_e = HGP_3D_Distance_Point_Point(m_p, s_e);
	double d_s_e = HGP_3D_Distance_Point_Point(s_s, s_e);

	if (d_m_s >= d_s_e)
		return s_e;
	if (d_m_e >= d_s_e)
		return s_s;
	return m_p;
}
/********************************************************
* Function name: HGP_3D_Distance_Point_Segment_Ref
* Description: Calculate the 3D Euclidean distance between a point and a line segment.
* Parameters:
* @v: The 3D point for which the distance is calculated.
* @s_0: The start point of the line segment.
* @s_1: The end point of the line segment.
* Return: The Euclidean distance between the point and the line segment.
*********************************************************/
extern "C" LIBHGP_EXPORT double HGP_3D_Distance_Point_Segment_Ref(const Vector3d & v, const Vector3d & s_0,
	const Vector3d & s_1) {
	return sqrt((double)CGAL::squared_distance(Point_3(v[0], v[1], v[2]),
		Segment_3(Point_3(s_0[0], s_0[1], s_0[2]),
			Point_3(s_1[0], s_1[1], s_1[2]))));
}
/********************************************************
* Function name: HGP_3D_Distance_Point_Polygon
* Description: Calculate the 3D Euclidean distance between a point and a 3D polygon.
* Parameters:
* @py: A vector of 3D points representing the vertices of the polygon.
* @p: The 3D point for which the distance is calculated.
* Return: The minimum Euclidean distance between the point and the polygon.
*********************************************************/
extern "C" LIBHGP_EXPORT double HGP_3D_Distance_Point_Polygon(const Vector3d1 & py, const Vector3d & p) {
	double distance = 1000000000000.0;
	for (int i = 0; i < py.size(); i++)
		distance = std::min(distance, HGP_3D_Distance_Point_Segment_Ref(p, py[i], py[(i + 1) % py.size()]));
	return distance;
}




//IO mesh
/***************************************************************************************************/

//void HGP_Load_Obj(const char* path, Vector1d1 &coords,Vector1i1 &tris) {
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
/********************************************************
* Function name: HGP_3D_Distance_Point_Line
* Description: Calculate the Euclidean distance between a 3D point and a 3D line defined by two points.
* Parameters:
* @p: The 3D point for which the distance is to be calculated.
* @l_s: The start point of the 3D line.
* @l_e: The end point of the 3D line.
* Return: The Euclidean distance between the point 'p' and the line defined by 'l_s' and 'l_e'.
*********************************************************/
extern "C" LIBHGP_EXPORT double HGP_3D_Distance_Point_Line(const Vector3d & p, const Vector3d & l_s, const Vector3d & l_e)
{
	return sqrt((double)CGAL::squared_distance(VectorPoint3d(p), Line_3(VectorPoint3d(l_s), VectorPoint3d(l_e))));
}
/********************************************************
* Function name: HGP_3D_Projection_Point_Line
* Description: Project a 3D point onto a 3D line defined by two points and return the projected point.
* Parameters:
* @p: The 3D point to be projected onto the line.
* @l_s: The start point of the 3D line.
* @l_e: The end point of the 3D line.
* Return: The 3D point representing the projection of 'p' onto the line defined by 'l_s' and 'l_e'.
*********************************************************/
extern "C" LIBHGP_EXPORT Vector3d HGP_3D_Projection_Point_Line(const Vector3d & p, const Vector3d & l_s, const Vector3d & l_e)
{
	Line_3 l(VectorPoint3d(l_s), VectorPoint3d(l_e));
	Point_3 o_p = l.projection(VectorPoint3d(p));
	return PointVector3d(o_p);
}
/********************************************************
* Function name: HGP_3D_Distance_Segment_Segment
* Description: Calculate the Euclidean distance between two 3D line segments defined by their start and end points.
* Parameters:
* @s_0_s: The start point of the first 3D line segment.
* @s_0_e: The end point of the first 3D line segment.
* @s_1_s: The start point of the second 3D line segment.
* @s_1_e: The end point of the second 3D line segment.
* Return: The Euclidean distance between the two 3D line segments.
*********************************************************/
extern "C" LIBHGP_EXPORT double HGP_3D_Distance_Segment_Segment(const Vector3d & s_0_s, const Vector3d & s_0_e, const Vector3d & s_1_s, const Vector3d & s_1_e)
{
	return sqrt((double)CGAL::squared_distance(Segment_3(VectorPoint3d(s_0_s), VectorPoint3d(s_0_e)), Segment_3(VectorPoint3d(s_1_s), VectorPoint3d(s_1_e))));
}
/********************************************************
* Function name: HGP_3D_Distance_Point_Plane
* Description: Calculate the Euclidean distance between a 3D point and a 3D plane defined by a point on the plane and its normal vector.
* Parameters:
* @v: The 3D point for which the distance to the plane is calculated.
* @plane_p: A point on the 3D plane.
* @plane_n: The normal vector of the 3D plane.
* Return: The Euclidean distance between the 3D point 'v' and the 3D plane defined by 'plane_p' and 'plane_n'.
*********************************************************/
extern "C" LIBHGP_EXPORT double HGP_3D_Distance_Point_Plane(const Vector3d & v, const Vector3d & plane_p, const Vector3d & plane_n)
{
	return sqrt((double)CGAL::squared_distance(VectorPoint3d(v), Plane_3(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]))));
}


/********************************************************
* Function name: HGP_3D_Intersection_Segment_Line
* Description: Check if a 3D line segment intersects a 3D line, and if so, calculate the intersection point.
* Parameters:
* @s_s: The starting point of the 3D line segment.
* @s_e: The ending point of the 3D line segment.
* @l_s: The starting point of the 3D line.
* @l_e: The ending point of the 3D line.
* @inter: A vector to store the intersection point if the segment and line intersect.
* Return: True if the segment and line intersect, false otherwise.
*********************************************************/
extern "C" LIBHGP_EXPORT bool HGP_3D_Intersection_Segment_Line(const Vector3d & s_s, const Vector3d & s_e, const Vector3d & l_s, const Vector3d & l_e, Vector3d & inter)
{
	double d = HGP_3D_Distance_Point_Point(s_s, s_e);

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
/********************************************************
* Function name: HGP_3D_Intersection_Segment_Segment
* Description: Check if two 3D line segments intersect, and if so, calculate the intersection point.
* Parameters:
* @s_0_s: The starting point of the first 3D line segment.
* @s_0_e: The ending point of the first 3D line segment.
* @s_1_s: The starting point of the second 3D line segment.
* @s_1_e: The ending point of the second 3D line segment.
* @iter: A vector to store the intersection point if the segments intersect.
* Return: True if the two segments intersect, false otherwise.
*********************************************************/
extern "C" LIBHGP_EXPORT bool HGP_3D_Intersection_Segment_Segment(const Vector3d & s_0_s, const Vector3d & s_0_e, const Vector3d & s_1_s, const Vector3d & s_1_e, Vector3d & iter)
{
	double d0 = HGP_3D_Distance_Point_Point(s_0_s, s_0_e);
	double d1 = HGP_3D_Distance_Point_Point(s_1_s, s_1_e);

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
/********************************************************
* Function name: HGP_3D_Intersection_Segment_Plane
* Description: Check if a 3D line segment intersects with a 3D plane and, if so, calculate the intersection point.
* Parameters:
* @s_s: The starting point of the 3D line segment.
* @s_e: The ending point of the 3D line segment.
* @plane_p: A point on the 3D plane.
* @plane_n: The normal vector of the 3D plane.
* @inter: A vector to store the intersection point if the segment and plane intersect.
* Return: True if the segment and plane intersect, false otherwise.
*********************************************************/
extern "C" LIBHGP_EXPORT bool HGP_3D_Intersection_Segment_Plane(const Vector3d & s_s, const Vector3d & s_e, const Vector3d & plane_p, const Vector3d & plane_n, Vector3d & inter)
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
/********************************************************
* Function name: HGP_3D_Intersection_Line_Plane
* Description: Check if a 3D line intersects with a 3D plane and, if so, calculate the intersection point.
* Parameters:
* @l_s: The starting point of the 3D line.
* @l_e: The ending point of the 3D line.
* @plane_p: A point on the 3D plane.
* @plane_n: The normal vector of the 3D plane.
* @inter: A vector to store the intersection point if the line and plane intersect.
* Return: True if the line and plane intersect, false otherwise.
*********************************************************/
extern "C" LIBHGP_EXPORT bool HGP_3D_Intersection_Line_Plane(const Vector3d & l_s, const Vector3d & l_e, const Vector3d & plane_p, const Vector3d & plane_n, Vector3d & inter)
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

/********************************************************
* Function name: HGP_3D_Projection_Point_Plane_C1
* Description: Project a 3D point onto a 3D plane and return the projection point.
* Parameters:
* @p: The 3D point to be projected.
* @plane_p: A point on the 3D plane.
* @plane_n: The normal vector of the 3D plane.
* Return: A 3D vector representing the projection of 'p' onto the plane.
*********************************************************/
extern "C" LIBHGP_EXPORT Vector3d HGP_3D_Projection_Point_Plane_C1(const Vector3d & p, const Vector3d & plane_p, const Vector3d & plane_n)
{
	Plane_3 plane(Point_3(plane_p[0], plane_p[1], plane_p[2]), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	Point_3 project = plane.projection(Point_3(p[0], p[1], p[1]));
	Vector3d result(project[0], project[1], project[2]);
	return result;
}
/********************************************************
* Function name: HGP_3D_Projection_Point_Plane_C2
* Description: Project a 3D point onto a 3D plane defined by three points on the plane.
* Parameters:
* @p: The 3D point to be projected.
* @plane_p_0: The first point on the 3D plane.
* @plane_p_1: The second point on the 3D plane.
* @plane_p_2: The third point on the 3D plane.
* Return: A 3D vector representing the projection of 'p' onto the plane.
*********************************************************/
extern "C" LIBHGP_EXPORT Vector3d HGP_3D_Projection_Point_Plane_C2(const Vector3d & p, const Vector3d & plane_p_0, const Vector3d & plane_p_1, const Vector3d & plane_p_2)
{
	Point_3 p0 = VectorPoint3d(plane_p_0);
	Point_3 p1 = VectorPoint3d(plane_p_1);
	Point_3 p2 = VectorPoint3d(plane_p_2);

	Plane_3 plane(p1, CGAL::cross_product(p2 - p1, p0 - p1));
	Point_3 project = plane.projection(VectorPoint3d(p));
	Vector3d result(project[0], project[1], project[2]);
	return result;

}
/********************************************************
* Function name: HGP_3D_Projection_3D_Point_Plane_2D_C1
* Description: Project a 3D point onto a 3D plane and represent the result as a 2D point.
* Parameters:
* @p: The 3D point to be projected.
* @plane_p: A point on the 3D plane.
* @plane_n: The normal vector of the 3D plane.
* Return: A 2D vector representing the projection of 'p' onto the plane.
*********************************************************/
extern "C" LIBHGP_EXPORT Vector2d HGP_3D_Projection_3D_Point_Plane_2D_C1(const Vector3d & p, const Vector3d & plane_p, const Vector3d & plane_n)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	Point_2 r = point_to_2d(VectorPoint3d(p), plane);
	return Vector2d(r[0], r[1]);
}
/********************************************************
* Function name: HGP_3D_Projection_3D_Point_Plane_2D_C2
* Description: Project a 3D point onto a 3D plane defined by three points and represent the result as a 2D point.
* Parameters:
* @p: The 3D point to be projected.
* @plane_p_0: A 3D point on the plane.
* @plane_p_1: Another 3D point on the plane.
* @plane_p_2: A third 3D point on the plane.
* Return: A 2D vector representing the projection of 'p' onto the plane.
*********************************************************/
extern "C" LIBHGP_EXPORT Vector2d HGP_3D_Projection_3D_Point_Plane_2D_C2(const Vector3d & p, const Vector3d & plane_p_0, const Vector3d & plane_p_1, const Vector3d & plane_p_2)
{
	Point_3 p0 = VectorPoint3d(plane_p_0);
	Point_3 p1 = VectorPoint3d(plane_p_1);
	Point_3 p2 = VectorPoint3d(plane_p_2);

	Plane_3 plane(p1, CGAL::cross_product(p2 - p1, p0 - p1));
	Point_2 r = point_to_2d(VectorPoint3d(p), plane);
	return Vector2d(r[0], r[1]);
}
/********************************************************
* Function name: HGP_3D_Plane_ABCD
* Description: Compute the coefficients of a 3D plane in the form Ax + By + Cz + D = 0.
* Parameters:
* @plane_p: A 3D point on the plane.
* @plane_n: The normal vector of the 3D plane.
* @a: Reference to a double variable to store the coefficient 'A'.
* @b: Reference to a double variable to store the coefficient 'B'.
* @c: Reference to a double variable to store the coefficient 'C'.
* @d: Reference to a double variable to store the coefficient 'D'.
* Return: void
*********************************************************/
extern "C" LIBHGP_EXPORT void HGP_3D_Plane_ABCD(const Vector3d & plane_p, const Vector3d & plane_n, double& a, double& b, double& c, double& d)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	a = plane.a();
	b = plane.b();
	c = plane.c();
	d = plane.d();
}
/********************************************************
* Function name: HGP_3D_Plane_Base_1
* Description: Compute the first base vector of a 3D plane defined by a point on the plane and its normal vector.
* Parameters:
* @plane_p: A 3D point on the plane.
* @plane_n: The normal vector of the 3D plane.
* Return: Vector3d representing the first base vector of the plane.
*********************************************************/
extern "C" LIBHGP_EXPORT Vector3d HGP_3D_Plane_Base_1(const Vector3d & plane_p, const Vector3d & plane_n)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	Vector_3 v = plane.base1();
	return Vector3d(v[0], v[1], v[2]);
}



/********************************************************
* Function name: HGP_Face_Normal
* Description: Compute the normal vector of a point within a triangle by interpolating the normals of its vertices.
* Parameters:
* @source: The 3D point for which to compute the interpolated normal.
* @tri_0: The first vertex of the triangle.
* @tri_1: The second vertex of the triangle.
* @tri_2: The third vertex of the triangle.
* @normal_0: The normal vector at the first vertex.
* @normal_1: The normal vector at the second vertex.
* @normal_2: The normal vector at the third vertex.
* Return: Vector3d representing the interpolated normal at the 'source' point.
*********************************************************/
extern "C" LIBHGP_EXPORT Vector3d HGP_Face_Normal(const Vector3d & source, const Vector3d & tri_0, const Vector3d & tri_1, const Vector3d & tri_2, Vector3d & normal_0, Vector3d & normal_1, Vector3d & normal_2)
{
	double u, v, w;
	Functs::Barycentric(source, tri_0, tri_1, tri_2, u, v, w);
	return (double)u * normal_0 + (double)v * normal_1 + (double)w * normal_2;
}

