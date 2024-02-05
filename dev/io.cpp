#include "geom.h"

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

// Test function
// Input: 3d vector from glm and char* to print
extern "C" LIBHGP_EXPORT void HGP_Test_PGL(const Vector3d & n, const char* str, const char* char_)
{
    Functs::CerrLine(char_);
    Functs::CerrLine(str);
    Functs::MAssert("Test PGL ...");
}

extern "C" LIBHGP_EXPORT void HGP_Vector_Base(const Vector3d & n, Vector3d & result) {
    Plane_3 plane(Point_3(0.0, 0.0, 0.0), Vector_3(n[0], n[1], n[2]));
    Vector_3 v = plane.base1();
    result = Vector3d(v[0], v[1], v[2]);
}

//3d plane related
Vector3d HGP_3D_Plane_Base_1(Vector3d plane_p, Vector3d plane_n) {
    Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
    Vector_3 v = plane.base1();
    return Vector3d(v[0], v[1], v[2]);
}