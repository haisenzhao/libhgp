#include "libhgp.h"

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

