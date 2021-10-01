#ifndef math_hpp
#define math_hpp


#include <vector>
#include <set>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <functional>
#include <queue>
#include <sstream>
#include <math.h>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>

typedef glm::highp_dvec2 Vector2d;
typedef glm::highp_dvec3 Vector3d;

namespace Math {
    static const double DOUBLE_EPSILON = 1.0E-05;
    static const double Math_PI = 3.14159265359;

    inline std::string IntString(int i) {
        std::stringstream ss;
        std::string str;
        ss << i;
        ss >> str;
        return str;
    }

    inline std::string DoubleString(double d) {
        std::stringstream ss;
        ss.precision(5);
        std::string str;
        ss << d;
        ss >> str;
        return str;
    }

    inline double GetLength(Vector3d v) {
        return length(v);
    }

    inline Vector3d RotationAxis(Vector3d p, double angle, Vector3d n) {
        glm::mat4 inputMatrix(0.0);
        inputMatrix[0][0] = p[0];
        inputMatrix[1][0] = p[1];
        inputMatrix[2][0] = p[2];
        inputMatrix[3][0] = 1.0;
        double u = n[0];
        double v = n[1];
        double w = n[2];

        glm::mat4 rotationMatrix;

        double L = (u * u + v * v + w * w);

        //angle = angle * M_PI / 180.0; //converting to radian value
        double u2 = u * u;
        double v2 = v * v;
        double w2 = w * w;

        rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
        rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
        rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
        rotationMatrix[0][3] = 0.0;

        rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
        rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
        rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
        rotationMatrix[1][3] = 0.0;

        rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
        rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
        rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
        rotationMatrix[2][3] = 0.0;

        rotationMatrix[3][0] = 0.0;
        rotationMatrix[3][1] = 0.0;
        rotationMatrix[3][2] = 0.0;
        rotationMatrix[3][3] = 1.0;

        double outputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 1; j++) {
                outputMatrix[i][j] = 0;
                for (int k = 0; k < 4; k++) {
                    outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
                }
            }
        }
        return Vector3d(outputMatrix[0][0], outputMatrix[0][1], outputMatrix[0][2]);
    }

    inline Vector3d GetCrossproduct(Vector3d &v1, Vector3d &v2) {
        return cross(v1, v2);
    }

    /*double GetDotProduct(Vector3d& v1, Vector3d& v2){
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }
*/
    inline Vector3d SetVectorLength(Vector3d &v, double length) {
        double l = GetLength(v);

        v[0] = v[0] / l * length;
        v[1] = v[1] / l * length;
        v[2] = v[2] / l * length;

        return v;
    }

    inline bool IsAlmostZero(double value) {
        return value < 10.0 * DOUBLE_EPSILON && value > -10.0 * DOUBLE_EPSILON;
    }

    inline bool IsAlmostZero_Double(double value, double EPSILON) {
        return value < 10.0 * DOUBLE_EPSILON && value > -10.0 * EPSILON;
    }

    inline double GetAngleBetween(Vector3d v1, Vector3d v2) {
        double d = dot(v1, v2) / (length(v1) * length(v2));

        if (IsAlmostZero(d - 1.0))
            return 0.0;

        if (IsAlmostZero(d + 1.0))
            return Math_PI;

        return glm::acos(dot(v1, v2) / (length(v1) * length(v2)));
    }

    inline double GetAngleBetween(const Vector2d &v1, const Vector2d &v2) {
        double d = dot(v1, v2) / (length(v1) * length(v2));

        if (IsAlmostZero(d - 1.0))
            return 0.0;

        if (IsAlmostZero(d + 1.0))
            return Math_PI;

        return glm::acos(dot(v1, v2) / (length(v1) * length(v2)));
    }

    inline Vector3d Vector3dBase(Vector3d v) {
        Vector3d n(1.0, 1.0, 1.0);
        if (!IsAlmostZero(v[0])) {
            n[0] = -(v[1] + v[2]) / v[0];
            return n;
        }
        if (!IsAlmostZero(v[1])) {
            n[1] = -(v[0] + v[2]) / v[1];
            return n;
        }
        if (!IsAlmostZero(v[2])) {
            n[2] = -(v[0] + v[1]) / v[2];
            return n;
        }
        return n;
    }
}
#endif
