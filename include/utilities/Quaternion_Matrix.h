#ifndef QUATERNION_MATRIX_H
#define QUATERNION_MATRIX_H

#include <yafray_constants.h>

#include <core_api/vector3d.h>

__BEGIN_YAFRAY


class Matrix_3x3
{
    public:
    vector3d_t operator* (vector3d_t const& v) const
    {
        vector3d_t res(0.0f);

        res.x = _m[0] * v.x + _m[1] * v.y + _m[2] * v.z;
        res.y = _m[3] * v.x + _m[4] * v.y + _m[5] * v.z;
        res.z = _m[6] * v.x + _m[7] * v.y + _m[8] * v.z;

        return res;
    }

    Matrix_3x3 operator* (Matrix_3x3 const& m) const
    {
        Matrix_3x3 res;

        res._m[0] = _m[0] * m._m[0] + _m[1] * m._m[3] + _m[2] * m._m[6];
        res._m[1] = _m[0] * m._m[1] + _m[1] * m._m[4] + _m[2] * m._m[7];
        res._m[2] = _m[0] * m._m[2] + _m[1] * m._m[5] + _m[2] * m._m[8];

        res._m[3] = _m[3] * m._m[0] + _m[4] * m._m[3] + _m[5] * m._m[6];
        res._m[4] = _m[3] * m._m[1] + _m[4] * m._m[4] + _m[5] * m._m[7];
        res._m[5] = _m[3] * m._m[2] + _m[4] * m._m[5] + _m[5] * m._m[8];

        res._m[6] = _m[6] * m._m[0] + _m[7] * m._m[3] + _m[8] * m._m[6];
        res._m[7] = _m[6] * m._m[1] + _m[7] * m._m[4] + _m[8] * m._m[7];
        res._m[8] = _m[6] * m._m[2] + _m[7] * m._m[5] + _m[8] * m._m[8];

        return res;
    }

    float _m[9]; // row order, i.e. 0, 1, 2 = first row etc.
    // private:
};

class Quaternion
{
    public:
    Quaternion(vector3d_t const& axis, float const angle)
    {
        float s = std::sin(angle / 2.0f);
        x = axis.x * s;
        y = axis.y * s;
        z = axis.z * s;
        w = std::cos(angle / 2.0f); // real part
    }

    Matrix_3x3 to_matrix() const
    {
        Matrix_3x3 res;
        res._m[0] = w * w + x * x - y * y - z * z;
        res._m[1] = 2 * x * y - 2 * w * z;
        res._m[2] = 2 * x * z + 2 * w * y;

        res._m[3] = 2 * x* y + 2 * w * z;
        res._m[4] = w * w - x * x + y * y - z * z;
        res._m[5] = 2 * y * z - 2 * w * x;

        res._m[6] = 2 * x * z - 2 * w * y;
        res._m[7] = 2 * y * z + 2 * w * x;
        res._m[8] = w * w - x * x - y * y + z * z;

        return res;
    }

    private:
    float x, y, z, w;
};

__END_YAFRAY

#endif // QUATERNION_MATRIX_H
