#pragma once
#include "matrix.h"

namespace lina
{
    template<typename T>
    struct quaternion
    {
        union
        {
            struct
            {
                T x, y, z, w;
            };

            struct
            {
                vector<T, 3> v;
            };

            T components[4];
        };

        static quaternion<T> rotate_euler(const T pitch, const T yaw, const T roll);

        static quaternion<T> rotate_angle_axis(const vector<T, 3>& axis, const T angle);

        static quaternion<T> rotate_matrix(const matrix<T, 3, 3>& m);

        quaternion(const T ix, const T iy, const T iz, const T iw);

        quaternion(const vector<T, 3>& v, const T w);

        quaternion(const quaternion<T>& o);

        quaternion(const T* o);

        quaternion();

        T length_squared() const;

        T length() const;

        quaternion<T>& normalize();

        quaternion<T>& safe_normalize(const T tolerance = 0.001f);

        quaternion<T>& conjugate();

        quaternion<T>& inverse();

        matrix<T, 3, 3> matrix3x3() const;

        matrix<T, 4, 3> matrix4x3() const;

        matrix<T, 4, 4> matrix4x4() const;

        vector<T, 3> euler() const;

        quaternion<T>& operator+=(const quaternion<T>& rhs);

        quaternion<T>& operator-=(const quaternion<T>& rhs);

        quaternion<T>& operator*=(const quaternion<T>& rhs);

        quaternion<T>& operator*=(const T rhs);

        quaternion<T>& operator/=(const T rhs);

        quaternion<T>& operator=(const quaternion<T>& rhs);

        quaternion<T> operator-() const;

        T& operator[](const unsigned int index);

        const T& operator[](const unsigned int index) const;
    };

    template<typename T>
    T dot(const quaternion<T>& lhs, const quaternion<T>& rhs);

    template<typename T>
    T length(const quaternion<T>& rhs);

    template<typename T>
    T length_squared(const quaternion<T>& val);

    template<typename T>
    quaternion<T> normalize(const quaternion<T>& rhs);

    template<typename T>
    quaternion<T> safe_normalize(const quaternion<T>& rhs, const T tolerance = 0.001f);

    template<typename T>
    quaternion<T> conjugate(const quaternion<T>& rhs);

    template<typename T>
    quaternion<T> inverse(const quaternion<T>& rhs);

    template<typename T>
    quaternion<T> slerp(quaternion<T> q0, quaternion<T> q1, const T t);

    template<typename T>
    quaternion<T> operator+(const quaternion<T>& lhs, const quaternion<T>& rhs);

    template<typename T>
    quaternion<T> operator-(const quaternion<T>& lhs, const quaternion<T>& rhs);

    template<typename T>
    quaternion<T> operator*(const quaternion<T>& lhs, const quaternion<T>& rhs);

    template<typename T>
    vector<T, 3> operator*(const vector<T, 3>& lhs, const quaternion<T>& rhs);

    template<typename T>
    vector<T, 3> operator*(const quaternion<T>& lhs, const vector<T, 3>& rhs);

    template<typename T>
    quaternion<T> operator*(const quaternion<T>& lhs, const T rhs);

    template<typename T>
    quaternion<T> operator*(const T lhs, const quaternion<T>& rhs);

    template<typename T>
    quaternion<T> operator/(const quaternion<T>& lhs, const T rhs);

    template<typename T>
    quaternion<T> quaternion<T>::rotate_euler(const T pitch, const T yaw, const T roll)
    {
        T halfRoll = roll * 0.5f;
        T halfPitch = pitch * 0.5f;
        T halfYaw = yaw * 0.5f;
        T cosHalfRoll = cos(halfRoll);
        T cosHalfPitch = cos(halfPitch);
        T cosHalfYaw = cos(halfYaw);
        T sinHalfRoll = sin(halfRoll);
        T sinHalfPitch = sin(halfPitch);
        T sinHalfYaw = sin(halfYaw);

        return
        {
            (cosHalfYaw * sinHalfPitch * cosHalfRoll) + (sinHalfYaw * cosHalfPitch * sinHalfRoll),
            (sinHalfYaw * cosHalfPitch * cosHalfRoll) - (cosHalfYaw * sinHalfPitch * sinHalfRoll),
            (cosHalfYaw * cosHalfPitch * sinHalfRoll) - (sinHalfYaw * sinHalfPitch * cosHalfRoll),
            (cosHalfYaw * cosHalfPitch * cosHalfRoll) + (sinHalfYaw * sinHalfPitch * sinHalfRoll)
        };
    }

    template<typename T>
    quaternion<T> quaternion<T>::rotate_angle_axis(const vector<T, 3>& axis, const T angle)
    {
        T a = angle * 0.5f;
        T s = sin(a);
        T c = cos(a);

        return
        {
            axis.x * s,
            axis.y * s,
            axis.z * s,
            c
        };
    }

    template<typename T>
    quaternion<T> quaternion<T>::rotate_matrix(const matrix<T, 3, 3>& m)
    {
        quaternion<T> o;
        const T m11 = m.p_data[0], m12 = m.p_data[1], m13 = m.p_data[2];
        const T m21 = m.p_data[3], m22 = m.p_data[4], m23 = m.p_data[5];
        const T m31 = m.p_data[6], m32 = m.p_data[7], m33 = m.p_data[8];

        // Determine which of w, x, y or z has the largest absolute value
        T fourWSquaredMinus1 = +m11 + m22 + m33;
        T fourXSquaredMinus1 = +m11 - m22 - m33;
        T fourYSquaredMinus1 = -m11 + m22 - m33;
        T fourZSquaredMinus1 = -m11 - m22 + m33;

        int biggestIndex = 0;
        T fourBiggestSquardeMinus1 = fourWSquaredMinus1;
        if (fourXSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourXSquaredMinus1;
            biggestIndex = 1;
        }
        if (fourYSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourYSquaredMinus1;
            biggestIndex = 2;
        }
        if (fourZSquaredMinus1 > fourBiggestSquardeMinus1)
        {
            fourBiggestSquardeMinus1 = fourZSquaredMinus1;
            biggestIndex = 3;
        }

        T biggestVal = sqrt(fourBiggestSquardeMinus1 + 1) * .5f;
        T mult = 0.25f / biggestVal;

        switch (biggestIndex)
        {
        case 0:
        {
            o.x = (m23 - m32) * mult;
            o.y = (m31 - m13) * mult;
            o.z = (m12 - m21) * mult;
            o.w = biggestVal;
            break;
        }
        case 1:
        {
            o.x = biggestVal;
            o.y = (m12 + m21) * mult;
            o.z = (m31 + m13) * mult;
            o.w = (m23 - m32) * mult;
            break;
        }
        case 2:
        {
            o.x = (m12 + m21) * mult;
            o.y = biggestVal;
            o.z = (m23 + m32) * mult;
            o.w = (m31 - m13) * mult;
            break;
        }
        case 3:
        {
            o.x = (m31 + m13) * mult;
            o.y = (m23 + m32) * mult;
            o.z = biggestVal;
            o.w = (m12 - m21) * mult;
            break;
        }
        default:
            o.x = 0.0f;
            o.y = 0.0f;
            o.z = 0.0f;
            o.w = 1.0f;
            break;
        }
        return o;
    }

    template<typename T>
    quaternion<T>::quaternion(const T ix, const T iy, const T iz, const T iw)
        : x(ix)
        , y(iy)
        , z(iz)
        , w(iw)
    {

    }

    template<typename T>
    quaternion<T>::quaternion(const vector<T, 3>& v, const T w)
        : x(v[0])
        , y(v[1])
        , z(v[2])
        , w(w)
    {

    }

    template<typename T>
    quaternion<T>::quaternion(const quaternion<T>& o)
        : x(o.x)
        , y(o.y)
        , z(o.z)
        , w(o.w)
    {

    }

    template<typename T>
    quaternion<T>::quaternion(const T* o)
        : x(o[0])
        , y(o[1])
        , z(o[2])
        , w(o[3])
    {

    }

    template<typename T>
    quaternion<T>::quaternion()
        : x(static_cast<T>(0))
        , y(static_cast<T>(0))
        , z(static_cast<T>(0))
        , w(static_cast<T>(1))
    {

    }

    template<typename T>
    T quaternion<T>::length_squared() const
    {
        return (x * x) + (y * y) + (z * z) + (w * w);
    }

    template<typename T>
    T quaternion<T>::length() const
    {
        return sqrtf((x * x) + (y * y) + (z * z) + (w * w));
    }

    template<typename T>
    quaternion<T>& quaternion<T>::normalize()
    {
        T il = 1.0f / sqrtf((x * x) + (y * y) + (z * z) + (w * w));
        x *= il;
        y *= il;
        z *= il;
        w *= il;
        return *this;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::safe_normalize(const T tolerance /*= 0.001f*/)
    {
        T l = sqrtf((x * x) + (y * y) + (z * z) + (w * w));
        if (l > tolerance)
        {
            T il = 1.0f / l;
            x *= il;
            y *= il;
            z *= il;
            w *= il;
        }

        return *this;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::conjugate()
    {
        x = -x;
        y = -y;
        z = -z;
        return *this;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::inverse()
    {
        T il = 1.0f / sqrtf((x * x) + (y * y) + (z * z) + (w * w));
        x = -x * il;
        y = -y * il;
        z = -z * il;
        w = w * il;
        return *this;
    }

    template<typename T>
    matrix<T, 3, 3> quaternion<T>::matrix3x3() const
    {
        T x2 = x * x;
        T y2 = y * y;
        T z2 = z * z;
        T wx = w * x;
        T wy = w * y;
        T wz = w * z;
        T xy = x * y;
        T xz = x * z;
        T yz = y * z;
        return
        {
            1.0f - (2.0f * y2) - (2.0f * z2),
            (2.0f * xy) + (2.0f * wz),
            (2.0f * xz) - (2.0f * wy),
            (2.0f * xy) - (2.0f * wz),
            1.0f - (2.0f * x2) - (2.0f * z2),
            (2.0f * yz) + (2.0f * wx),
            (2.0f * xz) + (2.0f * wy),
            (2.0f * yz) - (2.0f * wx),
            1.0f - (2.0f * x2) - (2.0f * y2)
        };
    }

    template<typename T>
    matrix<T, 4, 3> quaternion<T>::matrix4x3() const
    {
        T x2 = x * x;
        T y2 = y * y;
        T z2 = z * z;
        T wx = w * x;
        T wy = w * y;
        T wz = w * z;
        T xy = x * y;
        T xz = x * z;
        T yz = y * z;

        return
        {
            1.0f - (2.0f * y2) - (2.0f * z2),
            (2.0f * xy) + (2.0f * wz),
            (2.0f * xz) - (2.0f * wy),
            (2.0f * xy) - (2.0f * wz),
            1.0f - (2.0f * x2) - (2.0f * z2),
            (2.0f * yz) + (2.0f * wx),
            (2.0f * xz) + (2.0f * wy),
            (2.0f * yz) - (2.0f * wx),
            1.0f - (2.0f * x2) - (2.0f * y2),
            0.0f,
            0.0f,
            0.0f
        };
    }

    template<typename T>
    matrix<T, 4, 4> quaternion<T>::matrix4x4() const
    {
        T x2 = x * x;
        T y2 = y * y;
        T z2 = z * z;
        T wx = w * x;
        T wy = w * y;
        T wz = w * z;
        T xy = x * y;
        T xz = x * z;
        T yz = y * z;

        return
        {
            1.0f - (2.0f * y2) - (2.0f * z2),
            (2.0f * xy) + (2.0f * wz),
            (2.0f * xz) - (2.0f * wy),
            0.0f,
            (2.0f * xy) - (2.0f * wz),
            1.0f - (2.0f * x2) - (2.0f * z2),
            (2.0f * yz) + (2.0f * wx),
            0.0f,
            (2.0f * xz) + (2.0f * wy),
            (2.0f * yz) - (2.0f * wx),
            1.0f - (2.0f * x2) - (2.0f * y2),
            0.0f,
            0.0f,
            0.0f,
            0.0f,
            1.0f
        };
    }

    template<typename T>
    vector<T, 3> quaternion<T>::euler() const
    {
        vector<T, 3> o;
        auto sp = -2.0f * ((y * z) - (w * x));
        if (abs(sp) > 0.9999f)
        {
            o.x = 1.570796f * sp;
            o.y = atan2((-x * z) + (w * y), 0.5f - (y * y) - (z * z));
            o.z = 0.0f;
        }
        else
        {
            o.x = asin(sp);
            o.y = atan2((x * z) + (w * y), 0.5f - (x * x) - (y * y));
            o.z = atan2((x * y) + (w * z), 0.5f - (x * x) - (z * z));
        }
        return o;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::operator+=(const quaternion<T>& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::operator-=(const quaternion<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::operator*=(const quaternion<T>& rhs)
    {
        T ow = (w * rhs.w) - ((x * rhs.x) + (y * rhs.y) + (z * rhs.z));
        T ox = (w * rhs.x) + (rhs.w * x) + ((y * rhs.z) - (z * rhs.y));
        T oy = (w * rhs.y) + (rhs.w * y) + ((z * rhs.x) - (x * rhs.z));
        T oz = (w * rhs.z) + (rhs.w * z) + ((x * rhs.y) - (y * rhs.x));
        x = ox;
        y = oy;
        z = oz;
        w = ow;
        return *this;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::operator=(const quaternion<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::operator*=(const T rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;
        return *this;
    }

    template<typename T>
    quaternion<T>& quaternion<T>::operator/=(const T rhs)
    {
        T inv = 1.0f / rhs;
        x *= inv;
        y *= inv;
        z *= inv;
        w *= inv;
        return *this;
    }

    template<typename T>
    quaternion<T> quaternion<T>::operator-() const
    {
        quaternion<T> out;
        out.x = -x;
        out.y = -y;
        out.z = -z;
        out.w = -w;
        return out;
    }

    template<typename T>
    T& quaternion<T>::operator[](const unsigned int index)
    {
        return components[index];
    }

    template<typename T>
    const T& quaternion<T>::operator[](const unsigned int index) const
    {
        return components[index];
    }

    template<typename T, typename... Ts>
    quaternion<T> quaternion_set(Ts... values)
    {
        const T components[4] = { values... };
        return quaternion<T>(components);
    }

    template<typename T>
    T dot(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w);
    }

    template<typename T>
    T length_squared(const quaternion<T>& val)
    {
        return (val.x * val.x) + (val.y * val.y) + (val.z * val.z) + (val.w * val.w);
    }

    template<typename T>
    T length(const quaternion<T>& rhs)
    {
        return sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
    }

    template<typename T>
    quaternion<T> normalize(const quaternion<T>& rhs)
    {
        T il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    template<typename T>
    quaternion<T> safe_normalize(const quaternion<T>& rhs, const T tolerance /*= 0.001f*/)
    {
        T l = sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        T il = (l > tolerance) ? (1.0f / l) : 1.0f;
        return
        {
            rhs.x * il,
            rhs.y * il,
            rhs.z * il,
            rhs.w * il
        };
    }

    template<typename T>
    quaternion<T> conjugate(const quaternion<T>& rhs)
    {
        return { -rhs.x, -rhs.y, -rhs.z, rhs.w };
    }

    template<typename T>
    quaternion<T> inverse(const quaternion<T>& rhs)
    {
        T il = 1.0f / sqrt((rhs.x * rhs.x) + (rhs.y * rhs.y) + (rhs.z * rhs.z) + (rhs.w * rhs.w));
        return { -rhs.x * il, -rhs.y * il, -rhs.z * il, rhs.w * il };
    }

    template<typename T>
    quaternion<T> slerp(quaternion<T> q0, quaternion<T> q1, const T t)
    {
        T cosAngle = (q0.x * q1.x) + (q0.y * q1.y) + (q0.z * q1.z) + (q0.w * q1.w);
        if (cosAngle < 0.0f) {
            q1 = -q1;
            cosAngle = -cosAngle;
        }

        T k0, k1;

        // Check for divide by zero
        if (cosAngle > 0.9999f) {
            k0 = 1.0f - t;
            k1 = t;
        }
        else {
            T angle = acos(cosAngle);
            T oneOverSinAngle = 1.0f / sin(angle);

            k0 = ((sin(1.0f - t) * angle) * oneOverSinAngle);
            k1 = (sin(t * angle) * oneOverSinAngle);
        }

        q0 = q0 * k0;
        q1 = q1 * k1;

        return q0 + q1;
    }

    template<typename T>
    quaternion<T> operator+(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        return
        {
            lhs.x + rhs.x,
            lhs.y + rhs.y,
            lhs.z + rhs.z,
            lhs.w + rhs.w
        };
    }

    template<typename T>
    quaternion<T> operator-(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        return
        {
            lhs.x - rhs.x,
            lhs.y - rhs.y,
            lhs.z - rhs.z,
            lhs.w - rhs.w
        };
    }

    template<typename T>
    quaternion<T> operator*(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        return quaternion<T>(cross(lhs.v, rhs.v) + (lhs.w * rhs.v) + (rhs.w * lhs.v), lhs.w * rhs.w - dot(lhs.v, rhs.v));
    }

    template<typename T>
    vector<T, 3> operator*(const vector<T, 3>& lhs, const quaternion<T>& rhs)
    {
        return rhs * lhs;
    }

    template<typename T>
    vector<T, 3> operator*(const quaternion<T>& lhs, const vector<T, 3>& rhs)
    {
        vector<T, 3> vxp = cross(lhs.v, rhs);
        vector<T, 3> vxpxv = cross(lhs.v, vxp);
        return rhs + ((vxp * lhs.w) + vxpxv) * 2.0f;
    }

    template<typename T>
    quaternion<T> operator*(const quaternion<T>& lhs, const T rhs)
    {
        return
        {
            lhs.x * rhs,
            lhs.y * rhs,
            lhs.z * rhs,
            lhs.w * rhs
        };
    }

    template<typename T>
    quaternion<T> operator/(const quaternion<T>& lhs, const T rhs)
    {
        T inv = 1.0f / rhs;
        return
        {
            lhs.x * inv,
            lhs.y * inv,
            lhs.z * inv,
            lhs.w * inv
        };
    }

    template<typename T>
    quaternion<T> operator*(const T lhs, const quaternion<T>& rhs)
    {
        return
        {
            rhs.x * lhs,
            rhs.y * lhs,
            rhs.z * lhs,
            rhs.w * lhs
        };
    }
}