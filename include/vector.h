#pragma once
#include <cmath>

namespace lina
{
    template<typename T, unsigned int N>
    struct vector
    {
        T components[N];

        vector();

        vector(const T* p_components);

        vector(const T components);

        vector(const vector<T, N>& o);

        T length_manhattan() const;

        T length_squared() const;

        T length() const;

        vector<T, N>& vector<T, N>::operator=(const vector<T, N>& rhs);

        vector<T, N>& operator+=(const vector<T, N>& rhs);

        vector<T, N>& operator-=(const vector<T, N>& rhs);

        vector<T, N>& operator*=(const vector<T, N>& rhs);

        vector<T, N>& operator+=(const T rhs);

        vector<T, N>& operator-=(const T rhs);

        vector<T, N>& operator*=(const T rhs);

        vector<T, N>& operator/=(const T rhs);

        vector<T, N> operator-() const;

        T& operator[](const unsigned int index);

        const T& operator[](const unsigned int index) const;
    };

    template<typename T, unsigned int N, typename... Ts>
    vector<T, N> vector_set(Ts... values);

    template<typename T, unsigned int N>
    T dot(const vector<T, N>& lhs, const vector<T, N>& rhs);

    template<typename T, unsigned int N>
    T length_manhattan(const vector<T, N>& val);

    template<typename T, unsigned int N>
    T length_squared(const vector<T, N>& val);

    template<typename T, unsigned int N>
    T length(const vector<T, N>& val);

    template<typename T, unsigned int N>
    vector<T, N> project(const vector<T, N>& lhs, const vector<T, N>& rhs);

    template<typename T, unsigned int N>
    vector<T, N> reject(const vector<T, N>& lhs, const vector<T, N>& rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator+(const vector<T, N>& lhs, const vector<T, N>& rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator-(const vector<T, N>& lhs, const vector<T, N>& rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator*(const vector<T, N>& lhs, const vector<T, N>& rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator+(const vector<T, N>& lhs, const T rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator-(const vector<T, N>& lhs, const T rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator*(const vector<T, N>& lhs, const T rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator/(const vector<T, N>& lhs, const T rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator+(const T lhs, const vector<T, N>& rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator-(const T lhs, const vector<T, N>& rhs);

    template<typename T, unsigned int N>
    vector<T, N> operator*(const T lhs, const vector<T, N>& rhs);

    template<typename T>
    vector<T, 3> cross(const vector<T, 3>& lhs, const vector<T, 3>& rhs);

    template<typename T, unsigned int N>
    vector<T, N>::vector()
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] = T(0);
    }

    template<typename T, unsigned int N>
    vector<T, N>::vector(const T* p_components)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] = p_components[i];
    }

    template<typename T, unsigned int N>
    vector<T, N>::vector(const T component)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] = component;
    }

    template<typename T, unsigned int N>
    vector<T, N>::vector(const vector<T, N>& o)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] = o.components[i];
    }

    template<typename T, unsigned int N>
    T vector<T, N>::length_manhattan() const
    {
        T result = T(0);
        for (unsigned int i = 0; i < N; ++i)
            result += T(abs(components[i]));
        return result;
    }

    template<typename T, unsigned int N>
    T vector<T, N>::length_squared() const
    {
        T result = T(0);
        for (unsigned int i = 0; i < N; ++i)
            result += components[i] * components[i];
        return result;
    }

    template<typename T, unsigned int N>
    T vector<T, N>::length() const
    {
        T result = T(0);
        for (unsigned int i = 0; i < N; ++i)
            result += components[i] * components[i];
        return static_cast<T>(sqrt(result));
    }

    template<typename T, unsigned int N>
    vector<T, N>& vector<T, N>::operator=(const vector<T, N>& rhs)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] = rhs.components[i];
        return *this;
    }

    template<typename T, unsigned int N>
    vector<T, N>& vector<T, N>::operator+=(const vector<T, N>& rhs)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] += rhs.components[i];
        return *this;
    }

    template<typename T, unsigned int N>
    vector<T, N>& vector<T, N>::operator-=(const vector<T, N>& rhs)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] -= rhs.components[i];
        return *this;
    }

    template<typename T, unsigned int N>
    vector<T, N>& vector<T, N>::operator*=(const vector<T, N>& rhs)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] *= rhs.components[i];
        return *this;
    }

    template<typename T, unsigned int N>
    vector<T, N>& vector<T, N>::operator+=(const T rhs)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] += rhs;
        return *this;
    }

    template<typename T, unsigned int N>
    vector<T, N>& vector<T, N>::operator-=(const T rhs)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] -= rhs;
        return *this;
    }

    template<typename T, unsigned int N>
    vector<T, N>& vector<T, N>::operator*=(const T rhs)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] *= rhs;
        return *this;
    }

    template<typename T, unsigned int N>
    vector<T, N>& vector<T, N>::operator/=(const T rhs)
    {
        for (unsigned int i = 0; i < N; ++i)
            components[i] /= rhs;
        return *this;
    }

    template<typename T, unsigned int N>
    vector<T, N> vector<T, N>::operator-() const
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = -components[i];
        return result;
    }

    template<typename T, unsigned int N>
    T& vector<T, N>::operator[](const unsigned int index)
    {
        return components[index];
    }

    template<typename T, unsigned int N>
    const T& vector<T, N>::operator[](const unsigned int index) const
    {
        return components[index];
    }

    template<typename T, unsigned int N, typename... Ts>
    vector<T, N> vector_set(Ts... values)
    {
        const T components[N] = { values... };
        return vector<T, N>(components);
    }

    template<typename T, unsigned int N>
    T dot(const vector<T, N>& lhs, const vector<T, N>& rhs)
    {
        T result = static_cast<T>(0);
        for (unsigned int i = 0; i < N; ++i)
            result += lhs.components[i] * rhs.components[i];
        return result;
    }

    template<typename T, unsigned int N>
    T length_manhattan(const vector<T, N>& vec)
    {
        return vec.length_manhattan();
    }

    template<typename T, unsigned int N>
    T length_squared(const vector<T, N>& vec)
    {
        return vec.length_squared();
    }

    template<typename T, unsigned int N>
    T length(const vector<T, N>& vec)
    {
        return vec.length();
    }

    template<typename T, unsigned int N>
    vector<T, N> project(const vector<T, N>& lhs, const vector<T, N>& rhs)
    {
        return rhs * (dot(lhs, rhs) / length_squared(rhs));
    }

    template<typename T, unsigned int N>
    vector<T, N> reject(const vector<T, N>& lhs, const vector<T, N>& rhs)
    {
        return lhs - project(lhs, rhs);
    }

    template<typename T, unsigned int N>
    vector<T, N> operator+(const vector<T, N>& lhs, const vector<T, N>& rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs.components[i] + rhs.components[i];
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator-(const vector<T, N>& lhs, const vector<T, N>& rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs.components[i] - rhs.components[i];
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator*(const vector<T, N>& lhs, const vector<T, N>& rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs.components[i] * rhs.components[i];
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator+(const vector<T, N>& lhs, const T rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs.components[i] + rhs;
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator-(const vector<T, N>& lhs, const T rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs.components[i] - rhs;
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator*(const vector<T, N>& lhs, const T rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs.components[i] * rhs;
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator/(const vector<T, N>& lhs, const T rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs.components[i] / rhs;
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator+(const T lhs, const vector<T, N>& rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs + rhs.components[i];
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator-(const T lhs, const vector<T, N>& rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs - rhs.components[i];
        return result;
    }

    template<typename T, unsigned int N>
    vector<T, N> operator*(const T lhs, const vector<T, N>& rhs)
    {
        vector<T, N> result;
        for (unsigned int i = 0; i < N; ++i)
            result.components[i] = lhs * rhs.components[i];
        return result;
    }

    template<typename T>
    vector<T, 3> cross(const vector<T, 3>& lhs, const vector<T, 3>& rhs)
    {
        vector<T, 3> result;
        result[0] = (lhs[1] * rhs[2]) - (lhs[2] * rhs[1]);
        result[1] = (lhs[2] * rhs[0]) - (lhs[0] * rhs[2]);
        result[2] = (lhs[0] * rhs[1]) - (lhs[1] * rhs[0]);
        return result;
    }
}