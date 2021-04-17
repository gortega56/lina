#pragma once
#include "vector.h"

namespace lina
{
    template<typename T, unsigned int N, unsigned int M>
    struct matrix
    {
        static constexpr unsigned int s_num_entries = N * M;

        union
        {
            T entries1d[s_num_entries];
            T entries2d[N][M];
        };

        static matrix<T, N, M> identity();

        matrix();

        matrix(const T entry);

        matrix(const T* p_entries);

        matrix(const matrix<T, N, M>& o);

        T minor(unsigned int i, unsigned int j) const;

        T cofactor(unsigned int i, unsigned int j) const;

        T determinant() const;

        matrix<T, M, N> transpose() const;

        matrix<T, N, M> adjugate() const;

        matrix<T, N, M> inverse() const;

        matrix<T, N, M>& operator=(const matrix<T, N, M>& rhs);

        matrix<T, N, M>& operator+=(const matrix<T, N, M>& rhs);

        matrix<T, N, M>& operator-=(const matrix<T, N, M>& rhs);

        matrix<T, N, M>& operator+=(const T rhs);

        matrix<T, N, M>& operator-=(const T rhs);

        matrix<T, N, M>& operator*=(const T rhs);

        matrix<T, N, M>& operator/=(const T rhs);

        matrix<T, N, M> operator-() const;

        T* operator[](unsigned int i);
    };

    template<typename T, unsigned int N, unsigned int M, typename... Ts>
    matrix<T, N, M> matrix_set(Ts... values);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator+(const matrix<T, N, M>& lhs, const matrix<T, N, M>& rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator-(const matrix<T, N, M>& lhs, const matrix<T, N, M>& rhs);

    template<typename T, unsigned int N, unsigned int M, unsigned int P>
    matrix<T, N, P> operator*(const matrix<T, N, M>& lhs, const matrix<T, M, P>& rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator+(const matrix<T, N, M>& lhs, const T rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator-(const matrix<T, N, M>& lhs, const T rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator*(const matrix<T, N, M>& lhs, const T rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator/(const matrix<T, N, M>& lhs, const T rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator+(const T lhs, const matrix<T, N, M>& rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator-(const T lhs, const matrix<T, N, M>& rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator*(const T lhs, const matrix<T, N, M>& rhs);

    template<typename T, unsigned int N, unsigned int M>
    vector<T, N> operator*(const vector<T, N>& lhs, const matrix<T, N, M>& rhs);

    template<typename T, unsigned int N, unsigned int M>
    vector<T, M> operator*(const matrix<T, N, M>& lhs, const vector<T, M>& rhs);

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> matrix<T, N, M>::identity()
    {
        matrix<T, N, M> result = matrix<T, N, M>();
        for (unsigned int i = 0; i < N; ++i)
        {
            for (unsigned int j = 0; j < P; ++j)
            {
                if (i == j)
                    result.entries2d[i][j] = static_cast<T>(1);
            }
        }
        return result;
    }
    

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>::matrix()
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] = T(0);
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>::matrix(const T entry)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] = entry;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>::matrix(const T* p_entries)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] = p_entries[i];
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>::matrix(const matrix<T, N, M>& o)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] = o.entries1d[i];
    }

    template<typename T, unsigned int N, unsigned int M>
    T matrix<T, N, M>::minor(unsigned int i, unsigned int j) const
    {
        matrix<T, N - 1, M - 1> result;
        unsigned int row = 0;
        for (unsigned int ii = 0; ii < N; ++ii)
        {
            if (ii != i)
            {
                unsigned int col = 0;
                for (unsigned int jj = 0; jj < M; ++jj)
                {
                    if (jj != j)
                    {
                        result.entries2d[row][col] = entries2d[ii][jj];
                        ++col;
                    }
                }
                ++row;
            }
        }
        return result.determinant();
    }

    template<typename T, unsigned int N, unsigned int M>
    T matrix<T, N, M>::cofactor(unsigned int i, unsigned int j) const
    {
        T sign = (((i + j) & 1) == 0) ? T(1) : T(-1);
        return sign * minor(i, j);
    }

    template<typename T, unsigned int N, unsigned int M>
    T matrix<T, N, M>::determinant() const
    {
        T result = T(0);
        for (unsigned int i = 0; i < N; ++i)
            result += entries2d[i][M - 1] * cofactor(i, M - 1);
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, M, N> matrix<T, N, M>::transpose() const
    {
        matrix<T, M, N> result;
        for (unsigned int i = 0; i < N; ++i)
        {
            for (unsigned int j = 0; j < M; ++j)
                result.entries2d[i][j] = entries2d[j][i];
        }
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> matrix<T, N, M>::adjugate() const
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < N; ++i)
        {
            for (unsigned int j = 0; j < M; ++j)
                result.entries2d[i][j] = cofactor(i, j);
        }
        return result.transpose();
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> matrix<T, N, M>::inverse() const
    {
        matrix<T, N, M> adj = adjugate();
        T det = determinant();
        T inv = 1.0f / det;
        return adj * inv;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>& matrix<T, N, M>::operator=(const matrix<T, N, M>& rhs)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] = rhs.entries1d[i];
        return *this;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>& matrix<T, N, M>::operator+=(const matrix<T, N, M>& rhs)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] += rhs.entries1d[i];
        return *this;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>& matrix<T, N, M>::operator-=(const matrix<T, N, M>& rhs)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] -= rhs.entries1d[i];
        return *this;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>& matrix<T, N, M>::operator+=(const T rhs)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] += rhs;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>& matrix<T, N, M>::operator-=(const T rhs)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] -= rhs;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>& matrix<T, N, M>::operator*=(const T rhs)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] *= rhs;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M>& matrix<T, N, M>::operator/=(const T rhs)
    {
        for (unsigned int i = 0; i < s_num_entries; ++i)
            entries1d[i] /= rhs;
        return *this;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> matrix<T, N, M>::operator-() const
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < s_num_entries; ++i)
            result.entries1d[i] = -entries1d[i];
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    T* matrix<T, N, M>::operator[](unsigned int i)
    {
        return entries1d[i];
    }

    template<typename T, unsigned int N, unsigned int M, typename... Ts>
    matrix<T, N, M> matrix_set(Ts... values)
    {
        const T entries[matrix<T, N, M>::s_num_entries] = { values... };
        return matrix<T, N, M>(entries);
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator+(const matrix<T, N, M>& lhs, const matrix<T, N, M>& rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs.entries1d[i] + rhs.entries1d[i];
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator-(const matrix<T, N, M>& lhs, const matrix<T, N, M>& rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs.entries1d[i] - rhs.entries1d[i];
        return result;
    }

    template<typename T, unsigned int N, unsigned int M, unsigned int P>
    matrix<T, N, P> operator*(const matrix<T, N, M>& lhs, const matrix<T, M, P>& rhs)
    {
        matrix<T, N, P> result = matrix<T, N, P>();
        for (unsigned int i = 0; i < N; ++i)
        {
            for (unsigned int j = 0; j < P; ++j)
            {
                for (unsigned int k = 0; k < M; ++k)
                    result.entries2d[i][j] += lhs.entries2d[i][k] * rhs.entries2d[k][j];
            }
        }
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    vector<T, N> operator*(const vector<T, N>& lhs, const matrix<T, N, M>& rhs)
    {
        vector<T, N> result = vector<T, N>();
        for (unsigned int j = 0; j < M; ++j)
        {
            for (unsigned int k = 0; k < N; ++k)
                result.components[j] += lhs.components[k] * rhs.entries2d[k][j];
        }
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    vector<T, M> operator*(const matrix<T, N, M>& lhs, const vector<T, M>& rhs)
    {
        vector<T, M> result = vector<T, M>();
        for (unsigned int i = 0; i < N; ++i)
        {
            for (unsigned int k = 0; k < M; ++k)
                result.components[i] += lhs.entries2d[i] * rhs.components[i][k];
        }
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator+(const matrix<T, N, M>& lhs, const T rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs.entries1d[i] + rhs;
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator-(const matrix<T, N, M>& lhs, const T rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs.entries1d[i] - rhs;
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator*(const matrix<T, N, M>& lhs, const T rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs.entries1d[i] * rhs;
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator/(const matrix<T, N, M>& lhs, const T rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs.entries1d[i] / rhs;
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator+(const T lhs, const matrix<T, N, M>& rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs + rhs.entries1d[i];
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator-(const T lhs, const matrix<T, N, M>& rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs - rhs.entries1d[i];
        return result;
    }

    template<typename T, unsigned int N, unsigned int M>
    matrix<T, N, M> operator*(const T lhs, const matrix<T, N, M>& rhs)
    {
        matrix<T, N, M> result;
        for (unsigned int i = 0; i < matrix<T, N, M>::s_num_entries; ++i)
            result.entries1d[i] = lhs * rhs.entries1d[i];
        return result;
    }

    template<typename T>
    struct matrix<T, 2, 2>
    {
        static constexpr unsigned int s_num_entries = 4;

        union
        {
            T entries1d[s_num_entries];
            T entries2d[2][2];
        };

        T minor(unsigned int i, unsigned int j) const;

        T determinant() const;

        matrix<T, 2, 2> inverse() const;

        matrix<T, 2, 2> transpose() const;
    };

    template<typename T>
    T matrix<T, 2, 2>::minor(unsigned int i, unsigned int j) const
    {
        for (unsigned int ii = 0; ii < 2; ++ii)
        {
            if (ii != i)
            {
                for (unsigned int jj = 0; jj < 2; ++jj)
                {
                    if (jj != j)
                    {
                        return entries2d[ii][jj];
                    }
                }
            }
        }

        // Should never get here.
        return T(0);
    }

    template<typename T>
    T  matrix<T, 2, 2>::determinant() const
    {
        return (entries2d[0][0] * entries2d[1][1]) - (entries2d[0][1] * entries2d[1][0]);
    }

    template<typename T>
    matrix<T, 2, 2>  matrix<T, 2, 2>::inverse() const
    {
        T id = static_cast<T>(1.0 / determinant());
        return
        {
             entries2d[1][1] * id, -entries2d[0][1] * id,
            -entries2d[1][0] * id, entries2d[0][0] * id,
        };
    }

    template<typename T>
    matrix<T, 2, 2>  matrix<T, 2, 2>::transpose() const
    {
        return
        {
            entries2d[0][0], entries2d[1][0],
            entries2d[0][1], entries2d[1][1],
        };
    }
}
