#include "../include/vector.h"
#include "../include/matrix.h"
#include "../include/quaternion.h"

#include <gtest/gtest.h>

#define EXPECT_NEAR3(expected, actual, tolerance)\
 for (unsigned int i = 0; i < 3; ++i)\
        EXPECT_NEAR(expected[i], actual[i], tolerance)

#define EXPECT_NEAR4(expected, actual, tolerance)\
 for (unsigned int i = 0; i < 4; ++i)\
        EXPECT_NEAR(expected[i], actual[i], tolerance)

#define EXPECT_NEAR16(expected, actual, tolerance)\
 for (unsigned int i = 0; i < 16; ++i)\
        EXPECT_NEAR(expected.entries1d[i], actual.entries1d[i], tolerance)

namespace lina
{
    namespace unit_tests
    {
        TEST(float4, addition)
        {
            typedef vector<float, 4> operand_type;
            typedef vector<float, 4> result_type;

            operand_type operand0 = operand_type(-1.0f);
            operand_type operand1 = operand_type(1.0f);
            result_type expected = result_type(0.0f);
            result_type actual = operand0 + operand1;
            EXPECT_NEAR4(expected, actual, 0.001f);
        }

        TEST(float4, subtraction)
        {
            typedef vector<float, 4> operand_type;
            typedef vector<float, 4> result_type;

            operand_type operand0 = operand_type(1.0f);
            operand_type operand1 = operand_type(1.0f);
            result_type expected = result_type(0.0f);
            result_type actual = operand0 - operand1;
            EXPECT_NEAR4(expected, actual, 0.001f);
        }

        TEST(float4, multiplication)
        {
            typedef vector<float, 4> operand_type;
            typedef vector<float, 4> result_type;

            operand_type operand0 = operand_type(0.1f);
            operand_type operand1 = operand_type(10.0f);
            result_type expected = result_type(1.0f);
            result_type actual = operand0 * operand1;
            EXPECT_NEAR4(expected, actual, 0.001f);
        }

        TEST(float4, dot)
        {
            typedef vector<float, 4> operand_type;
            typedef float result_type;

            operand_type operand0 = vector_set<float, 4>(1.0f, 0.0f, 0.0f, 0.0f);
            operand_type operand1 = vector_set<float, 4>(0.0f, 0.0f, 1.0f, 0.0f);
            result_type expected = 0.0f;
            result_type actual = dot(operand0, operand1);
            EXPECT_NEAR(expected, actual, 0.001f);
        }

        TEST(float4, normalization)
        {
            typedef vector<float, 4> operand_type;
            typedef float result_type;
            
            operand_type operand0 = vector_set<float, 4>(1.0f, 3.0f, 5.0f, 1.0f);
            operand_type operand1 = operand0 / operand0.length();
            result_type expected = 1.0f;
            result_type actual = operand1.length();
            EXPECT_NEAR(expected, actual, 0.001f);
        }

        TEST(float4, projection)
        {
            typedef vector<float, 4> operand_type;
            typedef vector<float, 4> result_type;

            operand_type operand0 = vector_set<float, 4>(0.9f, 0.5f, 0.2f, 0.1f);
            operand_type operand1 = vector_set<float, 4>(0.4f, 0.3f, 0.1f, 0.0f);
            result_type expected = operand0;
            result_type actual = project(operand0, operand1) + reject(operand0, operand1);
            EXPECT_NEAR4(expected, actual, 0.001f);
        }

        TEST(matrix, addition)
        {
            typedef matrix<float, 4, 4> operand_type;
            typedef matrix<float, 4, 4> result_type;

            operand_type operand0 = operand_type(-1.0f);
            operand_type operand1 = operand_type(1.0f);
            result_type expected = result_type(0.0f);
            result_type actual = operand0 + operand1;
            EXPECT_NEAR16(expected, actual, 0.001f);
        }

        TEST(matrix, subtraction)
        {
            typedef matrix<float, 4, 4> operand_type;
            typedef matrix<float, 4, 4> result_type;

            operand_type operand0 = operand_type(1.0f);
            operand_type operand1 = operand_type(1.0f);
            result_type expected = result_type(0.0f);
            result_type actual = operand0 - operand1;
            EXPECT_NEAR16(expected, actual, 0.001f);
        }

        TEST(matrix, multiplication)
        {
            typedef matrix<float, 4, 4> operand_type;
            typedef matrix<float, 4, 4> result_type;

            operand_type operand0 = 
                matrix_set<float, 4, 4>(
                    10.f, 0.f, 0.f, 0.f,
                    0.f, 10.f, 0.f, 0.f,
                    0.f, 0.f, 10.f, 0.f,
                    0.f, 0.f, 0.f, 1.0f);
            operand_type operand1 = 
                matrix_set<float, 4, 4>(
                    1.0f, 0.f, 0.f, 0.f,
                    0.f, 1.0f, 0.f, 0.f,
                    0.f, 0.f, 1.0f, 0.f,
                    3.f, 4.f, 50.f, 1.0f);
            result_type expected =
                matrix_set<float, 4, 4>(
                    10.f, 0.f, 0.f, 0.f,
                    0.f, 10.f, 0.f, 0.f,
                    0.f, 0.f, 10.f, 0.f,
                    3.f, 4.f, 50.f, 1.0f);
            result_type actual = operand0 * operand1;
            EXPECT_NEAR16(expected, actual, 0.001f);
        }

        TEST(matrix, transformation)
        {
            typedef matrix<float, 4, 4> operand_type0;
            typedef vector<float, 4> operand_type1;
            typedef vector<float, 4> result_type;

            float c = cosf(0.524f);
            float s = sinf(0.524f);
            operand_type0 operand0 = 
                matrix_set<float, 4, 4>(
                    c, s, 0.f, 0.f,
                   -s, c, 0.f, 0.f,
                    0.f, 0.f, 1.f, 0.f,
                    0.f, 0.f, 0.f, 1.f);
            operand_type1 operand1= vector_set<float, 4>(1.0f, 0.0f, 0.0f, 1.0f);
            result_type expected = vector_set<float, 4>(0.866f, 0.5f, 0.0f, 1.0f);
            result_type actual = operand1 * operand0;
            EXPECT_NEAR4(expected, actual, 0.001f);
        }

        TEST(matrix, transpose)
        {
            typedef matrix<float, 4, 4> operand_type;
            typedef matrix<float, 4, 4> result_type;

            operand_type operand0 = 
                matrix_set<float, 4, 4>(
                    0.f, 1.f, 2.f, 3.f,
                    4.f, 5.f, 6.f, 7.f,
                    8.f, 9.f, 10.f, 11.f,
                    12.f, 13.f, 14.f, 15.f);
            result_type expected = 
                matrix_set<float, 4, 4>(
                    0.f, 4.f,  8.f, 12.f,
                    1.f, 5.f,  9.f, 13.f,
                    2.f, 6.f, 10.f,14.f,
                    3.f, 7.f, 11.f,15.f);
            result_type actual = operand0.transpose();
            EXPECT_NEAR16(expected, actual, 0.001f);
        }

        TEST(quaternion, multiplication)
        {
            typedef quaternion<float> operand_type0;
            typedef quaternion<float> result_type;

            operand_type0 operand0 = operand_type0::rotate_euler(0.0f, 0.0f, 0.26179939);
            operand_type0 operand1 = operand_type0::rotate_euler(0.0f, 0.0f, 0.26179939);
            result_type expected = operand_type0::rotate_euler(0.0f, 0.0f, 0.524f);
            result_type actual = operand1 * operand0;
            EXPECT_NEAR4(expected, actual, 0.001f);
        }

        TEST(quaternion, transformation)
        {
            typedef quaternion<float> operand_type0;
            typedef vector<float, 3> operand_type1;
            typedef vector<float, 3> result_type;

            operand_type0 operand0 = operand_type0::rotate_euler(0.0f, 0.0f, 0.524f);
            operand_type1 operand1 = vector_set<float,3>(1.0f, 0.0f, 0.0f);
            result_type expected = vector_set<float,3>(0.866f, 0.5f, 0.0f);
            result_type actual = operand0.transform_point(operand1);
            EXPECT_NEAR3(expected, actual, 0.001f);
        }

        TEST(quaternion, difference)
        {
            typedef quaternion<float> operand_type0;
            typedef quaternion<float> operand_type1;
            typedef quaternion<float> result_type;

            operand_type0 operand0 = operand_type0::rotate_euler(0.0f, 0.0f, 0.524f);
            operand_type1 operand1 = operand_type0::rotate_euler(0.0f, 0.524f, 0.0f);
            result_type difference = result_type::rotate_difference(operand0, operand1);
            result_type expected = difference * operand0;
            result_type actual = operand1;
            EXPECT_NEAR4(expected, actual, 0.001f);
        }
    }
}