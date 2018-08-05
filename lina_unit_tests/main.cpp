#include "lina.h"
#include <iostream>
#include <random>
#include <Windows.h>
#include <functional>
#include <assert.h>
#include <DirectXMath.h>

using namespace lina;
using namespace std;

template<typename T, size_t R, size_t C> void FillMatrix(mt19937& gen, T(&m)[R*C]);

template<typename T, size_t I, size_t J> bool Equals(const T* l, const T* r, const T& tolerance);

double TimeFunction(const char* name, std::function<void(void)> func);

void RunCorrectnessTests(mt19937& gen);

void RunSpeedTests(mt19937& gen);

template<typename T, size_t R, size_t C>
void FillMatrix(mt19937& gen, T(&m)[R*C])
{
    uniform_real_distribution<> dis(0.0, 100.0f);
    for (size_t r = 0; r < R; ++r)
    {
        for (size_t c = 0; c < C; ++c)
        {
            m[r*C + c] = (T)dis(gen);
        }
    }
}

template<typename T, size_t I, size_t J>
bool Equals(const T* l, const T* r, const T& tolerance)
{
    for (size_t i = 0; i < I; ++i)
    {
        for (size_t j = 0; j < J; ++j)
        {
            if (abs(l[i * J + j] - r[i * J + j]) > tolerance) return false;
        }
    }

    return true;
}

double TimeFunction(const char* name, std::function<void(void)> func)
{
    LARGE_INTEGER start, end;
    LARGE_INTEGER frequency;
    double elapsedTime;

    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start);

    func();

    QueryPerformanceCounter(&end);
    elapsedTime = (end.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;

    return elapsedTime;
}

void RunCorrectnessTests(mt19937& gen)
{
    bool correct = false;
    float eps = 0.01f;

    // float2
    {
        float v20[2], v21[2];
        FillMatrix<float, 1, 2>(gen, v20);
        FillMatrix<float, 1, 2>(gen, v21);

        DirectX::XMVECTOR xmv20, xmv21, xmv2sum, xmv2diff, xmv2prod, xmv2norm, xmv2dot, xmv2len, xmv2lensq;
        xmv20 = DirectX::XMVectorSet(v20[0], v20[1], 0, 0);
        xmv21 = DirectX::XMVectorSet(v21[0], v21[1], 0, 0);
        xmv2sum = DirectX::XMVectorAdd(xmv20, xmv21);
        xmv2diff = DirectX::XMVectorSubtract(xmv20, xmv21);
        xmv2prod = DirectX::XMVectorMultiply(xmv20, xmv21);
        xmv2norm = DirectX::XMVector2Normalize(xmv20);
        xmv2dot = DirectX::XMVector2Dot(xmv20, xmv21);
        xmv2len = DirectX::XMVector2Length(xmv20);
        xmv2lensq = DirectX::XMVector2LengthSq(xmv20);

        DirectX::XMFLOAT2 xmf2sum, xmf2diff, xmf2prod, xmf2norm;
        float xmf2dot, xmf2len, xmf2lensq;
        DirectX::XMStoreFloat2(&xmf2sum, xmv2sum);
        DirectX::XMStoreFloat2(&xmf2diff, xmv2diff);
        DirectX::XMStoreFloat2(&xmf2prod, xmv2prod);
        DirectX::XMStoreFloat2(&xmf2norm, xmv2norm);
        DirectX::XMStoreFloat(&xmf2dot, xmv2dot);
        DirectX::XMStoreFloat(&xmf2len, xmv2len);
        DirectX::XMStoreFloat(&xmf2lensq, xmv2lensq);

        float2 f20, f21, f2sum, f2diff, f2prod, f2norm;
        float f2dot, f2len, f2lensq;
        f20 = float2(v20);
        f21 = float2(v21);
        f2sum = f20 + f21;
        f2diff = f20 - f21;
        f2prod = f20 * f21;
        f2norm = normalize(f20);
        f2dot = dot(f20, f21);
        f2len = length(f20);
        f2lensq = length_squared(f20);

        correct = Equals<float, 1, 2>(xmv2sum.m128_f32, f2sum.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 2>(xmv2diff.m128_f32, f2diff.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 2>(xmv2prod.m128_f32, f2prod.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 2>(xmv2norm.m128_f32, f2norm.p_cols, eps);
        assert(correct);

        correct = xmf2dot == f2dot;
        assert(correct);

        correct = xmf2len == f2len;
        assert(correct);

        correct = xmf2lensq == f2lensq;
        assert(correct);
    }

    // float3
    {
        float v30[3], v31[3];
        FillMatrix<float, 1, 3>(gen, v30);
        FillMatrix<float, 1, 3>(gen, v31);

        DirectX::XMVECTOR xmv30, xmv31, xmv3sum, xmv3diff, xmv3prod, xmv3norm, xmv3dot, xmv3len, xmv3lensq;
        xmv30 = DirectX::XMVectorSet(v30[0], v30[1], v30[2], 0);
        xmv31 = DirectX::XMVectorSet(v31[0], v31[1], v31[2], 0);
        xmv3sum = DirectX::XMVectorAdd(xmv30, xmv31);
        xmv3diff = DirectX::XMVectorSubtract(xmv30, xmv31);
        xmv3prod = DirectX::XMVectorMultiply(xmv30, xmv31);
        xmv3norm = DirectX::XMVector3Normalize(xmv30);
        xmv3dot = DirectX::XMVector3Dot(xmv30, xmv31);
        xmv3len = DirectX::XMVector3Length(xmv30);
        xmv3lensq = DirectX::XMVector3LengthSq(xmv30);

        DirectX::XMFLOAT3 xmf3sum, xmf3diff, xmf3prod, xmf3norm;
        float xmf3dot, xmf3len, xmf3lensq;
        DirectX::XMStoreFloat3(&xmf3sum, xmv3sum);
        DirectX::XMStoreFloat3(&xmf3diff, xmv3diff);
        DirectX::XMStoreFloat3(&xmf3prod, xmv3prod);
        DirectX::XMStoreFloat3(&xmf3norm, xmv3norm);
        DirectX::XMStoreFloat(&xmf3dot, xmv3dot);
        DirectX::XMStoreFloat(&xmf3len, xmv3len);
        DirectX::XMStoreFloat(&xmf3lensq, xmv3lensq);

        float3 f30, f31, f3sum, f3diff, f3prod, f3norm;
        float f3dot, f3len, f3lensq;
        f30 = float3(v30);
        f31 = float3(v31);
        f3sum = f30 + f31;
        f3diff = f30 - f31;
        f3prod = f30 * f31;
        f3norm = normalize(f30);
        f3dot = dot(f30, f31);
        f3len = length(f30);
        f3lensq = length_squared(f30);

        correct = Equals<float, 1, 3>(xmv3sum.m128_f32, f3sum.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 3>(xmv3diff.m128_f32, f3diff.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 3>(xmv3prod.m128_f32, f3prod.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 3>(xmv3norm.m128_f32, f3norm.p_cols, eps);
        assert(correct);

        correct = fabsf(xmf3dot - f3dot) < 0.001f;
        assert(correct);

        correct = fabsf(xmf3len - f3len) < 0.001f;
        assert(correct);

        correct = fabsf(xmf3lensq - f3lensq) < 0.001f;
        assert(correct);
    }

    // float4
    {
        float v40[4], v41[4];
        FillMatrix<float, 1, 4>(gen, v40);
        FillMatrix<float, 1, 4>(gen, v41);

        DirectX::XMVECTOR xmv40, xmv41, xmv4sum, xmv4diff, xmv4prod, xmv4norm, xmv4dot, xmv4len, xmv4lensq;
        xmv40 = DirectX::XMVectorSet(v40[0], v40[1], v40[2], v40[3]);
        xmv41 = DirectX::XMVectorSet(v41[0], v41[1], v41[2], v41[3]);
        xmv4sum = DirectX::XMVectorAdd(xmv40, xmv41);
        xmv4diff = DirectX::XMVectorSubtract(xmv40, xmv41);
        xmv4prod = DirectX::XMVectorMultiply(xmv40, xmv41);
        xmv4norm = DirectX::XMVector4Normalize(xmv40);
        xmv4dot = DirectX::XMVector4Dot(xmv40, xmv41);
        xmv4len = DirectX::XMVector4Length(xmv40);
        xmv4lensq = DirectX::XMVector4LengthSq(xmv40);

        DirectX::XMFLOAT4 xmf4sum, xmf4diff, xmf4prod, xmf4norm;
        float xmf4dot, xmf4len, xmf4lensq;
        DirectX::XMStoreFloat4(&xmf4sum, xmv4sum);
        DirectX::XMStoreFloat4(&xmf4diff, xmv4diff);
        DirectX::XMStoreFloat4(&xmf4prod, xmv4prod);
        DirectX::XMStoreFloat4(&xmf4norm, xmv4norm);
        DirectX::XMStoreFloat(&xmf4dot, xmv4dot);
        DirectX::XMStoreFloat(&xmf4len, xmv4len);
        DirectX::XMStoreFloat(&xmf4lensq, xmv4lensq);

        float4 f40, f41, f4sum, f4diff, f4prod, f4norm;
        float f4dot, f4len, f4lensq;
        f40 = float4(v40);
        f41 = float4(v41);
        f4sum = f40 + f41;
        f4diff = f40 - f41;
        f4prod = f40 * f41;
        f4norm = normalize(f40);
        f4dot = dot(f40, f41);
        f4len = length(f40);
        f4lensq = length_squared(f40);

        correct = Equals<float, 1, 4>(xmv4sum.m128_f32, f4sum.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 4>(xmv4diff.m128_f32, f4diff.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 4>(xmv4prod.m128_f32, f4prod.p_cols, eps);
        assert(correct);

        correct = Equals<float, 1, 4>(xmv4norm.m128_f32, f4norm.p_cols, eps);
        assert(correct);

        correct = fabsf(xmf4dot - f4dot) < 0.001f;
        assert(correct);

        correct = fabsf(xmf4len - f4len) < 0.001f;
        assert(correct);

        correct = fabsf(xmf4lensq - f4lensq) < 0.001f;
        assert(correct);
    }

    // float2x2
    {
        float v220[4], v221[4];
        FillMatrix<float, 2, 2>(gen, v220);
        FillMatrix<float, 2, 2>(gen, v221);

        float2x2 f220, f221, f22sum, f22diff, f22prod, f22trans, f22inv;
        float f22det;
        f220 = float2x2(v220);
        f221 = float2x2(v221);

        f22sum = f220 + f221;
        f22diff = f220 - f221;
        f22prod = f220 * f221;
        f22trans = f220.transpose();
        f22inv = f220.inverse();
        f22det = f220.determinant();
    }

    // float3x3
    {
        float v330[9], v331[9];
        FillMatrix<float, 3, 3>(gen, v330);
        FillMatrix<float, 3, 3>(gen, v331);

        DirectX::XMFLOAT3X3 xmf330, xmf331, xmf33sum, xmf33diff, xmf33prod, xmf33trans, xmf33inv;
        float xmf33det;
        memcpy_s(xmf330.m, sizeof(float) * 9, v330, sizeof(float) * 9);
        memcpy_s(xmf331.m, sizeof(float) * 9, v331, sizeof(float) * 9);

        DirectX::XMMATRIX xmm330, xmm331, xmm33sum, xmm33diff, xmm33prod, xmm33trans, xmm33inv;
        xmm330 = DirectX::XMLoadFloat3x3(&xmf330);
        xmm331 = DirectX::XMLoadFloat3x3(&xmf331);
        xmm33sum = xmm330 + xmm331;
        xmm33diff = xmm330 - xmm331;
        xmm33prod = xmm330 * xmm331;
        xmm33trans = DirectX::XMMatrixTranspose(xmm330);
        auto det = DirectX::XMMatrixDeterminant(xmm330);
        xmm33inv = DirectX::XMMatrixInverse(&det, xmm330);

        DirectX::XMStoreFloat3x3(&xmf33sum, xmm33sum);
        DirectX::XMStoreFloat3x3(&xmf33diff, xmm33diff);
        DirectX::XMStoreFloat3x3(&xmf33prod, xmm33prod);
        DirectX::XMStoreFloat3x3(&xmf33trans, xmm33trans);
        DirectX::XMStoreFloat3x3(&xmf33inv, xmm33inv);
        DirectX::XMStoreFloat(&xmf33det, det);

        float3x3 f330, f331, f33sum, f33diff, f33prod, f33trans, f33inv;
        float f33det;
        f330 = float3x3(v330);
        f331 = float3x3(v331);

        f33sum = f330 + f331;
        f33diff = f330 - f331;
        f33prod = f330 * f331;
        f33trans = f330.transpose();
        f33inv = f330.inverse();
        f33det = f330.determinant();

        correct = Equals<float, 3, 3>(f33sum.p_data, &xmf33sum.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 3, 3>(f33diff.p_data, &xmf33diff.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 3, 3>(f33prod.p_data, &xmf33prod.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 3, 3>(f33trans.p_data, &xmf33trans.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 3, 3>(f33inv.p_data, &xmf33inv.m[0][0], eps);
        assert(correct);

        correct = abs(xmf33det - f33det) < 1.0;
        assert(correct);
    }

    // float4x3
    {
        float v430[12], v431[12];
        FillMatrix<float, 4, 3>(gen, v430);
        FillMatrix<float, 4, 3>(gen, v431);

        DirectX::XMFLOAT4X3 xmf430, xmf431, xmf43sum, xmf43diff, xmf43prod, xmf43trans, xmf43inv;
        float xmf43det;
        memcpy_s(xmf430.m, sizeof(float) * 12, v430, sizeof(float) * 12);
        memcpy_s(xmf431.m, sizeof(float) * 12, v431, sizeof(float) * 12);

        DirectX::XMMATRIX xmm430, xmm431, xmm43sum, xmm43diff, xmm43prod, xmm43trans, xmm43inv;
        xmm430 = DirectX::XMLoadFloat4x3(&xmf430);
        xmm431 = DirectX::XMLoadFloat4x3(&xmf431);
        xmm43sum = xmm430 + xmm431;
        xmm43diff = xmm430 - xmm431;
        xmm43prod = xmm430 * xmm431;
        xmm43trans = DirectX::XMMatrixTranspose(xmm430);
        auto det = DirectX::XMMatrixDeterminant(xmm430);
        xmm43inv = DirectX::XMMatrixInverse(&det, xmm430);

        DirectX::XMStoreFloat4x3(&xmf43sum, xmm43sum);
        DirectX::XMStoreFloat4x3(&xmf43diff, xmm43diff);
        DirectX::XMStoreFloat4x3(&xmf43prod, xmm43prod);
        DirectX::XMStoreFloat4x3(&xmf43trans, xmm43trans);
        DirectX::XMStoreFloat4x3(&xmf43inv, xmm43inv);
        DirectX::XMStoreFloat(&xmf43det, det);

        float4x3 f430, f431, f43sum, f43diff, f43prod, f43trans, f43inv;
        float f43det;
        f430 = float4x3(v430);
        f431 = float4x3(v431);

        f43sum = f430 + f431;
        f43diff = f430 - f431;
        f43prod = f430 * f431;
        f43trans = f430.transpose();
        f43inv = f430.inverse();
        f43det = f430.determinant();

        correct = Equals<float, 4, 3>(f43sum.p_data, &xmf43sum.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 4, 3>(f43diff.p_data, &xmf43diff.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 4, 3>(f43prod.p_data, &xmf43prod.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 4, 3>(f43trans.p_data, &xmf43trans.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 4, 3>(f43inv.p_data, &xmf43inv.m[0][0], eps);
        assert(correct);

        correct = abs(xmf43det - f43det) < 1.0;
        assert(correct);
    }

    // float4x4
    {
        float v440[16], v441[16];
        FillMatrix<float, 4, 4>(gen, v440);
        FillMatrix<float, 4, 4>(gen, v441);

        DirectX::XMFLOAT4X4 xmf440, xmf441, xmf44sum, xmf44diff, xmf44prod, xmf44trans, xmf44inv;
        float xmf44det;
        memcpy_s(xmf440.m, sizeof(float) * 16, v440, sizeof(float) * 16);
        memcpy_s(xmf441.m, sizeof(float) * 16, v441, sizeof(float) * 16);

        DirectX::XMMATRIX xmm440, xmm441, xmm44sum, xmm44diff, xmm44prod, xmm44trans, xmm44inv;
        xmm440 = DirectX::XMLoadFloat4x4(&xmf440);
        xmm441 = DirectX::XMLoadFloat4x4(&xmf441);
        xmm44sum = xmm440 + xmm441;
        xmm44diff = xmm440 - xmm441;
        xmm44prod = xmm440 * xmm441;
        xmm44trans = DirectX::XMMatrixTranspose(xmm440);
        auto det = DirectX::XMMatrixDeterminant(xmm440);
        xmm44inv = DirectX::XMMatrixInverse(&det, xmm440);

        DirectX::XMStoreFloat4x4(&xmf44sum, xmm44sum);
        DirectX::XMStoreFloat4x4(&xmf44diff, xmm44diff);
        DirectX::XMStoreFloat4x4(&xmf44prod, xmm44prod);
        DirectX::XMStoreFloat4x4(&xmf44trans, xmm44trans);
        DirectX::XMStoreFloat4x4(&xmf44inv, xmm44inv);
        DirectX::XMStoreFloat(&xmf44det, det);

        float4x4 f440, f441, f44sum, f44diff, f44prod, f44trans, f44inv;
        float f44det;
        f440 = float4x4(v440);
        f441 = float4x4(v441);

        f44sum = f440 + f441;
        f44diff = f440 - f441;
        f44prod = f440 * f441;
        f44trans = f440.transpose();
        f44inv = f440.inverse();
        f44det = f440.determinant();

        correct = Equals<float, 4, 4>(f44sum.p_data, &xmf44sum.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 4, 4>(f44diff.p_data, &xmf44diff.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 4, 4>(f44prod.p_data, &xmf44prod.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 4, 4>(f44trans.p_data, &xmf44trans.m[0][0], eps);
        assert(correct);

        correct = Equals<float, 4, 4>(f44inv.p_data, &xmf44inv.m[0][0], eps);
        assert(correct);

        auto d = abs(xmf44det - f44det);
        correct = d < 1.0;
        //assert(correct);
    }
}

void RunSpeedTests(mt19937& gen)
{
    uniform_real_distribution<> dis(0.0, 1000.0f);
    int iterations = 3000;

    double adds[2] = { 0, 0 };
    double subs[2] = { 0, 0 };
    double muls[2] = { 0, 0 };
    double nrms[2] = { 0, 0 };
    double crss[2] = { 0, 0 };
    double dots[2] = { 0, 0 };
    double lens[2] = { 0, 0 };
    double vmul[2] = { 0, 0 };
    double trns[2] = { 0, 0 };

    for (int i = 0; i < iterations; ++i)
    {
        float a0[3], a1[3], a2[12];
        FillMatrix<float, 1, 3>(gen, a0);
        FillMatrix<float, 1, 3>(gen, a1);
        FillMatrix<float, 4, 3>(gen, a2);

        float3 f0, f1;
        f0 = float3(a0);
        f1 = float3(a1);
        float4x3 f2(a2);

        DirectX::XMVECTOR xm0, xm1;
        DirectX::XMVectorSet(a0[0], a0[1], a0[2], 0.0f);
        DirectX::XMVectorSet(a1[0], a1[1], a1[2], 0.0f);
        DirectX::XMFLOAT4X3 xmf;
        memcpy_s(xmf.m, sizeof(float) * 12, a2, sizeof(float) * 12);
        DirectX::XMMATRIX xm2 = DirectX::XMLoadFloat4x3(&xmf);

        adds[0] += TimeFunction("float3::operator+", [&f0, &f1] { float3 fo = f0 + f1; });
        adds[1] += TimeFunction("XMVector::XMVectorAdd", [&xm0, &xm1] {DirectX::XMVECTOR xo = DirectX::XMVectorAdd(xm0, xm1); });

        subs[0] += TimeFunction("float3::operator-", [&f0, &f1] { float3 fo = f0 - f1; });
        subs[1] += TimeFunction("XMVector::XMVectorSubtract", [&xm0, &xm1] {DirectX::XMVECTOR xo = DirectX::XMVectorSubtract(xm0, xm1); });

        muls[0] += TimeFunction("float3::operator*", [&f0, &f1] { float3 fo = f0 * f1; });
        muls[1] += TimeFunction("XMVector::XMVectorMultiply", [&xm0, &xm1] {DirectX::XMVECTOR xo = DirectX::XMVectorMultiply(xm0, xm1); });

        nrms[0] += TimeFunction("float3::normalize", [&f0, &f1] { float3 fo = normalize(f0); });
        nrms[1] += TimeFunction("XMVector::XMVector3Normalize", [&xm0, &xm1] {DirectX::XMVECTOR xo = DirectX::XMVector3Normalize(xm0); });

        crss[0] += TimeFunction("float3::cross", [&f0, &f1] { float3 fo = cross(f0, f1); });
        crss[1] += TimeFunction("XMVector::XMVector3Cross", [&xm0, &xm1] {DirectX::XMVECTOR xo = DirectX::XMVector3Cross(xm0, xm1); });

        dots[0] += TimeFunction("float3::dot", [&f0, &f1] { float fo = dot(f0, f1); });
        dots[1] += TimeFunction("XMVector::XMVector3Dot", [&xm0, &xm1] {DirectX::XMVECTOR xo = DirectX::XMVector3Dot(xm0, xm1); });

        lens[0] += TimeFunction("float3::length", [&f0, &f1] { float fo = length(f0); });
        lens[1] += TimeFunction("XMVector::XMVector3Length", [&xm0, &xm1] {DirectX::XMVECTOR xo = DirectX::XMVector3Length(xm0); });

        vmul[0] += TimeFunction("float3x3::operator*", [&f0, &f2] { float3 fo = f0 * f2; });
        vmul[1] += TimeFunction("XMVector::XMVector3Transform", [&xm0, &xm2] {DirectX::XMVECTOR xo = DirectX::XMVector3Transform(xm0, xm2); });

        trns[0] += TimeFunction("float3x3::operator*", [&f2] { float4x3 fo = f2 * f2; });
        trns[1] += TimeFunction("XMVector::XMMatrixMultiply", [&xm0, &xm2] {DirectX::XMMATRIX xo = DirectX::XMMatrixMultiply(xm2, xm2); });
    }

    double inv = 1.0 / static_cast<double>(iterations);
    for (int i = 0; i < 2; ++i) adds[i] *= inv;
    for (int i = 0; i < 2; ++i) subs[i] *= inv;
    for (int i = 0; i < 2; ++i) muls[i] *= inv;
    for (int i = 0; i < 2; ++i) nrms[i] *= inv;
    for (int i = 0; i < 2; ++i) crss[i] *= inv;
    for (int i = 0; i < 2; ++i) dots[i] *= inv;
    for (int i = 0; i < 2; ++i) lens[i] *= inv;
    for (int i = 0; i < 2; ++i) vmul[i] *= inv;
    for (int i = 0; i < 2; ++i) trns[i] *= inv;

    printf("float3::operator+:\t\t%f\n", adds[0]);
    printf("XMVector::XMVectorAdd:\t\t%f\n", adds[1]);
    printf("\n");
    printf("float3::operator-:\t\t%f\n", subs[0]);
    printf("XMVector::XMVectorSubtract:\t%f\n", subs[1]);
    printf("\n");
    printf("float3::operator*:\t\t%f\n", muls[0]);
    printf("XMVector::XMVectorMultiply:\t%f\n", muls[1]);
    printf("\n");
    printf("float3::normalize:\t\t%f\n", nrms[0]);
    printf("XMVector::XMVector3Normalize:\t%f\n", nrms[1]);
    printf("\n");
    printf("float3::cross:\t\t\t%f\n", crss[0]);
    printf("XMVector::XMVector3Cross:\t%f\n", crss[1]);
    printf("\n");
    printf("float3::dot:\t\t\t%f\n", dots[0]);
    printf("XMVector::XMVector3Dot:\t\t%f\n", dots[1]);
    printf("\n");
    printf("float3::length:\t\t\t%f\n", lens[0]);
    printf("XMVector::XMVector3Length:\t%f\n", lens[1]);
    printf("\n");
    printf("float4x3::operator*:\t\t%f\n", vmul[0]);
    printf("XMVector::XMVector3Transform:\t%f\n", vmul[1]);
    printf("\n");
    printf("float4x3::operator*:\t\t%f\n", trns[0]);
    printf("XMVector::XMMatrixMultiply:\t%f\n", trns[1]);
}

int main(int argc, int* argv[])
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1000.0f);

    printf("Running correctness tests...\n");
    RunCorrectnessTests(gen);
    printf("Done\n");

    printf("Running speed tests... \n");
    RunSpeedTests(gen);
    printf("Done\n");

    getchar();

    return 1;
}