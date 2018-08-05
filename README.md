# lina
Lina includes templated types and methods for performing vector and matrix operations on n-dimensional arrays. Specializations exist for 2D, 3D, and 4D operations.

Usage:
```
float2 lf = { 1.0f, 3.0f };
float2 rf = { 4.0f, 5.0f };
float2 sum = lf + rf;

float2x2 rmf = { 2.0f, 0.0f, 0.0f, 2.0f };
float2 prodf = lf * rf;

float2x2 lmf = { 3.0f, 0.0f, 0.0f, 3.0f };
float2x3 pmf = lmf * rmf; 
```

Template types and methods exist for testing overlaps on areas and volumes.

Usage:
```
ray<2> r(0.0f, 0.0f, 1.0f, 1.0f); // origin, and direction
aabb<2> box(5.0f, 5.0f, 1.0f, 1.0f); // origin and extents

float poi[2], t;
if (overlap(r, box, poi, t))
{
  // overlap exists and poi contains point of intersection
}
```
