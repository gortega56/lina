# LINear Algebra Templates
This library is focused on implementing formulas found in linear alebra, specifically vector and matrix operations. It is not fast, just correct and useful for testing more optimized approaches. 

Usage:
```
float c = cosf(0.524f);
float s = sinf(0.524f);
matrix<float, 4, 4> rotation_axis_z = 
    matrix_set<float, 4, 4>(
          c,   s, 0.f, 0.f,
         -s,   c, 0.f, 0.f,
        0.f, 0.f, 1.f, 0.f,
        0.f, 0.f, 0.f, 1.f);
		
vector<float, 4> position = vector_set<float, 4>(1.0f, 0.0f, 0.0f, 1.0f);
vector<float, 4> position_rotated = position * rotation_axis_z;
```
