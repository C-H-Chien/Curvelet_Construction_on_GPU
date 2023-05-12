
#include "indices.hpp"
#include <stdio.h>
#include <math.h>

__device__ __inline__ float
To2Pi (float angle)
{
  #define M_PI 3.14159265358979323846
  float a;
  if (angle>=2*M_PI)
    a = fmodf(angle, (float)2*M_PI);
  else if (angle < 0)
    a = (2*M_PI + fmodf(angle, (float)2*M_PI));
  else 
    a= angle;

  // added by Nhon: these two lines of code is to fix the bug when
  // angle = -1.1721201390607859e-016
  // then after all the computation, we get
  // a = 6.2831853071795862 == 2*vnl_math::pi !!!!!!!
  // the only case this can happen is when a is very close to zero.
  if (!(a>=0 && a<2*M_PI)) {
    a = 0;
  }

  return a;
}

__device__ __inline__ float 
angle_from_pt_to_pt (float x1, float y1, float x2, float y2)
{
  return To2Pi (atan2f(y2 - y1, x2 - x1) );
}

__device__ __inline__ float 
dot (float v1, float v2)
{
  return cosf(v1)*cosf(v2) + sinf(v1)*sinf(v2);
}

