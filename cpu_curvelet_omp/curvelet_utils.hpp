#ifndef CURVELET_UTILS_HPP
#define CURVELET_UTILS_HPP

/***************************************************************************************
//file: curvelet_utils.hpp
//brief: utility functions for the curvelet grouping algorithm
//author: Chiang-Heng Chien
***************************************************************************************/

#define M_PI 3.14159265358979323846
#include <cmath>

template<typename T>
T sq_dist(T x1, T y1, T x2, T y2)
{
  return (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
}
/*
float sq_norm(const point_2d& A)
{
    return A.first*A.first+A.second*A.second;
}
*/
//: Convert an angle to [0, 2Pi) range
template<typename T>
T To2Pi(T angle)
{
  T a;
  if (angle>=T(2)*T(M_PI))
    a = std::fmod(angle, T(2)*T(M_PI));
  else if (angle < 0)
    a = (T(2)*T(M_PI)+std::fmod(angle, T(2)*T(M_PI)));
  else 
    a= angle;

  // added by Nhon: these two lines of code is to fix the bug when
  // angle = -1.1721201390607859e-016
  // then after all the computation, we get
  // a = 6.2831853071795862 == 2*vnl_math::pi !!!!!!!
  // the only case this can happen is when a is very close to zero.
  if (!(a>=0 && a<T(2)*T(M_PI))) {
    a = 0;
  }

  return a;
}

//: Convert an angle to [-Pi, Pi) range
template<typename T>
T ToPi(T angle)
{
    T a = angle+T(M_PI);
    a = To2Pi(a);
    
    return a-T(M_PI);
}

template<typename T>
T angle_from_pt_to_pt (T x1, T y1, T x2, T y2)
{
  return To2Pi (std::atan2(y2 - y1, x2 - x1) );
}

//: dot product between two angle
template<typename T>
T dot (T v1, T v2)
{
  return std::cos(v1)*std::cos(v2) + std::sin(v1)*std::sin(v2);
}
/*
//: rotate a vector by theta
point_2d rotate(const point_2d &pt, float theta)
{
    return point_2d(cos(theta)*pt.first-sin(theta)*pt.second, sin(theta)*pt.first+cos(theta)*pt.second);
}
*/

#endif  // CURVELET_UTILS_HPP
