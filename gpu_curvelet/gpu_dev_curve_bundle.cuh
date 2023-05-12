
#include "indices.hpp"
#include <utility> 
#include <iostream>
#include <cmath>

template<int _curves_num_in_bundle_pixel, int _curves_num_in_bundle_theta>
__device__ __inline__ bool
compute_curve_bundle( float _dx, float _dt, float _sx, float _st,
                      int le_idx, float* sBundle_min_ks, float* sBundle_max_ks, bool* sHyp_LookEdge,
                      float te_pt_x, float te_pt_y, float te_orient, float le_pt_x, float le_pt_y, float le_orient
) 
{  
    //> global indices to fill values in sBundle_min_ks and sBundle_max_ks fill
    unsigned bk_idx = 0;

    //> \Delta X, \Delta Y, and their squares
    const float DX = te_pt_x - le_pt_x;
    const float DY = te_pt_y - le_pt_y;
    const float DX2 = DX*DX;
    const float DY2 = DY*DY;

    //> some pre-calculated constants
    const float T0 = le_orient;
    const float T2 = te_orient;
    const float _dx2 = _dx * _dx;
    const float T0_p_dt = T0 + _dt;
    const float T0_m_dt = T0 - _dt;

    //> non-constant variables computed on the fly
    // TODO: change variable name to increase readability
    float k_dx_min, k_dx_max;
    float k_dt_min, k_dt_max;
    float delta_x2;
    float delta_t2;
    float T2t, tt;
    float eq_denom;
    float ttmdx2;
    float sin_t2_p_dt2_m_t0_m_dt;
    float sin_t2_p_dt2_m_t0_p_dt;
    float sin_T0pdt_times_DY;
    float cos_T0pdt_times_DX;
    float sin_T0mdt_times_DY;
    float cos_T0mdt_times_DX;

    //> if neighboring point is close to the target edgel, no static constraint can be computed of leave it as the default bundle
    //  and the max and min curvatures of the curve bundle remains as default
    if ( DX2 + DY2 < 4*_dx*_dx ) {
        sHyp_LookEdge[ le_idx ] = true;
        return true;
    }

    //> Proposition 4 in the paper: Circular Arc Bundle Transport
    if ( (DX * cos(T2) + DY * sin(T2)) < 0 ) {

        //> loop over all the possible curves to find the max and min curvature
        for (unsigned i = 0; i < _curves_num_in_bundle_pixel; i++) {
            for (unsigned j = 0; j < _curves_num_in_bundle_theta; j++) {

                //> calculate the relative difference in terms of pixel and radian
                delta_x2 = _sx * ((float)i - (float(_curves_num_in_bundle_pixel) - 1.0) / 2.0);
                delta_t2 = _st * ((float)j - (float(_curves_num_in_bundle_theta) - 1.0) / 2.0);

                //> ---------------- i) With Respect To delta x --------------
                //> absolute grid orientation of the curve
                T2t = T2 + delta_t2;

                //> intermediate calculations
                tt = sin(T2t) * DX - cos(T2t) * DY;
                ttmdx2 = tt - delta_x2;

                //> k_max / k_min equation denominator w.r.float. pixel
                eq_denom = (-2.0 * delta_x2 * tt + DX2 + DY2 - _dx2 + delta_x2 * delta_x2);

                //> max and min curvatures in terms of pixel
                k_dx_min = 2.0 * (ttmdx2 - _dx) / eq_denom;
                k_dx_max = 2.0 * (ttmdx2 + _dx) / eq_denom;

                //> compute the intersection of all surfaces
                //> 1) min
                if ( sBundle_min_ks(le_idx, bk_idx) < k_dx_min )
                    sBundle_min_ks(le_idx, bk_idx) = k_dx_min;

                //> 2) max
                if ( sBundle_max_ks(le_idx, bk_idx) > k_dx_max )
                    sBundle_max_ks(le_idx, bk_idx) = k_dx_max;
                
                //> ----------------- ii) With Respect To delta theta ------------------
                //> intermediate calculations
                sin_t2_p_dt2_m_t0_m_dt = sin(T2t - T0_m_dt);
                sin_t2_p_dt2_m_t0_p_dt = sin(T2t - T0_p_dt);
                sin_T0pdt_times_DY     = DY * sin(T0_p_dt);
                cos_T0pdt_times_DX     = DX * cos(T0_p_dt);
                sin_T0mdt_times_DY     = DY * sin(T0_m_dt);
                cos_T0mdt_times_DX     = DX * cos(T0_m_dt);

                //> max and min curvatures in terms of orientation
                k_dt_min = sin_t2_p_dt2_m_t0_m_dt / (-sin_t2_p_dt2_m_t0_m_dt * delta_x2 + cos_T0mdt_times_DX + sin_T0mdt_times_DY);
                k_dt_max = sin_t2_p_dt2_m_t0_p_dt / (-sin_t2_p_dt2_m_t0_p_dt * delta_x2 + cos_T0pdt_times_DX + sin_T0pdt_times_DY);

                /*if (te_id == 0 && le_id == 7) {
                    std::cout<<bk_idx+1<<" (k_dt_min, k_dt_max) = ("<<k_dt_min<<", "<<k_dt_max<<")"<<std::endl;      
                }*/

                //> compute the intersection of all surfaces
                //> 1) min
                if ( sBundle_min_ks(le_idx, bk_idx) < k_dt_min )
                    sBundle_min_ks(le_idx, bk_idx) = k_dt_min;
                
                //> 2) max
                if ( sBundle_max_ks(le_idx, bk_idx) > k_dt_max )
                    sBundle_max_ks(le_idx, bk_idx) = k_dt_max;

                bk_idx++;
            }
        }
    }
    else {
        //> loop over all the possible curves to find the max and min curvature
        for (unsigned i = 0; i < _curves_num_in_bundle_pixel; i++) {
            for (unsigned j = 0; j < _curves_num_in_bundle_theta; j++) {

                //> calculate the relative difference in terms of pixel and radian
                delta_x2 = _sx * ((float)i - (float(_curves_num_in_bundle_pixel) - 1.0) / 2.0);
                delta_t2 = _st * ((float)j - (float(_curves_num_in_bundle_theta) - 1.0) / 2.0);

                //> -------------- i) With Respect To delta x ----------------
                //> absolute grid orientation of the curve
                T2t = T2 + delta_t2;

                //> intermediate calculations
                tt = sin(T2t) * DX - cos(T2t) * DY;
                ttmdx2 = tt - delta_x2;

                //> k_max / k_min equation denominator w.r.float. pixel
                eq_denom = (-2.0 * delta_x2 * tt + DX2 + DY2 - _dx2 + delta_x2 * delta_x2);

                //> max and min curvatures in terms of pixel
                k_dx_min = 2.0 * (ttmdx2 - _dx) / eq_denom;
                k_dx_max = 2.0 * (ttmdx2 + _dx) / eq_denom;

                //> compute the intersection of all surfaces
                //> 1) min
                if ( sBundle_min_ks(le_idx, bk_idx) < k_dx_min )
                    sBundle_min_ks(le_idx, bk_idx) = k_dx_min;

                //> 2) max
                if ( sBundle_max_ks(le_idx, bk_idx) > k_dx_max )
                    sBundle_max_ks(le_idx, bk_idx) = k_dx_max;

                //> --------------- ii) With Respect To delta theta ---------------
                //> intermediate calculations
                sin_t2_p_dt2_m_t0_m_dt = sin(T2t - T0_m_dt);
                sin_t2_p_dt2_m_t0_p_dt = sin(T2t - T0_p_dt);
                sin_T0pdt_times_DY     = DY * sin(T0_p_dt);
                cos_T0pdt_times_DX     = DX * cos(T0_p_dt);
                sin_T0mdt_times_DY     = DY * sin(T0_m_dt);
                cos_T0mdt_times_DX     = DX * cos(T0_m_dt);
                
                //> max and min curvatures in terms of orientation
                k_dt_max = sin_t2_p_dt2_m_t0_m_dt/ (-sin_t2_p_dt2_m_t0_m_dt * delta_x2 + cos_T0mdt_times_DX + sin_T0mdt_times_DY);
                k_dt_min = sin_t2_p_dt2_m_t0_p_dt/ (-sin_t2_p_dt2_m_t0_p_dt * delta_x2 + cos_T0pdt_times_DX + sin_T0pdt_times_DY);

                //> compute the intersection of all surfaces
                //> 1) min
                if ( sBundle_min_ks(le_idx, bk_idx) < k_dt_min )
                    sBundle_min_ks(le_idx, bk_idx) = k_dt_min;
                
                //> 2) max
                if ( sBundle_max_ks(le_idx, bk_idx) > k_dt_max )
                    sBundle_max_ks(le_idx, bk_idx) = k_dt_max;

                bk_idx++;
            }
        }
    }

    //> check whether the built bundle is valid or not
    //bool valid = bundle_valid_check< _curves_num_in_bundle_pixel, _curves_num_in_bundle_theta >( le_idx, sBundle_min_ks, sBundle_max_ks );
    bool valid = false;
    //> check the curve bundle is valid or not
    for (int k_idx = 0; k_idx < (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta) ; k_idx++ ) {
        if ( sBundle_max_ks(le_idx, k_idx) > sBundle_min_ks(le_idx, k_idx) ) 
            valid = true;
    }

    //> this boolean list indicates which look edge is valid to be proceeded
    sHyp_LookEdge[ le_idx ] = (valid) ? true : false;

    return (valid) ? true : false;
}

//> keep curvelet bundle in a cmp bundle region
template<int _curves_num_in_bundle_pixel, int _curves_num_in_bundle_theta, int cmp_bundle_total_size, int intersect_bundle_total_size>
__device__ __inline__ void
move_to_cmp_bundle( 
    int cmp_idx, int le_idx, bool rep_by_intersection,
    float rCmp_bundle_min_ks[cmp_bundle_total_size], float rCmp_bundle_max_ks[cmp_bundle_total_size],
    float rIntersect_bundle_min_ks[intersect_bundle_total_size], float rIntersect_bundle_max_ks[intersect_bundle_total_size], 
    float* sBundle_min_ks, float* sBundle_max_ks
)
{
    //> keep the curvelet bundle
    //> cmp_idx is either 0 or 1
    if (!rep_by_intersection) {
        for (unsigned bidx = 0; bidx < (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta); bidx++) {
            rCmp_bundle_min_ks(cmp_idx, bidx) = sBundle_min_ks(le_idx, bidx);
            rCmp_bundle_max_ks(cmp_idx, bidx) = sBundle_max_ks(le_idx, bidx);
        }
    }
    else {
        for (unsigned bidx = 0; bidx < (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta); bidx++) {
            rCmp_bundle_min_ks(cmp_idx, bidx) = rIntersect_bundle_min_ks[ bidx ];
            rCmp_bundle_max_ks(cmp_idx, bidx) = rIntersect_bundle_max_ks[ bidx ];
        }
    }
}

//> compute bundle intersection device function
template<int _curves_num_in_bundle_pixel, int _curves_num_in_bundle_theta, int cmp_bundle_total_size, int intersect_bundle_total_size>
__device__ __inline__ void
bundle_intersection(float rCmp_bundle_min_ks[cmp_bundle_total_size], float rCmp_bundle_max_ks[cmp_bundle_total_size],
                    float rIntersect_bundle_min_ks[intersect_bundle_total_size], float rIntersect_bundle_max_ks[intersect_bundle_total_size]
)
{
    //> compare the curvatures in the rCmp_bundle_min_ks and rCmp_bundle_max_ks
    #pragma unroll
    for (unsigned j = 0; j < (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta) ; j++) {
        //> compare the min curvature, pick the larger one
        rIntersect_bundle_min_ks[ j ] = (rCmp_bundle_min_ks(0, j) < rCmp_bundle_min_ks(1, j)) ? (rCmp_bundle_min_ks(1, j)) : (rCmp_bundle_min_ks(0, j));

        //> compare the max curvature, pick the smaller one
        rIntersect_bundle_max_ks[ j ] = (rCmp_bundle_max_ks(0, j) > rCmp_bundle_max_ks(1, j)) ? (rCmp_bundle_max_ks(1, j)) : (rCmp_bundle_max_ks(0, j));
    }
}

//> check validity of bundle intersection device function
template<int _curves_num_in_bundle_pixel, int _curves_num_in_bundle_theta, int intersect_bundle_total_size>
__device__ __inline__ bool
bundle_intersection_valid_check(
    float rIntersect_bundle_min_ks[intersect_bundle_total_size], float rIntersect_bundle_max_ks[intersect_bundle_total_size]
)
{
    //> check the intersected curve bundle is valid or not
    for (int k_idx = 0; k_idx < (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta) ; k_idx++ ) {
        if ( rIntersect_bundle_max_ks[k_idx] > rIntersect_bundle_min_ks[k_idx] ) 
            return true;
    }
    return false;
}