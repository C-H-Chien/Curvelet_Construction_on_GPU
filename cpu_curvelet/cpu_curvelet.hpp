#ifndef CPU_CURVELET_HPP
#define CPU_CURVELET_HPP

#include <map>
#include <vector>
#include <deque>
#include <list>
#include <utility> 
#include <iostream>
#include <cmath>

#include "indices.hpp"
#include "curvelet_utils.hpp"

template<typename T>
class CurveletCPU
{
protected:
    int _sz_edge_data;
    int _num_edges;
    unsigned _group_max_sz;
    int _max_num_look_edges;

    // > target edge information
    T te_id;
    T te_pt_x;
    T te_pt_y;
    T te_orient;
    T te_grad_mag; 

    // > look edge information
    T le_id;
    T le_pt_x;
    T le_pt_y;
    T le_orient;
    T le_grad_mag;

    T _dx; //< variation of edge in pixel
    T _dt; //< variation of angle in radians
    T _sx; //< sample pixel
    T _st; //< sample angle in radians

    // > number of curves in a bundle in terms of location and orientation variations
    unsigned curves_num_in_bundle_pixel;
    unsigned curves_num_in_bundle_theta;

    // > used for building bundles
    T _max_k;

private:
    
    T *bundle_min_ks;
    T *bundle_max_ks;
    T *cmp_bundle_min_ks;
    T *cmp_bundle_max_ks;
    T *intersect_bundle_min_ks;
    T *intersect_bundle_max_ks;

    bool *hyp_LookEdge;
    unsigned *edge_chain_final;         //< store all final curvelets by edge ids
    unsigned *edge_chain_on_the_fly;    //< store a curvelet candidate on the fly
    unsigned *edge_chain_target;        //< store curvelets w.r.t. target edge

public:
    T *_edgeLookList;
    
    // > constructor
    CurveletCPU(int &num_edges, int &sz_edge_data, T *edgeLookList, int max_LookEdgeNum,
                T dx, T dt, T sx, T st, T max_k, unsigned group_max_sz):
                _num_edges(num_edges), _sz_edge_data(sz_edge_data), _max_num_look_edges(max_LookEdgeNum),
                _dx(dx), _dt(dt), _sx(sx), _st(st), _max_k(max_k), _group_max_sz(group_max_sz)
    {
        //_edgeLookList = edgeLookList;
        _edgeLookList = new T[(_sz_edge_data+1) * (_max_num_look_edges+1) * _num_edges ];
        //> assign edgeLookList to edgeLookList_fit
        for (unsigned i = 0; i < _num_edges; i++) {
            for (unsigned j = 0; j < (_max_num_look_edges+1)*5; j++) {
                _edgeLookList(i, j) = edgeLookList(i, j);
            }
        }

        //> 1) calculate number of curves in a bundle in terms of pixel and theta
        curves_num_in_bundle_pixel = 2*floor( _dx/_sx + 0.5 ) + 1;
        curves_num_in_bundle_theta = 2*floor( _dt/_st + 0.5 ) + 1;

        bundle_min_ks = new T[ _max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        bundle_max_ks = new T[ _max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        cmp_bundle_min_ks = new T[ 2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        cmp_bundle_max_ks = new T[ 2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        intersect_bundle_min_ks = new T[ curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        intersect_bundle_max_ks = new T[ curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        hyp_LookEdge = new bool[ _max_num_look_edges ];
        edge_chain_on_the_fly = new unsigned[ _group_max_sz ];
        edge_chain_target     = new unsigned[ (_group_max_sz+1) * _max_num_look_edges ];
        edge_chain_final      = new unsigned[ (_num_edges*_max_num_look_edges) * (_group_max_sz+1) ];
    }

    // > destructor
    ~CurveletCPU();

    void preprocessing();
    void build_curvelets_greedy();
    bool compute_curve_bundle( unsigned te_idx, unsigned le_idx );
    bool bundle_valid_check( unsigned le_idx );
    void retrieve_edge_data_from_edgeLookList(unsigned target_idx, unsigned look_idx, T &id, T &pt_x, T &pt_y, T &orient, T &strength);

    void move_to_cmp_bundle( unsigned cmp_idx, unsigned le_idx, bool rep_by_intersection );
    void bundle_intersection();
    bool bundle_intersection_valid_check();
    bool check_curvelet_exist( unsigned edge_chain_on_the_fly_sz );
};

template<typename T>
void CurveletCPU<T>::preprocessing() {
    
    for (unsigned i = 0; i < _max_num_look_edges; i++) {
        for (unsigned j = 0; j < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; j++) {
            bundle_min_ks(i, j) = -_max_k;
            bundle_max_ks(i, j) = _max_k;
        }

        hyp_LookEdge[i] = 0;
    }

    //std::cout<<curves_num_in_bundle_pixel<<", "<<curves_num_in_bundle_theta<<std::endl;

    //> initialize edge_chain_on_the_fly
    for (unsigned i = 0; i < _group_max_sz; i++) {
        edge_chain_on_the_fly[i] = 0;
    }

    //> initialize edge_chain_target
    for (unsigned i = 0; i < _max_num_look_edges; i++) {
        for (unsigned j = 0; j < (_group_max_sz+1); j++) {
            edge_chain_target(i, j) = 0;
        }
    }

    //> initialize edge_chain_final
    for (unsigned i = 0; i < _num_edges*_max_num_look_edges; i++) {
        for (unsigned j = 0; j < (_group_max_sz+1); j++) {
            edge_chain_final(i, j) = 0;
        }
    }
}

template<typename T>
void CurveletCPU<T>::retrieve_edge_data_from_edgeLookList(unsigned target_idx, unsigned look_idx, T &id, T &pt_x, T &pt_y, T &orient, T &strength)
{
    id       = _edgeLookList(target_idx, look_idx*5);
    pt_x     = _edgeLookList(target_idx, look_idx*5 + 1);
    pt_y     = _edgeLookList(target_idx, look_idx*5 + 2);
    orient   = _edgeLookList(target_idx, look_idx*5 + 3);
    strength = _edgeLookList(target_idx, look_idx*5 + 4);
}

template<typename T>
void CurveletCPU<T>::build_curvelets_greedy( )
{
    // > Process to build curvelets in a greedy way:
    // > 1) loop over all target edges, and form pairs of curves from target edges and look edges
    // > 2) examine all pairs of curves by intersection to gradually build a curve with maximal grouping size of 7 edges

    //> some variables
    bool valid_bundle_created = false;
    bool valid_bundle_intersection = false;
    int valid_edge_num = 0;
    unsigned edge_chain_target_idx = 0;
    unsigned last_edge_chain_target_idx = 0;
    unsigned edge_chain_lidx = 0;

    unsigned DEBUG_TE_ID = 20837;

    //> loop over all target edges
    for (unsigned te_idx = 0; te_idx < _num_edges; te_idx++) {

        //> retrieve target edge data
        retrieve_edge_data_from_edgeLookList(te_idx, 0, te_id, te_pt_x, te_pt_y, te_orient, te_grad_mag);

        //> for every target edge, try two times,
        //  each time goes different direction
        for (unsigned f_run = 0; f_run < 2; f_run++) {

            //> refresh data
            edge_chain_target_idx = 0;

            //> loop over all associated look edges (maximal _max_num_look_edges look edges), create pair-wise curvelet bundles
            for (unsigned le_idx = 0; le_idx < _max_num_look_edges; le_idx++) {

                //> retrieve look edge data
                retrieve_edge_data_from_edgeLookList(te_idx, le_idx+1, le_id, le_pt_x, le_pt_y, le_orient, le_grad_mag);
                if (le_id < 0)
                    break;

                //if ((f_run == 0) && te_id == DEBUG_TE_ID) {
                    //std::cout<<"("<<te_id<<", "<<le_id<<")"<<std::endl;
                //}
                //> DEBUG!!!!!
                /*if (te_id == DEBUG_TE_ID) {
                    std::cout<<"("<<le_id<<", "<<le_idx<<")"<<std::endl;
                }*/

                //> compute the angle direction from target edge to look edges
                const T _dir = angle_from_pt_to_pt(te_pt_x, te_pt_y, le_pt_x, le_pt_y);

                //> first try: forward == true && leading == true
                if ((f_run == 0) && (dot(_dir, le_orient) > 0)) {
                    //> create curve bundles between the target edge and all the look edges
                    valid_bundle_created = compute_curve_bundle( te_idx, le_idx );
                    valid_edge_num++;
                    
                }
                else if ((f_run == 1) && (dot(_dir, le_orient) < 0)) {
                    valid_bundle_created = compute_curve_bundle( te_idx, le_idx );
                    valid_edge_num++;
                }
                else {
                    continue;
                }

                //> DEBUG!!!!!
                /*if ((f_run == 0) && te_idx == DEBUG_TE_ID && le_id == 9) {
                    int count = 0;
                    std::cout<<"le_idx = "<<le_idx<<std::endl;
                    std::cout<<"valid_bundle_created = "<<valid_bundle_created<<std::endl;
                    for (unsigned i = 0; i < curves_num_in_bundle_pixel*curves_num_in_bundle_theta; i++) {
                        std::cout<<bundle_min_ks(0, i)<<","<<bundle_max_ks(0, i)<<"\t";
                        count++;
                        if (count == 7){
                            std::cout<<std::endl;
                            count = 0;
                        }
                    }
                }*/
                
            }   //> for loop over forming pairs of curve bundle

            //> DEBUG!!!!!!!
            #if 0
            if (te_id == DEBUG_TE_ID) {
                std::cout<<"look edge id:   ";
                for ( unsigned i = 0; i < _max_num_look_edges; i++ ) {
                    if (hyp_LookEdge[i]) {
                        retrieve_edge_data_from_edgeLookList(te_idx, i+1, le_id, le_pt_x, le_pt_y, le_orient, le_grad_mag);
                        std::cout<<le_id<<"  ";
                    }
                }
                std::cout<<std::endl;
                std::cout<<"look edge indx: ";
                for ( unsigned i = 0; i < _max_num_look_edges; i++ ) {
                    if (hyp_LookEdge[i]) {
                        retrieve_edge_data_from_edgeLookList(te_idx, i+1, le_id, le_pt_x, le_pt_y, le_orient, le_grad_mag);
                        std::cout<<i<<"  ";
                    }
                }
                std::cout<<std::endl;

                if (f_run == 0) {
                    for (unsigned l = 0; l < _max_num_look_edges; l++) {
                        std::cout<<hyp_LookEdge[l]<<"  ";
                    }
                    std::cout<<std::endl;
                }

                //> Print curvature bundle
                /*if (f_run == 0) {
                    for (unsigned l = 0; l < _max_num_look_edges; l++) {
                        if (hyp_LookEdge[l]) {
                            int count = 0;
                            std::cout<<"le_idx = "<<l<<std::endl;
                            //std::cout<<"le_idx = "<<le_idx<<std::endl;
                            //std::cout<<"valid_bundle_created = "<<valid_bundle_created<<std::endl;
                            for (unsigned i = 0; i < curves_num_in_bundle_pixel*curves_num_in_bundle_theta; i++) {
                                std::cout<<bundle_min_ks(l, i)<<","<<bundle_max_ks(l, i)<<"\t";
                                count++;
                                if (count == 7){
                                    std::cout<<std::endl;
                                    count = 0;
                                }
                            }
                            std::cout<<std::endl;
                        }
                    }
                }*/
                std::cout<<"==============================================="<<std::endl;
            }
            #endif

            //> now, for each pair-wise curvelet bundle hypothesis formed w.r.t. the target edge,
            //  fix one and examine curve bundle intersections with all the rest,
            //  and loop over all the look edges to do the same procedure
            for (unsigned le_idx = 0; le_idx < _max_num_look_edges; le_idx++) {

                //> refresh data
                edge_chain_lidx = 0;

                //> if the edge is valid
                if (hyp_LookEdge[le_idx]) {

                    //> put the target edge id in the edge chain list
                    //edge_chain(edge_chain_tidx, edge_chain_lidx) = te_id;
                    //edge_chain_lidx++;

                    //> keep the curvelet bundle
                    move_to_cmp_bundle( 0, le_idx, false );

                    //> loop over the rest of the hypothesis
                    for (unsigned le_remain_idx = 0; le_remain_idx < _max_num_look_edges; le_remain_idx++) {

                        //> if the look edge is valid
                        if (hyp_LookEdge[le_remain_idx]) {

                            //> retrieve all the valid look edges
                            retrieve_edge_data_from_edgeLookList(te_idx, le_remain_idx+1, le_id, le_pt_x, le_pt_y, le_orient, le_grad_mag);

                            //> put the reference edge id to the edge chain list
                            if (le_idx == le_remain_idx) {

                                //> put the look edge id to the edge chain list
                                //edge_chain(edge_chain_tidx, edge_chain_lidx) = le_id;
                                //edge_chain_lidx++;
                                edge_chain_on_the_fly[ edge_chain_lidx ] = le_id;
                                edge_chain_lidx++;

                                //> continue looping if the edge chain is not filled entirely by edge ids
                                if (edge_chain_lidx == (_group_max_sz-1)) {
                                    break;
                                }
                                else {
                                    continue;
                                }
                            }

                            //> keep the curvelet bundle
                            move_to_cmp_bundle( 1, le_remain_idx, false );

                            //> find bundle intersections: the absolute min and max curvatures
                            bundle_intersection();

                            //> check whether the intersected bundle is valid or not
                            valid_bundle_intersection = bundle_intersection_valid_check();

                            //> if the intersection is valid
                            if (valid_bundle_intersection) {

                                //> put the look edge id to the edge chain list
                                //edge_chain(edge_chain_tidx, edge_chain_lidx) = le_id;
                                //edge_chain_lidx++;
                                edge_chain_on_the_fly[ edge_chain_lidx ] = le_id;
                                edge_chain_lidx++;

                                //> change the cmp_bundle set by putting in the intersected bundle
                                move_to_cmp_bundle( 0, le_remain_idx, true );
                            }

                            //> continue looping if the edge chain is not filled entirely by edge ids
                            if (edge_chain_lidx == (_group_max_sz-1)) {
                                break;
                            }
                            else {
                                continue;
                            }

                        }
                        else {
                            continue;
                        }
                    }
                }
                else {
                    continue;
                }

                //> DEBUG!!!!! Print edge_chain_on_the_fly arrays
                #if 1
                if (te_id == DEBUG_TE_ID) {
                    std::cout<<le_idx<<": ";
                    for (unsigned di = 0; di < _group_max_sz; di++) {
                        std::cout<<edge_chain_on_the_fly[di]<<"  ";
                    }
                    std::cout<<std::endl;
                }
                #endif
                
                //> Should this edge chain be kept or not?
                //  This is examined by: i) edge chain size > 2, and
                //                      ii) curvelet does not exist before
                // EDGE_CHAIN_TARGET: (# of rows: _max_num_look_edges; # of cols: _group_max_sz + 1)
                // | le_id11 | le_id12 | ... | # of le_ids in this chain |
                // | le_id21 | le_id22 | ... | # of le_ids in this chain |
                // | ...     | ...     | ... | ...                       |
                if (edge_chain_lidx > 2 && !check_curvelet_exist(edge_chain_lidx)) {
                    
                    //> add the edge_chain_on_the_fly to edge_chain_target
                    for (unsigned chain_idx = 0; chain_idx < edge_chain_lidx; chain_idx++) {
                        edge_chain_target(edge_chain_target_idx, chain_idx) = edge_chain_on_the_fly[chain_idx];
                    }

                    //> keep the vliad edge chain size at the end
                    edge_chain_target(edge_chain_target_idx, _group_max_sz) = edge_chain_lidx;
                    edge_chain_target_idx++;
                }

                //> refresh the edge_chain_on_the_fly
                for (unsigned i = 0; i < _group_max_sz; i++) {
                    edge_chain_on_the_fly[i] = 0;
                }
                
            } //> for loop over all l_id

            //> move edge_chain_target to edge_chain_final
            // EDGE_CHAIN_FINAL: (# of rows: num_of_te * _max_num_look_edges; # of cols: _group_max_sz + 1)
            // | te_id1 | le_id11 | le_id12 | ...  |
            // | te_id2 | le_id21 | le_id22 | ...  |
            // | ...    | ...     | ...     | ...  |
            unsigned target_chain_idx = 0;
            unsigned edge_chain_size = 0;
            unsigned echain_final_ridx = 2;
            for (unsigned i = last_edge_chain_target_idx; i < (last_edge_chain_target_idx + edge_chain_target_idx); i++) {
                edge_chain_final(i, 0) = te_id;

                edge_chain_size = edge_chain_target(target_chain_idx, _group_max_sz);

                //> 1st run, go one direction
                if (f_run == 0) {

                    edge_chain_final(i, 1) = te_id;
                    echain_final_ridx = 2;

                    for (unsigned j = 0; j < _group_max_sz; j++) {
                        edge_chain_final(i, echain_final_ridx) = edge_chain_target(target_chain_idx, j);
                        echain_final_ridx++;
                    }
                }
                else if (f_run == 1) {
                    //> DEBUG!!!!!!!!!!!!!!!!
                    #if 0
                    if (te_id == DEBUG_TE_ID) {
                        std::cout<<"edge_chain_size = "<<edge_chain_size<<std::endl;
                    }
                    #endif

                    echain_final_ridx = 1;

                    for (int j = (edge_chain_size-1); j >= 0; j--) {
                        #if 0
                        if (te_id == DEBUG_TE_ID) {
                            std::cout<<edge_chain_target(target_chain_idx, j)<<"  "<<echain_final_ridx<<"  "<<j<<std::endl;
                        }
                        #endif
                        edge_chain_final(i, echain_final_ridx) = edge_chain_target(target_chain_idx, j);
                        echain_final_ridx++;
                    }

                    //std::cout<<"=== MAKE IT HERE? ======"<<std::endl;
                    edge_chain_final(i, echain_final_ridx) = te_id;
                }
                target_chain_idx++;
            }
            
            //> update last_edge_chain_target_idx
            last_edge_chain_target_idx += edge_chain_target_idx;

            //> DEBUG!!!!!!!!!!!!!!!!
            #if 0
            if (te_id == DEBUG_TE_ID) {
                std::cout<<"last_edge_chain_target_idx = "<<last_edge_chain_target_idx<<std::endl;
            }
            #endif

            //> refresh data for the use of new target edges
            for (unsigned i = 0; i < _max_num_look_edges; i++) {
                for (unsigned j = 0; j < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; j++) {
                    bundle_min_ks(i, j) = -_max_k;
                    bundle_max_ks(i, j) = _max_k;
                }
                hyp_LookEdge[i] = false;
            }
            for (unsigned i = 0; i < _max_num_look_edges; i++) {
                for (unsigned j = 0; j < (_group_max_sz+1); j++) {
                    edge_chain_target(i, j) = 0;
                }
            }

        }   //> for loop over 2 forward runs

        //> DEBUG!!!!!!!!!!!!!
        //  PRINT OUT THE FINAL EDGE CHAIN AFTER TWO RUNS ARE FINISHED
        #if 0
        if (te_id == DEBUG_TE_ID) {
            for (unsigned i = 0; i < last_edge_chain_target_idx; i++) {
                for (unsigned j = 0; j < (_group_max_sz+1); j++) {
                    std::cout<< edge_chain_final(i, j) << "  ";
                }
                std::cout<<std::endl;
            }
        }
        #endif
    } //> for loop over all t_id
}

//> compute the max and min curvature of the curve bundle from a pair of edges
template<typename T>
bool CurveletCPU<T>::compute_curve_bundle( unsigned te_idx, unsigned le_idx ) {
    
    //> global indices to fill values in bundle_min_ks and bundle_max_ks fill
    unsigned bk_idx = 0;

    //> \Delta X, \Delta Y, and their squares
    const T DX = te_pt_x - le_pt_x;
    const T DY = te_pt_y - le_pt_y;
    const T DX2 = DX*DX;
    const T DY2 = DY*DY;

    //> some pre-calculated constants
    const T T0 = le_orient;
    const T T2 = te_orient;
    const T _dx2 = _dx * _dx;
    const T T0_p_dt = T0 + _dt;
    const T T0_m_dt = T0 - _dt;

    //> final max and min curvature
    //T k_min = -_max_k;
    //T k_max = _max_k;

    //> non-constant variables computed on the fly
    // TODO: change variable name to increase readability
    T k_dx_min, k_dx_max;
    T k_dt_min, k_dt_max;
    T delta_x1;
    T delta_x2;
    T delta_t1;
    T delta_t2;
    T T2t, tt;
    T eq_denom;
    T ttmdx2;
    T sin_t2_p_dt2_m_t0_m_dt;
    T sin_t2_p_dt2_m_t0_p_dt;
    T sin_T0pdt_times_DY;
    T cos_T0pdt_times_DX;
    T sin_T0mdt_times_DY;
    T cos_T0mdt_times_DX;

    //> if neighboring point is close to the target edgel, no static constraint can be computed of leave it as the default bundle
    //  and the max and min curvatures of the curve bundle remains as default
    if ( DX2 + DY2 < 4*_dx*_dx ) {
        //if (te_id == 0 && le_id == 15)
        //    std::cout<<"failure here ..."<<std::endl;
        hyp_LookEdge[ le_idx ] = true;
        return true;
    }

    //> Proposition 4 in the paper: Circular Arc Bundle Transport
    if ( (DX * cos(T2) + DY * sin(T2)) < 0 ) {

        //> loop over all the possible curves to find the max and min curvature
        for (unsigned i = 0; i < curves_num_in_bundle_pixel; i++) {
            for (unsigned j = 0; j < curves_num_in_bundle_theta; j++) {

                //> calculate the relative difference in terms of pixel and radian
                delta_x2 = _sx * ((T)i - (T(curves_num_in_bundle_pixel) - 1.0) / 2.0);
                delta_t2 = _st * ((T)j - (T(curves_num_in_bundle_theta) - 1.0) / 2.0);

                //> ---------------- i) With Respect To delta x --------------
                //> absolute grid orientation of the curve
                T2t = T2 + delta_t2;

                //> intermediate calculations
                tt = sin(T2t) * DX - cos(T2t) * DY;
                ttmdx2 = tt - delta_x2;

                //const T ttmdx2 = tt - dx2;

                //> k_max / k_min equation denominator w.r.t. pixel
                eq_denom = (-2.0 * delta_x2 * tt + DX2 + DY2 - _dx2 + delta_x2 * delta_x2);

                //> max and min curvatures in terms of pixel
                k_dx_min = 2.0 * (ttmdx2 - _dx) / eq_denom;
                k_dx_max = 2.0 * (ttmdx2 + _dx) / eq_denom;

                //> compute the intersection of all surfaces
                //> 1) min
                if ( bundle_min_ks(le_idx, bk_idx) < k_dx_min )
                    bundle_min_ks(le_idx, bk_idx) = k_dx_min;

                //> 2) max
                if ( bundle_max_ks(le_idx, bk_idx) > k_dx_max )
                    bundle_max_ks(le_idx, bk_idx) = k_dx_max;
                
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

                //if (te_id == 0 && le_id == 7) {
                //    std::cout<<bk_idx+1<<" (k_dt_min, k_dt_max) = ("<<k_dt_min<<", "<<k_dt_max<<")"<<std::endl;      
                //}

                //> compute the intersection of all surfaces
                //> 1) min
                if ( bundle_min_ks(le_idx, bk_idx) < k_dt_min )
                    bundle_min_ks(le_idx, bk_idx) = k_dt_min;
                
                //> 2) max
                if ( bundle_max_ks(le_idx, bk_idx) > k_dt_max )
                    bundle_max_ks(le_idx, bk_idx) = k_dt_max;

                bk_idx++;
            }
        }
    }
    else {
        //> loop over all the possible curves to find the max and min curvature
        for (unsigned i = 0; i < curves_num_in_bundle_pixel; i++) {
            for (unsigned j = 0; j < curves_num_in_bundle_theta; j++) {

                //> calculate the relative difference in terms of pixel and radian
                delta_x2 = _sx * ((T)i - (T(curves_num_in_bundle_pixel) - 1.0) / 2.0);
                delta_t2 = _st * ((T)j - (T(curves_num_in_bundle_theta) - 1.0) / 2.0);

                //> -------------- i) With Respect To delta x ----------------
                //> absolute grid orientation of the curve
                T2t = T2 + delta_t2;

                //> intermediate calculations
                tt = sin(T2t) * DX - cos(T2t) * DY;
                ttmdx2 = tt - delta_x2;

                //> k_max / k_min equation denominator w.r.t. pixel
                eq_denom = (-2.0 * delta_x2 * tt + DX2 + DY2 - _dx2 + delta_x2 * delta_x2);

                //> max and min curvatures in terms of pixel
                k_dx_min = 2.0 * (ttmdx2 - _dx) / eq_denom;
                k_dx_max = 2.0 * (ttmdx2 + _dx) / eq_denom;

                //> compute the intersection of all surfaces
                //> 1) min
                if ( bundle_min_ks(le_idx, bk_idx) < k_dx_min )
                    bundle_min_ks(le_idx, bk_idx) = k_dx_min;

                //> 2) max
                if ( bundle_max_ks(le_idx, bk_idx) > k_dx_max )
                    bundle_max_ks(le_idx, bk_idx) = k_dx_max;

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
                if ( bundle_min_ks(le_idx, bk_idx) < k_dt_min )
                    bundle_min_ks(le_idx, bk_idx) = k_dt_min;
                
                //> 2) max
                if ( bundle_max_ks(le_idx, bk_idx) > k_dt_max )
                    bundle_max_ks(le_idx, bk_idx) = k_dt_max;

                bk_idx++;
            }
        }
    }

    //> check whether the built bundle is valid or not
    bool valid = bundle_valid_check( le_idx );

    //> this boolean list indicates which look edge is valid to be proceeded
    hyp_LookEdge[ le_idx ] = (valid) ? true : false;

    return (valid) ? true : false;
}

template<typename T>
bool CurveletCPU<T>::bundle_valid_check( unsigned le_idx )
{    
    //> check the curve bundle is valid or not
    for (int k_idx = 0; k_idx < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; k_idx++ ) {
        if ( bundle_max_ks(le_idx, k_idx) > bundle_min_ks(le_idx, k_idx) ) 
            return true;
    }
    return false;
}

template<typename T>
void CurveletCPU<T>::bundle_intersection()
{
    //> compare the curvatures in the cmp_bundle_min_ks and cmp_bundle_max_ks
    for (unsigned j = 0; j < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; j++) {
        //> compare the min curvature, pick the larger one
        intersect_bundle_min_ks[ j ] = (cmp_bundle_min_ks(0, j) < cmp_bundle_min_ks(1, j)) ? (cmp_bundle_min_ks(1, j)) : (cmp_bundle_min_ks(0, j));

        //> compare the max curvature, pick the smaller one
        intersect_bundle_max_ks[ j ] = (cmp_bundle_max_ks(0, j) > cmp_bundle_max_ks(1, j)) ? (cmp_bundle_max_ks(1, j)) : (cmp_bundle_max_ks(0, j));
    }
}

template<typename T>
bool CurveletCPU<T>::bundle_intersection_valid_check()
{
    //> check the intersected curve bundle is valid or not
    for (int k_idx = 0; k_idx < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; k_idx++ ) {
        if ( intersect_bundle_max_ks[k_idx] > intersect_bundle_min_ks[k_idx] ) 
            return true;
    }
    return false;
}

template<typename T>
void CurveletCPU<T>::move_to_cmp_bundle( unsigned cmp_idx, unsigned le_idx, bool rep_by_intersection )
{
    //> keep the curvelet bundle
    if (!rep_by_intersection) {
        for (unsigned bidx = 0; bidx < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta); bidx++) {
            cmp_bundle_min_ks(cmp_idx, bidx) = bundle_min_ks(le_idx, bidx);
            cmp_bundle_max_ks(cmp_idx, bidx) = bundle_max_ks(le_idx, bidx);
        }
    }
    else {
        for (unsigned bidx = 0; bidx < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta); bidx++) {
            cmp_bundle_min_ks(cmp_idx, bidx) = intersect_bundle_min_ks[ bidx ];
            cmp_bundle_max_ks(cmp_idx, bidx) = intersect_bundle_max_ks[ bidx ];
        }
    }
}

template<typename T>
bool CurveletCPU<T>::check_curvelet_exist( unsigned edge_chain_on_the_fly_sz ) 
{
    //> does the curvelet exist before?
    //  this is examined by comparing the edge_chain_on_the_fly 
    //  with edge_chain_target to see whether overlapping edge 
    //  chain exists

    unsigned edge_chain_sz;
    bool cvlet_exists = true;
    for (unsigned i = 0; i < _max_num_look_edges; i++) {
        //> retrieve the size
        edge_chain_sz = edge_chain_target(i, _group_max_sz);

        //> reset
        cvlet_exists = true;

        //> if reaching to the end of the edge_chain_target
        if (edge_chain_sz == 0) {
            return false;
        }

        //> if the size of the edge chain is the same
        if (edge_chain_sz == edge_chain_on_the_fly_sz) {

            for (unsigned chain_idx = 0; chain_idx < edge_chain_on_the_fly_sz; chain_idx++) {
                cvlet_exists = cvlet_exists && (edge_chain_on_the_fly[chain_idx] == edge_chain_target[chain_idx]);
            }

            //> cvlet_exists remains true only if all the edgel ids in edge chains match
            if (cvlet_exists) {
                return true;
            }
        }
        else {
            continue;
        }
    }

    return false;
}

template<typename T>
CurveletCPU<T>::~CurveletCPU() {
    delete[] _edgeLookList;
    delete[] bundle_min_ks;
    delete[] bundle_max_ks;
    delete[] cmp_bundle_min_ks;
    delete[] cmp_bundle_max_ks;
    delete[] intersect_bundle_min_ks;
    delete[] intersect_bundle_max_ks;
    delete[] hyp_LookEdge;

    delete[] edge_chain_on_the_fly;
    delete[] edge_chain_target;
    delete[] edge_chain_final;
}

#endif // CPU_CURVELET_HPP