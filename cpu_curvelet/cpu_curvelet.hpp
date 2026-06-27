#ifndef CPU_CURVELET_HPP
#define CPU_CURVELET_HPP

#include <map>
#include <vector>
#include <deque>
#include <list>
#include <utility> 
#include <iostream>
#include <cmath>

#include <omp.h>

#include "indices.hpp"
#include "curvelet_utils.hpp"
#include "preprocess.hpp"

template<typename T>
class CurveletCPU
{
protected:
    int _sz_edge_data;
    int _num_edges;
    unsigned _group_max_sz;
    int _max_num_look_edges;
    const T *_csr_edges;
    const int *_neighbor_offsets;
    const int *_neighbor_ids;

    T _dx; //< variation of edge in pixel
    T _dt; //< variation of angle in radians
    T _sx; //< sample pixel
    T _st; //< sample angle in radians

    //> number of curves in a bundle in terms of location and orientation variations
    unsigned curves_num_in_bundle_pixel;
    unsigned curves_num_in_bundle_theta;

    //> used for building bundles
    T _max_k;
    T _nrad;
    T *_to_edges;

    static const unsigned CURVELET_INFO_WIDTH = 10;

public:
    int omp_threads;
    unsigned _num_curvelets;
    unsigned _max_curvelets;
    unsigned _max_per_anchor;
    unsigned *_anchor_chain_count;
    unsigned *_edge_chain_final;
    double *_curvelet_info;

    void init_curvelet_storage(int max_LookEdgeNum)
    {
        _max_num_look_edges = max_LookEdgeNum;
        _num_curvelets = 0;
        _max_per_anchor = (unsigned)(max_LookEdgeNum + 1) * 2;
        _max_curvelets = (unsigned)_num_edges * (unsigned)(max_LookEdgeNum + 1) * 2;

        _edge_chain_final = new unsigned[ _max_curvelets * (_group_max_sz + 1) ];
        _curvelet_info = new double[ _max_curvelets * CURVELET_INFO_WIDTH ]();
        _anchor_chain_count = new unsigned[ _num_edges ];
        for (int i = 0; i < _num_edges; i++) {
            _anchor_chain_count[i] = 0;
        }

        curves_num_in_bundle_pixel = 2 * floor( _dx / _sx + 0.5 ) + 1;
        curves_num_in_bundle_theta = 2 * floor( _dt / _st + 0.5 ) + 1;
    }

    unsigned neighbor_count(unsigned te_idx) const
    {
        return static_cast<unsigned>( _neighbor_offsets[te_idx + 1] - _neighbor_offsets[te_idx] );
    }

    CurveletCPU(int &num_edges, int &sz_edge_data, const CPUNeighborGraph<T> &csr,
                T dx, T dt, T sx, T st, T max_k, unsigned group_max_sz, int nthreads,
                T *to_edges, T nrad):
                _num_edges(num_edges), _sz_edge_data(sz_edge_data),
                _dx(dx), _dt(dt), _sx(sx), _st(st), _max_k(max_k), _nrad(nrad), 
                _to_edges(to_edges), _group_max_sz(group_max_sz), omp_threads(nthreads),
                _csr_edges(csr.edges.data()),
                _neighbor_offsets(csr.neighbor_offsets.data()),
                _neighbor_ids(csr.neighbor_ids.data())
    {
        init_curvelet_storage(static_cast<int>(csr.max_num_of_neighbors));
    }

    unsigned num_curvelets() const { return _num_curvelets; }
    unsigned chain_width() const { return _group_max_sz + 1; }

    //> destructor
    ~CurveletCPU();

    void curvelet_preprocessing(T* bundle_min_ks, T* bundle_max_ks, bool* hyp_LookEdge, unsigned* edge_chain_on_the_fly, unsigned* edge_chain_target);
    void build_curvelets_greedy();
    void compact_curvelet_output();
    bool compute_curve_bundle( unsigned te_idx, unsigned le_idx, T* bundle_min_ks, T* bundle_max_ks, bool* hyp_LookEdg, 
                               T te_pt_x, T te_pt_y, T le_pt_x, T le_pt_y, T te_orient, T le_orient );
    bool bundle_valid_check( unsigned le_idx, T* bundle_min_ks, T* bundle_max_ks );
    void retrieve_edge_data(unsigned target_idx, unsigned look_idx, T &id, T &pt_x, T &pt_y, T &orient, T &strength);

    void move_to_cmp_bundle( unsigned cmp_idx, unsigned le_idx, bool rep_by_intersection, 
                             T* bundle_min_ks, T* bundle_max_ks, T* cmp_bundle_min_ks, T* cmp_bundle_max_ks,
                             T* intersect_bundle_min_ks, T* intersect_bundle_max_ks );
    void bundle_intersection( T* cmp_bundle_min_ks, T* cmp_bundle_max_ks, T* intersect_bundle_min_ks, T* intersect_bundle_max_ks );
    bool bundle_intersection_valid_check( T* intersect_bundle_min_ks, T* intersect_bundle_max_ks );
    bool check_curvelet_exist( unsigned edge_chain_on_the_fly_sz, unsigned* edge_chain_on_the_fly, unsigned* edge_chain_target );
    void chain_push_front( unsigned* chain, unsigned &lidx, unsigned id );
    void record_curvelet_chain( unsigned te_id, unsigned* chain, unsigned chain_sz, bool forward,
                                T ref_pt_x, T ref_pt_y, T ref_theta,
                                T* cmp_bundle_min_ks, T* cmp_bundle_max_ks );
    void fill_curvelet_info( unsigned row, bool forward, T ref_pt_x, T ref_pt_y, T ref_theta,
                             T* cmp_bundle_min_ks, T* cmp_bundle_max_ks,
                             unsigned* chain, unsigned chain_sz );
};

template<typename T>
void CurveletCPU<T>::curvelet_preprocessing( T* bundle_min_ks, T* bundle_max_ks, bool* hyp_LookEdge, 
                                    unsigned* edge_chain_on_the_fly, unsigned* edge_chain_target ) 
{
    
    for (unsigned i = 0; i < _max_num_look_edges; i++) {
        for (unsigned j = 0; j < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; j++) {
            bundle_min_ks(i, j) = -_max_k;
            bundle_max_ks(i, j) = _max_k;
        }

        hyp_LookEdge[i] = 0;
    }

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
    /*for (unsigned i = 0; i < _num_edges*_max_num_look_edges; i++) {
        for (unsigned j = 0; j < (_group_max_sz+1); j++) {
            edge_chain_final(i, j) = 0;
        }
    }*/
}

template<typename T>
void CurveletCPU<T>::retrieve_edge_data(unsigned target_idx, unsigned look_idx, T &id, T &pt_x, T &pt_y, T &orient, T &strength)
{
    int edge_id = 0;
    if (look_idx == 0) {
        edge_id = static_cast<int>(target_idx);
    } else {
        const unsigned nbr_slot = look_idx - 1;
        const int nb_begin = _neighbor_offsets[target_idx];
        const int nb_end = _neighbor_offsets[target_idx + 1];
        if (static_cast<int>(nbr_slot) >= nb_end - nb_begin) {
            id = T(-1);
            return;
        }
        edge_id = _neighbor_ids[nb_begin + static_cast<int>(nbr_slot)];
    }

    const int ebase = edge_id * _sz_edge_data;
    id       = T(edge_id);
    pt_x     = _csr_edges[ebase + 0];
    pt_y     = _csr_edges[ebase + 1];
    orient   = _csr_edges[ebase + 2];
    strength = _csr_edges[ebase + 3];
}

template<typename T>
void CurveletCPU<T>::chain_push_front( unsigned* chain, unsigned &lidx, unsigned id )
{
    if (lidx >= _group_max_sz) {
        return;
    }
    for (int k = (int)lidx; k > 0; k--) {
        chain[k] = chain[k - 1];
    }
    chain[0] = id;
    lidx++;
}

template<typename T>
void CurveletCPU<T>::fill_curvelet_info( unsigned row, bool forward, T ref_pt_x, T ref_pt_y, T ref_theta,
                                         T* cmp_bundle_min_ks, T* cmp_bundle_max_ks,
                                         unsigned* chain, unsigned chain_sz )
{
    T mind = T(100);
    unsigned mini = 0;
    unsigned minj = 0;
    for (unsigned ii = 0; ii < curves_num_in_bundle_pixel; ii++) {
        for (unsigned jj = 0; jj < curves_num_in_bundle_theta; jj++) {
            const unsigned bidx = ii * curves_num_in_bundle_theta + jj;
            if (cmp_bundle_max_ks(0, bidx) > cmp_bundle_min_ks(0, bidx)) {
                const T ddx = _sx * (T(ii) - (T(curves_num_in_bundle_pixel) - T(1)) / T(2));
                const T ddt = _st * (T(jj) - (T(curves_num_in_bundle_theta) - T(1)) / T(2));
                const T d = ddx * ddx + ddt * ddt;
                if (d < mind) {
                    mind = d;
                    mini = ii;
                    minj = jj;
                }
            }
        }
    }

    const unsigned best_bidx = mini * curves_num_in_bundle_theta + minj;
    const T k_max = cmp_bundle_max_ks(0, best_bidx);
    const T k_min = cmp_bundle_min_ks(0, best_bidx);
    const T k = (k_max + k_min) / T(2);
    const T dx = _sx * (T(mini) - (T(curves_num_in_bundle_pixel) - T(1)) / T(2));
    const T dt = _st * (T(minj) - (T(curves_num_in_bundle_theta) - T(1)) / T(2));
    const T theta = To2Pi(ref_theta + dt);
    const T pt_x = ref_pt_x - dx * std::sin(theta);
    const T pt_y = ref_pt_y + dx * std::cos(theta);

    T length = T(0);
    for (unsigned i = 0; i + 1 < chain_sz; i++) {
        const unsigned e1 = chain[i];
        const unsigned e2 = chain[i + 1];
        const T x1 = _to_edges[e1 * _sz_edge_data + 0];
        const T y1 = _to_edges[e1 * _sz_edge_data + 1];
        const T x2 = _to_edges[e2 * _sz_edge_data + 0];
        const T y2 = _to_edges[e2 * _sz_edge_data + 1];
        length += std::sqrt(sq_dist(x1, y1, x2, y2));
    }

    const T alpha3 = T(1);
    const T alpha4 = T(1);
    const T quality = (length > T(0) && chain_sz > 0) ? T(2) / (alpha3 * _nrad / length + alpha4 * length / T(chain_sz)) : T(0);

    _curvelet_info[0 * _max_curvelets + row] = forward ? 1.0 : 0.0;
    _curvelet_info[1 * _max_curvelets + row] = (double)k_max;
    _curvelet_info[2 * _max_curvelets + row] = (double)k_min;
    _curvelet_info[3 * _max_curvelets + row] = (double)ref_theta;
    _curvelet_info[4 * _max_curvelets + row] = (double)pt_x;
    _curvelet_info[5 * _max_curvelets + row] = (double)pt_y;
    _curvelet_info[6 * _max_curvelets + row] = (double)theta;
    _curvelet_info[7 * _max_curvelets + row] = (double)k;
    _curvelet_info[8 * _max_curvelets + row] = (double)length;
    _curvelet_info[9 * _max_curvelets + row] = (double)quality;
}

template<typename T>
void CurveletCPU<T>::record_curvelet_chain( unsigned te_id, unsigned* chain, unsigned chain_sz, bool forward,
                                            T ref_pt_x, T ref_pt_y, T ref_theta,
                                            T* cmp_bundle_min_ks, T* cmp_bundle_max_ks )
{
    // Per-anchor slots preserve greedy seed order within each anchor; compact by te_id after the parallel loop.
    unsigned &n = _anchor_chain_count[te_id];
    if (n >= _max_per_anchor) {
        return;
    }

    const unsigned row = te_id * _max_per_anchor + n++;
    for (unsigned j = 0; j < (_group_max_sz + 1); j++) {
        edge_chain_final(row, j) = 0;
    }
    edge_chain_final(row, 0) = te_id + 1;
    for (unsigned i = 0; i < chain_sz && (i + 1) < (_group_max_sz + 1); i++) {
        edge_chain_final(row, i + 1) = chain[i] + 1;
    }

    fill_curvelet_info(row, forward, ref_pt_x, ref_pt_y, ref_theta,
                       cmp_bundle_min_ks, cmp_bundle_max_ks, chain, chain_sz);
}

template<typename T>
void CurveletCPU<T>::compact_curvelet_output()
{
    const unsigned row_width = _group_max_sz + 1;
    unsigned dst = 0;
    for (unsigned te = 0; te < (unsigned)_num_edges; te++) {
        const unsigned n = _anchor_chain_count[te];
        for (unsigned k = 0; k < n; k++) {
            const unsigned src = te * _max_per_anchor + k;
            if (dst != src) {
                for (unsigned j = 0; j < row_width; j++) {
                    edge_chain_final(dst, j) = edge_chain_final(src, j);
                }
                for (unsigned col = 0; col < CURVELET_INFO_WIDTH; col++) {
                    _curvelet_info[col * _max_curvelets + dst] = _curvelet_info[col * _max_curvelets + src];
                }
            }
            dst++;
        }
    }
    _num_curvelets = dst;
}

template<typename T>
void CurveletCPU<T>::build_curvelets_greedy( )
{
    _num_curvelets = 0;
    for (int i = 0; i < _num_edges; i++) {
        _anchor_chain_count[i] = 0;
    }

    omp_set_num_threads(omp_threads);
    double start = omp_get_wtime();
    double time_direction_filter = 0.0;
    double time_bundle_transport = 0.0;
    double time_chain_growth = 0.0;
    double time_dedup = 0.0;
    #pragma omp parallel
    {
        double local_direction_filter = 0.0;
        double local_bundle_transport = 0.0;
        double local_chain_growth = 0.0;
        double local_dedup = 0.0;

        //> declare local arrays
        T *bundle_min_ks;
        T *bundle_max_ks;
        T *cmp_bundle_min_ks;
        T *cmp_bundle_max_ks;
        T *intersect_bundle_min_ks;
        T *intersect_bundle_max_ks;

        bool *hyp_LookEdge;
        unsigned *edge_chain_on_the_fly;    //< store a curvelet candidate on the fly
        unsigned *edge_chain_target;        //< store curvelets w.r.t. target edge

        //> dynamically allocated arrays
        bundle_min_ks = new T[ _max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        bundle_max_ks = new T[ _max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        cmp_bundle_min_ks = new T[ 2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        cmp_bundle_max_ks = new T[ 2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        intersect_bundle_min_ks = new T[ curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        intersect_bundle_max_ks = new T[ curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        hyp_LookEdge          = new bool[ _max_num_look_edges ];
        edge_chain_on_the_fly = new unsigned[ _group_max_sz ];
        edge_chain_target     = new unsigned[ (_group_max_sz+1) * _max_num_look_edges ];

        //> some variables
        bool valid_bundle_created = false;
        bool valid_bundle_intersection = false;
        int valid_edge_num = 0;
        unsigned edge_chain_target_idx = 0;
        unsigned edge_chain_lidx = 0;

        //> target edge information
        T te_id = 0;
        T te_pt_x = 0;
        T te_pt_y = 0;
        T te_orient = 0;
        T te_grad_mag = 0; 

        //> look edge information
        T le_id = 0;
        T le_pt_x = 0;
        T le_pt_y = 0;
        T le_orient = 0;
        T le_grad_mag = 0;

        //> initialize all local arrays of omp threads
        curvelet_preprocessing( bundle_min_ks, bundle_max_ks, hyp_LookEdge, edge_chain_on_the_fly, edge_chain_target);

        #pragma omp for schedule(dynamic)
        for (unsigned te_idx = 0; te_idx < (unsigned)_num_edges; te_idx++) {

            const unsigned n_look = neighbor_count(te_idx);

            //> retrieve target edge data
            retrieve_edge_data(te_idx, 0, te_id, te_pt_x, te_pt_y, te_orient, te_grad_mag);

            //> for every target edge, try two times, each goes in a different direction
            for (unsigned f_run = 0; f_run < 2; f_run++) {

                edge_chain_target_idx = 0;
                // cvlet_style=2 (anchor-leading bidirectional): ref_first is always true
                const bool ref_first = true;

                for (unsigned le_idx = 0; le_idx < n_look; le_idx++) {

                    double filter_t0 = omp_get_wtime();
                    //> retrieve look edge data
                    retrieve_edge_data(te_idx, le_idx+1, le_id, le_pt_x, le_pt_y, le_orient, le_grad_mag);

                    //> compute the angle direction from target edge to look edges
                    const T _dir = angle_from_pt_to_pt(te_pt_x, te_pt_y, le_pt_x, le_pt_y);
                    const bool forward_pass = (f_run == 0) && (dot(_dir, te_orient) > 0);
                    const bool backward_pass = (f_run == 1) && (dot(_dir, te_orient) < 0);
                    local_direction_filter += omp_get_wtime() - filter_t0;

                    //> first try: forward == true && leading == true
                    if (forward_pass || backward_pass) {
                        double bundle_t0 = omp_get_wtime();
                        //> create curve bundles between the target edge and all the look edges
                        valid_bundle_created = compute_curve_bundle( te_idx, le_idx, bundle_min_ks, bundle_max_ks, hyp_LookEdge, te_pt_x, te_pt_y, le_pt_x, le_pt_y, te_orient, le_orient );
                        valid_edge_num++;
                        local_bundle_transport += omp_get_wtime() - bundle_t0;
                        //T* bundle_min_ks, T* bundle_max_ks, bool* hyp_LookEdge
                    }                    
                }   //> for loop over forming pairs of curve bundle

                //> now, for each pair-wise curvelet bundle hypothesis formed w.r.t. the target edge,
                //  fix one and examine curve bundle intersections with all the rest,
                //  and loop over all the look edges to do the same procedure
                double growth_t0 = omp_get_wtime();
                for (unsigned le_idx = 0; le_idx < n_look; le_idx++) {

                    edge_chain_lidx = 0;

                    if (hyp_LookEdge[le_idx]) {

                        edge_chain_on_the_fly[ edge_chain_lidx++ ] = (unsigned)te_id;

                        move_to_cmp_bundle( 0, le_idx, false, bundle_min_ks, bundle_max_ks, cmp_bundle_min_ks, cmp_bundle_max_ks,
                                            intersect_bundle_min_ks, intersect_bundle_max_ks );

                        for (unsigned le_remain_idx = 0; le_remain_idx < n_look; le_remain_idx++) {

                            if (hyp_LookEdge[le_remain_idx]) {

                                retrieve_edge_data(te_idx, le_remain_idx+1, le_id, le_pt_x, le_pt_y, le_orient, le_grad_mag);

                                if (le_idx == le_remain_idx) {

                                    if (ref_first) {
                                        edge_chain_on_the_fly[ edge_chain_lidx++ ] = (unsigned)le_id;
                                    } 
                                    else {
                                        chain_push_front(edge_chain_on_the_fly, edge_chain_lidx, (unsigned)le_id);
                                    }

                                    if (edge_chain_lidx >= _group_max_sz) {
                                        break;
                                    }
                                    else {
                                        continue;
                                    }
                                }

                                move_to_cmp_bundle( 1, le_remain_idx, false, bundle_min_ks, bundle_max_ks, cmp_bundle_min_ks, cmp_bundle_max_ks,
                                                    intersect_bundle_min_ks, intersect_bundle_max_ks );

                                bundle_intersection(cmp_bundle_min_ks, cmp_bundle_max_ks, intersect_bundle_min_ks, intersect_bundle_max_ks);

                                valid_bundle_intersection = bundle_intersection_valid_check(intersect_bundle_min_ks, intersect_bundle_max_ks);

                                if (valid_bundle_intersection) {

                                    if (ref_first) {
                                        edge_chain_on_the_fly[ edge_chain_lidx++ ] = (unsigned)le_id;
                                    } else {
                                        chain_push_front(edge_chain_on_the_fly, edge_chain_lidx, (unsigned)le_id);
                                    }

                                    move_to_cmp_bundle( 0, le_remain_idx, true, bundle_min_ks, bundle_max_ks, cmp_bundle_min_ks, cmp_bundle_max_ks,
                                                        intersect_bundle_min_ks, intersect_bundle_max_ks );
                                }

                                if (edge_chain_lidx >= _group_max_sz) {
                                    break;
                                }
                            }
                        }
                    }
                    else {
                        local_chain_growth += omp_get_wtime() - growth_t0;
                        growth_t0 = omp_get_wtime();
                        continue;
                    }

                    local_chain_growth += omp_get_wtime() - growth_t0;

                    double dedup_t0 = omp_get_wtime();
                    bool accept_curvelet = false;
                    if (edge_chain_lidx > 2) {
                        accept_curvelet = !check_curvelet_exist(edge_chain_lidx, edge_chain_on_the_fly, edge_chain_target);
                    }
                    local_dedup += omp_get_wtime() - dedup_t0;

                    if (accept_curvelet) {
                        
                        for (unsigned chain_idx = 0; chain_idx < edge_chain_lidx; chain_idx++) {
                            edge_chain_target(edge_chain_target_idx, chain_idx) = edge_chain_on_the_fly[chain_idx];
                        }

                        edge_chain_target(edge_chain_target_idx, _group_max_sz) = edge_chain_lidx;
                        edge_chain_target_idx++;

                        record_curvelet_chain((unsigned)te_id, edge_chain_on_the_fly, edge_chain_lidx,
                                              f_run == 0, te_pt_x, te_pt_y, te_orient,
                                              cmp_bundle_min_ks, cmp_bundle_max_ks);
                    }

                    for (unsigned i = 0; i < _group_max_sz; i++) {
                        edge_chain_on_the_fly[i] = 0;
                    }

                    growth_t0 = omp_get_wtime();
                    
                }
                local_chain_growth += omp_get_wtime() - growth_t0;

                for (unsigned i = 0; i < (unsigned)_max_num_look_edges; i++) {
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
        } //> for loop over all t_id

        delete[] bundle_min_ks;
        delete[] bundle_max_ks;
        delete[] cmp_bundle_min_ks;
        delete[] cmp_bundle_max_ks;
        delete[] intersect_bundle_min_ks;
        delete[] intersect_bundle_max_ks;
        delete[] hyp_LookEdge;

        delete[] edge_chain_on_the_fly;
        delete[] edge_chain_target;

        #pragma omp atomic
        time_direction_filter += local_direction_filter;
        #pragma omp atomic
        time_bundle_transport += local_bundle_transport;
        #pragma omp atomic
        time_chain_growth += local_chain_growth;
        #pragma omp atomic
        time_dedup += local_dedup;
    }

    compact_curvelet_output();

    double curvelet_build_time = omp_get_wtime() - start;
    const double time_pairwise = time_direction_filter + time_bundle_transport;
    const double phase_cpu_total = time_pairwise + time_chain_growth + time_dedup;
    if (omp_threads > 1)
        std::cout<<"- Time of curvelet building (OpenMP, "<<omp_threads<<" threads): "<<curvelet_build_time*1000<<" (ms)"<<std::endl;
    else
        std::cout<<"- Time of curvelet building: "<<curvelet_build_time*1000<<" (ms)"<<std::endl;
    if (phase_cpu_total > 0.0) {
        std::cout<<"-   pairwise curve bundle formation: "<<time_pairwise*1000<<" (ms)"
                 <<"  ["<<(100.0*time_pairwise/phase_cpu_total)<<"% of compute]"<<std::endl;
        if (time_pairwise > 0.0) {
            std::cout<<"-     direction filtering: "<<time_direction_filter*1000<<" (ms)" << "  ["<<(100.0*time_direction_filter/time_pairwise)<<"% of pairwise]"<<std::endl;
            std::cout<<"-     bundle transport (compute_curve_bundle): "<<time_bundle_transport*1000<<" (ms)" <<"  ["<<(100.0*time_bundle_transport/time_pairwise)<<"% of pairwise]"<<std::endl;
        }
        std::cout<<"-   chain growth by bundle intersection: "<<time_chain_growth*1000<<" (ms)" << "  ["<<(100.0*time_chain_growth/phase_cpu_total)<<"% of compute]"<<std::endl;
        std::cout<<"-   curvelet dedup (check_curvelet_exist): "<<time_dedup*1000<<" (ms)" << "  [" << (100.0*time_dedup/phase_cpu_total)<<"% of compute]"<<std::endl;
    } 
    else {
        std::cout<<"-   pairwise curve bundle formation: 0 (ms)"<<std::endl;
        std::cout<<"-     direction filtering: 0 (ms)"<<std::endl;
        std::cout<<"-     bundle transport (compute_curve_bundle): 0 (ms)"<<std::endl;
        std::cout<<"-   chain growth by bundle intersection: 0 (ms)"<<std::endl;
        std::cout<<"-   curvelet dedup (check_curvelet_exist): 0 (ms)"<<std::endl;
    }
    std::cout<<"- Number of curvelets formed: "<<_num_curvelets<<std::endl;
}

//> compute the max and min curvature of the curve bundle from a pair of edges
template<typename T>
bool CurveletCPU<T>::compute_curve_bundle( unsigned te_idx, unsigned le_idx, T* bundle_min_ks, T* bundle_max_ks, bool* hyp_LookEdge,
                                           T te_pt_x, T te_pt_y, T le_pt_x, T le_pt_y, T te_orient, T le_orient ) 
{
    
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

                //> k_max / k_min equation denominator w.r.t. pixel (alternate transport branch)
                eq_denom = -(2.0 * delta_x2 * tt - DX2 - DY2 + _dx2 - delta_x2 * delta_x2) / 2.0;

                //> max and min curvatures in terms of pixel
                k_dx_min = (ttmdx2 - _dx) / eq_denom;
                k_dx_max = (ttmdx2 + _dx) / eq_denom;

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
                k_dt_max = sin_t2_p_dt2_m_t0_m_dt / (-sin_t2_p_dt2_m_t0_m_dt * delta_x2 + cos_T0mdt_times_DX + sin_T0mdt_times_DY);
                if ( bundle_max_ks(le_idx, bk_idx) > k_dt_max )
                    bundle_max_ks(le_idx, bk_idx) = k_dt_max;

                k_dt_min = sin_t2_p_dt2_m_t0_p_dt / (-sin_t2_p_dt2_m_t0_p_dt * delta_x2 + cos_T0pdt_times_DX + sin_T0pdt_times_DY);
                if ( bundle_min_ks(le_idx, bk_idx) < k_dt_min )
                    bundle_min_ks(le_idx, bk_idx) = k_dt_min;

                bk_idx++;
            }
        }
    }

    //> check whether the built bundle is valid or not
    bool valid = bundle_valid_check( le_idx, bundle_min_ks, bundle_max_ks );

    //> this boolean list indicates which look edge is valid to be proceeded
    hyp_LookEdge[ le_idx ] = (valid) ? true : false;

    return (valid) ? true : false;
}

template<typename T>
bool CurveletCPU<T>::bundle_valid_check( unsigned le_idx, T* bundle_min_ks, T* bundle_max_ks )
{    
    //> check the curve bundle is valid or not
    for (int k_idx = 0; k_idx < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; k_idx++ ) {
        if ( bundle_max_ks(le_idx, k_idx) > bundle_min_ks(le_idx, k_idx) ) 
            return true;
    }
    return false;
}

template<typename T>
void CurveletCPU<T>::bundle_intersection(T* cmp_bundle_min_ks, T* cmp_bundle_max_ks, T* intersect_bundle_min_ks, T* intersect_bundle_max_ks)
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
bool CurveletCPU<T>::bundle_intersection_valid_check( T* intersect_bundle_min_ks, T* intersect_bundle_max_ks )
{
    //> check the intersected curve bundle is valid or not
    for (int k_idx = 0; k_idx < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; k_idx++ ) {
        if ( intersect_bundle_max_ks[k_idx] > intersect_bundle_min_ks[k_idx] ) 
            return true;
    }
    return false;
}

template<typename T>
void CurveletCPU<T>::move_to_cmp_bundle( unsigned cmp_idx, unsigned le_idx, bool rep_by_intersection, 
                                         T* bundle_min_ks, T* bundle_max_ks, T* cmp_bundle_min_ks, T* cmp_bundle_max_ks,
                                         T* intersect_bundle_min_ks, T* intersect_bundle_max_ks )
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
bool CurveletCPU<T>::check_curvelet_exist( unsigned edge_chain_on_the_fly_sz, unsigned* edge_chain_on_the_fly, unsigned* edge_chain_target ) 
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
                cvlet_exists = cvlet_exists && (edge_chain_on_the_fly[chain_idx] == edge_chain_target(i, chain_idx));
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
    delete[] _anchor_chain_count;
    delete[] _edge_chain_final;
    delete[] _curvelet_info;
}

#endif // CPU_CURVELET_HPP