
#include <cmath>

//> Compute curvature bounds for one anchor-neighbor pair into slot-local bundle arrays.
__device__ __inline__ bool
get_circular_arc_bundle_transport(
    float dx, float dt, float sx, float st, float max_k,
    int curves_num_in_bundle_pixel, int curves_num_in_bundle_theta,
    float *slot_min_ks, float *slot_max_ks,
    float anchor_edge_x, float anchor_edge_y,
    float anchor_cos, float anchor_sin,
    float neighbor_edge_x, float neighbor_edge_y, float neighbor_edge_theta)
{
    unsigned bk_idx = 0;

    const float DX = anchor_edge_x - neighbor_edge_x;
    const float DY = anchor_edge_y - neighbor_edge_y;
    const float DX2 = DX * DX;
    const float DY2 = DY * DY;

    const float sin_T0_p_dt = sinf(neighbor_edge_theta + dt);
    const float cos_T0_p_dt = cosf(neighbor_edge_theta + dt);
    const float sin_T0_m_dt = sinf(neighbor_edge_theta - dt);
    const float cos_T0_m_dt = cosf(neighbor_edge_theta - dt);
    const float sin_T0pdt_times_DY = DY * sin_T0_p_dt;
    const float cos_T0pdt_times_DX = DX * cos_T0_p_dt;
    const float sin_T0mdt_times_DY = DY * sin_T0_m_dt;
    const float cos_T0mdt_times_DX = DX * cos_T0_m_dt;

    const float sin_T2_m_T0_m_dt = anchor_sin * cos_T0_m_dt - anchor_cos * sin_T0_m_dt;
    const float cos_T2_m_T0_m_dt = anchor_cos * cos_T0_m_dt + anchor_sin * sin_T0_m_dt;
    const float sin_T2_m_T0_p_dt = anchor_sin * cos_T0_p_dt - anchor_cos * sin_T0_p_dt;
    const float cos_T2_m_T0_p_dt = anchor_cos * cos_T0_p_dt + anchor_sin * sin_T0_p_dt;

    float k_dx_min, k_dx_max;
    float k_dt_min, k_dt_max;
    float delta_x2;
    float delta_t2;
    float tt;
    float eq_denom;
    float sin_t2_p_dt2_m_t0_m_dt;
    float sin_t2_p_dt2_m_t0_p_dt;
    float cos_dt2;
    float sin_dt2;

    //> If the anchor and neighbor edges are too close to each other, the curve bundle is essentially valid
    if (DX2 + DY2 < 4.f * dx * dx) {
        return true;
    }

    const bool alternate_transport = (DX * anchor_cos + DY * anchor_sin) >= 0.f;

    for (int i = 0; i < curves_num_in_bundle_pixel; i++) {
        for (int j = 0; j < curves_num_in_bundle_theta; j++) {
            delta_x2 = sx * (static_cast<float>(i) - (static_cast<float>(curves_num_in_bundle_pixel) - 1.f) / 2.f);
            delta_t2 = st * (static_cast<float>(j) - (static_cast<float>(curves_num_in_bundle_theta) - 1.f) / 2.f);

            cos_dt2 = cosf(delta_t2);
            sin_dt2 = sinf(delta_t2);
            tt = (anchor_sin * cos_dt2 + anchor_cos * sin_dt2) * DX - (anchor_cos * cos_dt2 - anchor_sin * sin_dt2) * DY;
            eq_denom = (-2.f * delta_x2 * tt + DX2 + DY2 - dx * dx + delta_x2 * delta_x2);

            k_dx_min = 2.f * (tt - delta_x2 - dx) / eq_denom;
            k_dx_max = 2.f * (tt - delta_x2 + dx) / eq_denom;

            //> Update the bundle min and max curvature bounds based on location perturbation
            if (slot_min_ks[bk_idx] < k_dx_min) {
                slot_min_ks[bk_idx] = k_dx_min;
            }
            if (slot_max_ks[bk_idx] > k_dx_max) {
                slot_max_ks[bk_idx] = k_dx_max;
            }

            sin_t2_p_dt2_m_t0_m_dt = sin_T2_m_T0_m_dt * cos_dt2 + cos_T2_m_T0_m_dt * sin_dt2;
            sin_t2_p_dt2_m_t0_p_dt = sin_T2_m_T0_p_dt * cos_dt2 + cos_T2_m_T0_p_dt * sin_dt2;

            if (alternate_transport) {
                k_dt_max = sin_t2_p_dt2_m_t0_m_dt / (-sin_t2_p_dt2_m_t0_m_dt * delta_x2 + cos_T0mdt_times_DX + sin_T0mdt_times_DY);
                k_dt_min = sin_t2_p_dt2_m_t0_p_dt / (-sin_t2_p_dt2_m_t0_p_dt * delta_x2 + cos_T0pdt_times_DX + sin_T0pdt_times_DY);
            } 
            else {
                k_dt_min = sin_t2_p_dt2_m_t0_m_dt / (-sin_t2_p_dt2_m_t0_m_dt * delta_x2 + cos_T0mdt_times_DX + sin_T0mdt_times_DY);
                k_dt_max = sin_t2_p_dt2_m_t0_p_dt / (-sin_t2_p_dt2_m_t0_p_dt * delta_x2 + cos_T0pdt_times_DX + sin_T0pdt_times_DY);
            }

            //> Update the bundle min and max curvature bounds based on orientation perturbation
            if (slot_min_ks[bk_idx] < k_dt_min) {
                slot_min_ks[bk_idx] = k_dt_min;
            }
            if (slot_max_ks[bk_idx] > k_dt_max) {
                slot_max_ks[bk_idx] = k_dt_max;
            }

            bk_idx++;
        }
    }

    const int bundle_cells = curves_num_in_bundle_pixel * curves_num_in_bundle_theta;
    for (int k_idx = 0; k_idx < bundle_cells; k_idx++) {
        //> Check if the curve bundle is valid
        if (slot_max_ks[k_idx] > slot_min_ks[k_idx]) {
            return true;
        }
    }
    return false;
}

//> Copy one slot bundle into the compare buffer (chain growth stage).
__device__ __inline__ void
move_to_cmp_bundle(
    int cmp_idx, bool rep_by_intersection, int bundle_cells,
    float *cmp_bundle_min_ks, float *cmp_bundle_max_ks,
    const float *intersect_bundle_min_ks, const float *intersect_bundle_max_ks,
    const float *slot_min_ks, const float *slot_max_ks)
{
    if (!rep_by_intersection) {
        for (int bidx = 0; bidx < bundle_cells; bidx++) {
            cmp_bundle_min_ks[cmp_idx * bundle_cells + bidx] = slot_min_ks[bidx];
            cmp_bundle_max_ks[cmp_idx * bundle_cells + bidx] = slot_max_ks[bidx];
        }
    } 
    else {
        for (int bidx = 0; bidx < bundle_cells; bidx++) {
            cmp_bundle_min_ks[cmp_idx * bundle_cells + bidx] = intersect_bundle_min_ks[bidx];
            cmp_bundle_max_ks[cmp_idx * bundle_cells + bidx] = intersect_bundle_max_ks[bidx];
        }
    }
}

__device__ __inline__ void
bundle_intersection(
    int bundle_cells,
    const float *cmp_bundle_min_ks, const float *cmp_bundle_max_ks,
    float *intersect_bundle_min_ks, float *intersect_bundle_max_ks)
{
    for (int j = 0; j < bundle_cells; j++) {
        const float min0 = cmp_bundle_min_ks[j];
        const float min1 = cmp_bundle_min_ks[bundle_cells + j];
        intersect_bundle_min_ks[j] = (min0 < min1) ? min1 : min0;

        const float max0 = cmp_bundle_max_ks[j];
        const float max1 = cmp_bundle_max_ks[bundle_cells + j];
        intersect_bundle_max_ks[j] = (max0 > max1) ? max1 : max0;
    }
}

__device__ __inline__ bool
bundle_intersection_valid_check(
    int bundle_cells,
    const float *intersect_bundle_min_ks, const float *intersect_bundle_max_ks)
{
    for (int k_idx = 0; k_idx < bundle_cells; k_idx++) {
        if (intersect_bundle_max_ks[k_idx] > intersect_bundle_min_ks[k_idx]) {
            return true;
        }
    }
    return false;
}
