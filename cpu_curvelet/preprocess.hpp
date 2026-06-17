#ifndef PREPROCESS_HPP
#define PREPROCESS_HPP

#include <map>
#include <vector>
#include <deque>
#include <list>
#include <utility> 
#include <iostream>
#include <cmath>
#include <algorithm>

#include "curvelet_utils.hpp"

// -- class structure storing information of each third-order edge --
template<typename T>
class edgel
{
public:
    T _pt_x;             //< subpix location x of the edgel
    T _pt_y;             //< subpix location y of the edgel
    T _orientation;      //< the orientation of the edgel
    T _strength;         //< the strength of the edgel (typically gradient magnitude)
    int _id;             //< unique id
    
    edgel(T x, T y, T tan, T edge_strength, unsigned id):
        _pt_x(x), _pt_y(y), _orientation(tan), _strength(edge_strength), _id((int)id) {}
    
    // -- copy constructor --
    edgel(const edgel &e):
        _pt_x(e._pt_x), _pt_y(e._pt_y), _orientation(e._orientation),
        _strength(e._strength), _id(e._id) {}
    
    // -- destructor --
    ~edgel() {}
};

// -- class structure building a "map of edges" from the third-order edge list --
template<typename T>
class edgeMap
{    
public:
    //: retinotopic map of edgels
    // std::map is a container that contains key-value pairs
    std::map< std::pair<int,int>, std::vector<edgel<T>*> > _map;
    
    //: local list of edgels for easier traversal
    std::map<unsigned, edgel<T>*> _list;
    
    //: constructor
    edgeMap(int num_edges, int sz_edge_data, T *to_edges)
    {
        // -- TODO?: assert if the TO_edges size is incorrect
        //if ( sz_edges < 4 ){
        //    std::cout << "Edge information size error" <<std::endl;
        //}
        for (unsigned i = 0; i < (unsigned)num_edges; i++)
        {
            //insert(to_edges.val(i,0), to_edges.val(i,1), to_edges.val(i,2), to_edges.val(i,3), i);

            //T& val(const unsigned &r, const unsigned &c) const { return _data[c*_h+r]; }

            const unsigned base = i * (unsigned)sz_edge_data;
            const T x = to_edges[base + 0];
            const T y = to_edges[base + 1];
            const T orient = to_edges[base + 2];
            const T strength = to_edges[base + 3];

            const std::pair<int, int> key((int)std::round(x), (int)std::round(y));
            edgel<T> *edge = new edgel<T>(x, y, orient, strength, i);

            // create an "edge map" _map
            if (_map.find(key) == _map.end() ) {
                std::vector<edgel<T>*> cell;
                cell.push_back(edge);
                _map.insert(std::make_pair(key,cell));
            }
            else {
                // CH: otherwise insert to the _map
                _map[key].push_back(edge);
            }
            _list.insert(std::make_pair(i ,edge));
        }
    }
    
    //: destructor
    ~edgeMap()
    {
        //go over each cell and delete the edgels
        for (typename std::map<unsigned, edgel<T>*>::iterator it = _list.begin(); it != _list.end(); ++it)
            delete it->second;
        
        _map.clear();
        //also clear the list of edgels
        _list.clear();
    }

    void print_map()
    {
        for (typename std::map<std::pair<int,int>, std::vector<edgel<T>*> >::iterator itm = _map.begin();
             itm != _map.end(); ++itm) {
            std::cout<<"Cell ("<<itm->first.first<<","<<itm->first.second<<") Edge ID:";
            for (typename std::vector<edgel<T>*>::iterator itv = itm->second.begin(); itv != itm->second.end(); ++itv) {
                std::cout<<(*itv)->_id<<"  ";
            }
            std::cout<<std::endl;
        }
    }
};

//> group look edges by its idx and squared distance
template<typename T>
struct LookEdgeData 
{
    T dist_to_te;
    edgel<T> *eN;
};

template<typename T>
class edgeNeighborList
{
protected:
    edgeMap<T> _edgeMap;
    unsigned nr;        //< neighbor region
    unsigned img_h;
    unsigned img_w;
    T _rad;
    int _sz_edge_data;
    int _num_edges;
    unsigned _look_slots;

public:
    //> constructor
    edgeNeighborList(int &img_width, int &img_height, int &num_edges, int &sz_edge_data,
                     T *to_edges, int group_mask_sz, const T &rad, unsigned look_slots = 64):
        _edgeMap(num_edges, sz_edge_data, to_edges),
        img_h((unsigned)img_height), img_w((unsigned)img_width),
        _rad(rad), _num_edges(num_edges), _sz_edge_data(sz_edge_data), _look_slots(look_slots)
    {
        // 7x7 spatial bucket (reg=3), matches original form_curvelet_process.cpp
        (void)group_mask_sz;
        nr = 3;
    }

    //> destructor
    ~edgeNeighborList() {}

    unsigned edge_look_stride() const
    {
        return (_sz_edge_data + 1) * _look_slots;
    }

    void init_edgeLookList(T *edgeLookList)
    {
        const unsigned stride = edge_look_stride();
        for (unsigned i = 0; i < (unsigned)_num_edges; i++) {
            for (unsigned j = 0; j < stride; j++) {
                edgeLookList[i * stride + j] = T(-1);
            }
        }
    }

    void create_edgeLookList(T *edgeLookList, unsigned &max_LookEdgeNum)
    {
        //_edgeMap.print_map();
        const unsigned stride = edge_look_stride();
        const T rad_sqr = _rad * _rad;
        unsigned LookEdgeIdx = 1;
        T dist_betwee_te_and_le;

        max_LookEdgeNum = 0;

        for (unsigned i = 0; i < _edgeMap._list.size(); i++) {
            edgel<T> *edgeTarget = _edgeMap._list[i];

            const unsigned x = (unsigned)std::round(edgeTarget->_pt_x);
            const unsigned y = (unsigned)std::round(edgeTarget->_pt_y);

            std::vector<LookEdgeData<T> > le_data;

            //> fill the target edge in the edgeLookList
            fill_edge_to_edgeLookList(edgeLookList, stride, i, edgeTarget, 0, true, le_data);

            LookEdgeIdx = 0;

            //> loop over the 7x7 neighbor
            for (unsigned p = x - nr; p <= (x + nr); p++) {
                for (unsigned q = y - nr; q <= (y + nr); q++) {

                    //> ignore if out of image boundary
                    if (p < 0 || p >= img_w || q < 0 || q >= img_h )
                        continue;

                    const std::pair<int, int> cell_key((int)p, (int)q);
                    typename std::map<std::pair<int, int>, std::vector<edgel<T>*> >::iterator cell_it =
                        _edgeMap._map.find(cell_key);
                    if (cell_it == _edgeMap._map.end())
                        continue;

                    for (unsigned k = 0; k < cell_it->second.size(); k++) 
                    {
                        //> neighbor edge
                        edgel<T> *edgeNeighbor = cell_it->second[k];

                        //> ignore the target edge itself
                        if (edgeNeighbor->_id == edgeTarget->_id) 
                            continue;
                    
                        //> compute the squared distance between the target and the look edge
                        dist_betwee_te_and_le = sq_dist<T>(edgeTarget->_pt_x, edgeTarget->_pt_y,
                                                           edgeNeighbor->_pt_x, edgeNeighbor->_pt_y);

                        //ordered_le_idx_by_dist[LookEdgeIdx-1] = dist_betwee_te_and_le;

                        //> do a better check of circular radius because the bucketing neighborhood is very coarse
                        if ( dist_betwee_te_and_le > rad_sqr) 
                            continue;

                        //> push back to the look edge vector structure
                        LookEdgeData<T> le_data_group = { dist_betwee_te_and_le, edgeNeighbor };
                        le_data.push_back(le_data_group);

                        //> fill the candidate neighbor edges to the edge LookList
                        //fill_edge_to_edgeLookList(edgeLookList, stride, i, LookEdgeIdx, edgeNeighbor, false);
                        LookEdgeIdx++;
                    }
                }
            }

            /*if (i == 0) {
                std::cout<<"Before: ";
                for (unsigned look_idx = 0; look_idx < LookEdgeIdx; look_idx++) {
                    std::cout<<le_data[look_idx].eN->_id<<"  ";
                }
                std::cout<<std::endl;
            }*/
       
            //> sort the look edges in the list by their distances in an ascending manner
            std::sort(std::begin(le_data), std::end(le_data),
                [](const LookEdgeData<T> &a, const LookEdgeData<T> &b) {
                    return a.dist_to_te < b.dist_to_te;
                });

            //> fill look edges to the edgeLookList
            fill_edge_to_edgeLookList(edgeLookList, stride, i, edgeTarget, LookEdgeIdx, false, le_data);

            //> find the maximal number of look edges among all target edges
            if (LookEdgeIdx > max_LookEdgeNum) {
                max_LookEdgeNum = LookEdgeIdx;
            }

            /*if (i == 0) {
                std::cout<<"After: ";
                for (unsigned look_idx = 0; look_idx < LookEdgeIdx; look_idx++) {
                    std::cout<<le_data[look_idx].eN->_id<<"  ";
                }
                std::cout<<std::endl;
            }*/
        }

        std::cout<<"maxmial number of look edges (from preprocess) = "<<max_LookEdgeNum<<std::endl;

        //print_edgeLookList(edgeLookList);
    }

    void fill_edge_to_edgeLookList(T *edgeLookList, unsigned stride, unsigned target_idx, edgel<T> *e,
                                   unsigned LookEdge_num, bool ToFinalEdgeLookList,
                                   std::vector<LookEdgeData<T> > &le) 
    {
        //> fill edge information to the edgeLookList
        if (ToFinalEdgeLookList) {
            //> target edge data
            const unsigned base = target_idx * stride;
            edgeLookList[base + 0] = T(e->_id);
            edgeLookList[base + 1] = e->_pt_x;
            edgeLookList[base + 2] = e->_pt_y;
            edgeLookList[base + 3] = e->_orientation;
            edgeLookList[base + 4] = e->_strength;
        }
        else {
        
            //> all look edge data associated to the target edge
            for (unsigned look_idx = 0; look_idx < LookEdge_num; look_idx++) {
                const unsigned base = target_idx * stride + (look_idx + 1) * 5;
                edgeLookList[base + 0] = T(le[look_idx].eN->_id);
                edgeLookList[base + 1] = le[look_idx].eN->_pt_x;
                edgeLookList[base + 2] = le[look_idx].eN->_pt_y;
                edgeLookList[base + 3] = le[look_idx].eN->_orientation;
                edgeLookList[base + 4] = le[look_idx].eN->_strength;
            }
        }
    }

    void print_edgeLookList(T *edgeLookList) 
    {
        const unsigned stride = edge_look_stride();
        /*for (unsigned i = 0; i < (unsigned)_num_edges; i++) {
            std::cout<<"Cell ID #"<<edgeLookList[i * stride + 0]<<" Look ID:";
            for (unsigned j = 0; j < _look_slots; j++) {
                if (edgeLookList[i * stride + (j+1)*5] < 0)
                    break;
                std::cout<<edgeLookList[i * stride + (j+1)*5]<<" ";
            }
            std::cout<<std::endl;
        }*/

        for (unsigned i = 20837; i < 20840; i++) {
            std::cout<<"Cell ID #"<<edgeLookList[i * stride + 0]<<" Look ID:";
            for (unsigned j = 0; j < _look_slots; j++) {
                if (edgeLookList[i * stride + (j+1)*5] < 0)
                    break;
                std::cout<<edgeLookList[i * stride + (j+1)*5]<<" ";
            }
            std::cout<<std::endl;
        }
    }
};

#endif // PREPROCESS_HPP
