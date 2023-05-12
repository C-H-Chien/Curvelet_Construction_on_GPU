#ifndef PREPROCESS_HPP
#define PREPROCESS_HPP

#include <map>
#include <vector>
#include <deque>
#include <list>
#include <utility> 
#include <iostream>
#include <math.h>
#include <algorithm>

#include "indices.hpp"
#include "curvelet_utils.hpp"

// -- class structure storing information of each third-order edge --
class edgel
{
public:
    float _pt_x;             //< subpix location x of the edgel
    float _pt_y;             //< subpix location y of the edgel
    float _orientation;      //< the orientation of the edgel
    float _strength;         //< the strength of the edgel (typically gradient magnitude)
    int _id;                  //< unique id
    
    edgel(float x, float y, float tan, float edge_strength, unsigned id):
    _pt_x(x), _pt_y(y), _orientation(tan), _strength(edge_strength), _id(id) {};
    
    // -- copy constructor --
    edgel(edgel & e):
    _pt_x(e._pt_x), _pt_y(e._pt_y), _orientation(e._orientation), _strength(e._strength), _id(e._id) {};
    
    // -- destructor --
    ~edgel(){};
};

// -- class structure building a "map of edges" from the third-order edge list --
class edgeMap
{    
public:
    //: retinotopic map of edgels
    // std::map is a container that contains key-value pairs
    std::map< std::pair<int,int>, std::vector<edgel*> > _map;
    
    //: local list of edgels for easier traversal
    std::map<unsigned, edgel*> _list;
    
    //: constructor
    edgeMap(int num_edges, int sz_edge_data, float *TO_edges)
    {
        // -- TODO?: assert if the TO_edges size is incorrect
        //if ( sz_edges < 4 ){
        //    std::cout << "Edge information size error" <<std::endl;
        //}
        for (unsigned i = 0; i < num_edges; i++)
        {
            //insert(TO_edges.val(i,0), TO_edges.val(i,1), TO_edges.val(i,2), TO_edges.val(i,3), i);

            //T& val(const unsigned &r, const unsigned &c) const { return _data[c*_h+r]; }

            std::pair<int, int> key(round(TO_edges(i, 0)),round(TO_edges(i, 1)));
            edgel* edge = new edgel(TO_edges(i, 0), TO_edges(i, 1), TO_edges(i, 2), TO_edges(i, 3), i);

            // create an "edge map" _map
            if (_map.find(key) == _map.end() ) {
                std::vector<edgel*> cell;
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
        //for (std::map<unsigned,edgel*>::iterator it=_list.begin(); it!=_list.end(); ++it)
        //    it->second->unref();
        
        _map.clear();
        //also clear the list of edgels
        _list.clear();
    }

    void print_map()
    {
        for (std::map<std::pair<int,int>, std::vector<edgel*> >::iterator itm=_map.begin(); itm!=_map.end(); ++itm) {
            std::cout<<"Cell ("<<itm->first.first<<","<<itm->first.second<<") Edge ID:";
            for(std::vector<edgel*>::iterator itv=itm->second.begin(); itv!=itm->second.end(); ++itv) {
                std::cout<<(*itv)->_id<<"  ";
            }
            std::cout<<std::endl;
        }
    }
};

//> group look edges by its idx and squared distance
struct LookEdgeData 
{
    float dist_to_te;
    edgel* eN;
};

class edgeNeighborList
{
protected:
    edgeMap _edgeMap;
    unsigned nr;        //< neighbor region
    unsigned img_h;
    unsigned img_w;
    float _rad;
    int _sz_edge_data;
    int _num_edges;

/*private:
    float *unsorted_edgeLookList;          //< temporarily stored unordered look edge list
    unsigned *ordered_le_idx_by_dist;    //< ordered list of look edge index sorted by the squared distance w.r.t. the target edge
*/
public:
    // > constructor
    edgeNeighborList(int &img_width, int &img_height, int &num_edges, int &sz_edge_data, float *TOED_edges, int group_mask_sz, const float &rad):
    _edgeMap(num_edges, sz_edge_data, TOED_edges), img_h(img_height), img_w(img_width), _rad(rad), _num_edges(num_edges), _sz_edge_data(sz_edge_data)
    {
        nr = (group_mask_sz-1) / 2;
        //_edgeMap.print_map();

        //unsorted_edgeLookList = new float[(_sz_edge_data+1) * 32 * _num_edges ];
        //ordered_le_idx_by_dist = new unsigned[ 32 ];
    }

    // > destructor
    ~edgeNeighborList();

    void init_edgeLookList(float *edgeLookList);
    void create_edgeLookList(float *edgeLookList, unsigned &max_LookEdgeNum);
    void fill_edge_to_edgeLookList(float *edgeLookList, unsigned target_idx, edgel* e, unsigned LookEdge_num, bool ToFinalEdgeLookList, std::vector<LookEdgeData> &le);
    void print_edgeLookList(float *edgeLookList);    
};

void edgeNeighborList::init_edgeLookList(float *edgeLookList)
{
    for (unsigned i = 0; i < _num_edges; i++) {
        for (unsigned j = 0; j < (_sz_edge_data+1) * 32; j++) {
            edgeLookList(i, j) = -1;
        }
    }
}

void edgeNeighborList::create_edgeLookList(float *edgeLookList, unsigned &max_LookEdgeNum)
{
    //_edgeMap.print_map();
    unsigned const rad_sqr = _rad*_rad;
    unsigned LookEdgeIdx = 1;
    //unsigned max_LookEdgeNum = 0;
    float dist_betwee_te_and_le;

    for (unsigned i = 0; i < _edgeMap._list.size(); i++) {
        edgel* edgeTarget = _edgeMap._list[i];

        const unsigned x = round(edgeTarget->_pt_x);
        const unsigned y = round(edgeTarget->_pt_y);

        std::vector<LookEdgeData> le_data;

        //> fill the target edge in the edgeLookList
        fill_edge_to_edgeLookList(edgeLookList, i, edgeTarget, 0, true, le_data);

        LookEdgeIdx = 0;

        // > loop over the 7x7 neighbor
        for (unsigned p = x-nr; p <= (x+nr); p++) {
            for (unsigned q = y-nr; q <= (y+nr); q++) {

                //> ignore if out of image boundary
                if (p < 0 || p >= img_w || q < 0 || q >= img_h )
                    continue;

                for (unsigned k = 0; k < _edgeMap._map[std::make_pair(p, q)].size(); k++) 
                {
                    //> neighbor edge
                    edgel* edgeNeighbor = _edgeMap._map[std::make_pair(p, q)][k];

                    //> ignore the target edge itself
                    if (edgeNeighbor->_id == edgeTarget->_id) 
                        continue;
                    
                    //> compute the squared distance between the target and the look edge
                    dist_betwee_te_and_le = sq_dist<float>(edgeTarget->_pt_x, edgeTarget->_pt_y, edgeNeighbor->_pt_x, edgeNeighbor->_pt_y);

                    //ordered_le_idx_by_dist[LookEdgeIdx-1] = dist_betwee_te_and_le;

                    //> do a better check of circular radius because the bucketing neighborhood is very coarse
                    if ( dist_betwee_te_and_le > rad_sqr) 
                        continue;

                    //> push back to the look edge vector structure
                    LookEdgeData le_data_group = { dist_betwee_te_and_le, edgeNeighbor };
                    le_data.push_back(le_data_group);

                    //> fill the candidate neighbor edges to the edge LookList
                    //fill_edge_to_edgeLookList(edgeLookList, i, LookEdgeIdx, edgeNeighbor, false);
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
        std::sort(std::begin(le_data), std::end(le_data), [](const LookEdgeData& a, const LookEdgeData& b) { return a.dist_to_te < b.dist_to_te; });

        //> fill look edges to the edgeLookList
        fill_edge_to_edgeLookList(edgeLookList, i, edgeTarget, LookEdgeIdx, false, le_data);

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

    print_edgeLookList(edgeLookList);
}

void edgeNeighborList::fill_edge_to_edgeLookList(float *edgeLookList, unsigned target_idx, edgel* e, unsigned LookEdge_num, bool ToFinalEdgeLookList, std::vector<LookEdgeData> &le) 
{
    //> fill edge information to the edgeLookList
    if (ToFinalEdgeLookList) {
        //> target edge data
        edgeLookList(target_idx, 0) = e->_id;
        edgeLookList(target_idx, 1) = e->_pt_x;
        edgeLookList(target_idx, 2) = e->_pt_y;
        edgeLookList(target_idx, 3) = e->_orientation;
        edgeLookList(target_idx, 4) = e->_strength;
    }
    else {
        
        //> all look edge data associated to the target edge
        for (unsigned look_idx = 0; look_idx < LookEdge_num; look_idx++) {
            edgeLookList(target_idx, (look_idx+1)*5    ) = le[look_idx].eN->_id;
            edgeLookList(target_idx, (look_idx+1)*5 + 1) = le[look_idx].eN->_pt_x;
            edgeLookList(target_idx, (look_idx+1)*5 + 2) = le[look_idx].eN->_pt_y;
            edgeLookList(target_idx, (look_idx+1)*5 + 3) = le[look_idx].eN->_orientation;
            edgeLookList(target_idx, (look_idx+1)*5 + 4) = le[look_idx].eN->_strength;
        }
    }
}

void edgeNeighborList::print_edgeLookList(float *edgeLookList) 
{
    /*for (unsigned i = 0; i < _num_edges; i++) {
        std::cout<<"Cell ID #"<<edgeLookList(i, 0)<<" Look ID:";
        for (unsigned j = 0; j < 32; j++) {
            if (edgeLookList(i, (j+1)*5) < 0)
                break;
            std::cout<<edgeLookList(i, (j+1)*5)<<" ";
        }
        std::cout<<std::endl;
    }*/

    for (unsigned i = 20837; i < 20840; i++) {
        std::cout<<"Cell ID #"<<edgeLookList(i, 0)<<" Look ID:";
        for (unsigned j = 0; j < 32; j++) {
            if (edgeLookList(i, (j+1)*5) < 0)
                break;
            std::cout<<edgeLookList(i, (j+1)*5)<<" ";
        }
        std::cout<<std::endl;
    }
}

edgeNeighborList::~edgeNeighborList() {
    //> free dynamic allocated memories
    //delete[] unsorted_edgeLookList;
    //delete[] ordered_le_idx_by_dist;
}

#endif // PREPROCESS_HPP