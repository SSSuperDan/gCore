//
// Created by ldd on 2021/9/15.
//

#ifndef CSM4GMG_GRAPH_H
#define CSM4GMG_GRAPH_H

#include "../Utils/ArrayUtils.h"
#include "../Header.h"

const int DEFAULT_EDGE_BUF_SIZE = 5000;

struct edge{
    int s{-1};
    int t{-1};
    edge() = default;
    edge(int s_, int t_):s(s_), t(t_){};
    bool operator < (const edge & e) const {
        return s < e.s || (s == e.s && t < e.t);
    }
};

class Graph {
public:
    Graph() = default;
    ~Graph();
    void BuildFromEdgeLst(edge* edge_buf, int num_of_vtx, int num_of_edge);
    void Enlarge(int new_num_of_vtx);
    static void LoadVtx2IdMap(string & map_file, unordered_map<long long, int> &vtx_to_id);
    static void LoadId2VtxMap(string & map_file, long long *id_to_vtx);

    [[nodiscard]] int GetN() const{
        return n;
    }

    [[nodiscard]] int GetM() const {
        return m ;
    }

    [[nodiscard]] int GetMaxDeg() const {
        return max_deg;
    }

    [[nodiscard]] int **GetAdjLst(){
        return adj_lst;
    }

protected:
    int *adj_lst_buf{nullptr};
    int **adj_lst{nullptr};

    int m{0};
    int n{0};
    int max_deg{0};
};


#endif //CSM4GMG_GRAPH_H
