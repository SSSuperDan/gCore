//
// Created by ldd on 2021/9/23.
//

#ifndef CSM4GMG_BIPARTITEGRAPH_H
#define CSM4GMG_BIPARTITEGRAPH_H

#include "UndirectedGraph.h"

class BipartiteGraph : public Graph {
public:
    void LoadEdges(const string &map_file_, const string &v2_map_file_, const string &file_name);
    void Enlarge(int new_num_of_vtx, int new_num_of_vtx2) ;
    void BuildInvertedGraph(BipartiteGraph * g);
    void BuildBijectionById(int n_);
    void PrintStatistics(bool print_detail);

    [[nodiscard]] int GetN2() const{
        return n2;
    }

private:
    string graph_file{};
    string map_file{};
    string v2_map_file{};
    int n2;
};


#endif //CSM4GMG_BIPARTITEGRAPH_H
