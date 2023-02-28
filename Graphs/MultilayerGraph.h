//
// Created by ldd on 2021/9/22.
//

#ifndef CSM4GMG_MULTILAYERGRAPH_H
#define CSM4GMG_MULTILAYERGRAPH_H

#include <utility>

#include "../Core/CoreDec.h"
#include "../Utils/StringUtils.h"
#include "BipartiteGraph.h"

///*
// *  RAND: random
// *  INTRA_DEN_INC: increasing order of intra-layer average density.
// *  INTRA_DEN_DEC: decreasing order of intra-layer average density.
// *  CROSS_DEN_INC: increasing order of cross-layer average density.
// *  CROSS_DEN_DEC: decreasing order of cross-layer average density.
// */
//
//enum ordering {
//    RAND, INTRA_DEN_INC, INTRA_DEN_DEC, CROSS_DEN_INC, CROSS_DEN_DEC
//};

// bipartite graph file
struct bg_file {
    int id1{};
    int id2{};
    string filename{};

    bg_file() = default;
    explicit bg_file(string filename_) : filename(move(filename_)) {}
    bg_file(int id1_, int id2_, string filename_) : id1(id1_), id2(id2_), filename(move(filename_)) {}

    bool operator<(const bg_file &file) const {
        return filename < file.filename;
    }
};

class MultilayerGraph {
public:
    MultilayerGraph() = default;
    ~MultilayerGraph();
    void LoadGraph(const string &input_path, const string &name_);
    void PrintSummary(bool print_detail = false);
    void PrintStatistics();

    UndirectedGraph *GetGraph(int gid) {
        return graph_layers[gid];
    }

    UndirectedGraph *GetPrimaryGraph() {
        return graph_layers[pid];
    }

    BipartiteGraph *GetClgFromPrime(int sid) {
        return ps_clg[sid];
    }

    BipartiteGraph *GetClgToPrime(int sid) {
        return sp_clg[sid];
    }

    string GetGraphName() {
        return name;
    }

    [[nodiscard]] int GetLayerNumber() const {
        return num_of_layers;
    }

    [[nodiscard]] int GetPrimaryGraphId() const {
        return pid;
    }

    [[nodiscard]] bool IsGeneral() const {
        return is_general;
    }

    static int ParseLayerNumber(const string &path);

private:
    UndirectedGraph **graph_layers{nullptr};
    BipartiteGraph **cross_layer_graphs{nullptr};

    BipartiteGraph **ps_clg{nullptr};
    BipartiteGraph **sp_clg{nullptr};

    int num_of_layers{1};
    int num_of_cross_layer_graphs{1};
    int pid{0};

    long num_of_total_vtx{0};
    long num_of_total_ie{0};
    long num_of_total_ce{0};
    
    bool is_general{};

    string path{};
    string name{};

    inline void Enlarge(int *new_num_of_vtx);

    void LoadIntraLayerEdges(vector<string> &intra_layer_file, const bool * is_bijective_layer);

    void LoadCrossLayerEdges(vector<bg_file> &cross_layer_file, const bool * is_bijective_layer = nullptr);

    static void
    ParseConfFile(const string &conf_file, vector<string> &intra_layer_file, vector<bg_file> &cross_layer_file);

};


#endif //CSM4GMG_MULTILAYERGRAPH_H
