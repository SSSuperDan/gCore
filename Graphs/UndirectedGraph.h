//
// Created by ldd on 2021/9/23.
//

#ifndef CSM4GMG_UNDIRECTEDGRAPH_H
#define CSM4GMG_UNDIRECTEDGRAPH_H

#include "Graph.h"

class UndirectedGraph : public Graph {
public:
    UndirectedGraph() = default;
    ~UndirectedGraph() = default;

    void LoadEdges(const string &path, const string &file_name, const string &map_file_ = "");
    void PrintStatistics(bool print_detail = false);

    string &GetMapFile() {
        return map_file;
    }

    string &GetGraphFile() {
        return graph_file;
    }

    void LoadVtx2IdMap(unordered_map<long long, int> &vtx_to_id) {
        Graph::LoadVtx2IdMap(map_file, vtx_to_id);
    }

    void LoadId2VtxMap(long long *id_to_vtx) {
        Graph::LoadId2VtxMap(map_file, id_to_vtx);
    }

    // Compute connected components using depth-first search.
    void CC(const int *vtx_subset, int length, vector<vector<int>> &ccs);

private:
    string graph_file{};
    string map_file{};

    inline void Dfs(bool *visit, bool *exists, int v, vector<int> &cc) {
        int u;
        visit[v] = true;
        for (int i = 1; i <= adj_lst[v][0]; i++) {
            u = adj_lst[v][i];
            if (exists[u] && !visit[u]) {
                cc.emplace_back(u);
                Dfs(visit, exists, u, cc);
            }
        }
    }

};


#endif //CSM4GMG_UNDIRECTEDGRAPH_H
