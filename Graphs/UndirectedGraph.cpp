//
// Created by ldd on 2021/9/23.
//

#include "UndirectedGraph.h"

void UndirectedGraph::LoadEdges(const string &path, const string &file_name, const string &map_file_) {
    int uid, vid, edge_buf_size, num_of_vtx, num_of_edge;
    long long u, v;
    unordered_map<long long, int> vtx_to_id;
    edge *edge_buf, *tmp_edge_buf;
    std::basic_ifstream<char> graph_in;
    std::basic_ofstream<char> map_out;

    graph_file = path + file_name;
    graph_in = ifstream(graph_file);

    if (map_file_.empty()) {
        map_file = path + "__map__" + file_name;
        map_out = ofstream(map_file);
        num_of_vtx = 0;
    }
    else {
        map_file = map_file_;
        LoadVtx2IdMap(vtx_to_id);
        map_out = ofstream(map_file, std::ios::app);
        num_of_vtx = (int) vtx_to_id.size();
    }

    // load graph
    num_of_edge = 0;
    edge_buf_size = DEFAULT_EDGE_BUF_SIZE;
    edge_buf = new edge[edge_buf_size];
    while (graph_in.peek() != EOF) {

        graph_in >> u >> v;

        // remove self-loop
        if (u != v) {
            auto iter1 = vtx_to_id.find(u);
            if (iter1 != vtx_to_id.end()) {
                uid = iter1->second;
            } else {
                uid = num_of_vtx++;
                vtx_to_id.emplace(u, uid);
                map_out << u << " " << uid << endl;
            }

            auto iter2 = vtx_to_id.find(v);
            if (iter2 != vtx_to_id.end()) {
                vid = iter2->second;
            } else {
                vid = num_of_vtx++;
                vtx_to_id.emplace(v, vid);
                map_out << v << " " << vid << endl;
            }

            if (num_of_edge + 2 > edge_buf_size) {
                edge_buf_size = edge_buf_size << 1;
                tmp_edge_buf = new edge[edge_buf_size];
                memcpy(tmp_edge_buf, edge_buf, num_of_edge * sizeof(edge));
                delete[] edge_buf;
                edge_buf = tmp_edge_buf;
            }
            edge_buf[num_of_edge++] = edge(uid, vid);
            edge_buf[num_of_edge++] = edge(vid, uid);
        }
    }

    graph_in.close();
    map_out.close();

    BuildFromEdgeLst(edge_buf, num_of_vtx, num_of_edge);
    delete[] edge_buf;
}

void UndirectedGraph::PrintStatistics(bool print_detail) {
    cout << "Graph path = \"" + graph_file + "\"" << ", ";
    cout << "|V| = " << n << ", |E| = " << (GetM() >> 1) << endl;
    if (print_detail) {
        long long id_to_vtx[n];
        LoadId2VtxMap(id_to_vtx);
        for (int i = 0; i < n; i++) {
            cout << "[" << id_to_vtx[i] << "]--";
            cout << Array2String(adj_lst[i] + 1, id_to_vtx, adj_lst[i][0]) << endl;
        }
    }
}

void UndirectedGraph::CC(const int *vtx_subset, int length, vector<vector<int>> &ccs) {
    int v;
    bool visit[n], exists[n];

    memset(visit, false, n * sizeof(bool));
    memset(exists, false, n * sizeof(bool));

    for (int i = 0; i < length; i++) exists[vtx_subset[i]] = true;
    for (int i = 0; i < length; i++) {
        v = vtx_subset[i];
        if (!visit[v]) {
            auto cc = vector<int>{v};
            Dfs(visit, exists, v, cc);
            ccs.emplace_back(cc);
        }
    }
}