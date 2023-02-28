//
// Created by ldd on 2021/9/23.
//

#include "BipartiteGraph.h"

void
BipartiteGraph::LoadEdges(const string &map_file_, const string &v2_map_file_, const string &file_name) {
    int u, uid, v, vid, num_of_vtx, num_of_vtx2, num_of_edge, edge_buf_size;
    unordered_map<long long, int> vtx_to_id, vtx2_to_id;
    edge *edge_buf, *tmp_edge_buf;

    graph_file = file_name;
    map_file = map_file_;
    v2_map_file = v2_map_file_;

    LoadVtx2IdMap(map_file, vtx_to_id);
    LoadVtx2IdMap(v2_map_file, vtx2_to_id);
    auto map_out = ofstream(map_file, std::ios::app);
    auto v2_map_out = ofstream(v2_map_file, std::ios::app);
    auto graph_in = ifstream(graph_file);

    num_of_vtx = (int) vtx_to_id.size();
    num_of_vtx2 = (int) vtx2_to_id.size();
    num_of_edge = 0;

    edge_buf_size = DEFAULT_EDGE_BUF_SIZE;
    edge_buf = new edge[edge_buf_size];

    while (graph_in.peek() != EOF) {
        graph_in >> u >> v;

        // dealing with vertices
        auto iter1 = vtx_to_id.find(u);
        if (iter1 != vtx_to_id.end()) {
            uid = iter1->second;
        } else {
            uid = num_of_vtx++;
            vtx_to_id.emplace(u, uid);
            map_out << u << " " << uid << endl;
        }

        auto iter2 = vtx2_to_id.find(v);
        if (iter2 != vtx2_to_id.end()) {
            vid = iter2->second;
        } else {
            vid = num_of_vtx2++;
            vtx2_to_id.emplace(v, vid);
            v2_map_out << v << " " << vid << endl;
        }

        if (num_of_edge + 1 > edge_buf_size) {
            edge_buf_size = edge_buf_size << 1;
            tmp_edge_buf = new edge[edge_buf_size];
            memcpy(tmp_edge_buf, edge_buf, num_of_edge * sizeof(edge));
            delete[] edge_buf;
            edge_buf = tmp_edge_buf;
        }
        edge_buf[num_of_edge++] = edge(uid, vid);
    }

    graph_in.close();
    map_out.close();
    v2_map_out.close();

    n2 = num_of_vtx2;
    BuildFromEdgeLst(edge_buf, num_of_vtx, num_of_edge);
    delete[] edge_buf;
}

void BipartiteGraph::BuildInvertedGraph(BipartiteGraph *g) {
    int v, **g_adj_lst = g->GetAdjLst();

    n = g->n2, n2 = g->n, m = g->m;
    graph_file = g->graph_file;
    map_file = g->v2_map_file;
    v2_map_file = g->map_file;

    adj_lst = new int *[n];
    adj_lst_buf = new int[m + n + 1];

    // bin sort
    for (int i = 0; i < n; i++) adj_lst_buf[i] = 1;

    for (int i = 0; i < n2; i++) {
        for (int j = 1; j <= g_adj_lst[i][0]; j++) {
            adj_lst_buf[g_adj_lst[i][j]]++;
        }
    }

    for (int i = 1; i < n; i++) {
        adj_lst_buf[i] += adj_lst_buf[i - 1];
    }

    // build
    for (int i = n - 1; i >= 1; i--) {
        adj_lst[i] = &adj_lst_buf[adj_lst_buf[i - 1]];
        adj_lst[i][0] = 0;
    }
    adj_lst[0] = adj_lst_buf;
    adj_lst[0][0] = 0;

    for (int i = 0; i < n2; i++) {
        for (int j = 1; j <= g_adj_lst[i][0]; j++) {
            v = g_adj_lst[i][j];
            adj_lst[v][++adj_lst[v][0]] = i;
        }
    }
    adj_lst_buf[m + n] = 0; // ended with 0, representing the intra-deg for all (newly added) vertices with no neighbors.
}


void BipartiteGraph::BuildBijectionById(int n_) {
    int i, j;

    n = n_;
    n2 = n_;
    m = n_;
    adj_lst = new int *[n];
    adj_lst_buf = new int[(n << 1) + 1];

    for (i = 0; i < n; i++) {
        j = i << 1;
        adj_lst[i] = &adj_lst_buf[j];
        adj_lst_buf[j] = 1;
        adj_lst_buf[j + 1] = i;
    }
    adj_lst_buf[i << 1] = 0;
}

void BipartiteGraph::Enlarge(int new_num_of_vtx1, int new_num_of_vtx2) {
    n2 = new_num_of_vtx2;
    Graph::Enlarge(new_num_of_vtx1);
}

void BipartiteGraph::PrintStatistics(bool print_detail) {
    cout << "Graph path = \"" + graph_file + "\"" << ", ";
    cout << "|V_1| = " << n << ", |V_2| = " << n2 << ", |E| = " << m << endl;

    if (print_detail) {
        long long id_to_vtx[n], id_to_vtx2[n2];
        LoadId2VtxMap(map_file, id_to_vtx);
        LoadId2VtxMap(v2_map_file, id_to_vtx2);
        for (int i = 0; i < n; i++) {
            cout << "[" << id_to_vtx[i] << "]--";
            cout << Array2String(adj_lst[i] + 1, id_to_vtx2, adj_lst[i][0]) << endl;
        }
    }
}
