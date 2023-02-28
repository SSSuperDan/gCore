//
// Created by ldd on 2021/9/15.
//

#include "Graph.h"

Graph::~Graph() {
    delete[] adj_lst;
    delete[] adj_lst_buf;
}

// Using CSR to store graph.
void Graph::BuildFromEdgeLst(edge *edge_buf, int num_of_vtx, int num_of_edge) {
    int i, j, pj;

    n = num_of_vtx;
    adj_lst = new int *[n];
    adj_lst_buf = new int[num_of_edge + n + 1];
    sort(edge_buf, edge_buf + num_of_edge);

    i = 0; // index of edge_buf
    j = 0; // index of adj_lst_buf
    for (int v = 0; v < n; v++) {
        adj_lst[v] = &adj_lst_buf[j];
        if (edge_buf[i].s > v || i >= num_of_edge) adj_lst_buf[j++] = 0;
        else {
            pj = j++; // index to store vtx degree
            adj_lst_buf[j++] = edge_buf[i++].t;
            while (i < num_of_edge && edge_buf[i].s == v) {
                if (edge_buf[i].t == edge_buf[i - 1].t) i++;  // duplicated edge
                else adj_lst_buf[j++] = edge_buf[i++].t;
            }
            adj_lst_buf[pj] = j - pj - 1;
            m += adj_lst_buf[pj];
            if (adj_lst_buf[pj] > max_deg) max_deg = adj_lst_buf[pj];
        }
    }
    adj_lst_buf[j] = 0; // ended with 0, representing the intra-deg for all (newly added) vertices with no neighbors.
}

void Graph::Enlarge(int new_num_of_vtx) {
    if (new_num_of_vtx <= n) return;

    int **tmp_adj_lst;

    tmp_adj_lst = new int *[new_num_of_vtx];
    memcpy(tmp_adj_lst, adj_lst, n * sizeof(int *));
    for (int i = n; i < new_num_of_vtx; i++) tmp_adj_lst[i] = &adj_lst_buf[m + n];
    n = new_num_of_vtx;

    delete[] adj_lst;
    adj_lst = tmp_adj_lst;
}

void Graph::LoadVtx2IdMap(string & map_file, unordered_map<long long, int> &vtx_to_id) {
    long long u;
    int uid;

    ifstream fin(map_file);
    while (fin.good() && !fin.eof()) {
        fin >> u >> uid;
        vtx_to_id.emplace(u, uid);
    }
    fin.close();
}

void Graph::LoadId2VtxMap(string & map_file, long long *id_to_vtx) {
    long long int u;
    int uid;

    ifstream fin(map_file);
    while (fin.good() && !fin.eof()) {
        fin >> u >> uid;
        id_to_vtx[uid] = u;
    }
    fin.close();
}

