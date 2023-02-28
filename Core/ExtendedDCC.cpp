//
// Created by ldd on 2022/6/23.
//

#include "ExtendedDCC.h"

ExtendedDCC::ExtendedDCC(MultilayerGraph &mg_) : mg(mg_) {
    ln = mg.GetLayerNumber();
    n = mg.GetPrimaryGraph()->GetN();

    k_vec = new int[ln];
    degs = new int *[ln];
    for (int i = 0; i < ln; i++) degs[i] = new int[n];
}

ExtendedDCC::~ExtendedDCC() {
    delete[] k_vec;
    if (degs) {
        for (int i = 0; i < ln; i++) delete[] degs[i];
        delete[] degs;
    }
}

void ExtendedDCC::Extract(vector<int> &k_vec_, int *dcc, int &length, dcc_impl impl) {
    int **adj_lst;

    // Initialize
    memcpy(k_vec, k_vec_.data(), ln * sizeof(int));
    for (int i = 0; i < ln; i++) {
        adj_lst = mg.GetGraph(i)->GetAdjLst();
        for (int j = 0; j < n; j++) {
            degs[i][j] = adj_lst[j][0];
        }
    }

    if (impl == Q) ExtractImplQueue(dcc, length);
    else ExtractImplCoreIndex(dcc, length);
}

void ExtendedDCC::ExtractImplQueue(int *dcc, int &length) {
    int k, v, u, **adj_lst;
    int queue[n], s, e, old_e;
    bool is_removal[n];

    memset(is_removal, false, n * sizeof(bool));

    s = 0, e = 0;
    for (int i = 0; i < ln; i++) {
        k = k_vec[i];
        for (int j = 0; j < n; j++) {
            if (!is_removal[j] && degs[i][j] < k) {
                queue[e++] = j;
                is_removal[j] = true;
            }
        }
    }

    while (s < e) {
        old_e = e;
        for (int i = 0; i < ln; i++) {
            k = k_vec[i];
            adj_lst = mg.GetGraph(i)->GetAdjLst();
            for (int j = s; j < old_e; j++) {
                v = queue[j];
                for (int l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (!is_removal[u]) {
                        degs[i][u]--;
                        if (degs[i][u] < k) {
                            queue[e++] = u;
                            is_removal[u] = true;
                        }
                    }
                }
            }
        }
        s = old_e;
    }

    length = 0;
    for (int i = 0; i < n; i++) {
        if (!is_removal[i]) dcc[length++] = i;
    }
}

void ExtendedDCC::ExtractImplCoreIndex(int *dcc, int &length) {
    int k, v, u, v_neighbor, **adj_lst, old_e;
    CoreIndex core_index;

    // Initialize
    core_index.Init(n);
    core_index.Set();
    auto &s = core_index.s;
    auto &e = core_index.e;

    for (int i = 0; i < ln; i++) {
        k = k_vec[i];
        for (int j = 0; j < n; j++) {
            if (core_index.pos[j] >= e && degs[i][j] < k) {
                core_index.Remove(j);
            }
        }
    }

    while (s + n < (e << 1)) {  // deletes more
        s = e;
        for (int i = 0; i < ln; i++) {
            k = k_vec[i];
            adj_lst = mg.GetGraph(i)->GetAdjLst();
            for (int j = e; j < n; j++) {
                v = core_index.vert[j];
                v_neighbor = 0;
                for (int l = 1; l <= adj_lst[v][0]; l++) {
                    if (core_index.pos[adj_lst[v][l]] >= s) v_neighbor++;
                }
                if (v_neighbor < k) core_index.Remove(v);
                else degs[i][v] = v_neighbor;
            }
        }
    }

    while (s < e) {
        old_e = e;
        for (int i = 0; i < ln; i++) {
            k = k_vec[i];
            adj_lst = mg.GetGraph(i)->GetAdjLst();
            for (int j = s; j < old_e; j++) {
                v = core_index.vert[j];
                for (int l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (core_index.pos[u] >= e) {
                        degs[i][u]--;
                        if (degs[i][u] < k) core_index.Remove(u);
                    }
                }
            }
        }
        s = old_e;
    }

    length = core_index.n - core_index.e;
    memcpy(dcc, core_index.vert + core_index.e, length * sizeof(int));
}
