//
// Created by ldd on 2021/10/30.
//

#include "RCD.h"

RCD::RCD(MultilayerGraph &mg_) : mg(mg_) {

    ln = mg.GetLayerNumber();
    pid = mg.GetPrimaryGraphId();

    k_vec = new int[ln];
    p_vec = new int[ln - 1];

    core_indices = new CoreIndex[ln];
    num_of_vtx = new int[ln];
    degs = new int *[ln];

    ps_degs = new int *[ln - 1];
    sp_degs = new int *[ln - 1];

    for (int i = 0; i < ln; i++) num_of_vtx[i] = mg.GetGraph(i)->GetN();
    for (int i = 0; i < ln; i++) degs[i] = new int[num_of_vtx[i]];
    for (int i = 0; i < ln; i++) core_indices[i].Init(num_of_vtx[i]);
    for (int i = 0; i < ln - 1; i++) ps_degs[i] = new int[num_of_vtx[pid]];
    for (int i = 0; i < ln - 1; i++) sp_degs[i] = new int[num_of_vtx[i]];
}

RCD::~RCD() {
    delete[] k_vec;
    delete[] p_vec;

    delete[] core_indices;
    delete[] num_of_vtx;
    if (degs) {
        for (int i = 0; i < ln; i++) delete[] degs[i];
        delete[] degs;
    }
    if (ps_degs) {
        for (int i = 0; i < ln - 1; i++) delete[] ps_degs[i];
        delete[] ps_degs;
    }
    if (sp_degs) {
        for (int i = 0; i < ln - 1; i++) delete[] sp_degs[i];
        delete[] sp_degs;
    }
}


int RCD::Extract(vector<int> &k_vec_, vector<int> &p_vec_, int *r_com, int &length) {
    int num_of_iterations;
    CoreIndex *p_core_index;

    // Init
    Init();
    memcpy(k_vec, k_vec_.data(), ln * sizeof(int));
    memcpy(p_vec, p_vec_.data(), (ln - 1) * sizeof(int));

    for (int i = 0; i < ln; i++) {
        InitKCore(i);
        Peel(i);
    }
    InitPartialCDegs();

    num_of_iterations = 0;
    p_core_index = &core_indices[pid];

    while (true) {
        num_of_iterations += 1;
        PeelPVtx();
        for (int i = 0; i < ln - 1; i++) {
            if (core_indices[i].s != core_indices[i].e) {
                PeelSVtx(i);
            }
        }
        // PrintCoreIndicesOffsets();
        if (p_core_index->s == p_core_index->e) break;
    }

    length = p_core_index->n - p_core_index->e;
    memcpy(r_com, p_core_index->vert + p_core_index->e, length * sizeof(int));
    return num_of_iterations;
}

void RCD::Init() {
    int n, **adj_lst;

    for (int i = 0; i < ln; i++) core_indices[i].Set();
    for (int i = 0; i < ln; i++) {
        n = num_of_vtx[i];
        adj_lst = mg.GetGraph(i)->GetAdjLst();
        for (int j = 0; j < n; j++) {
            degs[i][j] = adj_lst[j][0];
        }
    }
}

void RCD::InitKCore(int gid) {
    int k, n, *deg;
    CoreIndex *core_index;

    k = k_vec[gid];
    deg = degs[gid];
    core_index = &core_indices[gid];
    n = core_index->n;

    for (int i = 0; i < n; i++) {
        if (deg[i] < k) {
            core_index->Remove(i);
        }
    }
}

void RCD::Peel(int gid) {
    int k, v, u, *deg, **adj_lst;
    CoreIndex *core_index;

    k = k_vec[gid];
    deg = degs[gid];
    adj_lst = mg.GetGraph(gid)->GetAdjLst();
    core_index = &core_indices[gid];
    auto &s = core_index->s;
    auto &e = core_index->e;

    while (s + core_index->n < (e << 1)) {  // deletes more
        s = e;
        for (int i = s; i < core_index->n; i++) {
            v = core_index->vert[i];
            deg[v] = 0;
            for (int j = 1; j <= adj_lst[v][0]; j++) {
                if (core_index->pos[adj_lst[v][j]] >= s) deg[v]++;
            }
            if (deg[v] < k) core_index->Remove(v);
        }
    }
    for (; s < e; s++) {
        v = core_index->vert[s];
        for (int i = 1; i <= adj_lst[v][0]; i++) {
            u = adj_lst[v][i];
            if (core_index->pos[u] >= e) {
                deg[u]--;
                if (deg[u] < k) core_index->Remove(u);
            }
        }
    }
}

void RCD::InitPartialCDegs() {
    int v, u, **ps_adj_lst, **sp_adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    // init ps_degs
    p_core_index = &core_indices[pid];

    for (int i = 0; i < ln - 1; i++) {
        s_core_index = &core_indices[i];
        sp_adj_lst = mg.GetClgToPrime(i)->GetAdjLst();

        if (s_core_index->n < (s_core_index->e << 1)) {  // deletes more
            for (int j = p_core_index->e; j < p_core_index->n; j++) {
                ps_degs[i][p_core_index->vert[j]] = 0;
            }

            for (int j = s_core_index->e; j < s_core_index->n; j++) {
                v = s_core_index->vert[j];
                for (int l = 1; l <= sp_adj_lst[v][0]; l++) {
                    u = sp_adj_lst[v][l];
                    if (p_core_index->pos[u] >= p_core_index->e) ps_degs[i][u]++;
                }
            }
        } else {  // remains more
            ps_adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();

            for (int j = p_core_index->e; j < p_core_index->n; j++) {
                v = p_core_index->vert[j];
                ps_degs[i][v] = ps_adj_lst[v][0];
            }

            for (int j = 0; j < s_core_index->e; j++) {
                v = s_core_index->vert[j];
                for (int l = 1; l <= sp_adj_lst[v][0]; l++) {
                    u = sp_adj_lst[v][l];
                    if (p_core_index->pos[u] >= p_core_index->e) ps_degs[i][u]--;
                }
            }
        }
    }

    // init sp_degs
    for (int i = 0; i < ln - 1; i++) {
        s_core_index = &core_indices[i];
        ps_adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();

        if (p_core_index->n < (p_core_index->e << 1)) {  // deletes more
            for (int j = s_core_index->e; j < s_core_index->n; j++) {
                sp_degs[i][s_core_index->vert[j]] = 0;
            }

            for (int j = p_core_index->e; j < p_core_index->n; j++) {
                v = p_core_index->vert[j];
                for (int l = 1; l <= ps_adj_lst[v][0]; l++) {
                    u = ps_adj_lst[v][l];
                    if (s_core_index->pos[u] >= s_core_index->e) sp_degs[i][u]++;
                }
            }
        } else {  // remains more
            sp_adj_lst = mg.GetClgToPrime(i)->GetAdjLst();

            for (int j = s_core_index->e; j < s_core_index->n; j++) {
                v = s_core_index->vert[j];
                sp_degs[i][v] = sp_adj_lst[v][0];
            }

            for (int j = 0; j < p_core_index->e; j++) {
                v = p_core_index->vert[j];
                for (int l = 1; l <= ps_adj_lst[v][0]; l++) {
                    u = ps_adj_lst[v][l];
                    if (s_core_index->pos[u] >= s_core_index->e) sp_degs[i][u]--;
                }
            }
        }
    }

    ReduceVtx();
}

void RCD::ReduceVtx() {
    int v, p;
    CoreIndex *core_index;

    // Reduce primary layer vertices.
    core_index = &core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        p = p_vec[i];
        for (int j = core_index->e; j < core_index->n; j++) {
            v = core_index->vert[j];
            if (ps_degs[i][v] < p) {
                core_index->Remove(v);
            }
        }
    }

    // Reduce secondary layer vertices.
    for (int i = 0; i < ln - 1; i++) {
        core_index = &core_indices[i];
        for (int j = core_index->e; j < core_index->n; j++) {
            v = core_index->vert[j];
            if (sp_degs[i][v] == 0) core_index->Remove(v);
        }
    }
}

void RCD::PeelPVtx() {
    int v, u, k, old_s, *deg, *sp_deg, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    p_core_index = &core_indices[pid];
    auto &s = p_core_index->s;
    auto &e = p_core_index->e;
    old_s = s;

    k = k_vec[pid];
    deg = degs[pid];
    adj_lst = mg.GetPrimaryGraph()->GetAdjLst();

    while (s + p_core_index->n < (e << 1)) {  // deletes more
        s = e;
        for (int i = s; i < p_core_index->n; i++) {
            v = p_core_index->vert[i];
            deg[v] = 0;
            for (int j = 1; j <= adj_lst[v][0]; j++) {
                if (p_core_index->pos[adj_lst[v][j]] >= s) deg[v]++;
            }
            if (deg[v] < k) p_core_index->Remove(v);
        }
    }

    for (; s < e; s++) {  // remains more
        v = p_core_index->vert[s];
        for (int i = 1; i <= adj_lst[v][0]; i++) {
            u = adj_lst[v][i];
            if (p_core_index->pos[u] >= e) {
                deg[u]--;
                if (deg[u] < k) p_core_index->Remove(u);
            }
        }
    }

    for (int i = 0; i < ln - 1; i++) {
        sp_deg = sp_degs[i];
        s_core_index = &core_indices[i];
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();  // ps_adj_lst

        if (old_s + p_core_index->n < (e << 1)) {  // deletes more
            for (int j = s_core_index->e; j < s_core_index->n; j++) {
                sp_deg[s_core_index->vert[j]] = 0;
            }

            for (int j = e; j < p_core_index->n; j++) {
                v = p_core_index->vert[j];
                for (int l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (s_core_index->pos[u] >= s_core_index->e) {
                        sp_deg[u]++;
                    }
                }
            }

            for (int j = s_core_index->e; j < s_core_index->n; j++) {
                v = s_core_index->vert[j];
                if (sp_deg[v] == 0) s_core_index->Remove(v);
            }

        } else {  // remains more
            for (int j = old_s; j < e; j++) {
                v = p_core_index->vert[j];
                for (int l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (s_core_index->pos[u] >= s_core_index->e) {
                        sp_deg[u]--;
                        if (sp_deg[u] == 0) s_core_index->Remove(u);
                    }
                }
            }
        }
    }
}

void RCD::PeelSVtx(int sid) {
    int v, u, k, p, old_s, *deg, *ps_deg, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    s_core_index = &core_indices[sid];
    auto &s = s_core_index->s;
    auto &e = s_core_index->e;
    old_s = s;

    k = k_vec[sid];
    p = p_vec[sid];
    deg = degs[sid];
    adj_lst = mg.GetGraph(sid)->GetAdjLst();

    while (s + s_core_index->n < (e << 1)) {  // deletes more
        s = e;
        for (int i = s; i < s_core_index->n; i++) {
            v = s_core_index->vert[i];
            deg[v] = 0;
            for (int j = 1; j <= adj_lst[v][0]; j++) {
                if (s_core_index->pos[adj_lst[v][j]] >= s) deg[v]++;
            }
            if (deg[v] < k) s_core_index->Remove(v);
        }
    }

    for (; s < e; s++) {
        v = s_core_index->vert[s];
        for (int i = 1; i <= adj_lst[v][0]; i++) {
            u = adj_lst[v][i];
            if (s_core_index->pos[u] >= e) {
                deg[u]--;
                if (deg[u] < k) s_core_index->Remove(u);
            }
        }
    }

    ps_deg = ps_degs[sid];
    p_core_index = &core_indices[pid];
    adj_lst = mg.GetClgToPrime(sid)->GetAdjLst();  // sp_adj_lst

    if (old_s + s_core_index->n < (e << 1)) {  // deletes more
        for (int j = p_core_index->e; j < p_core_index->n; j++) {
            ps_deg[p_core_index->vert[j]] = 0;
        }

        for (int j = e; j < s_core_index->n; j++) {
            v = s_core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                if (p_core_index->pos[u] >= p_core_index->e) {
                    ps_deg[u]++;
                }
            }
        }

        for (int j = p_core_index->e; j < p_core_index->n; j++) {
            v = p_core_index->vert[j];
            if (ps_deg[v] < p) p_core_index->Remove(v);
        }

    } else {  // remains more
        for (int j = old_s; j < e; j++) {
            v = s_core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                if (p_core_index->pos[u] >= p_core_index->e) {
                    ps_deg[u]--;
                    if (ps_deg[u] < p) p_core_index->Remove(u);
                }
            }
        }
    }
}
