//
// Created by ldd on 2021/10/4.
//

#include "KPCore.h"

KPCore::KPCore(MultilayerGraph &mg_) : mg(mg_), ps_degs(nullptr), sp_degs(nullptr) {

    ln = mg.GetLayerNumber();
    pid = mg.GetPrimaryGraphId();

    k_vec = new int[ln];

    core_indices = new CoreIndex[ln];
    num_of_vtx = new int[ln];
    degs = new int *[ln];
    demands = new int *[ln - 1];

    for (int i = 0; i < ln; i++) num_of_vtx[i] = mg.GetGraph(i)->GetN();
    for (int i = 0; i < ln; i++) degs[i] = new int[num_of_vtx[i]];
    for (int i = 0; i < ln; i++) core_indices[i].Init(num_of_vtx[i]);
    for (int i = 0; i < ln - 1; i++) demands[i] = new int[num_of_vtx[pid]];
}

KPCore::~KPCore() {
    delete[] k_vec;
    delete[] core_indices;
    delete[] num_of_vtx;
    if (degs) {
        for (int i = 0; i < ln; i++) delete[] degs[i];
        delete[] degs;
    }
    if (demands) {
        for (int i = 0; i < ln - 1; i++) delete[] demands[i];
        delete[] demands;
    }
}


int KPCore::Extract(vector<int> &k_vec_, vector<float> &p_vec, int *core, int &length, kpcore_extractor extractor) {
    int n, ** adj_lst;

    Init();
    memcpy(k_vec, k_vec_.data(), ln * sizeof(int));
    n = num_of_vtx[pid];
    for (int i = 0; i < ln - 1; i++) {
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();  // ps_adj_lst
        for (int j = 0; j < n; j++) {
            demands[i][j] = int(ceil(float(std::max(1, adj_lst[j][0])) * p_vec[i]));
        }
    }

    if (extractor == WOCDEG) return ExtractImplWoCDeg(core, length);
    else if (extractor == WCDEG) return ExtractImplWCDeg(core, length);
    else if (extractor == MIX) return ExtractImplMix(core, length);
    else return ExtractImplOpt(core, length);
}

int KPCore::Extract(vector<int> &k_vec_, vector<Frac> &p_vec, int *core, int &length, kpcore_extractor extractor) {
    int n, ** adj_lst;

    Init();
    memcpy(k_vec, k_vec_.data(), ln * sizeof(int));
    n = num_of_vtx[pid];
    for (int i = 0; i < ln - 1; i++) {
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();  // ps_adj_lst
        for (int j = 0; j < n; j++) {
            demands[i][j] = int(ceil(float(std::max(1, adj_lst[j][0])) * float(p_vec[i].num) / float(p_vec[i].den)));
        }
    }

    if (extractor == WOCDEG) return ExtractImplWoCDeg(core, length);
    else if (extractor == WCDEG) return ExtractImplWCDeg(core, length);
    else if (extractor == MIX) return ExtractImplMix(core, length);
    else return ExtractImplOpt(core, length);
}

// (k,p)-core extraction without maintaining cross-layer degrees.
int KPCore::ExtractImplWoCDeg(int *core, int &length) {
    int num_of_iterations;
    CoreIndex *p_core_index;

    // Init
    for (int i = 0; i < ln; i++) InitKCore(i);

    num_of_iterations = 0;
    p_core_index = &core_indices[pid];
    while (true) {
        num_of_iterations += 1;

        Peel(pid);
        for (int i = 0; i < ln - 1; i++) {
            FilterSVtxByPNbrs(i);
            if (core_indices[i].s != core_indices[i].e) {
                Peel(i);
                FilterPVtxByFracDemand(i);
            }
        }
        if (p_core_index->s == p_core_index->e) break;
    }

    length = p_core_index->n - p_core_index->e;
    memcpy(core, p_core_index->vert + p_core_index->e, length * sizeof(int));
    return num_of_iterations;
}

// (k,p)-core extraction with maintaining cross-layer degrees.
int KPCore::ExtractImplWCDeg(int *core, int &length) {
    int num_of_iterations;
    CoreIndex *p_core_index;

    // Init
    for (int i = 0; i < ln; i++) InitKCore(i);
    InitClDegs();

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
        if (p_core_index->s == p_core_index->e) break;
    }

    length = p_core_index->n - p_core_index->e;
    memcpy(core, p_core_index->vert + p_core_index->e, length * sizeof(int));
    ReleaseClDegs();
    return num_of_iterations;
}

// Mix of ExtractImplWoCDeg and ExtractImplWCDeg
int KPCore::ExtractImplMix(int *core, int &length) {
    int num_of_iterations;
    CoreIndex *p_core_index;

    // Init
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
        if (p_core_index->s == p_core_index->e) break;
    }

    length = p_core_index->n - p_core_index->e;
    memcpy(core, p_core_index->vert + p_core_index->e, length * sizeof(int));
    ReleaseClDegs();
    return num_of_iterations;
}

// Efficient implementation of ExtractImplMix
int KPCore::ExtractImplOpt(int *core, int &length) {
    int num_of_iterations;
    CoreIndex *p_core_index;

    // Init
    for (int i = 0; i < ln; i++) {
        InitKCore(i);
        PeelOpt(i);
    }
    InitPartialCDegsOpt();

    num_of_iterations = 0;
    p_core_index = &core_indices[pid];

    while (true) {
        num_of_iterations += 1;
        PeelPVtxOpt();
        for (int i = 0; i < ln - 1; i++) {
            if (core_indices[i].s != core_indices[i].e) {
                PeelSVtxOpt(i);
            }
        }
        if (p_core_index->s == p_core_index->e) break;
    }

    length = p_core_index->n - p_core_index->e;
    memcpy(core, p_core_index->vert + p_core_index->e, length * sizeof(int));
    ReleaseClDegs();
    return num_of_iterations;
}

void KPCore::Init() {
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


// Inactivate vertices in layer gid with degree smaller than k_vec[gid].
void KPCore::InitKCore(int gid) {
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

// Init ps_degs and sp_degs, and remove useless vertices.
void KPCore::InitClDegs() {
    int n, **adj_lst;

    ps_degs = new int *[ln - 1];
    n = num_of_vtx[pid];
    for (int i = 0; i < ln - 1; i++) {
        ps_degs[i] = new int[n];
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();
        for (int j = 0; j < n; j++) ps_degs[i][j] = adj_lst[j][0];
    }

    sp_degs = new int *[ln - 1];
    for (int i = 0; i < ln - 1; i++) {
        n = num_of_vtx[i];
        sp_degs[i] = new int[n];
        adj_lst = mg.GetClgToPrime(i)->GetAdjLst();
        for (int j = 0; j < n; j++) sp_degs[i][j] = adj_lst[j][0];
    }

    ReduceVtx();
}

void KPCore::ReduceVtx() {
    int v;
    CoreIndex *core_index;

    // Reduce primary layer vertices.
    core_index = &core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        for (int j = core_index->e; j < core_index->n; j++) {
            v = core_index->vert[j];
            if (ps_degs[i][v] < demands[i][v]) {
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

// Init ps_degs and sp_degs using partially remained subgraph.
void KPCore::InitPartialCDegs() {
    int v, cnt, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    ps_degs = new int *[ln - 1];
    p_core_index = &core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        s_core_index = &core_indices[i];
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();  // ps_adj_lst

        ps_degs[i] = new int[p_core_index->n];
        for (int j = p_core_index->e; j < p_core_index->n; j++) {
            v = p_core_index->vert[j];
            cnt = 0;
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                if (s_core_index->pos[adj_lst[v][l]] >= s_core_index->e) cnt++;
            }
            ps_degs[i][v] = cnt;
        }
    }

    sp_degs = new int *[ln - 1];
    for (int i = 0; i < ln - 1; i++) {
        s_core_index = &core_indices[i];
        adj_lst = mg.GetClgToPrime(i)->GetAdjLst();  // sp_adj_lst

        sp_degs[i] = new int[s_core_index->n];
        for (int j = s_core_index->e; j < s_core_index->n; j++) {
            v = s_core_index->vert[j];
            cnt = 0;
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                if (p_core_index->pos[adj_lst[v][l]] >= p_core_index->e) cnt++;
            }
            sp_degs[i][v] = cnt;
        }
    }

    ReduceVtx();
}

// Release space for ps_degs and sp_degs.
void KPCore::ReleaseClDegs() {
    for (int i = 0; i < ln - 1; i++) {
        delete[] ps_degs[i];
    }
    delete[] ps_degs;

    for (int i = 0; i < ln - 1; i++) {
        delete[] sp_degs[i];
    }
    delete[] sp_degs;
}

// Peel layer gid.
void KPCore::Peel(int gid) {
    int k, v, u, *deg, **adj_lst;
    CoreIndex *core_index;

    k = k_vec[gid];
    deg = degs[gid];
    adj_lst = mg.GetGraph(gid)->GetAdjLst();
    core_index = &core_indices[gid];
    auto &s = core_index->s;
    auto &e = core_index->e;

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

// Remove primary layer vertices with neighbor coverage fraction w.r.t. layer sid smaller than p_vec[sid].
void KPCore::FilterPVtxByFracDemand(int sid) {
    int v, cnt, s_e, *s_pos, *demand, **adj_lst;
    CoreIndex *p_core_index;

    demand = demands[sid];
    s_e = core_indices[sid].e;
    s_pos = core_indices[sid].pos;
    adj_lst = mg.GetClgFromPrime(sid)->GetAdjLst();  // ps_adj_lst
    p_core_index = &core_indices[pid];

    for (int i = p_core_index->e; i < p_core_index->n; i++) {
        v = p_core_index->vert[i];
        cnt = 0;
        for (int j = 1; j <= adj_lst[v][0]; j++) {
            if (s_pos[adj_lst[v][j]] >= s_e) {
                cnt++;
            }
        }
        if (cnt < demand[v]) p_core_index->Remove(v);
    }
}

// Remove secondary layer vertices with no neighbors in primary layer.
void KPCore::FilterSVtxByPNbrs(int sid) {
    bool is_removable;
    int v, p_e, *p_pos, **adj_lst;
    CoreIndex *s_core_index;

    p_e = core_indices[pid].e;
    p_pos = core_indices[pid].pos;
    adj_lst = mg.GetClgToPrime(sid)->GetAdjLst();  // sp_adj_lst
    s_core_index = &core_indices[sid];

    for (int i = s_core_index->e; i < s_core_index->n; i++) {
        v = s_core_index->vert[i];
        is_removable = true;
        for (int j = 1; j <= adj_lst[v][0]; j++) {
            if (p_pos[adj_lst[v][j]] >= p_e) {
                is_removable = false;
                break;
            }
        }
        if (is_removable) s_core_index->Remove(v);
    }
}

// Peel on primary layer, remove affected secondary layer vertices.
void KPCore::PeelPVtx() {
    int v, u, k, old_s, *deg, *sp_deg, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    p_core_index = &core_indices[pid];
    auto &s = p_core_index->s;
    auto &e = p_core_index->e;
    old_s = s;

    k = k_vec[pid];
    deg = degs[pid];
    adj_lst = mg.GetPrimaryGraph()->GetAdjLst();
    for (; s < e; s++) {
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

// Peel on secondary layer, remove affected primary layer vertices.
void KPCore::PeelSVtx(int sid) {
    int v, u, k, old_s, *deg, *ps_deg, *demand, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    s_core_index = &core_indices[sid];
    auto &s = s_core_index->s;
    auto &e = s_core_index->e;
    old_s = s;

    k = k_vec[sid];
    deg = degs[sid];
    adj_lst = mg.GetGraph(sid)->GetAdjLst();
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
    demand = demands[sid];
    p_core_index = &core_indices[pid];
    adj_lst = mg.GetClgToPrime(sid)->GetAdjLst();  // sp_adj_lst
    for (int j = old_s; j < e; j++) {
        v = s_core_index->vert[j];
        for (int l = 1; l <= adj_lst[v][0]; l++) {
            u = adj_lst[v][l];
            if (p_core_index->pos[u] >= p_core_index->e) {
                ps_deg[u]--;
                if (ps_deg[u] < demand[u]) {
                    p_core_index->Remove(u);
                }
            }
        }
    }
}

// Efficient implementation of Peel.
void KPCore::PeelOpt(int gid) {
    int k, v, u, *deg, **adj_lst;
    CoreIndex *core_index;

    k = k_vec[gid];
    deg = degs[gid];
    adj_lst = mg.GetGraph(gid)->GetAdjLst();
    core_index = &core_indices[gid];
    auto &s = core_index->s;
    auto &e = core_index->e;

    while (s + core_index->n < (e << 1)) { // deletes more
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

// Efficient implementation of InitPartialCDegs.
void KPCore::InitPartialCDegsOpt() {
    int v, u, **ps_adj_lst, **sp_adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    // init ps_degs
    ps_degs = new int *[ln - 1];
    p_core_index = &core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        ps_degs[i] = new int[p_core_index->n];
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

    // init sp_deg
    sp_degs = new int *[ln - 1];
    for (int i = 0; i < ln - 1; i++) {
        s_core_index = &core_indices[i];
        sp_degs[i] = new int[s_core_index->n];
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

// Efficient implementation of PeelSVtx.
void KPCore::PeelSVtxOpt(int sid) {
    int v, u, k, old_s, *deg, *ps_deg, *demand, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    s_core_index = &core_indices[sid];
    auto &s = s_core_index->s;
    auto &e = s_core_index->e;
    old_s = s;

    k = k_vec[sid];
    deg = degs[sid];
    adj_lst = mg.GetGraph(sid)->GetAdjLst();

    // Peel secondary layer
    while (s + s_core_index->n < (e << 1)) { // deletes more
        s = e;
        for (int i = s; i < s_core_index->n; i++) {
            v = s_core_index->vert[i];
            deg[v] = 0;
            for (int j = 1; j <= adj_lst[v][0]; j++) {
                if (s_core_index->pos[adj_lst[v][j]] >= s) deg[v]++;  // == > e
            }
            if (deg[v] < k) s_core_index->Remove(v);
        }
    }

    for (; s < e; s++) { // remains more
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
    demand = demands[sid];
    p_core_index = &core_indices[pid];
    adj_lst = mg.GetClgToPrime(sid)->GetAdjLst();  // sp_adj_lst

    // Detect affected primary layer vertices.
    if (old_s + s_core_index->n < (e << 1)) { // deletes more
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
            if (ps_deg[v] < demand[v]) p_core_index->Remove(v);
        }

    } else { // remains more
        for (int j = old_s; j < e; j++) {
            v = s_core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                if (p_core_index->pos[u] >= p_core_index->e) {
                    ps_deg[u]--;
                    if (ps_deg[u] < demand[u]) p_core_index->Remove(u);
                }
            }
        }
    }
}

// Efficient implementation of PeelPVtx.
void KPCore::PeelPVtxOpt() {
    int v, u, k, old_s, *deg, *sp_deg, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    p_core_index = &core_indices[pid];
    auto &s = p_core_index->s;
    auto &e = p_core_index->e;
    old_s = s;

    k = k_vec[pid];
    deg = degs[pid];
    adj_lst = mg.GetPrimaryGraph()->GetAdjLst();

    // Peel primary layer.
    while (s + p_core_index->n < (e << 1)) { // deletes more
        s = e;
        for (int i = s; i < p_core_index->n; i++) {
            v = p_core_index->vert[i];
            deg[v] = 0;
            for (int j = 1; j <= adj_lst[v][0]; j++) {
                if (p_core_index->pos[adj_lst[v][j]] >= s) deg[v]++;   // == > e
            }
            if (deg[v] < k) p_core_index->Remove(v);
        }
    }

    for (; s < e; s++) { // remains more
        v = p_core_index->vert[s];
        for (int i = 1; i <= adj_lst[v][0]; i++) {
            u = adj_lst[v][i];
            if (p_core_index->pos[u] >= e) {
                deg[u]--;
                if (deg[u] < k) p_core_index->Remove(u);
            }
        }
    }

    // Detect affected secondary layer vertices.
    for (int i = 0; i < ln - 1; i++) {
        sp_deg = sp_degs[i];
        s_core_index = &core_indices[i];
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();  // ps_adj_lst

        if (old_s + p_core_index->n < (e << 1)) { // deletes more
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

        } else { // remains more
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


