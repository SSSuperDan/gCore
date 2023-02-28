//
// Created by ldd on 2022/6/11.
//

#include "KPTreeBuilder.h"

KPTreeBuilder::KPTreeBuilder(MultilayerGraph &mg_) : PTreeBuilder(mg_), step(0), degeneracy(nullptr) {}

void KPTreeBuilder::Execute(KPTree &kp_tree_, Frac2IntPri &frac2int, p_tree_builder opt, int step_,
                            const vector<int> &start_k_vec) {
    PNode *p_root;

    // Set parameters
    if (opt == NAIVE) BuildSubPTree = &KPTreeBuilder::BuildSubPTreeNaive;
    else if (opt == NE) BuildSubPTree = &KPTreeBuilder::BuildSubPTreeNe;
    else if (opt == SE) BuildSubPTree = &KPTreeBuilder::BuildSubPTreeSe;
    else BuildSubPTree = &KPTreeBuilder::BuildSubPTreeNeSe;
    step = step_;
    degeneracy = new int[ln];
    for (int i = 0; i < ln; i++) degeneracy[i] = CoreDec::GetDegeneracy(mg.GetGraph(i));

    kp_tree = &kp_tree_;
    if (start_k_vec.empty()) memset(k_vec, 0, ln * sizeof(int));
    else memcpy(k_vec, start_k_vec.data(), ln * sizeof(int));

    // Init
    InitKPCore();
    BuildKSubgraph();
    BuildFrac2IntPriMap(frac2int);
    BuildPValueHeap(frac2int);
    memset(p_vec, 0, (ln - 1) * sizeof(int));

    // Build root node
    p_tree = new PTree(ln - 1);
    p_root = p_tree->GetNewPTreeNode(p_vec, 0);
    p_tree->SetRootNode(p_root);
    kp_tree->Insert(k_vec, p_tree);

    (this->*BuildSubPTree)(0, p_root);

    // Recursively build.
    BuildSubKPTree(0);

    RestorePValueHeap();
    delete[] degeneracy;
}

void KPTreeBuilder::Execute(KPTree &kp_tree_, Frac2IntPri &frac2int, const vector<vector<int>> &k_vectors,
                            p_tree_builder opt) {
    PNode *p_root;
    int old_e[ln];

    // Set parameters
    if (opt == NAIVE) BuildSubPTree = &KPTreeBuilder::BuildSubPTreeNaive;
    else if (opt == NE) BuildSubPTree = &KPTreeBuilder::BuildSubPTreeNe;
    else if (opt == SE) BuildSubPTree = &KPTreeBuilder::BuildSubPTreeSe;
    else BuildSubPTree = &KPTreeBuilder::BuildSubPTreeNeSe;

    kp_tree = &kp_tree_;

    // compute minimum k
    for (int i = 0; i < ln; i++) k_vec[i] = INT_MAX;
    for (auto &k : k_vectors) {
        for (int i = 0; i < ln; i++) k_vec[i] = std::min(k_vec[i], k[i]);
    }

    InitKPCore();
    BuildKSubgraph();
    BuildFrac2IntPriMap(frac2int);
    BuildPValueHeap(frac2int);

    for (int i = 0; i < ln; i++) old_e[i] = k_core_indices[i].e;

    memset(p_vec, 0, (ln - 1) * sizeof(int));
    for (auto &k:k_vectors) {

        memcpy(k_vec, k.data(), ln * sizeof(int));
        InitIncKPCore();

        // Build p-tree root.
        p_tree = new PTree(ln - 1);
        p_root = p_tree->GetNewPTreeNode(p_vec, 0);
        p_tree->SetRootNode(p_root);
        kp_tree->Insert(k_vec, p_tree);

        // Recursively build p-tree.
        ShrinkSubgraph(old_e);
        (this->*BuildSubPTree)(0, p_root);

#ifdef MY_DEBUG
        cout << Array2String(k_vec, ln) << " " << p_tree->GetNumOfNodes() << endl;
#else
#endif

        RestoreSubgraph(old_e);
        Restore(old_e);
    }

    RestorePValueHeap();
}

void KPTreeBuilder::InitIncKPCore() {
    int v, k, *deg, old_e;
    CoreIndex *core_index;
    IntLinearHeap *p_value_heap;

    for (int i = 0; i < ln; i++) {
        k = k_vec[i];
        deg = degs[i];
        core_index = &k_core_indices[i];
        old_e = core_index->e;

        for (int j = old_e; j < core_index->n; j++) {
            v = core_index->vert[j];
            if (deg[v] < k) core_index->Remove(v);
        }
    }

    if (core_index->s != core_index->e) {  // p_core_index
        PeelPVtx();

        for (int i = 0; i < ln - 1; i++) {
            p_value_heap = &p_value_heaps[i];
            for (int j = old_e; j < core_index->e; j++) {
                p_value_heap->Remove(core_index->vert[j]);
            }
        }
    }
    for (int i = 0; i < ln - 1; i++) {
        if (k_core_indices[i].s != k_core_indices[i].e) {
            PeelSVtx(i);
        }
    }

}

// Compute (k, 0^{ln-1})-core for a newly generated coreness vector k, gid represents the dimensional changed.
void KPTreeBuilder::InitIncKPCore(int gid) {  // old_s/old/_e ???!!!
    int v, k, *deg, old_e;
    CoreIndex *core_index;
    IntLinearHeap *p_value_heap;

    k = k_vec[gid];
    deg = degs[gid];
    core_index = &k_core_indices[gid];
    old_e = core_index->e;

    for (int j = core_index->e; j < core_index->n; j++) {
        v = core_index->vert[j];
        if (deg[v] < k) core_index->Remove(v);
    }

    // fraction vector = 0^{ln-1},
    // so removing secondary layer vertices lead to no primary layer vertices removed.
    if (old_e < core_index->e) {  // has vertices deleted
        if (gid == pid) {
            PeelPVtx();

            for (int i = 0; i < ln - 1; i++) {
                p_value_heap = &p_value_heaps[i];
                for (int j = old_e; j < core_index->e; j++) {
                    p_value_heap->Remove(core_index->vert[j]);
                }
            }

            for (int i = 0; i < ln - 1; i++) {
                if (k_core_indices[i].s != k_core_indices[i].e) {
                    PeelSVtx(i);
                }
            }

        } else {
            PeelSVtx(gid);
        }
    }
}

// Recursively build sub-kp_tree given root r.
void KPTreeBuilder::BuildSubKPTree(int lk) {
    int old_e[ln];
    PNode *p_root;

    for (int i = 0; i < ln; i++) old_e[i] = k_core_indices[i].e;

    for (int i = lk; i < ln; i++) {

        k_vec[i] += step;

        //if (k_core_indices[i].e < k_core_indices[i].n) {  // has child
        if (k_vec[i] <= degeneracy[i]) {  // has child
            InitIncKPCore(i);

            // Build p-tree.
            p_tree = new PTree(ln - 1);
            p_root = p_tree->GetNewPTreeNode(p_vec, 0);
            p_tree->SetRootNode(p_root);
            kp_tree->Insert(k_vec, p_tree);

            ShrinkSubgraph(i, old_e);
            (this->*BuildSubPTree)(0, p_root);

#ifdef MY_DEBUG
            cout << Array2String(k_vec, ln) << " " << p_tree->GetNumOfNodes() << endl;
#else
#endif

            // Recursively build.
            BuildSubKPTree(i);

            RestoreSubgraph(i, old_e);
            Restore(old_e);
        }

        k_vec[i] -= step;
    }
}

// Generate reduced subgraph w.r.t. the current (k,0^{ln-1})-core.
void KPTreeBuilder::ShrinkSubgraph(int gid, const int *old_e) {

    if (old_e[gid] == k_core_indices[gid].e) return;

    // shrink intra-layer graph
    ShrinkLayer(&k_core_indices[gid], &k_core_indices[gid],
                old_e[gid], k_subgraph[gid], k_subgraph[gid]);

    if (gid == pid) {  // If the primary layer are reduced.
        for (int i = 0; i < ln - 1; i++) {
            // Reduce every affected secondary layers.
            ShrinkLayer(&k_core_indices[pid], &k_core_indices[i],
                        old_e[pid], ps_k_subgraph[i], sp_k_subgraph[i]);

            if (old_e[i] < k_core_indices[i].e) {
                ShrinkLayer(&k_core_indices[i], &k_core_indices[i],
                            old_e[i], k_subgraph[i], k_subgraph[i]);
                ShrinkLayer(&k_core_indices[i], &k_core_indices[pid],
                            old_e[i], sp_k_subgraph[i], ps_k_subgraph[i]);
            }
        }
    } else {
        // Reduce affected primary layers.
        ShrinkLayer(&k_core_indices[gid], &k_core_indices[pid],
                    old_e[gid], sp_k_subgraph[gid], ps_k_subgraph[gid]);
    }
}

void KPTreeBuilder::ShrinkSubgraph(const int *old_e) {

    if (old_e[pid] < k_core_indices[pid].e) {
        ShrinkLayer(&k_core_indices[pid], &k_core_indices[pid],
                    old_e[pid], k_subgraph[pid], k_subgraph[pid]);

        for (int i = 0; i < ln - 1; i++) {
            ShrinkLayer(&k_core_indices[pid], &k_core_indices[i],
                        old_e[pid], ps_k_subgraph[i], sp_k_subgraph[i]);
        }
    }

    for (int i = 0; i < ln - 1; i++) {
        if (old_e[i] < k_core_indices[i].e) {
            ShrinkLayer(&k_core_indices[i], &k_core_indices[i],
                        old_e[i], k_subgraph[i], k_subgraph[i]);
            ShrinkLayer(&k_core_indices[i], &k_core_indices[pid],
                        old_e[i], sp_k_subgraph[i], ps_k_subgraph[i]);
        }
    }
}

// Restore subgraph to the state before reduction.
void KPTreeBuilder::RestoreSubgraph(int gid, const int *old_e) {
    if (old_e[gid] == k_core_indices[gid].e) return;

    RestoreLayer(&k_core_indices[gid], &k_core_indices[gid], old_e[gid],
                 k_subgraph[gid], k_subgraph[gid]);
    if (gid == pid) {
        for (int i = 0; i < ln - 1; i++) {
            RestoreLayer(&k_core_indices[pid], &k_core_indices[i], old_e[pid],
                         ps_k_subgraph[i], sp_k_subgraph[i]);
            if (old_e[i] < k_core_indices[i].e) {
                RestoreLayer(&k_core_indices[i], &k_core_indices[i], old_e[i],
                             k_subgraph[i], k_subgraph[i]);
                RestoreLayer(&k_core_indices[i], &k_core_indices[pid], old_e[i],
                             sp_k_subgraph[i], ps_k_subgraph[i]);
            }
        }
    } else {
        RestoreLayer(&k_core_indices[gid], &k_core_indices[pid], old_e[gid],
                     sp_k_subgraph[gid], ps_k_subgraph[gid]);
    }
}

void KPTreeBuilder::RestoreSubgraph(const int *old_e) {

    if (old_e[pid] < k_core_indices[pid].e) {
        RestoreLayer(&k_core_indices[pid], &k_core_indices[pid], old_e[pid],
                     k_subgraph[pid], k_subgraph[pid]);
        for (int i = 0; i < ln - 1; i++) {
            RestoreLayer(&k_core_indices[pid], &k_core_indices[i], old_e[pid],
                         ps_k_subgraph[i], sp_k_subgraph[i]);
        }
    }

    for (int i = 0; i < ln - 1; i++) {
        if (old_e[i] < k_core_indices[i].e) {
            RestoreLayer(&k_core_indices[i], &k_core_indices[i], old_e[i],
                         k_subgraph[i], k_subgraph[i]);
            RestoreLayer(&k_core_indices[i], &k_core_indices[pid], old_e[i],
                         sp_k_subgraph[i], ps_k_subgraph[i]);
        }
    }
}

// Reduce target layer subgraph according to removed vertices in source layer.
void KPTreeBuilder::ShrinkLayer(CoreIndex *s_core_index, CoreIndex *t_core_index,
                                int s_old_e, int **st_adj_lst, int **ts_adj_lst) {
    int v, u, l, idx, tmp;

    // Identify affected target layer vertices.
    idx = 0;
    for (int i = s_old_e; i < s_core_index->e; i++) {
        v = s_core_index->vert[i];
        for (l = 1; l <= st_adj_lst[v][0]; l++) {
            u = st_adj_lst[v][l];
            if (t_core_index->pos[u] >= t_core_index->e && !aux_sign[u]) {
                aux_arr[idx++] = u;
                aux_sign[u] = true;
            }
        }
    }

    // Move the left source layers neighbors to the front of the adjacent list, and modify degree accordingly.
    for (int i = 0; i < idx; i++) {
        v = aux_arr[i];
        aux_sign[v] = false;

        l = ts_adj_lst[v][0];  // previous degree
        for (int j = 1; j <= l;) {
            if (s_core_index->pos[ts_adj_lst[v][j]] < s_core_index->e) {
                tmp = ts_adj_lst[v][j];
                ts_adj_lst[v][j] = ts_adj_lst[v][l];
                ts_adj_lst[v][l] = tmp;
                l--;
            } else j++;
        }
        ts_adj_lst[v][0] = l;
    }
}

// Restore target layer subgraph according to removed vertices in source layer.
void KPTreeBuilder::RestoreLayer(CoreIndex *s_core_index, CoreIndex *t_core_index, int s_old_e,
                                 int **st_adj_lst, int **ts_adj_lst) {
    int v, u, l;

    // adjacent lists of previously reduced vertices are not reduced,
    // which can therefore be used to restore that of all affected vertices.
    for (int i = s_old_e; i < s_core_index->e; i++) {
        v = s_core_index->vert[i];
        for (l = 1; l <= st_adj_lst[v][0]; l++) {
            u = st_adj_lst[v][l];
            if (t_core_index->pos[u] >= t_core_index->e) {
                ts_adj_lst[u][0]++;
            }
        }
    }
}
