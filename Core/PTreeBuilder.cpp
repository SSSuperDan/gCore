//
// Created by ldd on 2021/10/5.
//

#include "PTreeBuilder.h"

PTreeBuilder::PTreeBuilder(MultilayerGraph &mg_) : mg(mg_) {
    int pn, n, max_n;

    ln = mg.GetLayerNumber();
    pid = mg.GetPrimaryGraphId();

    k_core_indices = new CoreIndex[ln];
    p_value_heaps = new IntLinearHeap[ln - 1];

    ps_k_subgraph = new int **[ln - 1];
    sp_k_subgraph = new int **[ln - 1];
    k_subgraph = new int **[ln];

    init_ps_degs = new int *[ln - 1];
    ps_degs = new int *[ln - 1];
    sp_degs = new int *[ln - 1];
    degs = new int *[ln];

    num_of_vtx = new int[ln];
    k_vec = new int[ln];
    p_vec = new int[ln - 1];

    max_n = 0;
    for (int i = 0; i < ln; i++) {
        n = mg.GetGraph(i)->GetN();
        if (max_n < n) max_n = n;
        num_of_vtx[i] = n;
    }

    for (int i = 0; i < ln; i++) k_core_indices[i].Init(num_of_vtx[i]);
    for (int i = 0; i < ln; i++) degs[i] = new int[num_of_vtx[i]];
    for (int i = 0; i < ln; i++) k_subgraph[i] = new int *[num_of_vtx[i]];
    for (int i = 0; i < ln - 1; i++) sp_degs[i] = new int[num_of_vtx[i]];
    for (int i = 0; i < ln - 1; i++) sp_k_subgraph[i] = new int *[num_of_vtx[i]];

    pn = num_of_vtx[pid];
    for (int i = 0; i < ln - 1; i++) p_value_heaps[i].Init(pn);
    for (int i = 0; i < ln - 1; i++) ps_degs[i] = new int[pn];
    for (int i = 0; i < ln - 1; i++) init_ps_degs[i] = new int[pn];
    for (int i = 0; i < ln - 1; i++) ps_k_subgraph[i] = new int *[pn];

    aux_arr = new int[max_n];
    aux_sign = new bool[max_n];
    memset(aux_sign, false, max_n * sizeof(bool));
}

PTreeBuilder::~PTreeBuilder() {
    delete[] k_vec;
    delete[] p_vec;
    delete[] k_core_indices;
    delete[] p_value_heaps;

    if (ps_k_subgraph) {
        for (int i = 0; i < ln - 1; i++) delete[] ps_k_subgraph[i];
        delete[] ps_k_subgraph;
    }

    if (sp_k_subgraph) {
        for (int i = 0; i < ln - 1; i++) delete[] sp_k_subgraph[i];
        delete[] sp_k_subgraph;
    }

    if (k_subgraph) {
        for (int i = 0; i < ln; i++) delete[] k_subgraph[i];
        delete[] k_subgraph;
    }

    if (init_ps_degs) {
        for (int i = 0; i < ln - 1; i++) delete[] init_ps_degs[i];
        delete[] init_ps_degs;
    }

    if (ps_degs) {
        for (int i = 0; i < ln - 1; i++) delete[] ps_degs[i];
        delete[] ps_degs;
    }

    if (sp_degs) {
        for (int i = 0; i < ln - 1; i++) delete[] sp_degs[i];
        delete[] sp_degs;
    }

    if (degs) {
        for (int i = 0; i < ln; i++) delete[] degs[i];
        delete[] degs;
    }

    delete[] num_of_vtx;
    delete[] aux_arr;
    delete[] aux_sign;
}

void PTreeBuilder::Execute(vector<int> &k_vec_, PTree &p_tree_, Frac2IntPri &frac2int, p_tree_builder opt) {
    PNode *root;

    // Set parameters
    memcpy(k_vec, k_vec_.data(), ln * sizeof(int));
    this->p_tree = &p_tree_;

    // Init
    InitKPCore();
    BuildKSubgraph();
    BuildFrac2IntPriMap(frac2int);
    BuildPValueHeap(frac2int);
    memset(p_vec, 0, (ln - 1) * sizeof(int));

    // Build
    root = p_tree->GetNewPTreeNode(p_vec, 0);
    p_tree->SetRootNode(root);

    if (opt == NAIVE) BuildSubPTreeNaive(0, root);
    else if (opt == NE) BuildSubPTreeNe(0, root);
    else if (opt == SE) BuildSubPTreeSe(0, root, nullptr, nullptr);
    else if (opt == NESE) BuildSubPTreeNeSe(0, root, nullptr, nullptr);

    RestorePValueHeap();
}


// Compute the (k,p)-core represented by the root, i.e., p = 0^{ln-1},
// Initialize essential degree-related structures.
void PTreeBuilder::InitKPCore() {
    int k, n, *deg, **adj_lst;
    CoreIndex *core_index;

    for (int i = 0; i < ln; i++) k_core_indices[i].Set();

    for (int i = 0; i < ln; i++) {
        k = k_vec[i];
        deg = degs[i];
        n = num_of_vtx[i];
        adj_lst = mg.GetGraph(i)->GetAdjLst();

        if (k) {
            core_index = &k_core_indices[i];
            for (int j = 0; j < n; j++) {
                if (adj_lst[j][0] < k) core_index->Remove(j);
                else deg[j] = adj_lst[j][0];
            }
            Peel(i);
        } else for (int j = 0; j < n; j++) deg[j] = adj_lst[j][0];
    }

    InitPartialCDegs();

    for (int i = 0; i < ln - 1; i++) {
        if (k_core_indices[i].s != k_core_indices[i].e) {
            PeelS(i);
        }
    }
}

// Peel layer gid.
void PTreeBuilder::Peel(int gid) {
    int k, v, u, *deg, **adj_lst;
    CoreIndex *core_index;

    k = k_vec[gid];
    deg = degs[gid];
    adj_lst = mg.GetGraph(gid)->GetAdjLst();
    core_index = &k_core_indices[gid];
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

// Peel secondary vertices, modify ps_degs.
void PTreeBuilder::PeelS(int sid) {
    int v, u, old_s, *ps_deg, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    s_core_index = &k_core_indices[sid];
    old_s = s_core_index->s;

    Peel(sid);

    ps_deg = ps_degs[sid];
    p_core_index = &k_core_indices[pid];
    adj_lst = mg.GetClgToPrime(sid)->GetAdjLst();  // sp_adj_lst

    // Modify ps_degs
    if (old_s + s_core_index->n < (s_core_index->e << 1)) {  // deletes more
        for (int j = p_core_index->e; j < p_core_index->n; j++) {
            ps_deg[p_core_index->vert[j]] = 0;
        }

        for (int j = s_core_index->e; j < s_core_index->n; j++) {
            v = s_core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                if (p_core_index->pos[u] >= p_core_index->e) ps_deg[u]++;
            }
        }
    } else {  // remains more
        for (int j = old_s; j < s_core_index->e; j++) {
            v = s_core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                if (p_core_index->pos[u] >= p_core_index->e) ps_deg[u]--;
            }
        }
    }
}

// Init ps_degs and sp_degs.
void PTreeBuilder::InitPartialCDegs() {
    int v, u, *ps_deg, *sp_deg, **ps_adj_lst, **sp_adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    // init ps_degs
    p_core_index = &k_core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        ps_deg = ps_degs[i];
        s_core_index = &k_core_indices[i];
        sp_adj_lst = mg.GetClgToPrime(i)->GetAdjLst();

        if (s_core_index->n < (s_core_index->e << 1)) {  // deletes more
            for (int j = p_core_index->e; j < p_core_index->n; j++) {
                ps_deg[p_core_index->vert[j]] = 0;
            }

            for (int j = s_core_index->e; j < s_core_index->n; j++) {
                v = s_core_index->vert[j];
                for (int l = 1; l <= sp_adj_lst[v][0]; l++) {
                    u = sp_adj_lst[v][l];
                    if (p_core_index->pos[u] >= p_core_index->e) ps_deg[u]++;
                }
            }
        } else {  // remains more
            ps_adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();

            for (int j = p_core_index->e; j < p_core_index->n; j++) {
                v = p_core_index->vert[j];
                ps_deg[v] = ps_adj_lst[v][0];
            }

            for (int j = 0; j < s_core_index->e; j++) {
                v = s_core_index->vert[j];
                for (int l = 1; l <= sp_adj_lst[v][0]; l++) {
                    u = sp_adj_lst[v][l];
                    if (p_core_index->pos[u] >= p_core_index->e) ps_deg[u]--;
                }
            }
        }
    }

    // init sp_degs
    for (int i = 0; i < ln - 1; i++) {
        sp_deg = sp_degs[i];
        s_core_index = &k_core_indices[i];
        ps_adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();

        if (p_core_index->n < (p_core_index->e << 1)) {  // deletes more
            for (int j = s_core_index->e; j < s_core_index->n; j++) {
                sp_deg[s_core_index->vert[j]] = 0;
            }

            for (int j = p_core_index->e; j < p_core_index->n; j++) {
                v = p_core_index->vert[j];
                for (int l = 1; l <= ps_adj_lst[v][0]; l++) {
                    u = ps_adj_lst[v][l];
                    if (s_core_index->pos[u] >= s_core_index->e) {
                        sp_deg[u]++;
                    }
                }
            }
        } else {
            sp_adj_lst = mg.GetClgToPrime(i)->GetAdjLst();

            for (int j = s_core_index->e; j < s_core_index->n; j++) {
                v = s_core_index->vert[j];
                sp_deg[v] = sp_adj_lst[v][0];
            }

            for (int j = 0; j < p_core_index->e; j++) {
                v = p_core_index->vert[j];
                for (int l = 1; l <= ps_adj_lst[v][0]; l++) {
                    u = ps_adj_lst[v][l];
                    if (s_core_index->pos[u] >= s_core_index->e) sp_deg[u]--;
                }
            }
        }

        // Remove secondary layer vertices with no primary layer neighbors.
        for (int j = s_core_index->e; j < s_core_index->n; j++) {
            v = s_core_index->vert[j];
            if (sp_deg[v] == 0) s_core_index->Remove(v);
        }
    }
}

// Generate representations for multi-layer subgraph induced by k_vec[i]-core in each layer i.
void PTreeBuilder::BuildKSubgraph() {
    int v, u, id, *deg, *init_ps_deg, **adj_lst;
    CoreIndex *core_index, *p_core_index;

    // build intra-layer subgraph
    for (int i = 0; i < ln; i++) {
        deg = degs[i];
        core_index = &k_core_indices[i];
        adj_lst = mg.GetGraph(i)->GetAdjLst();

        // Construct sub edge lists.
        for (int j = core_index->e; j < core_index->n; j++) {
            v = core_index->vert[j];
            k_subgraph[i][v] = subgraph_buf.Allocate(deg[v] + 1);

            if (deg[v] == adj_lst[v][0]) {  // no changes in v's adjacent list
                memcpy(k_subgraph[i][v], adj_lst[v], (adj_lst[v][0] + 1) * sizeof(int));
            } else {
                k_subgraph[i][v][0] = deg[v];
                id = 1;
                for (int l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (core_index->pos[u] >= core_index->e) {
                        k_subgraph[i][v][id++] = u;
                    }
                }
            }
        }
    }

    // Build bipartite cross-layer subgraph, from primary layer to secondary layers.
    p_core_index = &k_core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        deg = ps_degs[i];  // ps_deg
        core_index = &k_core_indices[i];
        init_ps_deg = init_ps_degs[i];
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();  // ps_adj_lst

        // Construct sub edge lists.
        for (int j = p_core_index->e; j < p_core_index->n; j++) {
            v = p_core_index->vert[j];
            init_ps_deg[v] = adj_lst[v][0];  // initial ps_deg.

            ps_k_subgraph[i][v] = subgraph_buf.Allocate(deg[v] + 1);

            if (deg[v] == adj_lst[v][0]) {
                memcpy(ps_k_subgraph[i][v], adj_lst[v], (adj_lst[v][0] + 1) * sizeof(int));
            } else {
                ps_k_subgraph[i][v][0] = deg[v];
                id = 1;
                for (int l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (core_index->pos[u] >= core_index->e) {
                        ps_k_subgraph[i][v][id++] = u;
                    }
                }
            }
        }
    }

    // Build bipartite cross-layer subgraph, from secondary layer to primary layers.
    for (int i = 0; i < ln - 1; i++) {
        deg = sp_degs[i];  // sp_deg
        core_index = &k_core_indices[i];
        adj_lst = mg.GetClgToPrime(i)->GetAdjLst();  // sp_adj_lst

        // Construct sub edge lists.
        for (int j = core_index->e; j < core_index->n; j++) {
            v = core_index->vert[j];
            sp_k_subgraph[i][v] = subgraph_buf.Allocate(deg[v] + 1);

            if (deg[v] == adj_lst[v][0]) {
                memcpy(sp_k_subgraph[i][v], adj_lst[v], (adj_lst[v][0] + 1) * sizeof(int));
            } else {
                sp_k_subgraph[i][v][0] = deg[v];
                id = 1;
                for (int l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (p_core_index->pos[u] >= p_core_index->e) {
                        sp_k_subgraph[i][v][id++] = u;
                    }
                }
            }
        }
    }
}

// Init a Frac2IntPri map that maps possible fraction vectors to integer priority vectors.
void PTreeBuilder::BuildFrac2IntPriMap(Frac2IntPri &frac2int) {
    int d, l, max_deg, *init_ps_deg;
    CoreIndex *p_core_index;

    p_core_index = &k_core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        l = 0, max_deg = 0;
        init_ps_deg = init_ps_degs[i];

        for (int j = p_core_index->e; j < p_core_index->n; j++) {
            d = init_ps_deg[p_core_index->vert[j]];
            if (max_deg < d) max_deg = d;
            aux_arr[l++] = d;
        }

        frac2int.Build(i, aux_arr, l, max_deg);
    }
    priority = frac2int.priority;
}

// Init p_value_heaps.
void PTreeBuilder::BuildPValueHeap(Frac2IntPri &frac2int) {
    int v, *ps_deg, *init_ps_deg, **pri;
    CoreIndex *p_core_index;
    IntLinearHeap *p_value_heap;

    p_core_index = &k_core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        ps_deg = ps_degs[i];
        init_ps_deg = init_ps_degs[i];
        pri = frac2int.priority[i];
        p_value_heap = &p_value_heaps[i];

        p_value_heap->SetBin(frac2int.frac_size[i] - 1);
        for (int j = p_core_index->e; j < p_core_index->n; j++) {
            v = p_core_index->vert[j];
            p_value_heap->Insert(v, pri[init_ps_deg[v]][ps_deg[v]]);
        }
    }
}

// Recursively build sub-p_tree given root r.
void PTreeBuilder::BuildSubPTreeNaive(int lp, PNode *r) {
    bool has_child;
    int old_e[ln];
    PNode *child;
    Branch *branch;

    for (int i = 0; i < ln; i++) old_e[i] = k_core_indices[i].e;

    // Build child node in terms of every dimension from lp to ln-2.
    for (int i = ln - 2; i >= lp; i--) {

        branch = p_tree->get_relative_branch(r, i);
        has_child = BuildBranch(i, branch);

        if (has_child) {
            child = p_tree->GetNewPTreeNode(p_vec, i);
            branch->SetChild(child);

            BuildSubPTreeNaive(i, child);
        }

        // restore
        p_vec[i]--;
        Restore(old_e);
    }
}

// Recursively build sub-p_tree given root r using node elimination.
void PTreeBuilder::BuildSubPTreeNe(int lp, PNode *r) {
    bool has_child;
    int old_e[ln], old_p, *rp, *hp, *chp;
    PNode *child;
    Branch *branch;

    for (int i = 0; i < ln; i++) old_e[i] = k_core_indices[i].e;

    rp = PTree::get_p(r);
    hp = p_tree->get_hp(r);

    // Init the maximal fraction vector w.r.t. the subtree rooted by r.
    for (int i = 0; i <= lp; i++) {
        hp[i] = p_value_heaps[i].GetMinValue();
    }

    // Build child node in terms of every dimension from ln-2 to lp.
    for (int i = ln - 2; i >= lp; i--) {

        // record the ith component of the current p_vec
        old_p = p_vec[i];

        // node elimination
        if (i == lp && hp[i] > rp[i]) {
            rp[i] = hp[i];
            p_vec[i] = hp[i];
        }

        branch = p_tree->get_relative_branch(r, i);
        has_child = BuildBranch(i, branch);
        if (has_child) {
            child = p_tree->GetNewPTreeNode(p_vec, i);
            branch->SetChild(child);

            BuildSubPTreeNe(i, child);

            // update the maximal fraction vector w.r.t. ST_r.
            chp = p_tree->get_hp(child);
            for (int j = 0; j <= lp; j++) {
                hp[j] = std::min(hp[j], chp[j]);
            }
        }

        // restoreNaive
        p_vec[i] = old_p;
        Restore(old_e);
    }
}

// Recursively build sub-p_tree given root r using subgraph elimination.
void PTreeBuilder::BuildSubPTreeSe(int lp, PNode *r, PNode *anc, PNode *pr) {
    bool has_child, is_sub_eli;
    int old_e[ln], old_p, *hp, *chp;
    PNode *cpr, *child;
    Branch *branch;

    for (int i = 0; i < ln; i++) old_e[i] = k_core_indices[i].e;

    // Init the maximal fraction vector w.r.t. the subtree rooted by r.
    hp = p_tree->get_hp(r);
    for (int i = 0; i <= lp; i++) {
        hp[i] = p_value_heaps[i].GetMinValue();
    }

    // Build child node in terms of every dimension from ln-2 to lp.
    for (int i = ln - 2; i >= lp; i--) {

        old_p = p_vec[i];

        branch = p_tree->get_relative_branch(r, i);
        has_child = BuildBranch(i, branch);

        if (has_child) {

            // subtree elimination
            if (i == lp) cpr = pr;
            else if (anc) cpr = p_tree->get_relative_branch(anc, i)->child;
            else cpr = nullptr;

            while (cpr && PTree::get_p(cpr)[i] < p_vec[i]) {
                cpr = p_tree->get_relative_branch(cpr, i)->child;
            }
            is_sub_eli = cpr && ArrGe(p_tree->get_hp(cpr), p_vec, i);

            if (is_sub_eli) {
                child = cpr;
                branch->SetChild(child);
            }
            else {
                child = p_tree->GetNewPTreeNode(p_vec, i);
                branch->SetChild(child);

                BuildSubPTreeSe(i, child, r, cpr);
            }

            // update the maximal fraction vector w.r.t. ST_r.
            chp = p_tree->get_hp(child);
            for (int j = 0; j <= lp; j++) {
                hp[j] = std::min(hp[j], chp[j]);
            }
        }

        // restoreNaive
        p_vec[i] = old_p;
        Restore(old_e);
    }
}

// Recursively build sub-p_tree given root r using both node elimination and subgraph elimination.
void PTreeBuilder::BuildSubPTreeNeSe(int lp, PNode *r, PNode *anc, PNode *pr) {
    bool has_child, is_sub_eli;
    int old_e[ln], old_p, *rp, *hp, *chp;
    PNode *cpr, *child;
    Branch *branch;

    for (int i = 0; i < ln; i++) old_e[i] = k_core_indices[i].e;

    rp = PTree::get_p(r);
    hp = p_tree->get_hp(r);

    // Init the maximal fraction vector w.r.t. the subtree rooted by r.
    for (int i = 0; i <= lp; i++) {
        hp[i] = p_value_heaps[i].GetMinValue();
    }

    // Build child node in terms of every dimension from ln-2 to lp.
    for (int i = ln - 2; i >= lp; i--) {

        old_p = p_vec[i];

        // node elimination
        if (i == lp && hp[i] > rp[i]) {
            rp[i] = hp[i];
            p_vec[i] = hp[i];
        }

        branch = p_tree->get_relative_branch(r, i);
        has_child = BuildBranch(i, branch);

        if (has_child) {

            // subtree elimination
            if (i == lp) cpr = pr;
            else if (anc) cpr = p_tree->get_relative_branch(anc, i)->child;
            else cpr = nullptr;

            while (cpr && PTree::get_p(cpr)[i] < p_vec[i]) {
                cpr = p_tree->get_relative_branch(cpr, i)->child;
            }
            is_sub_eli = cpr && ArrGe(p_tree->get_hp(cpr), p_vec, i);

            if (is_sub_eli) {
                child = cpr;
                branch->SetChild(child);
            }
            else {
                child = p_tree->GetNewPTreeNode(p_vec, i);
                branch->SetChild(child);

                BuildSubPTreeNeSe(i, child, r, cpr);
            }

            // update the maximal fraction vector w.r.t. ST_r.
            chp = p_tree->get_hp(child);
            for (int j = 0; j <= lp; j++) {
                hp[j] = std::min(hp[j], chp[j]);
            }
        }

        // restoreNaive
        p_vec[i] = old_p;
        Restore(old_e);
    }
}

// Build child node of branch in terms of dimension i.
bool PTreeBuilder::BuildBranch(int i, Branch *branch) {
    int old_e, num_of_removed_vtx;
    CoreIndex *p_core_index;
    IntLinearHeap *p_value_heap;

    p_core_index = &k_core_indices[pid];
    p_value_heap = &p_value_heaps[i];

    // Try to get vertices with fraction vector p satisfying p[i] = p_vec[i] as the first batch to delete.
    num_of_removed_vtx = p_value_heap->GetElements(p_vec[i], aux_arr);
    p_vec[i] += 1;

    if (num_of_removed_vtx) {  // if such vertices exist
        old_e = p_core_index->e;
        for (int j = 0; j < num_of_removed_vtx; j++) {
            p_core_index->Remove(aux_arr[j]);
        }
        ExtractKPCore();

        // Store difference along the right-most paths.
        if (i == ln - 2) p_tree->SetBranch(branch, p_core_index->vert + old_e, p_core_index->e - old_e);

    } else {
        branch->num_of_vtx = 0;
    }
    return p_core_index->e < p_core_index->n;
}

// Compute the (k,p_vec)-core.
void PTreeBuilder::ExtractKPCore() {
    int old_s;
    CoreIndex *p_core_index;
    IntLinearHeap *p_value_heap;

    p_core_index = &k_core_indices[pid];
    old_s = p_core_index->s;

    while (true) {
        PeelPVtx();
        for (int i = 0; i < ln - 1; i++) {
            if (k_core_indices[i].s != k_core_indices[i].e) {
                PeelSVtx(i);
            }
        }
        if (p_core_index->s == p_core_index->e) break;
    }

    // Remove inactive vertices from p_value_heap.
    for (int i = 0; i < ln - 1; i++) {
        p_value_heap = &p_value_heaps[i];
        for (int j = old_s; j < p_core_index->e; j++) {
            p_value_heap->Remove(p_core_index->vert[j]);
        }
    }
}

// Peel on primary layer, remove affected secondary layer vertices.
void PTreeBuilder::PeelPVtx() {
    int v, u, k, old_s, *deg, *sp_deg, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;

    p_core_index = &k_core_indices[pid];
    auto &s = p_core_index->s;
    auto &e = p_core_index->e;
    old_s = s;

    k = k_vec[pid];
    deg = degs[pid];
    adj_lst = k_subgraph[pid];

    // Peel
    for (; s < e; s++) {
        v = p_core_index->vert[s];
        for (int i = 1; i <= adj_lst[v][0]; i++) {
            u = adj_lst[v][i];
            deg[u]--;
            if (p_core_index->pos[u] >= e && deg[u] < k) {
                p_core_index->Remove(u);
            }
        }
    }

    // Remove affected vertices
    for (int i = 0; i < ln - 1; i++) {
        sp_deg = sp_degs[i];
        s_core_index = &k_core_indices[i];
        adj_lst = ps_k_subgraph[i];

        for (int j = old_s; j < e; j++) {
            v = p_core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                sp_deg[u]--;
                if (s_core_index->pos[u] >= s_core_index->e && sp_deg[u] == 0) {
                    s_core_index->Remove(u);
                }
            }
        }
    }
}

// Peel on secondary layer, remove affected primary layer vertices.
void PTreeBuilder::PeelSVtx(int sid) {
    int v, u, k, p, len, old_s, *deg, *ps_deg, *init_ps_deg, **pri, **adj_lst;
    CoreIndex *p_core_index, *s_core_index;
    IntLinearHeap *p_value_heap;

    s_core_index = &k_core_indices[sid];
    auto &s = s_core_index->s;
    auto &e = s_core_index->e;
    old_s = s;

    k = k_vec[sid];
    deg = degs[sid];
    adj_lst = k_subgraph[sid];

    // Peel
    for (; s < e; s++) {
        v = s_core_index->vert[s];
        for (int i = 1; i <= adj_lst[v][0]; i++) {
            u = adj_lst[v][i];
            deg[u]--;
            if (s_core_index->pos[u] >= e && deg[u] < k) {
                s_core_index->Remove(u);
            }
        }
    }

    // Identify affected primary layer vertices, modify ps_deg and neighbor coverage fraction in p_value_heap.
    p = p_vec[sid];
    pri = priority[sid];
    ps_deg = ps_degs[sid];
    init_ps_deg = init_ps_degs[sid];
    adj_lst = sp_k_subgraph[sid];  // sp_adj_lst
    p_core_index = &k_core_indices[pid];

    len = 0;
    for (int j = old_s; j < e; j++) {
        v = s_core_index->vert[j];
        for (int l = 1; l <= adj_lst[v][0]; l++) {
            u = adj_lst[v][l];
            ps_deg[u]--;

            // Identify vertices with neighbor coverage fraction changed but not removed.
            // Only currently active vertices are considered updated.
            if (p_core_index->pos[u] >= p_core_index->e) {
                if (pri[init_ps_deg[u]][ps_deg[u]] < p) {
                    p_core_index->Remove(u);
                    aux_sign[u] = false;
                } else if (!aux_sign[u]) {
                    aux_sign[u] = true;
                    aux_arr[len++] = u;
                }
            }
        }
    }

    // Update p_value_heap.
    p_value_heap = &p_value_heaps[sid];
    for (int i = 0; i < len; i++) {
        v = aux_arr[i];
        if (aux_sign[v]) {
            p_value_heap->Update(v, pri[init_ps_deg[v]][ps_deg[v]]);
            aux_sign[v] = false;
        }
    }
}

// Restore all auxiliary structures to the state before generating current branch.
void PTreeBuilder::Restore(const int *old_e) {
    int v, u, len, *deg, *sp_deg, *ps_deg, *init_ps_deg, **pri, **adj_lst;
    CoreIndex *core_index, *p_core_index;
    IntLinearHeap *p_value_heap;

    // restore deg
    for (int i = 0; i < ln; i++) {
        deg = degs[i];
        adj_lst = k_subgraph[i];
        core_index = &k_core_indices[i];
        for (int j = old_e[i]; j < core_index->e; j++) {
            v = core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                deg[adj_lst[v][l]]++;
            }
        }
    }

    // restore sp_deg
    p_core_index = &k_core_indices[pid];
    for (int i = 0; i < ln - 1; i++) {
        sp_deg = sp_degs[i];
        adj_lst = ps_k_subgraph[i];  // ps_adj_lst
        for (int j = old_e[pid]; j < p_core_index->e; j++) {
            v = p_core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                sp_deg[adj_lst[v][l]]++;
            }
        }
    }

    // restore ps_deg and p_value_heap
    for (int i = 0; i < ln - 1; i++) {
        ps_deg = ps_degs[i];
        adj_lst = sp_k_subgraph[i];  // sp_adj_lst
        core_index = &k_core_indices[i];

        len = 0;
        for (int j = old_e[i]; j < core_index->e; j++) {
            v = core_index->vert[j];
            for (int l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                ps_deg[u]++;

                // Identify vertices with neighbor coverage fraction changed but not removed.
                if (p_core_index->pos[u] >= p_core_index->e && !aux_sign[u]) {
                    aux_sign[u] = true;
                    aux_arr[len++] = u;
                }
            }
        }

        pri = priority[i];
        init_ps_deg = init_ps_degs[i];
        p_value_heap = &p_value_heaps[i];

        // Update vertices in p_value_heap with neighbor coverage fraction changed.
        for (int j = 0; j < len; j++) {
            v = aux_arr[j];
            p_value_heap->Update(v, pri[init_ps_deg[v]][ps_deg[v]]);
            aux_sign[v] = false;
        }

        // Insert all removed vertices back to the p_value_heap.
        for (int j = old_e[pid]; j < p_core_index->e; j++) {
            v = p_core_index->vert[j];
            p_value_heap->Insert(v, pri[init_ps_deg[v]][ps_deg[v]]);
        }
    }

    // restore k_core_indices
    for (int i = 0; i < ln; i++) {
        k_core_indices[i].e = old_e[i];
        k_core_indices[i].s = old_e[i];
    }
}

void PTreeBuilder::RestorePValueHeap() {
    for (int i = 0; i < ln - 1; i++) {
        p_value_heaps[i].ReleaseBin();
    }
}