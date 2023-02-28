//
// Created by ldd on 2021/10/5.
//

#ifndef CSM4GMG_PTREEBUILDER_H
#define CSM4GMG_PTREEBUILDER_H

#include "../Structures/IntLinearHeap.h"
#include "../Graphs/MultilayerGraph.h"

#include "Frac2IntPri.h"
#include "CoreIndex.h"
#include "PTree.h"

enum p_tree_builder {
    NAIVE, NE, SE, NESE
};

class PTreeBuilder {
public:
    explicit PTreeBuilder(MultilayerGraph &mg_);
    ~PTreeBuilder();
    void Execute(vector<int> &k_vec_, PTree & p_tree_, Frac2IntPri & frac2int, p_tree_builder opt = NESE);

protected:
    MultilayerGraph &mg;
    int *k_vec;
    int *p_vec;

    PTree *p_tree{};
    int *** priority{};

    CoreIndex *k_core_indices;
    IntLinearHeap *p_value_heaps;
    DataBuf<int> subgraph_buf;

    int ***ps_k_subgraph;
    int ***sp_k_subgraph;
    int ***k_subgraph;

    int **init_ps_degs;
    int **ps_degs;
    int **sp_degs;
    int **degs;

    int *num_of_vtx;
    int ln;
    int pid;

    int *aux_arr;
    bool *aux_sign;

    // Init
    void InitKPCore();
    void InitPartialCDegs();
    void Peel(int gid);
    void PeelS(int sid);
    void BuildKSubgraph();
    void BuildFrac2IntPriMap(Frac2IntPri & frac2int);
    void BuildPValueHeap(Frac2IntPri &frac2int);

    // Recursively build.
    void BuildSubPTreeNaive(int lp, PNode *r);
    void BuildSubPTreeNe(int lp, PNode *r);
    void BuildSubPTreeSe(int lp, PNode *r, PNode *anc, PNode *pr);
    void BuildSubPTreeNeSe(int lp, PNode *r, PNode *anc, PNode *pr);
    void ExtractKPCore();
    void PeelPVtx();
    void PeelSVtx(int sid);
    bool BuildBranch(int i, Branch *branch);
    void Restore(const int *old_e);
    void RestorePValueHeap();

};


#endif //CSM4GMG_PTREEBUILDER_H
