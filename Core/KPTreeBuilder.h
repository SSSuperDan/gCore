//
// Created by ldd on 2022/6/11.
//

#ifndef CSM4GMG_KPTREEBUILDER_H
#define CSM4GMG_KPTREEBUILDER_H

#include "PTreeBuilder.h"
#include "KPTree.h"

class KPTreeBuilder : PTreeBuilder {
public:
    explicit KPTreeBuilder(MultilayerGraph &mg_);
    ~KPTreeBuilder() = default;
    void Execute(KPTree &kp_tree_, Frac2IntPri &frac2int, p_tree_builder opt = NESE, int step_ = 1, const vector<int> &start_k_vec = {});
    void Execute(KPTree &kp_tree_, Frac2IntPri &frac2int, const vector<vector<int>> &k_vectors, p_tree_builder opt = NESE);

private:
    KPTree *kp_tree{};
    void (KPTreeBuilder::*BuildSubPTree)(int, PNode *){};

    int step;
    int *degeneracy;

    void InitIncKPCore();
    void InitIncKPCore(int gid);
    void BuildSubKPTree(int lk);
    void ShrinkSubgraph(const int *old_e);
    void ShrinkSubgraph(int gid, const int *old_e);
    void RestoreSubgraph(const int *old_e);
    void RestoreSubgraph(int gid, const int *old_e);
    void ShrinkLayer(CoreIndex *s_core_index, CoreIndex *t_core_index, int s_old_e, int **st_adj_lst, int **ts_adj_lst);
    static void RestoreLayer(CoreIndex *s_core_index, CoreIndex *t_core_index, int s_old_e, int **st_adj_lst, int **ts_adj_lst);

    inline void BuildSubPTreeSe(int lp, PNode *r) {
        PTreeBuilder::BuildSubPTreeSe(lp, r, nullptr, nullptr);
    }

    inline void BuildSubPTreeNeSe(int lp, PNode *r) {
        PTreeBuilder::BuildSubPTreeNeSe(lp, r, nullptr, nullptr);
    }
};


#endif //CSM4GMG_KPTREEBUILDER_H
