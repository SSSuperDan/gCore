//
// Created by ldd on 2021/10/30.
//

#ifndef CSM4GMG_RCD_H
#define CSM4GMG_RCD_H

#include "KPCore.h"

class RCD {
public:
    explicit RCD(MultilayerGraph &mg_);
    ~RCD();

    int Extract(vector<int> &k_vec_, vector<int> &p_vec, int *r_com, int &length);

private:
    MultilayerGraph &mg;
    int* k_vec;
    int* p_vec;

    CoreIndex *core_indices;
    int *num_of_vtx;
    int **degs;
    int **ps_degs;
    int **sp_degs;

    int ln;
    int pid;

    void Init();
    void InitKCore(int gid);
    void InitPartialCDegs();
    void ReduceVtx();

    void Peel(int gid);
    void PeelSVtx(int sid);
    void PeelPVtx();
};


#endif //CSM4GMG_RCD_H
