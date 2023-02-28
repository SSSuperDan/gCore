//
// Created by ldd on 2021/10/4.
//

#ifndef CSM4GMG_KPCORE_H
#define CSM4GMG_KPCORE_H

#include "../Graphs/MultilayerGraph.h"
#include "../Utils/ArrayUtils.h"
#include "CoreIndex.h"

enum kpcore_extractor {
    WOCDEG, WCDEG, MIX, OPT
};


class KPCore {
public:

    explicit KPCore(MultilayerGraph &mg_);
    ~KPCore();

    int Extract(vector<int> &k_vec_, vector<float> &p_vec, int *core, int &length, kpcore_extractor extractor = OPT);
    int Extract(vector<int> &k_vec_, vector<Frac> &p_vec, int *core, int &length, kpcore_extractor extractor = OPT);

private:
    // Parameters
    MultilayerGraph &mg;
    int *k_vec;

    CoreIndex *core_indices;
    int *num_of_vtx;
    int **degs;
    int **demands;
    int **ps_degs;
    int **sp_degs;

    int ln;
    int pid;

    int ExtractImplWoCDeg(int *core, int &length);  // Without maintaining cross-layer degree arrays.
    int ExtractImplWCDeg(int *core, int &length);  // Maintaining cross-layer degree arrays.
    int ExtractImplMix(int *core, int &length);  // Mix of the above two.
    int ExtractImplOpt(int *core, int &length);

    void Init();
    void InitKCore(int gid);
    void InitClDegs();
    void InitPartialCDegs();
    void ReleaseClDegs();
    void ReduceVtx();
    void Peel(int gid);
    void FilterSVtxByPNbrs(int sid);
    void FilterPVtxByFracDemand(int sid);
    void PeelSVtx(int sid);
    void PeelPVtx();
    void PeelOpt(int gid);
    void InitPartialCDegsOpt();
    void PeelSVtxOpt(int sid);
    void PeelPVtxOpt();

};

#endif //CSM4GMG_KPCORE_H
