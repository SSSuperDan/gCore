//
// Created by ldd on 2022/6/23.
//

#ifndef CSM4GMG_EXTENDEDDCC_H
#define CSM4GMG_EXTENDEDDCC_H

#include "../Utils/CollectionUtils.h"
#include "../Graphs/MultilayerGraph.h"
#include "KPCore.h"

enum dcc_impl {Q, CI};

class ExtendedDCC {
public:
    explicit ExtendedDCC(MultilayerGraph &mg_);
    ~ExtendedDCC();

    void Extract(vector<int> &k_vec_, int *dcc, int &length, dcc_impl impl = CI);
    void ExtractImplQueue(int *dcc, int &length);
    void ExtractImplCoreIndex(int *dcc, int &length);

private:
    MultilayerGraph &mg;
    int * k_vec;

    int **degs;
    int ln;
    int n;

};


#endif //CSM4GMG_EXTENDEDDCC_H
