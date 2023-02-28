//
// Created by ldd on 2022/7/9.
//

#ifndef CSM4GMG_KPCEVALTEST_H
#define CSM4GMG_KPCEVALTEST_H

#include "Core/KPTreeBuilder.h"
#include "Core/KPCore.h"
#include "Core/RCD.h"
#include "Core/CoreDec.h"

#include "Utils/FileUtils.h"

class KPCEvalTest {
public:
    // introduced metrics to measure the closeness of a vertex and other vertices reflected by indirect connections
    static void ComputeKValue(MultilayerGraph &mg, vector<int> &k, vector<int> &ck, vector<float> &p, const string &path);
    static void ComputePValue(MultilayerGraph &mg, vector<int> &k, vector<int> &ck, vector<float> &p, const string &path);

    // size of (k,p)-core when varying k and p
    static void ComputeSizeMatrix(MultilayerGraph &mg, int pk, float step, const string &path);
    static void ComputeSizeMatrixFixK(MultilayerGraph &mg, vector<int> &k, const string &path);

private:
    static void GetCrossLayerNeighbor(MultilayerGraph &mg, const int *core, int length, int layer, int *neighbors, int &num_of_neighbors);

    static double GetPValue(const int *adj_lst, const bool *in_s_core);
    static int GetKValue(const int *adj_lst, const int *nbr_pos, const int *core_number, float p);

    static string DumpTestSettings(vector<int> &k, vector<int> &ck, vector<float> &p, int layer);


    static inline int GetPCeilingPos(int n, Frac &p) {
        return ceil(((float) (n * p.num)) / (float) p.den);
    }

    static inline int GetPCeilingPos(int n, float p) {
        return ceil((float) n * p);
    }
};

#endif //CSM4GMG_KPCEVALTEST_H
