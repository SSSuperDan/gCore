//
// Created by ldd on 2022/10/14.
//

#ifndef CSM4GMG_KPVECGENERATOR_H
#define CSM4GMG_KPVECGENERATOR_H


#include "Graphs/MultilayerGraph.h"
#include "Core/CoreDec.h"

class KPVecGenerator {
public:
    static int LoadKVectors(int dim, const string &file, vector<vector<int>> &k_vectors);
    static int LoadPVectors(int dim, const string &file, vector<vector<float>> &p_vectors);
    static int LoadKPVectors(int dim, const string &file, vector<vector<int>> &k_vectors, vector<vector<float>> &p_vectors);

    static void GenRandomPVectors(int dim, int size, const string &output);
    static void GenRandomKVectors(MultilayerGraph &mg, int size, const string &output, const vector<int> &end_k = {});
    static void GenRandomKPVectors(MultilayerGraph &mg, int size, const string &output, const vector<int> &end_k = {});
    static void GenRandomKPVectors(MultilayerGraph &mg, const string &sampled_k_vec_file, int size, const string &output);

private:
    static void GetDegeneracy(MultilayerGraph &mg, int *degeneracy);
};


#endif //CSM4GMG_KPVECGENERATOR_H
