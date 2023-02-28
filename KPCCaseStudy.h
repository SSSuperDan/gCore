//
// Created by ldd on 2022/8/24.
//

#ifndef CSM4GMG_KPCCASESTUDY_H
#define CSM4GMG_KPCCASESTUDY_H

#include "Core/KPTreeBuilder.h"
#include "Core/CoreDec.h"
#include "Core/KPCore.h"

#include "Utils/FileUtils.h"

#define DEFAULT_TOP_FEATURE_CNT 5
#define MAXIMUM_LAYER_NUMBER 10

// default Settings
struct Settings {
    const char *dataset;
    int k[MAXIMUM_LAYER_NUMBER];
    float p[MAXIMUM_LAYER_NUMBER-1];
};

static Settings default_settings[] = {
        {"twitter", {5,  10, 5}, {0.5, 0.5}},
        {"dblp", {10, 10}, {-1}},

        {"", {}, {}}
};

class KPCCaseStudy {
public:
    static void TwitterCaseAnalyze(const string &output, const vector<int> & k = {}, const vector<float> & p = {});
    static void DBLPCaseAnalyze(const string &output, const vector<int> &k = {}, const vector<float> &p = {});

private:
    static string ComputeCore(MultilayerGraph & mg, const vector<int> &k_, const vector<float> &p_, vector<vector<int>> & kp_core_ccs, vector<vector<int>> & k_core_ccs);
    static void Compare(MultilayerGraph &mg, vector<vector<int>> &kp_core_ccs, vector<vector<int>> &k_core_ccs, unordered_map<long long, string> *vtx2name, const string &file_pref);

    static void PrintDegeneracy(MultilayerGraph &mg);
    static bool GetDefaultParameters(const string &dataset, int ln, vector<int> &k, vector<float> &p);

    static string GetFilePref(const string &dataset, const vector<int> &k, const vector<float> &p);
    static void GetCrossLayerNeighbor(BipartiteGraph *bg, const vector<int> &cc, int *neighbors, int &num_of_neighbors, int * cnt);
    static void MatchCC(vector<vector<int>> &ccs, vector<vector<int>> &k_core_ccs, int *mapped_cc_index);
    static void DumpCore(int id, const int *cc, int length, long long *id2vtx, std::unordered_map<long long, string> vtx2name, std::basic_ofstream<char> &out);
    static void DumpCoreCNbr(int id, const int *nbr, int length, long long *id2vtx, std::unordered_map<long long, string> vtx2name, int *cnt, std::basic_ofstream<char> &out);

};


#endif //CSM4GMG_KPCCASESTUDY_H