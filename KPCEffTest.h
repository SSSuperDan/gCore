//
// Created by ldd on 2022/6/28.
//

#ifndef CSM4GMG_KPCEFFTEST_H
#define CSM4GMG_KPCEFFTEST_H

#include "Core/KPTreeBuilder.h"
#include "Core/ExtendedDCC.h"
#include "Core/CoreDec.h"
#include "Core/RCD.h"

#include "Utils/FileUtils.h"
#include "Utils/Timer.h"

#include "KPVecGenerator.h"
#include "KPTreeUtils.h"

struct TS {
    double time;
    int visits;

    TS() = default;
    TS& operator = (const TS &ts_) = default;
};

struct P_runtime {
    double k_core_runtime{0};
    double dCC_runtime{0};
    double kp_core_runtime{0};

    [[nodiscard]] string to_string() const {
        return "k_core_runtime = " + std::to_string(k_core_runtime) +
        "\n dCC_runtime = " + std::to_string(dCC_runtime) +
        "\n kp_core_runtime = " + std::to_string(kp_core_runtime);
    }
};

struct Runtime {
    double k_core_runtime{0};
    double dCC_runtime{0};
    double rcd_runtime{0};
    double kp_core_runtime{0};
    TS kp_core_naive_search_ts;
    TS kp_core_ne_search_ts;
    TS kp_core_se_search_ts;
    TS kp_core_ne_se_search_ts;

    [[nodiscard]] string to_string() const {
        return "k_core_runtime = " + std::to_string(k_core_runtime) +
        "\n dCC_runtime = " + std::to_string(dCC_runtime) +
        "\n rcd_runtime = " + std::to_string(rcd_runtime) +
        "\n kp_core_runtime = " + std::to_string(kp_core_runtime) +
        "\n kp_core_naive_search_runtime = " + std::to_string(kp_core_naive_search_ts.time) + ", visits = " + std::to_string(kp_core_naive_search_ts.visits) +
        "\n kp_core_ne_search_runtime = " + std::to_string(kp_core_ne_search_ts.time) + ", visits = " + std::to_string(kp_core_ne_search_ts.visits) +
        "\n kp_core_se_search_runtime = " + std::to_string(kp_core_se_search_ts.time) + ", visits = " + std::to_string(kp_core_se_search_ts.visits) +
        "\n kp_core_ne_se_search_runtime = " + std::to_string(kp_core_ne_se_search_ts.time) + ", visits = " + std::to_string(kp_core_ne_se_search_ts.visits);
    }
};

class KPCEffTest {

public:

    static void P_GCS(MultilayerGraph &mg, vector<int> &k, const string &path = "");
    static void P_GCS(MultilayerGraph &mg, const string &sampled_k_vec_file, const string &path = "");

    static void GCS(MultilayerGraph &mg, vector<KP_tree_f2i> &kp_tree_f2i, vector<int> &k, vector<float> &p, const string &path = "");
    static void GCS(MultilayerGraph &mg, vector<KP_tree_f2i> &kp_tree_f2i, const string &sampled_kp_file, const string &path = "");
    static void GCS(MultilayerGraph &mg, vector<P_tree_f2i> &p_tree_f2i, vector<int> &k, vector<float> &p, const string &path = "");

    static void GCIStatistics(P_tree_f2i &p_tree_f2i, const string &path = "");
    static void GCIStatistics(KP_tree_f2i &kp_tree_f2i, const string &path = "");

private:

    static string DumpPTreeInfo(string &p_tree_name, long long number_of_nodes, long double p_tree_mem, long double f2i_mem);
    static string DumpKPTreeInfo(string &kp_tree_name, long long number_of_nodes, long double kp_tree_mem,
                                 long double f2i_mem);

    static double kCore(MultilayerGraph &mg, vector<int> &k);
    static double kCore(MultilayerGraph &mg, vector<vector<int>> &ks);

    static double dCC(MultilayerGraph &mg, vector<int> &k);
    static double dCC(MultilayerGraph &mg, vector<vector<int>> &ks);

    static double rcd(MultilayerGraph &mg, vector<int> &k, vector<int> &ck);
    static double rcd(MultilayerGraph &mg, vector<vector<int>> &ks, vector<int> &ck);

    static double kpCore(MultilayerGraph &mg, vector<int> &k, vector<float> &p);
    static double kpCore(MultilayerGraph &mg, vector<vector<int>> &ks, vector<float> &p);
    static double kpCore(MultilayerGraph &mg, vector<int> &k, vector<vector<float>> &ps);
    static double kpCore(MultilayerGraph &mg, vector<vector<int>> &ks, vector<vector<float>> &ps);

    static TS kpCoreSearch(MultilayerGraph &mg, PTree &p_tree, Frac2IntPri &frac2int, vector<float> &p);
    static TS kpCoreSearch(MultilayerGraph &mg, KPTree &kp_tree, Frac2IntPri &frac2int, vector<int> &k,vector<float> &p);
    static TS kpCoreSearch(MultilayerGraph &mg, KPTree &kp_tree, Frac2IntPri &frac2int, vector<vector<int>> &ks,
                           vector<vector<float>> &ps);

};


#endif //CSM4GMG_KPCEFFTEST_H
