//
// Created by ldd on 2022/6/28.
//

#include "KPCEffTest.h"


void KPCEffTest::P_GCS(MultilayerGraph &mg, vector<int> &k, const string &path) {
    string filename, runtime_info;
    vector<float> p(mg.GetLayerNumber() - 1, 1.0);

    P_runtime runtime;
    runtime.k_core_runtime = kCore(mg, k);
    runtime.dCC_runtime = dCC(mg, k);
    runtime.kp_core_runtime = kpCore(mg, k, p);

    runtime_info = runtime.to_string();

    if (!path.empty()) {
        filename = GetNextValidFile(path + mg.GetGraphName() + "_pillar_gcs_results.txt");
        auto f = std::ofstream(filename);
        f << "k = " << Vec2String(k) << endl;
        f << runtime_info << endl;
        f.close();
    }
    cout << runtime_info << endl;
}

void KPCEffTest::P_GCS(MultilayerGraph &mg, const string &sampled_k_vec_file, const string &path) {
    string filename, runtime_info;
    vector<float> p(mg.GetLayerNumber() - 1, 1.0);

    vector<vector<int>> sampled_k_vec;
    KPVecGenerator::LoadKVectors(mg.GetLayerNumber(), sampled_k_vec_file, sampled_k_vec);

    P_runtime runtime;
    runtime.k_core_runtime = kCore(mg, sampled_k_vec);
    runtime.dCC_runtime = dCC(mg, sampled_k_vec);
    runtime.kp_core_runtime = kpCore(mg, sampled_k_vec, p);

    runtime_info = runtime.to_string();

    if (!path.empty()) {
        filename = GetNextValidFile(path + mg.GetGraphName() + "_pillar_gcs_results.txt");
        auto f = std::ofstream(filename);
        f << "sampled_k_file = " << sampled_k_vec_file << endl;
        f << runtime_info << endl;
        f.close();
    }
    cout << runtime_info << endl;
}


void
KPCEffTest::GCS(MultilayerGraph &mg, vector<KP_tree_f2i> &kp_tree_f2i, vector<int> &k, vector<float> &p, const string &path) {

    string filename, runtime_info;

    Runtime runtime;
    runtime.k_core_runtime = kCore(mg, k);
    if (!mg.IsGeneral()) runtime.dCC_runtime = dCC(mg, k);
    else {
        vector<int> ck(mg.GetLayerNumber() - 1, 1);
        runtime.rcd_runtime = rcd(mg, k, ck);
    }
    runtime.kp_core_runtime = kpCore(mg, k, p);

    for (auto &ktf:kp_tree_f2i) {
        if (ktf.builder == NAIVE) {
            runtime.kp_core_naive_search_ts = kpCoreSearch(mg, *ktf.kp_tree, *ktf.f2i, k, p);
        } else if (ktf.builder == NE) {
            runtime.kp_core_ne_search_ts = kpCoreSearch(mg, *ktf.kp_tree,*ktf.f2i, k, p);
        } else if (ktf.builder == SE) {
            runtime.kp_core_se_search_ts = kpCoreSearch(mg, *ktf.kp_tree, *ktf.f2i, k, p);
        } else {
            runtime.kp_core_ne_se_search_ts = kpCoreSearch(mg, *ktf.kp_tree, *ktf.f2i, k, p);
        }
    }

    runtime_info = runtime.to_string();

    if (!path.empty()) {
        filename = GetNextValidFile(path + mg.GetGraphName() + "_gcs_results.txt");
        auto f = std::ofstream(filename);
        f << "k = " << Vec2String(k) << ", p = " << Vec2String(p) << endl;
        f << runtime_info << endl;
        f.close();
    }
    cout << runtime_info << endl;
}


void KPCEffTest::GCS(MultilayerGraph &mg, vector<KP_tree_f2i> &kp_tree_f2i, const string &sampled_kp_file,
                     const string &path) {
    string filename, runtime_info;

    vector<vector<int>> sampled_k_vec;
    vector<vector<float>> sampled_p_vec;
    KPVecGenerator::LoadKPVectors(mg.GetLayerNumber(), sampled_kp_file, sampled_k_vec, sampled_p_vec);

    Runtime runtime;
    runtime.k_core_runtime = kCore(mg, sampled_k_vec);
    if (!mg.IsGeneral()) runtime.dCC_runtime = dCC(mg, sampled_k_vec);
    else {
        vector<int> ck(mg.GetLayerNumber() - 1, 1);
        runtime.rcd_runtime = rcd(mg, sampled_k_vec, ck);
    }
    runtime.kp_core_runtime = kpCore(mg, sampled_k_vec, sampled_p_vec);

    for (auto &ktf:kp_tree_f2i) {
        if (ktf.builder == NAIVE) {
            runtime.kp_core_naive_search_ts = kpCoreSearch(mg, *ktf.kp_tree, *ktf.f2i, sampled_k_vec, sampled_p_vec);
        } else if (ktf.builder == NE) {
            runtime.kp_core_ne_search_ts = kpCoreSearch(mg, *ktf.kp_tree,*ktf.f2i, sampled_k_vec, sampled_p_vec);
        } else if (ktf.builder == SE) {
            runtime.kp_core_se_search_ts = kpCoreSearch(mg, *ktf.kp_tree, *ktf.f2i, sampled_k_vec, sampled_p_vec);
        } else {
            runtime.kp_core_ne_se_search_ts = kpCoreSearch(mg, *ktf.kp_tree, *ktf.f2i, sampled_k_vec, sampled_p_vec);
        }
    }

    runtime_info = runtime.to_string();

    if (!path.empty()) {
        filename = GetNextValidFile(path + mg.GetGraphName() + "_gcs_results.txt");
        auto f = std::ofstream(filename);
        f << "sampled_kp_file = " << sampled_kp_file << endl;
        f << runtime_info << endl;
        f.close();
    }
    cout << runtime_info << endl;
}

void KPCEffTest::GCS(MultilayerGraph &mg, vector<P_tree_f2i> &p_tree_f2i, vector<int> &k, vector<float> &p, const string &path) {

    string filename, runtime_info;

    Runtime runtime;
    runtime.k_core_runtime = kCore(mg, k);
    if (!mg.IsGeneral()) runtime.dCC_runtime = dCC(mg, k);
    else {
        vector<int> ck(mg.GetLayerNumber() - 1, 1);
        runtime.rcd_runtime = rcd(mg, k, ck);
    }
    runtime.kp_core_runtime = kpCore(mg, k, p);

    for (auto &ptf:p_tree_f2i) {
        if (ptf.builder == NAIVE) {
            runtime.kp_core_naive_search_ts = kpCoreSearch(mg, *ptf.p_tree, *ptf.f2i, p);
        } else if (ptf.builder == NE) {
            runtime.kp_core_ne_search_ts = kpCoreSearch(mg, *ptf.p_tree,*ptf.f2i, p);
        } else if (ptf.builder == SE) {
            runtime.kp_core_se_search_ts = kpCoreSearch(mg, *ptf.p_tree, *ptf.f2i, p);
        } else {
            runtime.kp_core_ne_se_search_ts = kpCoreSearch(mg, *ptf.p_tree, *ptf.f2i, p);
        }
    }

    runtime_info = runtime.to_string();

    if (!path.empty()) {
        filename = GetNextValidFile(path + mg.GetGraphName() + "_gcs_results.txt");
        auto f = std::ofstream(filename);
        f << "k = " << Vec2String(k) << ", p = " << Vec2String(p) << endl;
        f << runtime_info << endl;
        f.close();
    }
    cout << runtime_info << endl;

}

void KPCEffTest::GCIStatistics(P_tree_f2i &p_tree_f2i, const string &path) {

#ifdef USE_MEMTRACK
    long long number_of_nodes;
    long double p_tree_mem, f2i_mem, mem_usage1, mem_usage2;

    number_of_nodes = p_tree_f2i.p_tree->GetNumOfNodes();

    mem_usage1 = GetTotalMemoryUsageInMB();
    p_tree_f2i.ReleasePTree();
    mem_usage2 = GetTotalMemoryUsageInMB();

    p_tree_mem = mem_usage1 - mem_usage2;

    p_tree_f2i.ReleaseF2i();
    mem_usage1 = GetTotalMemoryUsageInMB();
    f2i_mem = mem_usage2 - mem_usage1;

    string p_tree_name = KPTreeUtils::GetPTreeFilePref(p_tree_f2i) + "_p_tree";
    auto runtime_info = DumpPTreeInfo(p_tree_name, number_of_nodes, p_tree_mem, f2i_mem);

    if (!path.empty()) {
        auto f = std::ofstream(path + p_tree_f2i.dataset + "_tree_mem_info.txt", std::ios::app);
        f << runtime_info << endl;
        f.close();
    }
    cout << runtime_info << endl;

#endif
}

void KPCEffTest::GCIStatistics(KP_tree_f2i &kp_tree_f2i, const string &path) {
#ifdef USE_MEMTRACK

    long long number_of_nodes;
    long double kp_tree_mem, f2i_mem, mem_usage1, mem_usage2;

    number_of_nodes = kp_tree_f2i.kp_tree->GetNumOfPNodes();

    mem_usage1 = GetTotalMemoryUsageInMB();
    kp_tree_f2i.ReleaseKPTree();
    mem_usage2 = GetTotalMemoryUsageInMB();

    kp_tree_mem = mem_usage1 - mem_usage2;

    kp_tree_f2i.ReleaseF2i();
    mem_usage1 = GetTotalMemoryUsageInMB();
    f2i_mem = mem_usage2 - mem_usage1;

    // output
    string kp_tree_name = KPTreeUtils::GetKPTreeFilePref(kp_tree_f2i) + "_kp_tree";
    auto runtime_info = DumpKPTreeInfo(kp_tree_name, number_of_nodes, kp_tree_mem, f2i_mem);

    if (!path.empty()) {
        auto f = std::ofstream(path + kp_tree_f2i.dataset + "_tree_mem_info.txt", std::ios::app);
        f << runtime_info << endl;
        f.close();
    }
    cout << runtime_info << endl;
#endif
}


string KPCEffTest::DumpPTreeInfo(string &p_tree_name, long long number_of_nodes, long double p_tree_mem,
                                 long double f2i_mem) {
    return "p_tree = " + p_tree_name + ", #nodes = " + std::to_string(number_of_nodes) + ", p_tree_mem = " +
           std::to_string(p_tree_mem) + "MB, f2i_mem = " + std::to_string(f2i_mem) + "MB";
}

string KPCEffTest::DumpKPTreeInfo(string &kp_tree_name, long long number_of_nodes, long double kp_tree_mem,
                                  long double f2i_mem) {
    return "kp_tree = " + kp_tree_name + ", #nodes = " + std::to_string(number_of_nodes) + ", kp_tree_mem = " +
           std::to_string(kp_tree_mem) + "MB, f2i_mem = " + std::to_string(f2i_mem) + "MB";
}

double KPCEffTest::kCore(MultilayerGraph &mg, vector<int> &k) {
    UndirectedGraph *g = mg.GetPrimaryGraph();
    int pid = mg.GetPrimaryGraphId(), n = g->GetN(), core[n], length;

    Timer timer;
    length = CoreDec::KCore(g, k[pid], core);
    timer.Stop();
    return timer.GetTimeInMs();
}

double KPCEffTest::kCore(MultilayerGraph &mg, vector<vector<int>> &ks) {
    UndirectedGraph *g = mg.GetPrimaryGraph();
    int pid = mg.GetPrimaryGraphId(), n = g->GetN(), core[n], length;

    Timer timer;
    for (auto &k:ks) {
        length = CoreDec::KCore(g, k[pid], core);
    }
    timer.Stop();
    return timer.GetTimeInMs();
}

double KPCEffTest::dCC(MultilayerGraph &mg, vector<int> &k) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length;

    Timer timer;

    ExtendedDCC dcc_extractor(mg);
    dcc_extractor.Extract(k, core, length);
    timer.Stop();
    return timer.GetTimeInMs();
}

double KPCEffTest::dCC(MultilayerGraph &mg, vector<vector<int>> &ks) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length;

    Timer timer;

    ExtendedDCC dcc_extractor(mg);
    for (auto &k : ks) {
        dcc_extractor.Extract(k, core, length);
    }

    timer.Stop();
    return timer.GetTimeInMs();
}

double KPCEffTest::rcd(MultilayerGraph &mg, vector<int> &k, vector<int> &ck) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length;

    Timer timer;

    RCD rc_extractor(mg);
    rc_extractor.Extract(k, ck, core, length);
    timer.Stop();

    return timer.GetTimeInMs();
}

double KPCEffTest::rcd(MultilayerGraph &mg, vector<vector<int>> &ks, vector<int> &ck) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length;

    Timer timer;

    RCD rc_extractor(mg);
    for (auto &k : ks) {
        rc_extractor.Extract(k, ck, core, length);
    }
    timer.Stop();

    return timer.GetTimeInMs();
}

double KPCEffTest::kpCore(MultilayerGraph &mg, vector<int> &k, vector<float> &p) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length;

    Timer timer;
    KPCore kp_core_extractor(mg);
    kp_core_extractor.Extract(k, p, core, length);

    timer.Stop();
    return timer.GetTimeInMs();
}

double KPCEffTest::kpCore(MultilayerGraph &mg, vector<vector<int>> &ks, vector<float> &p) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length;

    Timer timer;

    KPCore kp_core_extractor(mg);
    for (auto &k : ks) {
        kp_core_extractor.Extract(k, p, core, length);
    }
    timer.Stop();

    return timer.GetTimeInMs();
}

double KPCEffTest::kpCore(MultilayerGraph &mg, vector<int> &k, vector<vector<float>> &ps) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length;

    Timer timer;

    KPCore kp_core_extractor(mg);
    for (auto &p : ps) {
        kp_core_extractor.Extract(k, p, core, length);
    }
    timer.Stop();

    return timer.GetTimeInMs();
}

double KPCEffTest::kpCore(MultilayerGraph &mg, vector<vector<int>> &ks, vector<vector<float>> &ps) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length, id;

    Timer timer;
    KPCore kp_core_extractor(mg);

    id = 0;
    for (auto &k : ks) {
        kp_core_extractor.Extract(k, ps[id], core, length);
        id++;
    }

    timer.Stop();
    return timer.GetTimeInMs();
}

TS KPCEffTest::kpCoreSearch(MultilayerGraph &mg, PTree &p_tree, Frac2IntPri &frac2int, vector<float> &p) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length, num_of_visits;

    Timer timer;

    vector<int> int_p(mg.GetLayerNumber() - 1);
    frac2int.Convert(p, int_p);
    num_of_visits = p_tree.Search(int_p, core, length);

    timer.Stop();

    return {timer.GetTimeInMs(), num_of_visits};
}

TS KPCEffTest::kpCoreSearch(MultilayerGraph &mg, KPTree &kp_tree, Frac2IntPri &frac2int, vector<int> &k,
                            vector<float> &p) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length, num_of_visits;

    Timer timer;

    vector<int> int_p(mg.GetLayerNumber() - 1);
    frac2int.Convert(p, int_p);
    num_of_visits = kp_tree.Search(k, int_p, core, length);

    timer.Stop();
    return {timer.GetTimeInMs(), num_of_visits};
}

TS KPCEffTest::kpCoreSearch(MultilayerGraph &mg, KPTree &kp_tree, Frac2IntPri &frac2int, vector<vector<int>> &ks,
                            vector<vector<float>> &ps) {
    int n = mg.GetPrimaryGraph()->GetN(), core[n], length, id, num_of_visits = 0;

    Timer timer;

    vector<int> int_p(mg.GetLayerNumber() - 1);
    id = 0;
    for (auto &k : ks) {
        frac2int.Convert(ps[id], int_p);
        num_of_visits += kp_tree.Search(k, int_p, core, length);
        id++;
    }

    timer.Stop();
    return {timer.GetTimeInMs(), num_of_visits};
}