//
// Created by ldd on 2022/10/15.
//

#include "KPCCaseStudy.h"

void KPCCaseStudy::TwitterCaseAnalyze(const string &output, const vector<int> &k, const vector<float> &p) {

    string graph_path = "../dataset/twitter/";
    string dataset = "twitter";
    string setting;
    int ln;

    MultilayerGraph mg;
    mg.LoadGraph(graph_path, dataset);
    ln = mg.GetLayerNumber();

    vector<vector<int>> kp_core_ccs, k_core_ccs;
    setting = ComputeCore(mg, k, p, kp_core_ccs, k_core_ccs);

    unordered_map<long long, string> vtx2name[ln];
    LoadColumn(graph_path + "hashtag_labels.txt", 2, vtx2name[0], ' ');
    LoadColumn(graph_path + "keyword_labels.txt", 2, vtx2name[1], ' ');
    LoadColumn(graph_path + "mention_labels.txt", 2, vtx2name[2], ' ');
    Compare(mg, kp_core_ccs, k_core_ccs, vtx2name, output + setting);
}

void KPCCaseStudy::DBLPCaseAnalyze(const string &output, const vector<int> &k,
                                   const vector<float> &p) {

    string graph_path = "../dataset/dblp/";
    string dataset = "dblp";
    string setting;
    int ln;

    MultilayerGraph mg;
    mg.LoadGraph(graph_path, dataset);
    ln = mg.GetLayerNumber();

    vector<vector<int>> kp_core_ccs, k_core_ccs;
    setting = ComputeCore(mg, k, p, kp_core_ccs, k_core_ccs);

    unordered_map<long long, string> vtx2name[ln];
    LoadColumn(graph_path + "term_name.txt", 2, vtx2name[0], '\t');
    LoadColumn(graph_path + "author_name.txt", 2, vtx2name[1], '\t');
    Compare(mg, kp_core_ccs, k_core_ccs, vtx2name, output + setting);
}

string KPCCaseStudy::GetFilePref(const string &dataset, const vector<int> &k, const vector<float> &p) {
    return dataset + "_" + Vec2String(k) + "_" + Vec2String(p);
}

string KPCCaseStudy::ComputeCore(MultilayerGraph &mg, const vector<int> &k_, const vector<float> &p_,
                                 vector<vector<int>> &kp_core_ccs, vector<vector<int>> &k_core_ccs) {

    UndirectedGraph *pg = mg.GetPrimaryGraph();
    int ln = mg.GetLayerNumber(), n = pg->GetN(), core[n], length;
    vector<int> k;
    vector<float> p;

#ifdef MY_DEBUG
    PrintDegeneracy(mg);
#endif

    if (k_.empty() || p_.empty()) {
        if (GetDefaultParameters(mg.GetGraphName(), ln, k, p)) {
            if (p[0] == -1) {
                PTreeBuilder p_tree_builder(mg);
                PTree p_tree(ln - 1);
                Frac2IntPri f2i(ln - 1);

                p_tree_builder.Execute(k, p_tree, f2i);
                auto pb = p_tree.GetBoundary(0);
                p[0] = (float) f2i.fractions[0][pb].num / (float) f2i.fractions[0][pb].den;

            }
        } else {
            cerr << "Coreness vector (k) and fraction vector (p) must be provided." << endl;
            exit(-1);
        }
    } else {
        k.assign(k_.begin(), k_.end());
        p.assign(p_.begin(), p_.end());
    }

#ifdef MY_DEBUG
    cout << "Settings: k = " << Vec2String(k) << ", p = " << Vec2String(p) << "." << endl;
#endif

    KPCore kp_core_extractor(mg);
    kp_core_extractor.Extract(k, p, core, length);
    pg->CC(core, length, kp_core_ccs);

#ifdef MY_DEBUG
    cout << "Size of the (k,p)-core : " << length << ", with " << kp_core_ccs.size() << " connected components."
         << endl;
    cout << "Size of connected components = ";
    for (auto &cc:kp_core_ccs) cout << cc.size() << " ";
    cout << endl;
#endif

    length = CoreDec::KCore(pg, k[ln - 1], core);
    pg->CC(core, length, k_core_ccs);

#ifdef MY_DEBUG
    cout << "Size of the k_{" << ln - 1 << "}-core : " << length << ", with " << k_core_ccs.size()
         << " connected components." << endl;
    cout << "Size of connected components = ";
    for (auto &cc:k_core_ccs) cout << cc.size() << " ";
    cout << endl;
#endif

    return GetFilePref(mg.GetGraphName(), k, p);
}

bool KPCCaseStudy::GetDefaultParameters(const string &dataset, int ln, vector<int> &k, vector<float> &p) {

    for (int i = 0; default_settings[i].dataset; i++) {
        if (!strcmp(dataset.c_str(), default_settings[i].dataset)) {
            k.assign(default_settings[i].k, default_settings[i].k + ln);
            p.assign(default_settings[i].p, default_settings[i].p + ln - 1);
            return true;
        }
    }
    return false;
}

void KPCCaseStudy::PrintDegeneracy(MultilayerGraph &mg) {
    int ln = mg.GetLayerNumber();

    for (int i = 0; i < ln; i++) {
        cout << "Degeneracy of layer " << i << " = " << CoreDec::GetDegeneracy(mg.GetGraph(i)) << endl;
    }
}

void KPCCaseStudy::Compare(MultilayerGraph &mg, vector<vector<int>> &kp_core_ccs, vector<vector<int>> &k_core_ccs,
                           unordered_map<long long, string> *vtx2name, const string &file_pref) {

    int id, ln = mg.GetLayerNumber(), num_of_vtx[ln], s_max_n = 0;
    for (int j = 0; j < ln; j++) num_of_vtx[j] = mg.GetGraph(j)->GetN();
    for (int j = 0; j < ln - 1; j++) s_max_n = std::max(s_max_n, num_of_vtx[j]);
    int neighbors[s_max_n], num_of_neighbors, term_cnt[s_max_n], matched_cc_index[kp_core_ccs.size()];
    long long *id2vtx[ln];

    for (int j = 0; j < ln; j++) {
        id2vtx[j] = new long long[num_of_vtx[j]];
        mg.GetGraph(j)->LoadId2VtxMap(id2vtx[j]);
    }

    auto file_suf = {"_kp_core.txt", "_k_core.txt"};
    auto f_iter = file_suf.begin();

    for (auto &core_cc:{kp_core_ccs, k_core_ccs}) {
        auto out = ofstream(file_pref + *f_iter);

        id = 0;
        for (auto &cc:core_cc) {

            DumpCore(id, cc.data(), (int) cc.size(), id2vtx[ln - 1], vtx2name[ln - 1], out);

            for (int i = 0; i < ln - 1; i++) {

                auto bg = mg.GetClgFromPrime(i);
                GetCrossLayerNeighbor(bg, cc, neighbors, num_of_neighbors, term_cnt);

                // sort neighbors
                std::sort(neighbors, neighbors + num_of_neighbors,
                          [&term_cnt](int x, int y) { return term_cnt[x] > term_cnt[y]; });

                DumpCoreCNbr(i, neighbors, num_of_neighbors, id2vtx[i], vtx2name[i], term_cnt, out);
            }

            id++;
            out << endl;
        }

        out.close();
        f_iter++;
    }

    MatchCC(kp_core_ccs, k_core_ccs, matched_cc_index);
    auto out = std::ofstream(file_pref + "_match.txt");
    for (id = 0; id < kp_core_ccs.size(); id++) {
        out << id << " " << matched_cc_index[id] << endl;
    }
    out.close();
}

void
KPCCaseStudy::GetCrossLayerNeighbor(BipartiteGraph *bg, const vector<int> &cc, int *neighbors, int &num_of_neighbors,
                                    int *cnt) {
    int u, **adj_lst = bg->GetAdjLst();

    memset(cnt, 0, bg->GetN2() * sizeof(int));

    num_of_neighbors = 0;
    for (auto v : cc) {
        for (int j = 1; j <= adj_lst[v][0]; j++) {
            u = adj_lst[v][j];
            if (!cnt[u]) neighbors[num_of_neighbors++] = u;
            cnt[u]++;
        }
    }
}


void KPCCaseStudy::MatchCC(vector<vector<int>> &ccs, vector<vector<int>> &k_core_ccs, int *mapped_cc_index) {
    int id = 0;

    for (auto &cc:ccs) {
        for (int j = 0; j < k_core_ccs.size(); j++) {

            auto &kcc = k_core_ccs[j];
            auto iter = std::find(kcc.begin(), kcc.end(), cc[0]);
            if (iter != kcc.end()) {
                mapped_cc_index[id++] = j;

#ifdef MY_DEBUG
                cout << cc.size() << "->" << kcc.size() << endl;
#endif
                break;
            }
        }
    }
}

void KPCCaseStudy::DumpCore(int id, const int *cc, int length, long long *id2vtx,
                            std::unordered_map<long long, string> vtx2name, std::basic_ofstream<char> &out) {
    int i;

    out << "C" << id << ": ";
    for (i = 0; i < length - 1; i++) {
        out << vtx2name[id2vtx[cc[i]]] + ", ";
    }
    out << vtx2name[id2vtx[cc[i]]] << " (" << length << ")" << endl;
}

void KPCCaseStudy::DumpCoreCNbr(int id, const int *nbr, int length, long long *id2vtx,
                                std::unordered_map<long long, string> vtx2name,
                                int *cnt, std::basic_ofstream<char> &out) {
    int i;

    length = std::min(length, DEFAULT_TOP_FEATURE_CNT);

    out << "N" << id << ": ";
    for (i = 0; i < length - 1; i++) {
        out << vtx2name[id2vtx[nbr[i]]] << "(" << cnt[nbr[i]] << "), ";
    }
    out << vtx2name[id2vtx[nbr[i]]] << "(" << cnt[nbr[i]] << ")" << endl;
}

