//
// Created by ldd on 2022/7/15.
//

#include "KPCEvalTest.h"


void KPCEvalTest::ComputeKValue(MultilayerGraph &mg, vector<int> &k, vector<int> &ck, vector<float> &p,
                                const string &path) {

    int ln = mg.GetLayerNumber(), n = mg.GetPrimaryGraph()->GetN(), s_max_n = 0, **adj_lst;
    for (int i = 0; i < ln - 1; i++) s_max_n = std::max(s_max_n, mg.GetGraph(i)->GetN());

    int kp_core[n], rc[n], k_core[n], length1, length2, length3;
    int core_number[s_max_n], s_nbr[s_max_n], s_nbr_pos[s_max_n], s_length;
    int k_value_1[n], k_value_2[n], k_value_3[n];

    // compute (k,p)-core
    KPCore kp_core_extractor(mg);
    kp_core_extractor.Extract(k, p, kp_core, length1);

    // compute rc
    RCD rc_extractor(mg);
    rc_extractor.Extract(k, ck, rc, length2);

    // compute k-core
    UndirectedGraph *g = mg.GetPrimaryGraph();
    length3 = CoreDec::KCore(g, k[ln - 1], k_core);

    auto k_number_file = path + mg.GetGraphName() + "_k_number.txt";
    auto out = std::ofstream(k_number_file);
    for (int i = 0; i < ln - 1; i++) {
        g = mg.GetGraph(i);
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();

        // (k,p)-core
        GetCrossLayerNeighbor(mg, kp_core, length1, i, s_nbr, s_length);
        for (int j = 0; j < s_length; j++) s_nbr_pos[s_nbr[j]] = j;

        CoreDec::ArrayBasedCoreDec(g, s_nbr, s_length, core_number);
        for (int j = 0; j < length1; j++) {
            k_value_1[j] = GetKValue(adj_lst[kp_core[j]], s_nbr_pos, core_number, p[i]);
        }

        // rcd
        GetCrossLayerNeighbor(mg, rc, length2, i, s_nbr, s_length);
        for (int j = 0; j < s_length; j++) s_nbr_pos[s_nbr[j]] = j;

        CoreDec::ArrayBasedCoreDec(g, s_nbr, s_length, core_number);
        for (int j = 0; j < length2; j++) {
            k_value_2[j] = GetKValue(adj_lst[rc[j]], s_nbr_pos, core_number, p[i]);
        }

        // k-core
        GetCrossLayerNeighbor(mg, k_core, length3, i, s_nbr, s_length);
        for (int j = 0; j < s_length; j++) s_nbr_pos[s_nbr[j]] = j;

        CoreDec::ArrayBasedCoreDec(g, s_nbr, s_length, core_number);
        for (int j = 0; j < length3; j++) {
            k_value_3[j] = GetKValue(adj_lst[k_core[j]], s_nbr_pos, core_number, p[i]);
        }

        std::sort(k_value_1, k_value_1 + length1);
        std::sort(k_value_2, k_value_2 + length2);
        std::sort(k_value_3, k_value_3 + length3);

        out << DumpTestSettings(k, ck, p, i) << endl;

        for (int j = 0; j < length1; j++) out << k_value_1[j] << " ";
        out << endl;

        for (int j = 0; j < length2; j++) out << k_value_2[j] << " ";
        out << endl;

        for (int j = 0; j < length3; j++) out << k_value_3[j] << " ";
        out << endl;
    }

    out.close();
}

void KPCEvalTest::ComputePValue(MultilayerGraph &mg, vector<int> &k, vector<int> &ck, vector<float> &p,
                                const string &path) {

    int ln = mg.GetLayerNumber(), n = mg.GetPrimaryGraph()->GetN(), s_max_n = 0, **adj_lst;
    for (int i = 0; i < ln - 1; i++) s_max_n = std::max(s_max_n, mg.GetGraph(i)->GetN());

    int kp_core[n], rc[n], k_core[n], length1, length2, length3;
    int s_nbr[s_max_n], s_length, s_kp_core[s_max_n], s_rc[s_max_n], s_k_core[s_max_n];
    bool in_s_core[s_max_n];
    double p_value_1[n], p_value_2[n], p_value_3[n];

    // compute (k,p)-core
    KPCore kp_core_extractor(mg);
    kp_core_extractor.Extract(k, p, kp_core, length1);

    // compute rc
    RCD rc_extractor(mg);
    rc_extractor.Extract(k, ck, rc, length2);

    // compute k-core
    UndirectedGraph *g = mg.GetPrimaryGraph();
    length3 = CoreDec::KCore(g, k[ln - 1], k_core);

    auto p_number_file = path + mg.GetGraphName() + "_p_number.txt";
    auto out = std::ofstream(p_number_file);

    for (int i = 0; i < ln - 1; i++) {

        g = mg.GetGraph(i);
        adj_lst = mg.GetClgFromPrime(i)->GetAdjLst();

        // (k,p)-core
        GetCrossLayerNeighbor(mg, kp_core, length1, i, s_nbr, s_length);
        s_length = CoreDec::KCore(g, s_nbr, s_length, k[i], s_kp_core);

        memset(in_s_core, false, s_max_n * sizeof(bool));
        for (int j = 0; j < s_length; j++) in_s_core[s_kp_core[j]] = true;

        for (int j = 0; j < length1; j++) {
            p_value_1[j] = GetPValue(adj_lst[kp_core[j]], in_s_core);
        }

        // rcd
        GetCrossLayerNeighbor(mg, rc, length2, i, s_nbr, s_length);
        s_length = CoreDec::KCore(g, s_nbr, s_length, k[i], s_rc);

        memset(in_s_core, false, s_max_n * sizeof(bool));
        for (int j = 0; j < s_length; j++) in_s_core[s_rc[j]] = true;

        for (int j = 0; j < length2; j++) {
            p_value_2[j] = GetPValue(adj_lst[rc[j]], in_s_core);
        }

        // k-core
        GetCrossLayerNeighbor(mg, k_core, length3, i, s_nbr, s_length);
        s_length = CoreDec::KCore(g, s_nbr, s_length, k[i], s_k_core);

        memset(in_s_core, false, s_max_n * sizeof(bool));
        for (int j = 0; j < s_length; j++) in_s_core[s_k_core[j]] = true;

        for (int j = 0; j < length3; j++) {
            p_value_3[j] = GetPValue(adj_lst[k_core[j]], in_s_core);
        }

        std::sort(p_value_1, p_value_1 + length1);
        std::sort(p_value_2, p_value_2 + length2);
        std::sort(p_value_3, p_value_3 + length3);

        out << DumpTestSettings(k, ck, p, i) << endl;

        for (int j = 0; j < length1; j++) out << p_value_1[j] << " ";
        out << endl;

        for (int j = 0; j < length2; j++) out << p_value_2[j] << " ";
        out << endl;

        for (int j = 0; j < length3; j++) out << p_value_3[j] << " ";
        out << endl;
    }

    out.close();

}

void KPCEvalTest::ComputeSizeMatrix(MultilayerGraph &mg, int pk, float step, const string &path) {

    int ln = mg.GetLayerNumber(), n = mg.GetPrimaryGraph()->GetN(), core[n], length, mcn;

    vector<int> k(ln), int_p(ln - 1);
    vector<float> p(ln - 1);

    PTreeBuilder p_tree_builder(mg);
    k[ln - 1] = pk;

    for (int i = 0; i < ln - 1; i++) {
        auto size_matrix_file = path + mg.GetGraphName() + "_size_matrix_pk" + to_string(pk) + "_l" + to_string(i) + ".txt";
        auto out = std::ofstream(size_matrix_file);

        mcn = CoreDec::GetDegeneracy(mg.GetGraph(i));
        out << mcn << ", " << step << endl;  // output "mcn, step"

        memset(k.data(), 0, (ln - 1) * sizeof(int));
        for (int j = 0; j <= mcn; j++) {

            k[i] = j;
            PTree p_tree(ln - 1);
            Frac2IntPri f2i(ln - 1);
            p_tree_builder.Execute(k, p_tree, f2i, NESE);

            memset(p.data(), 0, (ln - 1) * sizeof(float));

            while (p[i] - 1.0 < 1e-6) {
                f2i.Convert(p, int_p);
                p_tree.Search(int_p, core, length);
                out << k[i] << ", " << p[i] << ", " << length << endl;
                p[i] += step;
            }
        }
        out.close();
    }
}

void KPCEvalTest::ComputeSizeMatrixFixK(MultilayerGraph &mg, vector<int> &k, const string &path) {

    int ln = mg.GetLayerNumber(), n = mg.GetPrimaryGraph()->GetN(), core[n], length;
    Frac *fractions;
    vector<int> p(ln - 1);

    PTreeBuilder p_tree_builder(mg);
    PTree p_tree(ln - 1);
    Frac2IntPri f2i(ln - 1);
    p_tree_builder.Execute(k, p_tree, f2i, NESE);

    for (int i = 0; i < ln - 1; i++) {
        auto size_matrix_file = path + mg.GetGraphName() + "_size_matrix_k" + Vec2String(k) + "_l" + to_string(i) + ".txt";
        auto out = std::ofstream(size_matrix_file);

        fractions = f2i.fractions[i];
        memset(p.data(), 0, (ln - 1) * sizeof(int));

        for (int j = 0; j < f2i.frac_size[i]; j++) {
            p[i] = j;
            p_tree.Search(p, core, length);
            out << (fractions[j].den) << ", " << (fractions[j].num) << ", " << length << endl;
        }

        out.close();
    }
}

void
KPCEvalTest::GetCrossLayerNeighbor(MultilayerGraph &mg, const int *core, int length, int layer, int *neighbors,
                                   int &num_of_neighbors) {
    int v, u, n = mg.GetGraph(layer)->GetN();
    int **adj_lst = mg.GetClgFromPrime(layer)->GetAdjLst();
    bool visit[n];

    num_of_neighbors = 0;
    memset(visit, false, n * sizeof(bool));
    for (int i = 0; i < length; i++) {
        v = core[i];
        for (int j = 1; j <= adj_lst[v][0]; j++) {
            u = adj_lst[v][j];
            if (!visit[u]) {
                visit[u] = true;
                neighbors[num_of_neighbors++] = u;
            }
        }
    }
}

int KPCEvalTest::GetKValue(const int *adj_lst, const int *nbr_pos, const int *core_number, float p) {
    if (!adj_lst[0]) return 0;

    int size_of_nbr = adj_lst[0], c[size_of_nbr], pos = size_of_nbr - GetPCeilingPos(size_of_nbr, p);
    for (int i = 1; i <= size_of_nbr; i++) {
        c[i - 1] = core_number[nbr_pos[adj_lst[i]]];
    }
    std::nth_element(c, c + pos, c + size_of_nbr);

    return c[pos];
}

double KPCEvalTest::GetPValue(const int *adj_lst, const bool *in_s_core) {
    if (!adj_lst[0]) return 0;

    int cnt = 0;
    for (int i = 1; i <= adj_lst[0]; i++) {
        if (in_s_core[adj_lst[i]]) cnt++;
    }
    return (double) cnt / (double) adj_lst[0];
}


string KPCEvalTest::DumpTestSettings(vector<int> &k, vector<int> &ck, vector<float> &p, int layer) {
    return "k = " + Vec2String(k) + ", ck = " + Vec2String(ck) + ", p = " + Vec2String(p) + ", l = " + std::to_string(layer);
}
