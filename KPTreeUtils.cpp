//
// Created by ldd on 2022/10/14.
//

#include "KPTreeUtils.h"


string KPTreeUtils::GetPTreeFilePref(P_tree_f2i &p_tree_f2i) {
    return GetPTreeFilePref(p_tree_f2i.dataset, p_tree_f2i.k, p_tree_f2i.builder);
}

string KPTreeUtils::GetKPTreeFilePref(KP_tree_f2i &kp_tree_f2i) {
    if (kp_tree_f2i.sampled_k_size) {
        return GetKPTreeFilePref(kp_tree_f2i.dataset, kp_tree_f2i.sampled_k_size, kp_tree_f2i.builder);
    } else {
        return GetKPTreeFilePref(kp_tree_f2i.dataset, kp_tree_f2i.builder, kp_tree_f2i.step, kp_tree_f2i.ini_k);
    }
}

void KPTreeUtils::BuildPTree(MultilayerGraph &mg, vector<int> &k, p_tree_builder builder, const string &path) {

    long long int num_of_nodes;
    string p_tree_file_prefix = path + GetPTreeFilePref(mg.GetGraphName(), k, builder);

    Timer timer;
    int dim = mg.GetLayerNumber() - 1;

    PTree p_tree(dim);
    Frac2IntPri f2i(dim);
    PTreeBuilder p_tree_builder(mg);
    p_tree_builder.Execute(k, p_tree, f2i, builder);

    timer.Stop();

    // write to file
    p_tree.Flush(p_tree_file_prefix + "_p_tree");
    f2i.Flush(p_tree_file_prefix + "_f2i");

    // output
    num_of_nodes = p_tree.GetNumOfNodes();

    auto f = std::ofstream(p_tree_file_prefix + "_info");
    f << "construction_time = " << timer.GetTimeInMs() << endl;
    f << "#nodes = " << num_of_nodes << endl;
    f.close();

    cout << "PTree has been constructed in " << timer.GetTimeInMs() << "ms, with total " << num_of_nodes
         << " nodes generated." << endl;
}

void KPTreeUtils::BuildKPTree(MultilayerGraph &mg, p_tree_builder builder, const string &path, int step,
                              const vector<int> &start_k) {

    long long int num_of_nodes;
    string file_name_prefix = path + GetKPTreeFilePref(mg.GetGraphName(), builder, step, start_k);

    Timer timer;
    int ln = mg.GetLayerNumber();

    KPTree kp_tree(ln);
    Frac2IntPri f2i(ln - 1);
    KPTreeBuilder kp_tree_builder(mg);
    kp_tree_builder.Execute(kp_tree, f2i, builder, step, start_k);

    timer.Stop();

    // write to file
    kp_tree.Flush(file_name_prefix + "_kp_tree");
    f2i.Flush(file_name_prefix + "_f2i");

    // output
    num_of_nodes = kp_tree.GetNumOfPNodes();

    auto f = std::ofstream(file_name_prefix + "_info");
    f << "construction_time = " << timer.GetTimeInMs() << endl;
    f << "#nodes = " << num_of_nodes << endl;
    f.close();

    cout << "KPTree has been constructed in " << timer.GetTimeInMs() << "ms, with total " << num_of_nodes
         << " nodes generated." << endl;

}

void KPTreeUtils::BuildKPTree(MultilayerGraph &mg, const string &sampled_k_vec_file, p_tree_builder builder,
                              const string &path) {

    long long int num_of_nodes;

    vector<vector<int>> sampled_k_vec;
    KPVecGenerator::LoadKVectors(mg.GetLayerNumber(), sampled_k_vec_file, sampled_k_vec);

    string file_name_prefix = path + GetKPTreeFilePref(mg.GetGraphName(), (int) sampled_k_vec.size(), builder);

    Timer timer;
    int ln = mg.GetLayerNumber();

    KPTree kp_tree(ln);
    Frac2IntPri f2i(ln - 1);
    KPTreeBuilder kp_tree_builder(mg);
    kp_tree_builder.Execute(kp_tree, f2i, sampled_k_vec, builder);

    timer.Stop();

    // write to file
    kp_tree.Flush(file_name_prefix + "_kp_tree");
    f2i.Flush(file_name_prefix + "_f2i");

    // output
    num_of_nodes = kp_tree.GetNumOfPNodes();

    auto f = std::ofstream(file_name_prefix + "_info");
    f << "construction_time = " << timer.GetTimeInMs() << endl;
    f << "#nodes = " << num_of_nodes << endl;
    f.close();

    cout << "KPTree has been constructed in " << timer.GetTimeInMs() << "ms, with total " << num_of_nodes
         << " nodes generated." << endl;
}


void KPTreeUtils::LoadPTree(string &p_tree_file, string &f2i_file, P_tree_f2i &p_tree_f2i) {

    Timer timer;

    p_tree_f2i.p_tree = PTree::Load(p_tree_file);
    p_tree_f2i.f2i = Frac2IntPri::Load(f2i_file);

    ParseTreeInfo(p_tree_file, p_tree_f2i);

    timer.Stop();
    cout << "PTree has been loaded in " << timer.GetTimeInMs() << "ms." << endl;
}

void KPTreeUtils::LoadKPTree(string &kp_tree_file, string &f2i_file, KP_tree_f2i &kp_tree_f2i) {

    Timer timer;

    kp_tree_f2i.kp_tree = KPTree::Load(kp_tree_file);
    kp_tree_f2i.f2i = Frac2IntPri::Load(f2i_file);

    ParseTreeInfo(kp_tree_file, kp_tree_f2i);

    timer.Stop();
    cout << "KPTree has been loaded in " << timer.GetTimeInMs() << "ms." << endl;
}

string KPTreeUtils::GetPTreeFilePref(const string &dataset, vector<int> &k, p_tree_builder builder) {
    return dataset + "_" + Vec2String(k) + "_" + GetBuilderName(builder);
}


string KPTreeUtils::GetKPTreeFilePref(const string &dataset, p_tree_builder builder, int step, const vector<int> &start_k) {
    return start_k.empty() ? (dataset + "_" + GetBuilderName(builder)) + "_" + to_string(step) :
           (dataset + "_" + GetBuilderName(builder) + "_" + to_string(step) + "_" + Vec2String(start_k));
}

string KPTreeUtils::GetKPTreeFilePref(const string &dataset, int sampled_k_vector_size, p_tree_builder builder) {
    return dataset + "_" + to_string(sampled_k_vector_size) + "_" + GetBuilderName(builder);
}

void KPTreeUtils::ParseTreeInfo(string &p_tree_file, P_tree_f2i &p_tree_f2i) {
    vector<string> tokens;

    Split(p_tree_file.substr((int) p_tree_file.find_last_of('/') + 1), tokens, false, '_');

    p_tree_f2i.dataset = tokens[0];
    String2Vec(tokens[1], p_tree_f2i.k);
    p_tree_f2i.builder = GetBuilderMode(tokens[2]);
}

void KPTreeUtils::ParseTreeInfo(string &kp_tree_file, KP_tree_f2i &kp_tree_f2i) {
    vector<string> tokens;

    Split(kp_tree_file.substr((int) kp_tree_file.find_last_of('/') + 1), tokens, false, '_');

    kp_tree_f2i.dataset = tokens[0];
    try {
        kp_tree_f2i.builder = GetBuilderMode(tokens[1]);
        kp_tree_f2i.step = stoi(tokens[2]);
        if (tokens[3].substr(0, 1) == "[") String2Vec(tokens[3], kp_tree_f2i.ini_k);

    } catch (std::exception &e) {
        kp_tree_f2i.sampled_k_size = stoi(tokens[1]);
        kp_tree_f2i.builder = GetBuilderMode(tokens[2]);
    }
}

p_tree_builder KPTreeUtils::GetBuilderMode(string &builder) {
    if (builder == "naive") return NAIVE;
    else if (builder == "ne") return NE;
    else if (builder == "se") return SE;
    else if (builder == "nese") return NESE;
    else throw std::exception();
}

string KPTreeUtils::GetBuilderName(p_tree_builder builder) {

    if (builder == NAIVE) return "naive";
    else if (builder == NE) return "ne";
    else if (builder == SE) return "se";
    else if (builder == NESE) return "nese";
    else exit(-1);
}
