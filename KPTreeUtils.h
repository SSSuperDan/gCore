//
// Created by ldd on 2022/10/14.
//

#ifndef CSM4GMG_KPTREEUTILS_H
#define CSM4GMG_KPTREEUTILS_H

#include "Core/KPTreeBuilder.h"
#include "Utils/Timer.h"

#include "KPVecGenerator.h"

struct KP_tree_f2i {
    KPTree *kp_tree;
    Frac2IntPri *f2i;

    p_tree_builder builder;
    vector<int> ini_k;
    int step = 1;
    int sampled_k_size = 0;
    string dataset;

    void ReleaseKPTree() {
        delete kp_tree;
        kp_tree = nullptr;
    }

    void ReleaseF2i() {
        delete f2i;
        f2i = nullptr;
    }

    ~KP_tree_f2i() {
        if (kp_tree) ReleaseKPTree();
        if (f2i) ReleaseF2i();
    }
};

struct P_tree_f2i {
    PTree *p_tree;
    Frac2IntPri *f2i;

    p_tree_builder builder;
    vector<int> k;
    string dataset;

    void ReleasePTree() {
        delete p_tree;
        p_tree = nullptr;
    }

    void ReleaseF2i() {
        delete f2i;
        f2i = nullptr;
    }

    ~P_tree_f2i() {
        if (p_tree) ReleasePTree();
        if (f2i) ReleaseF2i();
    }
};


/*
 * This class is used to build and load P/KP-trees.
 *
 * >> For P-trees
 *      all files has a prefix of 'path/dataset_k_builder'
 *      the tree is generated and flush into file 'prefix_p_tree'
 *      the f2i map is generated and flush into file 'prefix_f2i'
 *
 * >> For KP-trees
 *      all files has a prefix of 'path/dataset_builder_step' or
 *                                'path/dataset_builder_step_initK' or
 *                                'path/dataset_sampledN_builder'
 *      the tree is generated and flush into file 'prefix_kp_tree'
 *      the f2i map is generated and flush into file 'prefix_f2i'
 *
 */

class KPTreeUtils {
public:

    static string GetPTreeFilePref(P_tree_f2i & p_tree_f2i);
    static string GetKPTreeFilePref(KP_tree_f2i &kp_tree_f2i);

    static void ParseTreeInfo(string & p_tree_file, P_tree_f2i &p_tree_f2i);
    static void ParseTreeInfo(string & kp_tree_file, KP_tree_f2i &kp_tree_f2i);

    static void BuildPTree(MultilayerGraph &mg, vector<int> &k, p_tree_builder builder, const string &path);
    static void BuildKPTree(MultilayerGraph &mg, p_tree_builder builder, const string &path, int step = 1, const vector<int> &start_k = {});
    static void BuildKPTree(MultilayerGraph &mg, const string &sampled_k_vec_file, p_tree_builder builder, const string &path);

    static void LoadPTree(string &p_tree_file, string &f2i_file, P_tree_f2i &p_tree_f2i);
    static void LoadKPTree(string &kp_tree_file, string &f2i_file, KP_tree_f2i &kp_tree_f2i);


private:
    static p_tree_builder GetBuilderMode(string &builder);
    static string GetBuilderName(p_tree_builder builder);

    static string GetPTreeFilePref(const string &dataset, vector<int> &k, p_tree_builder builder);
    static string GetKPTreeFilePref(const string &dataset, p_tree_builder builder, int step = 1, const vector<int> &start_k = {});
    static string GetKPTreeFilePref(const string &dataset, int sampled_k_vector_size, p_tree_builder builder);
};


#endif //CSM4GMG_KPTREEUTILS_H
