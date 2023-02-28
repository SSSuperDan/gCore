//
// Created by ldd on 2022/10/12.
//

#ifndef CSM4GMG_PARAMS_H
#define CSM4GMG_PARAMS_H

#include "Core/PTreeBuilder.h"

struct Params {

    // execution mode
    int mode;

    // arguments
    int num_of_testcase = -1;
    int step = -1;
    int pk = -1;
    int dim = -1;
    float size_step = -1;
    p_tree_builder builder = NESE;
    vector<int> ik;
    vector<int> ck;
    vector<int> ek;
    vector<int> k;
    vector<float> p;
    vector<string> p_tree_file;
    vector<string> kp_tree_file;
    vector<string> f2i_file;
    string sampled_k_file;
    string sampled_p_file;
    string sampled_kp_file;
    string path;
    string dataset;
    string output;
};

#endif //CSM4GMG_PARAMS_H
