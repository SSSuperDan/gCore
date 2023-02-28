//
// Created by ldd on 2022/10/12.
//

#ifndef CSM4GMG_CMDPARSER_HPP
#define CSM4GMG_CMDPARSER_HPP

#include "Utils/ArrayUtils.h"
#include "Params.h"
#include <getopt.h>

// (k,p)-core related running mode
#define RUN_CASE_STUDY 10
#define GEN_RANDOM_K 11
#define GEN_RANDOM_P  12
#define GEN_RANDOM_KP 13
#define CMP_PILLAR_GCS 14
#define BUILD_P_TREE 15
#define BUILD_KP_TREE 16
#define CMP_GCS 17
#define CMP_GCI_INFO 18

#define CMP_K_VALUE 23
#define CMP_P_VALUE 24
#define CMP_SIZE_MATRIX 25
#define CMP_SIZE_MATRIX_K 26
#define PRT_GRAPH_INFO 27

// P-tree builder
#define NAIVE_BUILD 40
#define NE_BUILD 41
#define SE_BUILD 42
#define NESE_BUILD 43


// Parameter
#define PROVIDE_TESTCASE_NUM 50
#define PROVIDE_SAMPLED_K_FILE 51
#define PROVIDE_SAMPLED_P_FILE 52
#define PROVIDE_SAMPLED_KP_FILE 53
#define PROVIDE_K_VECTOR 54
#define PROVIDE_P_VECTOR 55
#define PROVIDE_P_TREE_BUILDER 56
#define PROVIDE_KP_TREE_STEP 57
#define PROVIDE_OUTPUT_PATH 58
#define PROVIDE_INIT_K_VECTOR 59
#define PROVIDE_CK_VECTOR 60
#define PROVIDE_PK 61
#define PROVIDE_P_TREE_FILE 62
#define PROVIDE_KP_TREE_FILE 63
#define PROVIDE_F2I_FILE 64
#define PROVIDE_P_VALUE_STEP 65
#define PROVIDE_END_K_VECTOR 66

#define PROVIDE_DIM 68
#define SPECIFY_GRAPH 69

#define DIR_CREATE_MODE 0777

struct String2Mode {
    const char *name;
    int id;
};

static option long_options[] = {

        {"ntc",  required_argument, nullptr, PROVIDE_TESTCASE_NUM},
        {"skf",  required_argument, nullptr, PROVIDE_SAMPLED_K_FILE},
        {"spf",  required_argument, nullptr, PROVIDE_SAMPLED_P_FILE},
        {"skpf", required_argument, nullptr, PROVIDE_SAMPLED_KP_FILE},
        {"k",    required_argument, nullptr, PROVIDE_K_VECTOR},
        {"p",    required_argument, nullptr, PROVIDE_P_VECTOR},
        {"b",    required_argument, nullptr, PROVIDE_P_TREE_BUILDER},
        {"pb",   required_argument, nullptr, PROVIDE_P_TREE_BUILDER},
        {"s",    required_argument, nullptr, PROVIDE_KP_TREE_STEP},
        {"o",    required_argument, nullptr, PROVIDE_OUTPUT_PATH},
        {"ik",   required_argument, nullptr, PROVIDE_INIT_K_VECTOR},
        {"ck",   required_argument, nullptr, PROVIDE_CK_VECTOR},
        {"pk",   required_argument, nullptr, PROVIDE_PK},
        {"ptf",  required_argument, nullptr, PROVIDE_P_TREE_FILE},
        {"kptf", required_argument, nullptr, PROVIDE_KP_TREE_FILE},
        {"f2if", required_argument, nullptr, PROVIDE_F2I_FILE},
        {"ps",   required_argument, nullptr, PROVIDE_P_VALUE_STEP},
        {"ek",   required_argument, nullptr, PROVIDE_END_K_VECTOR},
        {"d",    required_argument, nullptr, PROVIDE_DIM},
        {"dim",  required_argument, nullptr, PROVIDE_DIM},
        {"g",    required_argument, nullptr, SPECIFY_GRAPH},
};

const String2Mode execution_mode_options[] = {
        {"cs",    RUN_CASE_STUDY},
        {"grk",   GEN_RANDOM_K},
        {"grp",   GEN_RANDOM_P},
        {"grkp",  GEN_RANDOM_KP},
        {"pgcs",  CMP_PILLAR_GCS},
        {"bp",    BUILD_P_TREE},
        {"bkp",   BUILD_KP_TREE},
        {"gcs",   CMP_GCS},
        {"gcii",  CMP_GCI_INFO},
        {"kv",    CMP_K_VALUE},
        {"pv",    CMP_P_VALUE},
        {"sm",    CMP_SIZE_MATRIX},
        {"smk",   CMP_SIZE_MATRIX_K},
        {"info",  PRT_GRAPH_INFO},

        {nullptr, 0}
};

const String2Mode p_tree_builder_options[] = {
        {"naive", NAIVE_BUILD},
        {"ne",    NE_BUILD},
        {"se",    SE_BUILD},
        {"nese",  NESE_BUILD},

        {nullptr, 0}
};


static int GetMode(const String2Mode *options, const char *str_mode) {
    for (int i = 0; options[i].name; i++) {
        if (!strcasecmp(options[i].name, str_mode)) {
            return options[i].id;
        }
    }
    return -1;
}


void ParseCmd(Params &params, int argc, char *argv[]) {
    int c, option_index;

    // read options
    while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case PROVIDE_TESTCASE_NUM:
                if (optarg) params.num_of_testcase = std::stoi(optarg);
                break;
            case PROVIDE_SAMPLED_K_FILE:
                if (optarg) params.sampled_k_file = optarg;
                break;
            case PROVIDE_SAMPLED_P_FILE:
                if (optarg) params.sampled_p_file = optarg;
                break;
            case PROVIDE_SAMPLED_KP_FILE:
                if (optarg) params.sampled_kp_file = optarg;
                break;
            case PROVIDE_K_VECTOR:
                if (optarg) String2Vec(optarg, params.k);
                break;
            case PROVIDE_P_VECTOR:
                if (optarg) String2Vec(optarg, params.p);
                break;
            case PROVIDE_P_TREE_BUILDER:
                if (optarg) {
                    auto bid = GetMode(p_tree_builder_options, optarg);
                    if (bid == NAIVE_BUILD) params.builder = NAIVE;
                    else if (bid == NE_BUILD) params.builder = NE;
                    else if (bid == SE_BUILD) params.builder = SE;
                    else if (bid == NESE_BUILD) params.builder = NESE;
                    else {
                        cerr << "Invalid P-tree builder." << endl;
                        exit(-1);
                    }
                }
                break;
            case PROVIDE_KP_TREE_STEP:
                if (optarg) params.step = std::stoi(optarg);
                break;
            case PROVIDE_OUTPUT_PATH:
                if (optarg) params.output = optarg;
                break;
            case PROVIDE_INIT_K_VECTOR:
                if (optarg) String2Vec(optarg, params.ik);
                break;
            case PROVIDE_CK_VECTOR:
                if (optarg) String2Vec(optarg, params.ck);
                break;
            case PROVIDE_PK:
                if (optarg) params.pk = std::stoi(optarg);
                break;
            case PROVIDE_P_TREE_FILE:
                if (optarg) params.p_tree_file.emplace_back(optarg);
                break;
            case PROVIDE_KP_TREE_FILE:
                if (optarg) params.kp_tree_file.emplace_back(optarg);
                break;
            case PROVIDE_F2I_FILE:
                if (optarg) params.f2i_file.emplace_back(optarg);
                break;
            case PROVIDE_P_VALUE_STEP:
                if (optarg) params.size_step = std::stof(optarg);
                break;
            case PROVIDE_END_K_VECTOR:
                if (optarg) String2Vec(optarg, params.ek);
                break;
            case PROVIDE_DIM:
                if (optarg) params.dim = std::stoi(optarg);
                break;
            case SPECIFY_GRAPH:
                if (optarg) params.dataset = optarg;
                break;

            default:
                cerr << "Invalid option(s)." << endl;
                exit(-1);
        }
    }

    // read operands
    if (argc > optind) {
        params.mode = GetMode(execution_mode_options, argv[optind]);
        if (params.mode == -1) {
            cerr << "Invalid execution mode." << endl;
            exit(-1);
        }
    } else {
        cerr << "Execution mode is missing." << endl;
        exit(-1);
    }

    if (params.mode == RUN_CASE_STUDY) return;

    optind++;
    if (argc > optind) {
        params.path = argv[optind];

        optind++;
        if (argc > optind) {
            params.dataset = argv[optind];
        }
    }
}

static bool is_invalid_path(const string &file) {
    struct stat buffer{};

    return stat(file.c_str(), &buffer);
}

static void assert_non_empty_graph_path(Params &params) {

    if (params.path.empty()) {
        cerr << "Graph path is missing." << endl;
        exit(-1);
    }

    if (is_invalid_path(params.path)) {
        cerr << "Provided graph path doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_non_empty_graph_name(Params &params) {
    if (params.dataset.empty()) {
        cerr << "Graph name is missing." << endl;
        exit(-1);
    }
}

static void assert_case_study_graph_specified(Params &params) {
    if (params.dataset.empty()) {
        cerr << "Graph for case study is not specified." << endl;
        exit(-1);
    }
}


static void assert_non_empty_graph_info(Params &params) {
    assert_non_empty_graph_path(params);
    assert_non_empty_graph_name(params);
}

static void warn_empty_output_path(Params &params) {

    if (params.output.empty()) {
        params.output = "output/";
        if (is_invalid_path(params.output)) {
            mkdir(params.output.c_str(), DIR_CREATE_MODE);
            cerr << "Output path is not provided, and is set to \"output/\" by default." << endl;
        }
    } else if (is_invalid_path(params.output)) {
        cerr << "Provided output path doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_non_empty_k(Params &params) {
    if (params.k.empty()) {
        cerr << "A coreness vector (k) must be provided." << endl;
        exit(-1);
    }
}

static void assert_non_empty_p(Params &params) {
    if (params.p.empty()) {
        cerr << "A fraction vector (p) must be provided." << endl;
        exit(-1);
    }
}

static void assert_non_empty_kp(Params &params) {
    if (params.k.empty() || params.p.empty()) {
        cerr << "A coreness vector (k) and a fraction vector (p) must be provided." << endl;
        exit(-1);
    }
}

static void warn_empty_step(Params &params) {
    if (params.step <= 0) {
        params.step = 1;
        cerr << "Incremental step of coreness vectors is not provided, and is set to 1 by default." << endl;
    }
}

static void warn_empty_ck(Params &params) {
    if (params.ck.empty()) {

        int dim = (int) params.k.size() - 1;
        params.ck.assign(dim, 1);
        cerr << "Cross-layer degree threshold vector (ck) is not provided, and is set to [1]^" << dim << " by default."
             << endl;
    }
}

static void assert_non_empty_num_of_testcase(Params &params) {
    if (params.num_of_testcase <= 0) {
        cerr << "Number of testcases must be specified." << endl;
        exit(-1);
    }
}

static void assert_valid_sampled_k_file(Params &params) {
    if (is_invalid_path(params.sampled_k_file)) {
        cerr << "File \"" << params.sampled_k_file << "\" doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_valid_sampled_p_file(Params &params) {
    if (is_invalid_path(params.sampled_p_file)) {
        cerr << "File \"" << params.sampled_p_file << "\" doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_valid_sampled_kp_file(Params &params) {
    if (is_invalid_path(params.sampled_kp_file)) {
        cerr << "File \"" << params.sampled_kp_file << "\" doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_non_empty_k_sampled_k_file(Params &params) {
    if (params.k.empty()) {
        if (params.sampled_k_file.empty()) {
            cerr << "Either a coreness vector (k) or a file of sampled coreness vectors has to be provided." << endl;
            exit(-1);
        }
        assert_valid_sampled_k_file(params);
    }
}

static void assert_non_empty_kp_sampled_kp_file(Params &params) {
    if (params.k.empty() || params.p.empty()) {
        if (params.sampled_kp_file.empty()) {
            cerr << "Either a coreness vector (k) and a fraction vector (p),"
                    " or a file of sampled (k,p) pairs have to be provided." << endl;
            exit(-1);

        }
        assert_valid_sampled_kp_file(params);
    }
}

static void assert_valid_files(vector<string> &files, const string &type, bool check_last) {
    vector<string>::iterator iter;

    if (check_last) {
        iter = files.end() - 1;
        if (files.size() > 1) {
            cerr << "Multiple " << type << " files are provided, only the last one will be considered." << endl;
        }
    } else {
        iter = files.begin();
    }

    while (iter != files.end()) {
        if (is_invalid_path(*iter)) {
            cerr << "File \"" << *iter << "\" doesn't exist." << endl;
            exit(-1);
        }
        iter++;
    }
}

static void assert_valid_p_tree_file(Params &params, bool check_last) {
    assert_valid_files(params.p_tree_file, "P-tree", check_last);
}

static void assert_valid_kp_tree_file(Params &params, bool check_last) {
    assert_valid_files(params.kp_tree_file, "KP-tree", check_last);
}

static void assert_non_empty_p_tree_kp_tree_file(Params &params, bool check_last) {
    if (params.p_tree_file.empty() && params.kp_tree_file.empty()) {
        if (check_last) {
            cerr << "Either a P-tree file or a KP-tree file must be provided." << endl;
        } else {
            cerr << "Either P-tree files or KP-tree files must be provided." << endl;
        }

        exit(-1);
    }
}

static void assert_non_empty_f2i_file(Params &params, bool check_last) {

    if (params.f2i_file.empty()) {
        if (check_last) cerr << "A fraction to integer map file has to be provided." << endl;
        else cerr << "Fraction to integer map files have to be provided." << endl;
        exit(-1);
    }

    assert_valid_files(params.f2i_file, "fraction to integer map", check_last);
}

static void assert_p_tree_file_matched(Params &params) {
    if (params.f2i_file.size() != params.p_tree_file.size()) {
        cerr << "The number of provided P-tree files and fraction to integer map files are inconsistent." << endl;
        exit(-1);
    }
}

static void assert_kp_tree_file_matched(Params &params) {
    if (params.f2i_file.size() != params.kp_tree_file.size()) {
        cerr << "The number of provided KP-tree files and fraction to integer map files are inconsistent." << endl;
        exit(-1);
    }
}

static void assert_non_empty_pk(Params &params) {
    if (params.pk <= 0) {
        cerr << "A valid coreness threshold for primary layer must be provided." << endl;
        exit(-1);
    }
}

static void assert_non_empty_p_size_step(Params &params) {
    if (params.size_step <= 0) {
        cerr << "A valid incremental step of fraction vectors for computing size matrices must be provided." << endl;
        exit(-1);
    }
}

static void assert_mlg_pillar(MultilayerGraph &mg) {
    if (mg.IsGeneral()) {
        cerr << "Provided multi-layer graph should be pillar, i.e., multiplex." << endl;
        exit(-1);
    }
}

#endif //CSM4GMG_CMDPARSER_HPP