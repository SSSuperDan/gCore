//
// Created by ldd on 2021/9/24.
//

#include "CmdParser.hpp"

#include "KPTreeUtils.h"

#include "KPCEvalTest.h"
#include "KPCEffTest.h"
#include "KPCCaseStudy.h"


void LoadGraph(MultilayerGraph &mg, Params &params) {

    Timer timer;
    mg.LoadGraph(params.path, params.dataset);
    timer.Stop();

    cout << "Graph " << params.dataset << " has been loaded in " << timer.GetTimeInMs() << "ms." << endl;
}

int main(int argc, char *argv[]) {

    rlimit limit{};
    limit.rlim_cur = 500 * 1024 * 1024;  // 500MB
    limit.rlim_max = RLIM_INFINITY;
    if (setrlimit(RLIMIT_STACK, &limit)) {
        printf("set limit failed\n");
        return -1;
    }

    bool release_mode = true, enable_mem_track = false;

#ifdef MY_DEBUG
    release_mode = false;
#endif

#ifdef USE_MEMTRACK
    enable_mem_track = true;
#endif

    // print running mode.
    if (release_mode) {
        if (enable_mem_track) cout << "Release mode + Mem track." << endl;
        else cout << "Release mode." << endl;
    } else {
        if (enable_mem_track) cout << "Debug mode + Mem track." << endl;
        else cout << "Debug mode." << endl;
    }


    Params params;
    ParseCmd(params, argc, argv);

    if (params.mode == RUN_CASE_STUDY) {
        assert_case_study_graph_specified(params);
        warn_empty_output_path(params);

        if (!strcasecmp(params.dataset.c_str(), "dblp")) {
            KPCCaseStudy::DBLPCaseAnalyze(params.output, params.k, params.p);
        } else if (!strcasecmp(params.dataset.c_str(), "twitter")) {
            KPCCaseStudy::TwitterCaseAnalyze(params.output, params.k, params.p);
        } else {
            cerr << "Only case studies for the DBLP and Twitter datasets are provided." << endl;
            return -1;
        }
        return 0;
    }

    if (params.mode == CMP_GCI_INFO) {

#ifdef USE_MEMTRACK
        assert_non_empty_p_tree_kp_tree_file(params, true);

        if (params.kp_tree_file.empty()) {
            assert_valid_p_tree_file(params, true);
            assert_non_empty_f2i_file(params, true);

            P_tree_f2i p_tree_f2i;
            KPTreeUtils::LoadPTree(params.p_tree_file.back(), params.f2i_file.back(), p_tree_f2i);
            KPCEffTest::GCIStatistics(p_tree_f2i, params.output);
        } else {
            assert_valid_kp_tree_file(params, true);
            assert_non_empty_f2i_file(params, true);

            KP_tree_f2i kp_tree_f2i;
            KPTreeUtils::LoadKPTree(params.kp_tree_file.back(), params.f2i_file.back(), kp_tree_f2i);
            KPCEffTest::GCIStatistics(kp_tree_f2i, params.output);
        }

        return 0;
#else
        cerr << "To compute memory consumption, please recompile the program with option -DUSE_MEMTRACK=on." << endl;
        exit(-1);
#endif

    }

    if (params.mode == GEN_RANDOM_P) {
        assert_non_empty_num_of_testcase(params);
        warn_empty_output_path(params);

        if (params.dim) KPVecGenerator::GenRandomPVectors(params.dim, params.num_of_testcase, params.output);
        else {
            assert_non_empty_graph_path(params);
            KPVecGenerator::GenRandomPVectors(MultilayerGraph::ParseLayerNumber(params.path) - 1,
                                              params.num_of_testcase, params.output);
        }

        return 0;
    }

    MultilayerGraph mg;
    assert_non_empty_graph_info(params);

    if (params.mode == BUILD_P_TREE) {
        assert_non_empty_k(params);
        warn_empty_output_path(params);

        LoadGraph(mg, params);
        KPTreeUtils::BuildPTree(mg, params.k, params.builder, params.output);

    } else if (params.mode == BUILD_KP_TREE) {
        warn_empty_output_path(params);

        LoadGraph(mg, params);
        if (params.sampled_k_file.empty()) {
            warn_empty_step(params);
            KPTreeUtils::BuildKPTree(mg, params.builder, params.output, params.step, params.ik);
        } else {
            assert_valid_sampled_k_file(params);
            KPTreeUtils::BuildKPTree(mg, params.sampled_k_file, params.builder, params.output);
        }

    } else if (params.mode == GEN_RANDOM_K) {
        assert_non_empty_num_of_testcase(params);
        warn_empty_output_path(params);

        LoadGraph(mg, params);
        KPVecGenerator::GenRandomKVectors(mg, params.num_of_testcase, params.output, params.ek);

    } else if (params.mode == GEN_RANDOM_KP) {
        assert_non_empty_num_of_testcase(params);
        warn_empty_output_path(params);

        LoadGraph(mg, params);
        if (params.sampled_k_file.empty()) {
            KPVecGenerator::GenRandomKPVectors(mg, params.num_of_testcase, params.output, params.ek);
        } else {
            assert_valid_sampled_k_file(params);
            KPVecGenerator::GenRandomKPVectors(mg, params.sampled_k_file, params.num_of_testcase, params.output);
        }
    } else if (params.mode == CMP_PILLAR_GCS) {
        assert_non_empty_k_sampled_k_file(params);

        LoadGraph(mg, params);
        assert_mlg_pillar(mg);

        if (params.sampled_k_file.empty()) {
            KPCEffTest::P_GCS(mg, params.k, params.output);
        } else KPCEffTest::P_GCS(mg, params.sampled_k_file, params.output);

    } else if (params.mode == CMP_GCS) {

        if (!params.p_tree_file.empty()) {
            assert_valid_p_tree_file(params, false);
            assert_non_empty_f2i_file(params, false);
            assert_p_tree_file_matched(params);

            assert_non_empty_kp(params);

            int num_of_p_tree = (int) params.p_tree_file.size();
            vector<P_tree_f2i> p_tree_f2is(num_of_p_tree);
            for (int i = 0; i < num_of_p_tree; i++) {
                KPTreeUtils::LoadPTree(params.p_tree_file[i], params.f2i_file[i], p_tree_f2is[i]);
            }

            LoadGraph(mg, params);
            KPCEffTest::GCS(mg, p_tree_f2is, params.k, params.p, params.output);
        } else {
            assert_non_empty_kp_sampled_kp_file(params);

            vector<KP_tree_f2i> kp_tree_f2is;
            if (!params.kp_tree_file.empty()) {

                assert_valid_kp_tree_file(params, false);
                assert_non_empty_f2i_file(params, false);
                assert_kp_tree_file_matched(params);

                int num_of_kp_tree = (int) params.kp_tree_file.size();
                kp_tree_f2is.resize(num_of_kp_tree);
                for (int i = 0; i < num_of_kp_tree; i++) {
                    KPTreeUtils::LoadKPTree(params.kp_tree_file[i], params.f2i_file[i], kp_tree_f2is[i]);
                }
            }

            LoadGraph(mg, params);
            if (params.sampled_kp_file.empty()) {
                KPCEffTest::GCS(mg, kp_tree_f2is, params.k, params.p, params.output);
            } else {
                KPCEffTest::GCS(mg, kp_tree_f2is, params.sampled_kp_file, params.output);
            }
        }
    } else if (params.mode == CMP_K_VALUE) {
        assert_non_empty_k(params);
        assert_non_empty_p(params);
        warn_empty_ck(params);

        warn_empty_output_path(params);

        LoadGraph(mg, params);
        KPCEvalTest::ComputeKValue(mg, params.k, params.ck, params.p, params.output);
    } else if (params.mode == CMP_P_VALUE) {
        assert_non_empty_k(params);
        assert_non_empty_p(params);
        warn_empty_ck(params);

        warn_empty_output_path(params);

        LoadGraph(mg, params);
        KPCEvalTest::ComputePValue(mg, params.k, params.ck, params.p, params.output);
    } else if (params.mode == CMP_SIZE_MATRIX) {
        assert_non_empty_pk(params);
        assert_non_empty_p_size_step(params);
        warn_empty_output_path(params);

        LoadGraph(mg, params);
        KPCEvalTest::ComputeSizeMatrix(mg, params.pk, params.size_step, params.output);
    } else if (params.mode == CMP_SIZE_MATRIX_K) {
        assert_non_empty_k(params);
        warn_empty_output_path(params);

        LoadGraph(mg, params);
        KPCEvalTest::ComputeSizeMatrixFixK(mg, params.k, params.output);
    } else if (params.mode == PRT_GRAPH_INFO) {
        LoadGraph(mg, params);
        mg.PrintStatistics();
    }

    return 0;
}

