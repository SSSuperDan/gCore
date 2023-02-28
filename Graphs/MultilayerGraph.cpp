//
// Created by ldd on 2021/9/22.
//

#include "MultilayerGraph.h"


MultilayerGraph::~MultilayerGraph() {
    if (graph_layers) {
        for (int i = 0; i < num_of_layers; i++) {
            delete graph_layers[i];
        }
        delete[] graph_layers;
    }

    if (cross_layer_graphs) {
        for (int i = 0; i < num_of_cross_layer_graphs; i++) {
            delete cross_layer_graphs[i];
        }
        delete[] cross_layer_graphs;
    }

    delete[] ps_clg;
    delete[] sp_clg;
}


/*
 * Example "mlg.conf":
 * line 1: intra_layer_graph_file_1
 * line 2: intra_layer_graph_file_2
 * line 3: intra_layer_graph_file_3  <== primary graph
 * line 4:
 * line 5: 1, 3, cross_layer_graph_file_13
 * line 6: 2, 3, cross_layer_graph_file_23
 */

void MultilayerGraph::LoadGraph(const string &input_path, const string &name_) {
    int num_of_clf;
    bool has_bijection_layer, *is_bijective_layer;

    struct stat buffer{};
    string &&conf_path = input_path + "mlg.conf";

    if (stat((conf_path).c_str(), &buffer)) {
        cerr << "File \"mlg.conf\" is needed for building multi-layer graphs." << endl;
        exit(-1);
    }

    path = input_path;
    name = name_;

    vector<string> intra_layer_file;
    vector<bg_file> cross_layer_file;

    ParseConfFile(conf_path, intra_layer_file, cross_layer_file);

    num_of_layers = (int) intra_layer_file.size();
    num_of_clf = (int) cross_layer_file.size();
    has_bijection_layer = num_of_clf + 1 < num_of_layers;
    num_of_cross_layer_graphs = (num_of_clf << 1) + (has_bijection_layer ? 1 : 0);
    is_general = num_of_clf;

    is_bijective_layer = new bool[num_of_layers - 1];
    memset(is_bijective_layer, true, (num_of_layers - 1) * sizeof(bool));
    for (auto &clf:cross_layer_file) {
        if (clf.id1 == num_of_layers) {
            is_bijective_layer[clf.id2 - 1] = false;
        } else if (clf.id2 == num_of_layers) {
            is_bijective_layer[clf.id1 - 1] = false;
        }
    }

    LoadIntraLayerEdges(intra_layer_file, is_bijective_layer);
    if (has_bijection_layer) LoadCrossLayerEdges(cross_layer_file, is_bijective_layer);
    else LoadCrossLayerEdges(cross_layer_file);

    num_of_total_vtx += graph_layers[pid]->GetN();
    num_of_total_ie += (graph_layers[pid]->GetM() >> 1);
    for (int i = 0; i < num_of_layers - 1; i++) {
        if (!is_bijective_layer[i]) num_of_total_vtx += graph_layers[i]->GetN();
        num_of_total_ie += (graph_layers[i]->GetM() >> 1);
        num_of_total_ce += ps_clg[i]->GetM();
    }

    delete[] is_bijective_layer;
}

void MultilayerGraph::LoadIntraLayerEdges(vector<string> &intra_layer_file, const bool *is_bijective_layer) {
    pid = num_of_layers - 1;

    graph_layers = new UndirectedGraph *[num_of_layers];
    graph_layers[pid] = new UndirectedGraph();
    graph_layers[pid]->LoadEdges(path, intra_layer_file[pid]);
    auto &pl_map_file = graph_layers[pid]->GetMapFile();

    for (int i = 0; i < num_of_layers - 1; i++) {
        graph_layers[i] = new UndirectedGraph();
        if (is_bijective_layer[i]) graph_layers[i]->LoadEdges(path, intra_layer_file[i], pl_map_file);
        else graph_layers[i]->LoadEdges(path, intra_layer_file[i]);

    }
}

void MultilayerGraph::LoadCrossLayerEdges(vector<bg_file> &cross_layer_file, const bool *is_bijective_layer) {
    int id, id1, id2, num_of_vtx[num_of_layers], max_n;
    BipartiteGraph *bg1, *bg2;

    cross_layer_graphs = new BipartiteGraph *[num_of_cross_layer_graphs];
    ps_clg = new BipartiteGraph *[num_of_layers - 1];
    sp_clg = new BipartiteGraph *[num_of_layers - 1];
    for (int i = 0; i < num_of_layers; i++) num_of_vtx[i] = graph_layers[i]->GetN();

    // many-to-many mapping
    id = 0;
    for (auto &bgf :cross_layer_file) {
        id1 = bgf.id1 - 1, id2 = bgf.id2 - 1;

        bg1 = new BipartiteGraph();
        bg1->LoadEdges(graph_layers[id1]->GetMapFile(), graph_layers[id2]->GetMapFile(), path + bgf.filename);

        bg2 = new BipartiteGraph();
        bg2->BuildInvertedGraph(bg1);

        if (id1 == pid) {
            ps_clg[id2] = bg1;
            sp_clg[id2] = bg2;
        } else {
            ps_clg[id1] = bg2;
            sp_clg[id1] = bg1;
        }
        cross_layer_graphs[id++] = bg1;
        cross_layer_graphs[id++] = bg2;

        num_of_vtx[id1] = std::max(num_of_vtx[id1], bg1->GetN());
        num_of_vtx[id2] = std::max(num_of_vtx[id2], bg2->GetN());

    }


    // one-to-one mapping
    if (is_bijective_layer) {
        max_n = num_of_vtx[pid];
        for (int i = 0; i < num_of_layers - 1; i++) {
            if (is_bijective_layer[i]) max_n = std::max(max_n, num_of_vtx[i]);
        }

        bg1 = new BipartiteGraph();
        bg1->BuildBijectionById(max_n);
        cross_layer_graphs[id] = bg1;

        num_of_vtx[pid] = max_n;
        for (int i = 0; i < num_of_layers - 1; i++) {
            if (is_bijective_layer[i]) {
                ps_clg[i] = bg1;
                sp_clg[i] = bg1;

                num_of_vtx[i] = max_n;
            }
        }
    }

    Enlarge(num_of_vtx);
}


void MultilayerGraph::PrintSummary(bool print_detail) {

    cout << "Intra-layer:" << endl;
    for (int i = 0; i < num_of_layers; i++) {
        graph_layers[i]->PrintStatistics(print_detail);
    }

    cout << "Cross-layer (P==>S):" << endl;
    for (int i = 0; i < num_of_layers - 1; i++) {
        ps_clg[i]->PrintStatistics(print_detail);
    }

    cout << "Cross-layer (S==>P):" << endl;
    for (int i = 0; i < num_of_layers - 1; i++) {
        sp_clg[i]->PrintStatistics(print_detail);
    }

    cout << "Summary statistics: " << endl;
    PrintStatistics();
}

void MultilayerGraph::PrintStatistics() {
    int degeneracy[num_of_layers];
    cout << "|V| = " << num_of_total_vtx <<
    ", |IE| = " << num_of_total_ie << ", |CE|" << num_of_total_ce << ", degeneracy = ";
    for (int i = 0; i < num_of_layers; i++) {
        degeneracy[i] = CoreDec::GetDegeneracy(graph_layers[i]);
    }
    cout << Array2String(degeneracy, num_of_layers) << endl;
}

inline void MultilayerGraph::Enlarge(int *new_num_of_vtx) {

    for (int i = 0; i < num_of_layers; i++) {
        graph_layers[i]->Enlarge(new_num_of_vtx[i]);
    }

    for (int i = 0; i < num_of_layers - 1; i++) {
        ps_clg[i]->Enlarge(new_num_of_vtx[pid], new_num_of_vtx[i]);
        sp_clg[i]->Enlarge(new_num_of_vtx[i], new_num_of_vtx[pid]);
    }
}

int MultilayerGraph::ParseLayerNumber(const string &path) {
    int ln = 0;
    struct stat buffer{};
    string line, &&conf_path = move(path + "mlg.conf");

    if (stat((conf_path).c_str(), &buffer)) {
        cerr << "File \"mlg.conf\" is needed for building multi-layer graphs." << endl;
        exit(-1);
    }

    auto fin = ifstream(conf_path);
    while (fin.peek() != EOF) {
        getline(fin, line);
        if (line.empty()) break;
        ln += 1;
    }

    fin.close();
    return ln;
}


void MultilayerGraph::ParseConfFile(const string &conf_file, vector<string> &intra_layer_file,
                                    vector<bg_file> &cross_layer_file) {
    int id1, id2;
    size_t sep1, sep2;
    string line;
    bool read_intra_layer = true;

    auto fin = ifstream(conf_file);
    while (fin.peek() != EOF) {

        getline(fin, line);

        if (line.empty()) {
            if (!read_intra_layer) break;
            read_intra_layer = false;
            continue;
        }

        if (read_intra_layer) {
            intra_layer_file.emplace_back(line);
        } else {
            sep1 = line.find_first_of(',', 0);
            id1 = std::stoi(line.substr(0, sep1));

            sep2 = line.find_first_of(',', sep1 + 1);
            id2 = std::stoi(line.substr(sep1 + 1, sep2 - sep1 - 1));

            cross_layer_file.emplace_back(id1, id2, line.substr(sep2 + 1, line.size() - sep2 - 1));
        }
    }

    fin.close();
}


//void MultilayerGraph::SetGraphOrdering(ordering order) {
//
//    if (order == RAND) return;
//
//    // Re-ordering graph layers.
//    int num_of_slayers = num_of_layers - 1;
//    float avg_density[num_of_slayers];
//
//    if (order == INTRA_DEN_INC or order == INTRA_DEN_DEC) {
//        for (int i = 0; i < num_of_slayers; i++)
//            avg_density[i] = (float) graph_layers[i]->GetM() / (float) graph_layers[i]->GetN();
//
//        if (order == INTRA_DEN_INC)
//            sort(graph_layer_indices, graph_layer_indices + num_of_slayers,
//                 [&avg_density](int id1, int id2) { return avg_density[id1] < avg_density[id2]; });
//        else
//            sort(graph_layer_indices, graph_layer_indices + num_of_slayers,
//                 [&avg_density](int id1, int id2) { return avg_density[id1] > avg_density[id2]; });
//
//    } else if (order == CROSS_DEN_INC or order == CROSS_DEN_DEC) {
//        for (int i = 0; i < num_of_slayers; i++) {
//            // divide the number of primary graph layer vertices
//            avg_density[i] = (float) cross_layer_graphs[num_of_slayers][i]->GetM() /
//                    (float) cross_layer_graphs[num_of_slayers][i]->GetN();
//        }
//
//        if (order == CROSS_DEN_INC)
//            sort(graph_layer_indices, graph_layer_indices + num_of_slayers,
//                 [&avg_density](int id1, int id2) { return avg_density[id1] < avg_density[id2]; });
//        else
//            sort(graph_layer_indices, graph_layer_indices + num_of_slayers,
//                 [&avg_density](int id1, int id2) { return avg_density[id1] > avg_density[id2]; });
//    }
//}
//
//void MultilayerGraph::PrintGraphOrdering() {
//    cout << "Graph ordering:" << endl;
//    for (int i = 0; i < num_of_layers; i++) {
//        cout << graph_layers[graph_layer_indices[i]]->GetGraphFile() << " ==> " << "Graph " << i << endl;
//    }
//}