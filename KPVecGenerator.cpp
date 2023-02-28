//
// Created by ldd on 2022/10/14.
//


#include "KPVecGenerator.h"

int KPVecGenerator::LoadKVectors(int dim, const string &file, vector<vector<int>> &k_vectors) {
    int size;

    auto f = ifstream(file);
    f.read(reinterpret_cast<char *> (&size), sizeof(int));

    k_vectors.resize(size);
    for (int i = 0; i < size; i++) {
        k_vectors[i].resize(dim);
        f.read(reinterpret_cast<char *>(k_vectors[i].data()), long(dim * sizeof(int)));
    }

    return size;
}

int KPVecGenerator::LoadPVectors(int dim, const string &file, vector<vector<float>> &p_vectors) {
    int size;

    auto f = ifstream(file);
    f.read(reinterpret_cast<char *> (&size), sizeof(int));

    p_vectors.resize(size);
    for (int i = 0; i < size; i++) {
        p_vectors[i].resize(dim);
        f.read(reinterpret_cast<char *>(p_vectors[i].data()), long(dim * sizeof(float)));
    }

    return size;
}

int KPVecGenerator::LoadKPVectors(int dim, const string &file, vector<vector<int>> &k_vectors,
                                  vector<vector<float>> &p_vectors) {
    int size, p_dim = dim - 1;

    auto f = ifstream(file);
    f.read(reinterpret_cast<char *> (&size), sizeof(int));

    k_vectors.resize(size);
    p_vectors.resize(size);
    for (int i = 0; i < size; i++) {
        k_vectors[i].resize(dim);
        p_vectors[i].resize(p_dim);

        f.read(reinterpret_cast<char *>(k_vectors[i].data()), long(dim * sizeof(int)));
        f.read(reinterpret_cast<char *>(p_vectors[i].data()), long(p_dim * sizeof(float)));
    }

    return size;
}

void KPVecGenerator::GenRandomKVectors(MultilayerGraph &mg, int size, const string &output, const vector<int> &end_k) {

    int ln = mg.GetLayerNumber(), end_k_vec[ln], k_vec[ln];

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis[ln];

    GetDegeneracy(mg, end_k_vec);
    if (!end_k.empty()) for (int i = 0; i < ln; i++) end_k_vec[i] = std::min(end_k[i], end_k_vec[i]);
    for (int i = 0; i < ln; i++) dis[i] = std::uniform_int_distribution<int>{0, end_k_vec[i]};

    auto file = output + mg.GetGraphName() + "_sampled_k_vec_" + to_string(size);
    auto f = ofstream(file);
    f.write(reinterpret_cast<char *>(&size), sizeof(int));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < ln; j++) k_vec[j] = dis[j](gen);
        f.write(reinterpret_cast<char *>(k_vec), long(ln * sizeof(int)));
    }
    f.close();

    cout << "Random coreness vectors (k) have been generated, size = " << size << "." << endl;
}

void KPVecGenerator::GenRandomPVectors(int dim, int size, const string &output) {

    float p_vec[dim];

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<float> dis[dim];

    for (int i = 0; i < dim; i++) dis[i] = std::uniform_real_distribution<float>{0, 1};

    auto file = output + "dim" + to_string(dim) + "_sampled_p_vec_" + to_string(size);
    auto f = ofstream(file);
    f.write(reinterpret_cast<char *>(&size), sizeof(int));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < dim; j++) p_vec[j] = dis[j](gen);
        f.write(reinterpret_cast<char *>(p_vec), long(dim * sizeof(float)));
    }
    f.close();

    cout << "Random fraction vectors (p) have been generated, size = " << size << "." << endl;
}

void KPVecGenerator::GenRandomKPVectors(MultilayerGraph &mg, int size, const string &output, const vector<int> &end_k) {
    int ln = mg.GetLayerNumber(), dim = ln - 1, end_k_vec[ln], k_vec[ln];
    float p_vec[dim];

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> k_dis[ln];
    uniform_real_distribution<float> p_dis[dim];

    GetDegeneracy(mg, end_k_vec);
    if (!end_k.empty()) for (int i = 0; i < ln; i++) end_k_vec[i] = std::min(end_k[i], end_k_vec[i]);
    for (int i = 0; i < ln; i++) k_dis[i] = std::uniform_int_distribution<int>{0, end_k_vec[i]};
    for (int i = 0; i < dim; i++) p_dis[i] = std::uniform_real_distribution<float>{0, 1};

    auto file = output + mg.GetGraphName() + "_sampled_kp_vec_" + to_string(size);
    auto f = ofstream(file);
    f.write(reinterpret_cast<char *>(&size), sizeof(int));
    for (int i = 0; i < size; i++) {

        for (int j = 0; j < ln; j++) k_vec[j] = k_dis[j](gen);
        for (int j = 0; j < dim; j++) p_vec[j] = p_dis[j](gen);

        f.write(reinterpret_cast<char *>(k_vec), long(ln * sizeof(int)));
        f.write(reinterpret_cast<char *>(p_vec), long(dim * sizeof(float)));
    }
    f.close();

    cout << "Random (k,p) pairs have been generated, size = " << size << "." << endl;
}

void KPVecGenerator::GenRandomKPVectors(MultilayerGraph &mg, const string &sampled_k_vec_file, int size,
                                        const string &output) {
    int ln = mg.GetLayerNumber(), dim = ln - 1, id;
    float p_vec[dim];

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> k_dis;
    uniform_real_distribution<float> p_dis[dim];

    vector<vector<int>> sampled_k_vec;
    LoadKVectors(ln, sampled_k_vec_file, sampled_k_vec);

    k_dis = std::uniform_int_distribution<int>{0, (int) sampled_k_vec.size() - 1};
    for (int i = 0; i < dim; i++) p_dis[i] = std::uniform_real_distribution<float>{0, 1};

    auto file = output + mg.GetGraphName() + "_sampled_kp_vec_" + to_string(size);
    auto f = ofstream(file);
    f.write(reinterpret_cast<char *>(&size), sizeof(int));
    for (int i = 0; i < size; i++) {

        id = k_dis(gen);
        for (int j = 0; j < dim; j++) p_vec[j] = p_dis[j](gen);

        f.write(reinterpret_cast<char *>(sampled_k_vec[id].data()), long(ln * sizeof(int)));
        f.write(reinterpret_cast<char *>(p_vec), long(dim * sizeof(float)));
    }
    f.close();

    cout << "Random (k,p) pairs have been generated, size = " << size << "." << endl;
}

void KPVecGenerator::GetDegeneracy(MultilayerGraph &mg, int *degeneracy) {
    int ln = mg.GetLayerNumber();
    for (int i = 0; i < ln; i++) {
        degeneracy[i] = CoreDec::GetDegeneracy(mg.GetGraph(i));
    }
}