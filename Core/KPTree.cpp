//
// Created by ldd on 2022/6/9.
//

#include "KPTree.h"


KPTree::KPTree(int dim_, int bucket_size_) : dim(dim_), bucket_size(bucket_size_), num_of_k_nodes(0) {
    k_size = dim * sizeof(int);
    h_size = k_size + sizeof(KNode);  // size of the whole hash node
    bucket = new HashNode *[bucket_size];

    for (int i = 0; i < bucket_size; i++) bucket[i] = nullptr;
}

KPTree::~KPTree() {
    KNode *k_node;
    HashNode *hash_node;

    for (int i = 0; i < bucket_size; i++) {
        hash_node = bucket[i];
        while (hash_node) {
            k_node = get_k_node_addr(hash_node);
            delete k_node->p_tree;
            hash_node = k_node->next_node;
        }
    }

    delete[] bucket;
}

int KPTree::Search(const vector<int> &k_vec, const vector<int> &p_vec, int *core, int &length) {
    PTree *p_tree = Find(k_vec.data());
    if (p_tree) return p_tree->Search(p_vec, core, length);
    else return 0;
}

inline PTree *KPTree::Find(const int *k) {
    KNode *k_node;
    HashNode *hash_node;

    hash_node = bucket[range_hash_k(k)];
    while (hash_node) {
        k_node = get_k_node_addr(hash_node);
        if (!memcmp(hash_node, k, k_size)) {
            return k_node->p_tree;
        } else hash_node = k_node->next_node;
    }
    return nullptr;
}

long long int KPTree::GetNumOfPNodes() {
    long long int cnt = 0;
    KNode *k_node;
    HashNode *hash_node;

    for (int i = 0; i < bucket_size; i++) {
        hash_node = bucket[i];
        while (hash_node) {
            k_node = get_k_node_addr(hash_node);
            cnt += k_node->p_tree->GetNumOfNodes();
            hash_node = k_node->next_node;
        }
    }

    return cnt;
}

void KPTree::Flush(const string &file) {
    KNode *k_node;
    HashNode *hash_node;

    auto f = std::ofstream(file);

    f.write(reinterpret_cast<char *> (&dim), sizeof(int));
    f.write(reinterpret_cast<char *> (&bucket_size), sizeof(int));
    f.write(reinterpret_cast<char *> (&num_of_k_nodes), sizeof(int));

    for (int i = 0; i < bucket_size; i++) {
        hash_node = bucket[i];
        while (hash_node) {
            f.write(reinterpret_cast<char *> (hash_node), (long) k_size);
            k_node = get_k_node_addr(hash_node);
            k_node->p_tree->Flush(f);
            hash_node = k_node->next_node;
        }
    }

    f.close();
}

KPTree *KPTree::Load(const string &file) {
    int dimension, buk_size;

    auto f = std::ifstream(file);

    f.read(reinterpret_cast<char *> (&dimension), sizeof(int));
    f.read(reinterpret_cast<char *> (&buk_size), sizeof(int));

    KPTree *tree = new KPTree(dimension, buk_size);
    tree->Load(f);
    f.close();

    return tree;

}

void KPTree::Load(std::basic_ifstream<char> &f) {
    int k_vec[dim], size;

    f.read(reinterpret_cast<char *> (&size), sizeof(int));

    for (int i = 0; i < size; i++) {
        f.read(reinterpret_cast<char *> (&k_vec[0]), (long) k_size);
        PTree *p_tree = new PTree(dim - 1);
        p_tree->Load(f);
        Insert(k_vec, p_tree);
    }
}

bool KPTree::Equals(KPTree &kp_tree) {
    int size;
    KNode *k_node;
    HashNode *hash_node;

    auto k_vec_cmp = [this](pair<int *, PTree *> &kpp1, pair<int *, PTree *> &kpp2) {
        int *k1 = kpp1.first, *k2 = kpp2.first;
        for (int i = 0; i < dim; i++) {
            if (k1[i] < k2[i]) return true;
            else if (k1[i] > k2[i]) return false;
        }
        return false;
    };

    // Equals basic elements
    if (dim != kp_tree.dim || bucket_size != kp_tree.bucket_size || num_of_k_nodes != kp_tree.num_of_k_nodes) {
        return false;
    }

    for (int i = 0; i < bucket_size; i++) {
        vector<pair<int *, PTree *>> kp_pair, kp_tree_kp_pair;

        hash_node = bucket[i];
        while (hash_node) {
            k_node = get_k_node_addr(hash_node);
            kp_pair.emplace_back(get_k_vec_addr(hash_node), k_node->p_tree);
            hash_node = k_node->next_node;
        }

        hash_node = kp_tree.bucket[i];
        while (hash_node) {
            k_node = get_k_node_addr(hash_node);
            kp_pair.emplace_back(get_k_vec_addr(hash_node), k_node->p_tree);
            hash_node = k_node->next_node;
        }

        if (kp_pair.size() != kp_tree_kp_pair.size()) return false;
        if (kp_pair.empty()) continue;

        // sort and compare
        sort(kp_pair.begin(), kp_pair.end(), k_vec_cmp);
        sort(kp_tree_kp_pair.begin(), kp_tree_kp_pair.end(), k_vec_cmp);

        size = (int) kp_pair.size();
        for (int j = 0; j < size; j++) {
            if (memcmp(kp_pair[i].first, kp_tree_kp_pair[i].first, k_size) != 0 ||
                !kp_pair[i].second->Equals(*kp_tree_kp_pair[i].second)) {
                return false;
            }
        }
    }
    return true;
}

void KPTree::GetKFrontier(const vector<int> &p_vec, vector<vector<int>> &k_frontiers) {
    int k_vec[dim];
    memset(k_vec, 0, k_size);
    ExploreFrontiers(k_vec, 0, p_vec, k_frontiers);
}

bool KPTree::ExploreFrontiers(int *k_vec, int lk, const vector<int> &p_vec, vector<vector<int>> &k_frontiers) {
    bool has_frontier;
    PTree *p_tree = Find(k_vec);

    if (!p_tree || !p_tree->Search(p_vec)) return false;
    else {
        has_frontier = false;
        for (int i = lk; i < dim; i++) {
            k_vec[i]++;
            has_frontier |= ExploreFrontiers(k_vec, i, p_vec, k_frontiers);
            k_vec[i]--;
        }

        if (!has_frontier) k_frontiers.emplace_back(k_vec, k_vec + dim);
        return true;
    }
}

void KPTree::PrintBucket() {
    KNode *k_node;
    HashNode *hash_node;

    cout << "KPTree hashtable bucket layout:" << endl;
    for (int i = 0; i < bucket_size; i++) {
        if (bucket[i]) {
            cout << "Bucket " << i << ":";
            hash_node = bucket[i];
            while (hash_node) {
                k_node = get_k_node_addr(hash_node);
                std::cout << Array2String(get_k_vec_addr(hash_node), dim) << "->";
                hash_node = k_node->next_node;
            }
            cout << endl;
        }
    }
}
