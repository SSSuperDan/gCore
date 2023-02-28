//
// Created by ldd on 2021/10/11.
//

#ifndef CSM4GMG_KPTREE_H
#define CSM4GMG_KPTREE_H


#include "PTree.h"

const int HASH_BUCKET_SIZE = (1 << 10); // 1024

/*
 * hash node alignment: | === k === | === KNode === |
 */
struct HashNode;

struct KNode {
    PTree *p_tree{};
    HashNode *next_node{nullptr};
};


struct KPTree {
public:
    explicit KPTree(int dim_, int bucket_size_ = HASH_BUCKET_SIZE);
    ~KPTree();

    [[nodiscard]] size_t KHashNodeSize() const {
        return h_size;
    }

    [[nodiscard]] int GetNumOfNodes() const {
        return num_of_k_nodes;
    }

    inline void Insert(int *k, PTree *p_tree) {
        int buk_id;
        KNode *k_node;
        HashNode *new_node;

        buk_id = range_hash_k(k);
        new_node = Allocate();
        memcpy(new_node, k, k_size);

        k_node = get_k_node_addr(new_node);
        k_node->p_tree = p_tree;
        k_node->next_node = bucket[buk_id];
        bucket[buk_id] = new_node;

        num_of_k_nodes++;
    }

    int Search(const vector<int> &k_vec, const vector<int> &p_vec, int *core, int &length);
    long long int GetNumOfPNodes();

    void Flush(const string & file);
    static KPTree* Load(const string & file);

    bool Equals(KPTree & kp_tree);

    void GetKFrontier(const vector<int> &p_vec, vector<vector<int>> &k_frontiers);
    void PrintBucket();

private:
    MDataBuf<char> hash_node_buf;
    HashNode **bucket;

    int dim;
    int bucket_size;
    int num_of_k_nodes;

    size_t k_size; // size of k_vec, offset of k_node
    size_t h_size; // size of hash node

    inline HashNode *Allocate() {
        return reinterpret_cast<HashNode *>(hash_node_buf.Allocate(h_size));
    }

    // boost/functional/hash/hash.hpp
    inline int range_hash_k(const int *k) const {
        std::size_t seed = 0;
        hash<int> hash;

        for (int i = 0; i < dim; i++) {
            seed ^= hash(k[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        return (int) (seed % bucket_size);
    }

    static inline int *get_k_vec_addr(HashNode *hash_node) {
        return reinterpret_cast<int *> (hash_node);
    }

    inline KNode *get_k_node_addr(HashNode *hash_node) const {
        return reinterpret_cast<KNode *> (reinterpret_cast<char *> (hash_node) + k_size);
    }

    inline PTree *Find(const int *k);
    void Load(std::basic_ifstream<char> &f);
    bool ExploreFrontiers(int *k_vec, int lk, const vector<int> &p_vec, vector<vector<int>> &k_frontiers);
};

#endif //CSM4GMG_KPTREE_H
