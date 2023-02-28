//
// Created by ldd on 2021/10/5.
//


#ifndef CSM4GMG_PTREE_H
#define CSM4GMG_PTREE_H

#include "../Structures/DataBuf.h"
#include "../Utils/ArrayUtils.h"
#include "../Utils/Viz.h"

/*
 *  p-node alignment: | === p === | === hp === | last(p) |       branch 1      |       branch 2      | ...
 *  p-node alignment: | === p === | === hp === | last(p) | nov | p_vtx | p_chl | nov | p_vtx | p_chl | ...
 */

struct PNode;

struct Branch {
    int num_of_vtx{};  // size of difference
    int *vtx{nullptr};  // difference
    PNode *child{nullptr};

    Branch() = default;
    ~Branch() = default;

    inline void SetChild(PNode *child_) {
        child = child_;
    }
};

class PTree {
public:
    explicit PTree(int dim_);
    ~PTree() = default;

    [[nodiscard]] size_t PNodeSize(int lp) const {
        return GetPNodeSize(lp);
    }

    [[nodiscard]] long long int GetNumOfNodes() const {
        return num_of_nodes;
    }

    inline PNode *GetNewPTreeNode(int *p_, int lp) {
        PNode *new_node;

        new_node = Allocate(GetPNodeSize(lp));
        memcpy(new_node, p_, dim * sizeof(int));
        *get_lp_addr(new_node) = lp;

        num_of_nodes++;
        return new_node;
    }

    static inline int *get_p(PNode *p_node) {
        return reinterpret_cast<int *>(p_node);
    }

    inline int *get_hp(PNode *p_node) const {
        return reinterpret_cast<int *>(reinterpret_cast<char *> (p_node) + hp_offset);
    }

    inline Branch *get_relative_branch(PNode *p_node, int offset) const {
        offset -= *get_lp_addr(p_node);
        return reinterpret_cast<Branch *>(reinterpret_cast<char *> (p_node) + bc_offset + offset * sizeof(Branch));
    }

    inline int *get_lp_addr(PNode *p_node) const {
        return reinterpret_cast<int *>(reinterpret_cast<char *> (p_node) + lp_offset);
    }

    void SetRootNode(PNode *r) {
        root = r;
    }

    inline void SetBranch(Branch *branch, int *vtx, int length) {
        branch->num_of_vtx = length;
        branch->vtx = difference_buf.Allocate(length);
        memcpy(branch->vtx, vtx, length * sizeof(int));
    }

    PNode *Search(const vector<int> &p_vec);
    int Search(const vector<int> &p_vec, int *core, int &length);

    void Flush(const string &file);
    void Flush(std::basic_ofstream<char> &f);
    void Load(std::basic_ifstream<char> &f);
    static PTree* Load(const string &file);

    bool Equals(PTree & p_tree);

    int GetBoundary(int layer);
    void GetBoundary(vector<vector<int>> &boundaries);

    void PrintSearchPath(const vector<int> &p_vec);
    void Visualize(bool vis_detail, const vector<int> &p_vec = {});
    void DumpPTree();

private:

    MDataBuf<char> tree_node_buf;
    DataBuf<int> difference_buf;

    PNode *root;
    int dim;
    long long int num_of_nodes;

    size_t hp_offset;  // offset of hp
    size_t lp_offset;  // offset of the value of last(p)
    size_t bc_offset;  // offset of branches

    [[nodiscard]] inline size_t GetPNodeSize(int lp) const {
        return bc_offset + (dim - lp) * sizeof(Branch);
    }

    inline PNode *Allocate(size_t p_size) {
        return reinterpret_cast<PNode *>(tree_node_buf.Allocate(p_size));
    }

    inline Branch *get_leftmost_branch(PNode *p_node) const {
        return reinterpret_cast<Branch *>(reinterpret_cast<char *> (p_node) + bc_offset);
    }

    inline Branch *get_absolute_branch(PNode *p_node, int offset) const {
        return reinterpret_cast<Branch *>(reinterpret_cast<char *> (p_node) + bc_offset + offset * sizeof(Branch));
    }

    inline Branch *get_rightmost_branch(PNode *p_node) const {
        int lp = *get_lp_addr(p_node);
        return reinterpret_cast<Branch *>(reinterpret_cast<char *> (p_node) + bc_offset +
                                          (dim - lp - 1) * sizeof(Branch));
    }

    bool IsEqual(PNode *node, PNode * p_tree_node);
    bool IsEqual(Branch *branch, Branch * p_tree_branch);
    void FlushPNode(PNode *node, unordered_map<PNode *, int> &p_node2id, std::basic_ofstream<char> &f) const;
    void LoadPNode(PNode *node, unordered_map<int, PNode *> &id2p_node, std::basic_ifstream<char> &f);
    void ExploreBoundaries(PNode *node, int *p_vec, vector<vector<int>> &boundaries) const;
    void Traverse(PNode *node, unordered_map<PNode *, int> &p_node2id);
    void PTree2DotFile(unordered_map<PNode *, int> &p_node2id, bool vis_detail) const;
    void PTree2EdgeFile(unordered_map<PNode *, int> &p_node2id) const;
};

#endif //CSM4GMG_PTREE_H

