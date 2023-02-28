//
// Created by ldd on 2022/6/6.
//

#include "PTree.h"

PTree::PTree(int dim_) : root(nullptr), dim(dim_), num_of_nodes(0) {
    hp_offset = dim_ * sizeof(int);
    lp_offset = hp_offset << 1;
    bc_offset = lp_offset + sizeof(int);
}

PNode *PTree::Search(const vector<int> &p_vec) {
    int i = 0;
    PNode *node = root;

    while (true) {
        while (node && get_p(node)[i] < p_vec[i]) node = get_leftmost_branch(node)->child;
        if (!node) break;
        else {
            i++;
            while (i < dim && get_p(node)[i] >= p_vec[i]) i++;
            if (i == dim) break;
            else node = get_relative_branch(node, i)->child;
        }
    }
    return node;
}

// Retrieve (k,p)-core given fraction vector p_vec.
int PTree::Search(const vector<int> &p_vec, int *core, int &length) {
    int i, num_of_visits;
    PNode *node;
    Branch *branch;

    i = 0;
    node = root;
    num_of_visits = 1;  // visiting root

    // locate (k,p)-node
    while (true) {
        while (node && get_p(node)[i] < p_vec[i]) {
            node = get_leftmost_branch(node)->child;
            num_of_visits++;
        }
        if (!node) break;
        else {
            i++;
            while (i < dim && get_p(node)[i] >= p_vec[i]) i++;
            if (i == dim) break;
            else {
                node = get_relative_branch(node, i)->child;
                num_of_visits++;
            }
        }
    }

    length = 0;
    while (node) {
        branch = get_rightmost_branch(node);  // collect vertices along the right-most path.
        memcpy(core + length, branch->vtx, branch->num_of_vtx * sizeof(int));
        length += branch->num_of_vtx;
        node = branch->child;
    }
    return num_of_visits;
}

void PTree::Flush(const string &file) {
    auto f = std::ofstream(file);

    f.write(reinterpret_cast<char *> (&dim), sizeof(int));
    Flush(f);
    f.close();
}

void PTree::Flush(std::basic_ofstream<char> &f) {
    unordered_map<PNode *, int> p_node2id;

    p_node2id.emplace(root, 0);
    FlushPNode(root, p_node2id, f);
}

void PTree::FlushPNode(PNode *node, unordered_map<PNode *, int> &p_node2id, std::basic_ofstream<char> &f) const {
    Branch *branch;
    PNode *child;
    int child_id;

    f.write(reinterpret_cast<char *> (node), (long) bc_offset);

    for (int i = *get_lp_addr(node); i < dim; i++) {
        branch = get_relative_branch(node, i);
        child = branch->child;

        f.write(reinterpret_cast<char *>(&branch->num_of_vtx), sizeof(int));
        if (branch->num_of_vtx) {
            f.write(reinterpret_cast<char *> (branch->vtx), long(branch->num_of_vtx * sizeof(int)));
        }

        child_id = -1;
        if (child) {
            auto iter = p_node2id.find(child);
            if (iter == p_node2id.end()) {
                child_id = (int) p_node2id.size();
                f.write(reinterpret_cast<char *>(&child_id), sizeof(int));
                p_node2id.emplace(child, child_id);
                FlushPNode(child, p_node2id, f);
            } else {
                f.write(reinterpret_cast<char *>(&iter->second), sizeof(int));
            }
        } else f.write(reinterpret_cast<char *>(&child_id), sizeof(int));;
    }
}


PTree *PTree::Load(const string &file) {
    int dimension;
    auto f = std::ifstream(file);

    f.read(reinterpret_cast<char *> (&dimension), sizeof(int));
    PTree *tree = new PTree(dimension);
    tree->Load(f);
    f.close();
    return tree;
}

void PTree::Load(std::basic_ifstream<char> &f) {
    unordered_map<int, PNode *> id2p_node;

    root = Allocate(GetPNodeSize(0));
    num_of_nodes++;
    id2p_node.emplace(0, root);

    LoadPNode(root, id2p_node, f);
}

void PTree::LoadPNode(PNode *node, unordered_map<int, PNode *> &id2p_node, std::basic_ifstream<char> &f) {
    Branch *branch;
    PNode *child;
    int child_id;

    f.read(reinterpret_cast<char *> (node), (long) bc_offset);

    for (int i = *get_lp_addr(node); i < dim; i++) {
        // load each branch

        branch = get_relative_branch(node, i);
        f.read(reinterpret_cast<char *>(&branch->num_of_vtx), sizeof(int));
        if (branch->num_of_vtx) {
            branch->vtx = difference_buf.Allocate(branch->num_of_vtx);
            f.read(reinterpret_cast<char *> (branch->vtx), long(branch->num_of_vtx * sizeof(int)));
        }

        f.read(reinterpret_cast<char *>(&child_id), sizeof(int));

        if (child_id != -1) {
            auto iter = id2p_node.find(child_id);
            if (iter == id2p_node.end()) {
                child = Allocate(GetPNodeSize(i));
                num_of_nodes++;

                branch->child = child;
                id2p_node.emplace(child_id, child);

                LoadPNode(child, id2p_node, f);
            } else branch->child = iter->second;
        }
    }
}

bool PTree::Equals(PTree &p_tree) {
    // Equals basic elements
    if (dim != p_tree.dim || num_of_nodes != p_tree.num_of_nodes) {
        return false;
    }

    return IsEqual(root, p_tree.root);
}

bool PTree::IsEqual(PNode *node, PNode *p_tree_node) {
    int num_of_branches;

    if (node == nullptr && p_tree_node == nullptr) return true;
    else if (node == nullptr || p_tree_node == nullptr) return false;

    // Equals p
    if (!ArrEq(get_p(node), get_p(p_tree_node), dim)) {
        return false;
    }

    // Equals hp
    if (!ArrEq(get_hp(node), get_hp(p_tree_node), dim)) {
        return false;
    }

    // Equals last(p)
    if ((*get_lp_addr(node)) != (*get_lp_addr(p_tree_node))) {
        return false;
    }

    num_of_branches = dim - (*get_lp_addr(node));
    for (int i = 0; i < num_of_branches; i++) {
        if (!IsEqual(get_absolute_branch(node, i), get_absolute_branch(p_tree_node, i))) {
            return false;
        }
    }
    return true;
}

bool PTree::IsEqual(Branch *branch, Branch *p_tree_branch) {
    // Equals nov
    if (branch->num_of_vtx != p_tree_branch->num_of_vtx) {
        return false;
    }

    // Equals p_vtx
    if (!ArrEq(branch->vtx, p_tree_branch->vtx, branch->num_of_vtx)) {
        return false;
    }

    // Equals child
    return IsEqual(branch->child, p_tree_branch->child);
}


int PTree::GetBoundary(int layer) {
    PNode *pre_node = nullptr, *node = root;

    while (node) {
        pre_node = node;
        node = get_relative_branch(node, layer)->child;
    }

    if (pre_node) return get_p(pre_node)[layer];
    return -1;
}

void PTree::GetBoundary(vector<vector<int>> &boundaries) {
    int p_vec[dim];
    memset(p_vec, 0, dim * sizeof(int));
    ExploreBoundaries(root, p_vec, boundaries);
}

void PTree::ExploreBoundaries(PNode *node, int *p_vec, vector<vector<int>> &boundaries) const {
    PNode *child;
    int lp = *get_lp_addr(node), num_of_branch = dim - lp, old_p;
    bool has_child = false;

    old_p = p_vec[lp];
    p_vec[lp] = get_p(node)[lp];

    for (int i = 0; i < num_of_branch; i++) {
        child = get_absolute_branch(node, i)->child;
        if (child) {
            has_child = true;
            ExploreBoundaries(child, p_vec, boundaries);
        }
    }

    if (!has_child) {
        boundaries.emplace_back(p_vec, p_vec + dim);
    }

    p_vec[lp] = old_p;
}

// Print search path of locating the node representing fraction vector p_vec.
void PTree::PrintSearchPath(const vector<int> &p_vec) {
    int i;
    PNode *node;

    i = 0;
    node = root;

    // locate (k,p)-node
    while (true) {
        while (node && get_p(node)[i] < p_vec[i]) {
            cout << Array2String(get_p(node), dim) << ", match prefix of length" << i << "." << endl;
            node = get_leftmost_branch(node)->child;
        }
        if (!node) break;
        else {
            i++;
            while (i < dim && get_p(node)[i] >= p_vec[i]) i++;
            if (i == dim) break;
            else {
                cout << Array2String(get_p(node), dim) << ", match prefix of length" << i << "." << endl;
                node = get_relative_branch(node, i)->child;
            }
        }
    }
    if (node) cout << "Target node = " << Array2String(get_p(node), dim) << "." << endl;
    else cout << "No target node found." << endl;
}

// Visualize PTree with GraphViz.
// If p_vec is not null, only visualize the subtree rooted by node with fraction vector p_vec.
void PTree::Visualize(bool vis_detail, const vector<int> &p_vec) {
    PNode *node;
    unordered_map<PNode *, int> p_node2id;

    node = p_vec.empty() ? root : Search(p_vec);
    if (node) {
        Traverse(node, p_node2id);
        PTree2DotFile(p_node2id, vis_detail);
        RunViz("p_tree");
    } else {
        cerr << "No valid root found." << endl;
        exit(-1);
    }
}

// Save PTree to disk.
void PTree::DumpPTree() {
    unordered_map<PNode *, int> p_node2id;
    Traverse(root, p_node2id);
    PTree2EdgeFile(p_node2id);
}

// Traverse PTree, collect all distinct nodes.
void PTree::Traverse(PNode *node, unordered_map<PNode *, int> &p_node2id) {

    PNode *child;
    p_node2id[node] = (int) p_node2id.size() + 1;

    for (int i = *get_lp_addr(node); i < dim; i++) {
        child = get_relative_branch(node, i)->child;
        if (child && p_node2id.find(child) == p_node2id.end()) {
            Traverse(child, p_node2id);
        }
    }
}

// Generate .dot file for visualizing PTree with GraphViz.
void PTree::PTree2DotFile(unordered_map<PNode *, int> &p_node2id, bool vis_detail) const {
    PNode *node;
    Branch *branch;
    int lp, num_of_leaves = 1;

    auto p_tree_out = ofstream("TmpFile/p_tree.dot");
    p_tree_out << "digraph {" << endl;

    // write node labels
    for (const auto &node_id : p_node2id) {
        p_tree_out << node_id.second << "[ label = \"(" << node_id.second << ") " << Array2String(get_p(node_id.first), dim)
                   << "\" ]; \n";
    }

    // write edge labels
    if (vis_detail) {
        for (const auto &node_id : p_node2id) {
            node = node_id.first;
            lp = *get_lp_addr(node);
            for (int i = lp; i < dim; i++) {
                branch = get_relative_branch(node, i);

                if (branch->child) {  // not leaf
                    p_tree_out << node_id.second << " -> " << p_node2id[branch->child] <<
                               "[ label = \" B" << i << ":"
                               << Array2String(branch->vtx, branch->num_of_vtx)
                               << "\" ]; \n";
                } else {
                    p_tree_out << node_id.second << " -> " << "leaf_" << num_of_leaves++ <<
                               "[ label = \" B" << i << ":"
                              << Array2String(branch->vtx, branch->num_of_vtx)
                              << "\" ]; \n";
                }
            }
        }
    } else {
        for (const auto &node_id : p_node2id) {
            node = node_id.first;
            lp = *get_lp_addr(node);
            for (int i = lp; i < dim; i++) {
                branch = get_relative_branch(node, i);

                if (branch->child) {  // not leaf
                    p_tree_out << node_id.second << " -> " << p_node2id[branch->child] <<
                               "[ label = \" B" << i << ":"
                               << branch->num_of_vtx << "\" ]; \n";
                } else {
                    p_tree_out << node_id.second << " -> " << "leaf_" << num_of_leaves++ <<
                               "[ label = \" B" << i << ":"
                              << branch->num_of_vtx << "\" ]; \n";
                }
            }
        }
    }

    p_tree_out << "}" << endl;
    p_tree_out.close();
}

// Generate .txt file saving all PTree edges.
void PTree::PTree2EdgeFile(unordered_map<PNode *, int> &p_node2id) const {
    PNode *node;
    Branch *branch;

    auto p_tree_out = ofstream("TmpFile/p_tree.txt");

    for (const auto &node_id : p_node2id) {
        node = node_id.first;
        for (int i = *get_lp_addr(node); i < dim; i++) {
            branch = get_relative_branch(node, i);

            if (branch->child) {  // not leaf
                p_tree_out << Array2String(get_p(node), dim) << " -> "
                           << Array2String(get_p(branch->child), dim) << endl;
            } else {
                p_tree_out << Array2String(get_p(node), dim) << " -> Leaf " << endl;
            }
        }
    }

    p_tree_out.close();
}
