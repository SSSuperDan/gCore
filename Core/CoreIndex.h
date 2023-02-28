//
// Created by ldd on 2022/7/2.
//

#ifndef CSM4GMG_COREINDEX_H
#define CSM4GMG_COREINDEX_H

#include "../Header.h"

/*
 * ============= CoreIndex =============
 * Integer n represents the total number of vertices.
 * Array vert contains all vertices.
 * Array pos contains the positions of vertices in array vert.
 * Offsets s and t divide vert into 3 parts:
 * |***********|****************|**************************|
 * 0-----------s----------------e--------------------------n
 * Part 1: vert[0:s) contains all discarded vertices.
 * Part 2: vert[s:e) contains all inactive vertices.
 * Part 3: vert[e:n) contains all active vertices.
 */

struct CoreIndex {
    int n{0};
    int s{0};
    int e{0};
    int *vert{nullptr};
    int *pos{nullptr};

    CoreIndex() = default;

    ~CoreIndex() {
        delete[] vert;
        delete[] pos;
    }

    void Init(int n_) {
        n = n_;
        vert = new int[n];
        pos = new int[n];
    }

    // Init CoreIndex with n vertices.
    // All vertices are active initially.
    void Set() {
        s = 0;
        e = 0;
        for (int i = 0; i < n; i++) {
            vert[i] = i;
            pos[i] = i;
        }
    }

    // Move v from Part 3 to Part 2, i.e., convert v from active state to inactive state.
    inline void Remove(int v) {
        int v_pos;

        v_pos = pos[v];
        vert[v_pos] = vert[e];
        vert[e] = v;
        pos[vert[v_pos]] = v_pos;
        pos[v] = e++;
    }
};

#endif //CSM4GMG_COREINDEX_H
