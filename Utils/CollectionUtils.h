//
// Created by ldd on 2021/10/8.
//

#ifndef CSM4GMG_COLLECTIONUTILS_H
#define CSM4GMG_COLLECTIONUTILS_H

#include "../Header.h"

static void GetCombinations(int idx, int *vec, int len, const int *bnds, int step, vector<vector<int>> &comb) {
    comb.emplace_back(vec, vec + len);
    for (int i = idx; i < len; i++) {
        vec[i] += step;
        if (vec[i] < bnds[i]) GetCombinations(i, vec, len, bnds, step, comb);
        vec[i] -= step;
    }
}


#endif //CSM4GMG_COLLECTIONUTILS_H
