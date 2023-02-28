//
// Created by ldd on 2021/10/5.
//

#ifndef CSM4GMG_FRACUTILS_H
#define CSM4GMG_FRACUTILS_H

#include "../Header.h"

static bool IsNotEqual(int d1, int n1, int d2, int n2) {
    return n2 * d1 != d2 * n1;
}

static void FracBinSort(Frac *&fractions, int max_den, int frac_size) {
    int num, den, loc, cnt[max_den + 1], index[frac_size], keys[frac_size][2];
    Frac new_fractions[frac_size];
    long long l, den_square;

    // Build keys
    den_square = max_den * max_den;
    keys[0][0] = 0;
    keys[0][1] = 0;

    for (int i = 1; i < frac_size; i++) {
        num = fractions[i].num;
        den = fractions[i].den;
        l = den_square * num / den;
        keys[i][0] = int(l / max_den);
        keys[i][1] = int(l - keys[i][0] * max_den);
    }

    // BinSort
    // sort by dimension 1
    for (int j = 0; j <= max_den; j++) cnt[j] = 0;
    for (int j = 0; j < frac_size; j++) cnt[keys[j][1]]++;
    for (int j = 1; j <= max_den; j++) cnt[j] += cnt[j - 1];
    for (int j = frac_size - 1; j >= 0; j--) {
        loc = --cnt[keys[j][1]];
        new_fractions[loc].den = fractions[j].den;
        new_fractions[loc].num = fractions[j].num;
        index[loc] = j;
    }

    // sort by dimension 0
    for (int j = 0; j <= max_den; j++) cnt[j] = 0;
    for (int j = 0; j < frac_size; j++) cnt[keys[j][0]]++;
    for (int j = 1; j <= max_den; j++) cnt[j] += cnt[j - 1];
    for (int j = frac_size - 1; j >= 0; j--) {
        loc = --cnt[keys[index[j]][0]];
        fractions[loc].den = new_fractions[j].den;
        fractions[loc].num = new_fractions[j].num;
    }
}

static void FracBinSort_(Frac *&fractions, int max_den, int frac_size) {
    int num, den, l1, l2, loc, cnt[max_den + 1], ***keys;
    long long l, den_square;
    Frac *new_fractions;

    den_square = max_den * max_den;
    keys = new int **[max_den + 1];
    for (int i = 0; i <= max_den; i++) keys[i] = nullptr;

    // Build keys
    keys[0] = new int *[1];
    keys[0][0] = new int[2]{0, 0};

    for (int i = 1; i < frac_size; i++) {
        num = fractions[i].num;
        den = fractions[i].den;
        l = den_square * num / den;
        l1 = int(l / max_den);
        l2 = int(l - l1 * max_den);

        if (keys[den] == nullptr) {
            keys[den] = new int *[den + 1];
            for (int j = 0; j <= den; j++) keys[den][j] = nullptr;
        }
        keys[den][num] = new int[2]{l1, l2};
    }

    // binSort
    for (int i = 1; i >= 0; i--) {
        new_fractions = new Frac[frac_size];
        for (int j = 0; j <= max_den; j++) cnt[j] = 0;
        for (int j = 0; j < frac_size; j++) {
            cnt[keys[fractions[j].den][fractions[j].num][i]]++;
        }
        for (int j = 1; j <= max_den; j++) cnt[j] += cnt[j - 1];
        for (int j = frac_size - 1; j >= 0; j--) {
            loc = --cnt[keys[fractions[j].den][fractions[j].num][i]];
            new_fractions[loc].den = fractions[j].den;
            new_fractions[loc].num = fractions[j].num;
        }

        delete[] fractions;
        fractions = new_fractions;
    }

    // release keys
    for (den = 0; den <= max_den; den++) {
        if (keys[den] != nullptr) {
            for (num = 0; num <= den; num++) {
                if (keys[den][num] != nullptr) delete[] keys[den][num];
            }
            delete[] keys[den];
        }
    }
    delete[] keys;
}


#endif //CSM4GMG_FRACUTILS_H
