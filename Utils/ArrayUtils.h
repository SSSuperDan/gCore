//
// Created by ldd on 2021/9/22.
//

#ifndef CSM4GMG_ARRAYUTILS_H
#define CSM4GMG_ARRAYUTILS_H

#include <cmath>

#include "../Header.h"
#include "StringUtils.h"

const float epsilon = 1e-6;

static std::basic_string<char> fto_string(const Frac & frac){
    return to_string(frac.num) + "/" + to_string(frac.den);
}

template<class T>
inline bool ArrGe(const T *arr1, const T *arr2, int length) {
    for (int i = 0; i < length; i++) {
        if (arr1[i] < arr2[i]) {
            return false;
        }
    }
    return true;
}

template<class T>
inline bool ArrEq(const T *arr1, const T *arr2, int length) {
    for (int i = 0; i < length; i++) {
        if (arr1[i] != arr2[i]) {
            return false;
        }
    }
    return true;
}


template<class T>
static string Array2String(const T *arr, int length) {
    int i;
    string s = "[";

    for (i = 0; i < length - 1; i++) {
        s += to_string(arr[i]) + ",";
    }
    s += to_string(arr[i]) + "]";

    return s;
}

static string Array2String(const Frac *arr, int length) {
    int i;
    string s = "[";

    for (i = 0; i < length - 1; i++) {
        s += fto_string(arr[i]) + ",";
    }
    s += fto_string(arr[i]) + "]";

    return s;
}


template<class T>
static string Array2String(const int *index, const T *arr, int length) {
    int i;
    string s = "[";

    for (i = 0; i < length - 1; i++) {
        s += to_string(arr[index[i]]) + ",";
    }
    s += to_string(arr[index[i]]) + "]";

    return s;
}

static string FracVec2String(const int *index, Frac **frac_arr, int dim) {
    int i;
    string s = "[";

    for (i = 0; i < dim - 1; i++) {
        s += fto_string(frac_arr[i][index[i]]) + ",";
    }
    s += fto_string(frac_arr[i][index[i]]) +"]";

    return s;
}

template<class T>
static string Vec2String(const vector<T> &arr) {
    int i, length = (int) arr.size();

    string s = "[";
    for (i = 0; i < length - 1; i++) {
        s += to_string(arr[i]) + ",";
    }
    s += to_string(arr[i]) + "]";

    return s;
}


static string Vec2String(const vector<Frac> &arr) {
    int i, length = (int) arr.size();

    string s = "[";
    for (i = 0; i < length - 1; i++) {
        s += fto_string(arr[i]) + ",";
    }
    s += fto_string(arr[i]) + "]";

    return s;
}

static void String2Vec(string s, vector<int> &arr) {
    vector<string> str_arr;
    s = s.substr(1, s.size() - 2);  // remove '[' and ']'
    Split(s, str_arr, false, ',');
    for (auto &ss:str_arr) arr.emplace_back(stoi(ss));
}

static void String2Vec(string s, vector<float> &arr) {
    vector<string> str_arr;
    s = s.substr(1, s.size() - 2);  // remove '[' and ']'
    Split(s, str_arr, false, ',');
    for (auto &ss:str_arr) arr.emplace_back(stof(ss));
}

// Exact linear search.
template<class T>
inline int SequentialSearch(T *array, int length, const T &ele) {
    for (int i = 0; i < length; i++) {
        if (array[i] == ele) {
            return i;
        }
    }
    return -1;  // not exist.
}

template<class T>
inline int SequentialSearch(vector<T> array, const T &ele) {
    int length = (int) array.size();
    for (int i = 0; i < length; i++) {
        if (array[i] == ele) {
            return i;
        }
    }
    return -1;  // not exist.
}

// Exponential search over int array.
inline int BiSearch(const int *arr, int length, int v) {
    int i, j, m;

    if (arr[0] == v) return 0;
    else {
        i = 1;
        while (i < length) {
            if (arr[i] == v) return i;
            else if (arr[i] > v) break;
            else i *= 2;
        }

        // binary search
        j = std::min(length, i);
        i = i >> 1;
        while (i < j) {
            m = (i + j) >> 1;
            if (arr[m] < v) i = m + 1;
            else if (arr[m] > v) j = m;
            else break;
        }
        return (i + j) >> 1;
    }
}

// Exponential search over fraction array.
inline int BiSearch(Frac *fractions, int length, const Frac &frac) {
    int i, j, l, r, m;

    if (frac.num == 0) return 0;
    else {
        i = 1;
        while (i < length) {
            l = fractions[i].num * frac.den;
            r = frac.num * fractions[i].den;
            if (l == r) return i;
            else if (l > r) break;
            else i *= 2;
        }

        // binary search
        j = std::min(length, i);
        i = i >> 1;
        while (i < j) {
            m = (i + j) >> 1;
            l = fractions[m].num * frac.den;
            r = frac.num * fractions[m].den;
            if (l < r) i = m + 1;
            else if (l > r) j = m;
            else break;
        }
        return (i + j) >> 1;
    }
}


// Exponential search over fraction array.
inline int BiSearch(Frac *fractions, int length, float f) {
    int i, j, m;
    float fm;

    if (fabs(f - 0) < epsilon) return 0;
    else {
        i = 1;
        while (i < length) {
            fm = (float) fractions[i].num / (float) fractions[i].den;
            if (fabs(fm - f) < epsilon) return i;
            else if (fm > f) break;
            else i *= 2;
        }

        // binary search
        j = std::min(length, i);
        i = i >> 1;
        while (i < j) {
            m = (i + j) >> 1;
            fm = (float) fractions[m].num / (float) fractions[m].den;
            if (fabs(fm - f) < epsilon) break;
            else if (fm < f) i = m + 1;
            else j = m;
        }
        return (i + j) >> 1;
    }
}

#endif //CSM4GMG_ARRAYUTILS_H
