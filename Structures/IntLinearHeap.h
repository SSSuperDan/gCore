//
// Created by ldd on 2021/9/24.
//

/*
 * The linear heap data structure is proposed by Chang et al. in:
 *
 * " Lijun Chang and Lu Qin. Cohesive Subgraph Computation over
 * Large Sparse Graphs. Springer Series in the Data Sciences, 2018 "
 *
 * Below is the array-based implementation, referencing the source code in
 * https://github.com/LijunChang/Cohesive_subgraph_book
 *
 */

#ifndef CSM4GMG_INTLINEARHEAP_H
#define CSM4GMG_INTLINEARHEAP_H

#include "../Header.h"

class IntLinearHeap {

public:
    IntLinearHeap() = default;
    ~IntLinearHeap();
    void Clear();
    void Init(int n_);
    void SetBin(int bin_size_);
    void ReleaseBin();

    void Insert(int v, int value);
    void Remove(int v);
    void Update(int v, int new_value);
    void Decrease(int v);
    void PrintHeap();
    bool SelectMin(int &v, int &value);
    int GetElements(int value, int * arr);
    int GetHeapSize();

    [[nodiscard]] int GetMinValue() {
        Tighten();
        return min_value;
    }

    [[nodiscard]] int GetMaxValue() {
        Tighten();
        return max_value;
    }

    int GetValue(int v) {
        return heap[v];
    }

    int* GetValue() {
        return heap;
    }

    bool Empty() {
        Tighten();
        return min_value > max_value;
    }

private:
    inline void Tighten() {
        while (min_value <= max_value && bin[min_value] == -1) ++min_value;
        while (min_value <= max_value && bin[max_value] == -1) --max_value;
    }

    int *pre{nullptr};
    int *next{nullptr};
    int *bin{nullptr};
    int *heap{nullptr};
    int n{0};
    int bin_size{0};
    int max_value{0}, min_value{0};
};


#endif //CSM4GMG_INTLINEARHEAP_H
