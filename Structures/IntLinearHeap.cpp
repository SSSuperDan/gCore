//
// Created by ldd on 2021/9/24.
//

#include "IntLinearHeap.h"

void IntLinearHeap::Init(int n_) {
    n = n_;

    pre = new int[n];
    next = new int[n];
    heap = new int[n];
}

IntLinearHeap::~IntLinearHeap() {
    delete[] pre;
    delete[] next;
    delete[] bin;
    delete[] heap;
}

inline void IntLinearHeap::Clear() {
    for (int i = 0; i <= bin_size; i++) bin[i] = -1;
    min_value = bin_size;
    max_value = 0;
}

void IntLinearHeap::SetBin(int bin_size_) {
    bin_size = bin_size_;
    bin = new int[bin_size + 1];
    Clear();
}

void IntLinearHeap::ReleaseBin() {
    delete[] bin;
    bin = nullptr;
}

void IntLinearHeap::Insert(int v, int value) {
    int head;

    head = bin[value];
    heap[v] = value;
    pre[v] = -1;
    next[v] = head;
    bin[value] = v;
    if (head != -1) pre[head] = v;

    if (value < min_value) min_value = value;
    if (value > max_value) max_value = value;
}

void IntLinearHeap::Remove(int v) {
    if (pre[v] == -1) { // v is the first element of the bin
        bin[heap[v]] = next[v];
        if (next[v] != -1) pre[next[v]] = -1;
    } else { // v is not the first element
        int pv = pre[v];
        next[pv] = next[v];
        if (next[v] != -1) pre[next[v]] = pv;
    }
}

bool IntLinearHeap::SelectMin(int &v, int &value) {
    if (Empty()) return false;
    value = min_value;
    v = bin[min_value];
    return true;
}

void IntLinearHeap::Update(int v, int new_value) {
    Remove(v);
    Insert(v, new_value);
}

void IntLinearHeap::Decrease(int v) {
    Remove(v);
    heap[v]--;
    Insert(v, heap[v]);
}

int IntLinearHeap::GetElements(int value, int *arr) {
    int curr, i;

    i = 0;
    if (value >= min_value && value <= max_value) {
        curr = bin[value];
        while (curr != -1) {
            arr[i++] = curr;
            curr = next[curr];
        }
    }
    return i;
}

void IntLinearHeap::PrintHeap() {
    int curr, cnt;

    for (int i = min_value; i <= max_value; i++) {
        if (bin[i] != -1) {
            cout << "[" << i << "] : ";

            cnt = 0;
            curr = bin[i];
            while (curr != -1) {
                cnt ++;
                cout << curr << "-->";
                curr = next[curr];
            }
            cout << "-1 (" << cnt << ")" << endl;
        }
    }
    cout << "=========" << endl;
}

int IntLinearHeap::GetHeapSize() {
    int curr, cnt;

    cnt = 0;
    for (int i = min_value; i <= max_value; i++) {
        curr = bin[i];
        while (curr != -1) {
            cnt++;
            curr = next[curr];
        }
    }
    return cnt;
}