//
// Created by ldd on 2021/9/24.
//

#include "CoreDec.h"

int CoreDec::KCore(UndirectedGraph *g, int k, int *core) {
    int n = g->GetN(), **adj_lst, deg[n], cnt = 0;
    int v, u, queue[n], start = 0, end = 0;
    bool in_queue[n];

    memset(in_queue, false, n * sizeof(bool));
    adj_lst = g->GetAdjLst();
    for (v = 0; v < n; v++) {
        deg[v] = adj_lst[v][0];
        if (deg[v] < k) {
            queue[end++] = v;
            in_queue[v] = true;
        }
    }

    while (start < end) {
        v = queue[start++];
        for (int i = 1; i <= adj_lst[v][0]; i++) {
            u = adj_lst[v][i];
            deg[u]--;
            if (deg[u] < k && !in_queue[u]) {
                queue[end++] = u;
                in_queue[u] = true;
            }
        }
    }

    for (v = 0; v < n; v++) {
        if (!in_queue[v]) core[cnt++] = v;
    }

    return cnt;
}

int CoreDec::KCore(UndirectedGraph *g, const int *com, int length, int k, int *core) {
    int n = g->GetN(), **adj_lst, deg[length], cnt = 0;
    int v, u, queue[length], start = 0, end = 0, in_com[n];
    bool in_queue[length];

    memset(in_com, -1, n * sizeof(int));
    for (int i = 0; i < length; i++) in_com[com[i]] = i;

    memset(in_queue, false, length * sizeof(bool));
    adj_lst = g->GetAdjLst();
    for (int i = 0; i < length; i++) {
        deg[i] = GetDeg(adj_lst[com[i]], in_com);
        if (deg[i] < k) {
            queue[end++] = i;
            in_queue[i] = true;
        }
    }

    while (start < end) {
        v = queue[start++];
        for (int i = 1; i <= adj_lst[com[v]][0]; i++) {
            u = in_com[adj_lst[com[v]][i]];
            if (u >= 0) {
                deg[u]--;
                if (deg[u] < k && !in_queue[u]) {
                    queue[end++] = u;
                    in_queue[u] = true;
                }
            }
        }
    }

    for (int i = 0; i < length; i++) {
        if (!in_queue[i]) core[cnt++] = com[i];
    }

    return cnt;

}

int CoreDec::GetDegeneracy(UndirectedGraph *g) {
    int core_number[g->GetN()];
    return ArrayBasedCoreDec(g, core_number);
}

int CoreDec::ArrayBasedCoreDec(UndirectedGraph *g, int *core_number) {
    int n = g->GetN(), max_deg = g->GetMaxDeg(), start, num, v, pu, pw, w, u, **adj_lst;
    int pos[n], vert[n], bin[max_deg + 1];

    // Init
    adj_lst = g->GetAdjLst();
    for (int i = 0; i < n; i++) {
        core_number[i] = adj_lst[i][0];
    }
    memset(bin, 0, (max_deg + 1) * sizeof(int));

    // Bin sort
    for (int i = 0; i < n; i++) {
        bin[core_number[i]]++;
    }

    start = 0;
    for (int d = 0; d <= max_deg; d++) {
        num = bin[d];
        bin[d] = start;
        start += num;
    }

    for (int i = 0; i < n; i++) {
        pos[i] = bin[core_number[i]];
        vert[pos[i]] = i;
        bin[core_number[i]]++;
    }

    for (int d = max_deg; d >= 1; d--) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;

    // Compute core numbers.
    for (int i = 0; i < n; i++) {
        v = vert[i];
        for (int j = 1; j <= adj_lst[v][0]; j++) {
            u = adj_lst[v][j];
            if (core_number[u] > core_number[v]) {
                pu = pos[u], pw = bin[core_number[u]], w = vert[pw];

                vert[pw] = u, vert[pu] = w;
                pos[u] = pw, pos[w] = pu;

                bin[core_number[u]]++;
                core_number[u]--;
            }
        }
    }

    // Return the maximum core number.
    return core_number[v];
}

int CoreDec::ArrayBasedCoreDec(UndirectedGraph *g, const int *com, int length, int *core_number) {
    int n = g->GetN(), max_deg = std::min(g->GetMaxDeg(), length - 1), start, num, v, pu, pw, w, u, **adj_lst;
    int pos[length], vert[length], bin[max_deg + 1], in_com[n];

    memset(in_com, -1, n * sizeof(int));
    for (int i = 0; i < length; i++) in_com[com[i]] = i;

    // Init
    adj_lst = g->GetAdjLst();
    for (int i = 0; i < length; i++) {
        core_number[i] = GetDeg(adj_lst[com[i]], in_com);
    }

    memset(bin, 0, (max_deg + 1) * sizeof(int));

    // Bin sort
    for (int i = 0; i < length; i++) {
        bin[core_number[i]] ++;
    }

    start = 0;
    for (int d = 0; d <= max_deg; d++) {
        num = bin[d];
        bin[d] = start;
        start += num;
    }

    for (int i = 0; i < length; i++) {
        pos[i] = bin[core_number[i]];
        vert[pos[i]] = i;
        bin[core_number[i]] ++;
    }

    for (int d = max_deg; d >= 1; d--) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;

    // Compute core numbers.
    for (int i = 0; i < length; i++) {
        v = vert[i];
        for (int j = 1; j <= adj_lst[com[v]][0]; j++) {
            u = in_com[adj_lst[com[v]][j]];

            if (u >= 0 && core_number[u] > core_number[v]) {
                pu = pos[u], pw = bin[core_number[u]], w = vert[pw];

                vert[pw] = u, vert[pu] = w;
                pos[u] = pw, pos[w] = pu;

                bin[core_number[u]]++;
                core_number[u]--;
            }
        }
    }

    // Return the maximum core number.
    return core_number[v];
}


int CoreDec::ILHBasedCoreDec(UndirectedGraph *g, int *core_number) {
    int n, v, u, value, **adj_lst, *coreness;
    IntLinearHeap heap;

    // Init
    n = g->GetN();
    adj_lst = g->GetAdjLst();
    heap.Init(n);
    heap.SetBin(g->GetMaxDeg());
    coreness = heap.GetValue();

    for (int i = 0; i < n; i++) {
        heap.Insert(i, adj_lst[i][0]);
    }

    // Compute core numbers.
    for (int i = 0; i < n; i++) {
        heap.SelectMin(v, value);
        heap.Remove(v);
        for (int j = 1; j <= adj_lst[v][0]; j++) {
            u = adj_lst[v][j];
            if (coreness[u] > value) heap.Decrease(u);
        }
    }

    memcpy(core_number, coreness, n * sizeof(int));
    return value;  // Return the maximum core number.
}

void CoreDec::HIndexBasedCoreDec(UndirectedGraph *g, int *core_number) {
    int n = g->GetN(), v, u, old_cn, idx1, idx2, **adj_lst;
    int cnt[n], queue[n << 1], counter[g->GetMaxDeg() + 1], *q1, *q2;
    bool in_q1[n];

    // Init
    q1 = queue, q2 = &queue[n];
    adj_lst = g->GetAdjLst();
    for (int i = 0; i < n; i++) {
        core_number[i] = adj_lst[i][0];
        in_q1[i] = false;
    }

    idx1 = 0;
    for (int i = 0; i < n; i++) {
        cnt[i] = 0;
        for (int j = 1; j <= adj_lst[i][0]; j++) if (core_number[adj_lst[i][j]] >= core_number[i]) cnt[i]++;
        if (cnt[i] < core_number[i]) q1[idx1++] = i;
    }

    while (idx1) {
        idx2 = 0;
        for (int i = 0; i < idx1; i++) in_q1[q1[i]] = true;
        for (int i = 0; i < idx1; i++) {
            v = q1[i];
            in_q1[v] = false;
            old_cn = core_number[v];

            // Compute the h-index of v.
            memset(counter, 0, (old_cn + 1) * sizeof(int));
            for (int j = 1; j <= adj_lst[v][0]; j++) {
                u = adj_lst[v][j];
                if (core_number[u] >= old_cn) counter[old_cn] += 1;
                else counter[core_number[u]] += 1;
            }
            for (int j = old_cn; j >= 1; j--) {
                if (counter[j] >= j) {
                    core_number[v] = j;
                    break;
                }
                counter[j - 1] += counter[j];
            }

            cnt[v] = 0;
            for (int j = 1; j <= adj_lst[v][0]; j++) {
                u = adj_lst[v][j];
                if (core_number[u] >= core_number[v]) cnt[v]++;
                if (!in_q1[u] && core_number[u] > core_number[v] && core_number[u] <= old_cn) {
                    if (cnt[u] == core_number[u]) q2[idx2++] = u;
                    cnt[u]--;
                }
            }
        }
        idx1 = idx2;
        std::swap(q1, q2);
    }
}
