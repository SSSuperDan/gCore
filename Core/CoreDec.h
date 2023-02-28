//
// Created by ldd on 2021/9/24.
//

#ifndef CSM4GMG_COREDEC_H
#define CSM4GMG_COREDEC_H

#include "../Graphs/UndirectedGraph.h"
#include "../Structures/IntLinearHeap.h"

class CoreDec {
public:
    static int KCore(UndirectedGraph *g, int k, int *core);
    static int KCore(UndirectedGraph *g, const int *com, int length, int k, int *core);
    static int GetDegeneracy(UndirectedGraph *g);
    static int ArrayBasedCoreDec(UndirectedGraph *g, int *core_number);
    static int ArrayBasedCoreDec(UndirectedGraph *g, const int *com, int length, int *core_number);
    static int ILHBasedCoreDec(UndirectedGraph *g, int *core_number);
    static void HIndexBasedCoreDec(UndirectedGraph *g, int *core_number);
private:
    static inline int GetDeg(const int *adj_lst, const int *valid) { // valid[i] >= 0
        int deg = 0;

        for (int j = 1; j <= adj_lst[0]; j++) {
            if (valid[adj_lst[j]] >= 0) deg++;
        }

        return deg;
    }
};


#endif //CSM4GMG_COREDEC_H
