//
// Created by ldd on 2022/6/6.
//

#ifndef CSM4GMG_FRAC2INTPRI_H
#define CSM4GMG_FRAC2INTPRI_H

#include "../Header.h"
#include "../Utils/FracUtils.h"
#include "../Utils/ArrayUtils.h"

// Map from fraction vector to integer vector.
class Frac2IntPri {
public:
    explicit Frac2IntPri(int dim_);
    ~Frac2IntPri();
    void Build(int i, const int *den_arr, int den_size, int md) const;
    void Convert(vector<Frac> & frac_vec, vector<int> & pri_vec) const;  // Convert(source, target);
    void Convert(vector<float> & f_vec, vector<int> & pri_vec) const;
    bool Convert(vector<int> & pri_vec, vector<Frac> & frac_vec) const;
    void Flush(const string & file);
    static Frac2IntPri* Load(const string & file);
    bool Equals(Frac2IntPri & f2i) const;

    Frac **fractions;
    int ***priority;
    int *frac_size;
    int *max_den;
    int dim;

    int GenFractions(int i, const int *den_arr, int den_size, int md) const;
    void SetPriority(int i) const;

private:
    void Load(std::basic_ifstream<char> &f) const;
};


#endif //CSM4GMG_FRAC2INTPRI_H
