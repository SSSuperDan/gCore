//
// Created by ldd on 2022/6/6.
//

#include "Frac2IntPri.h"

Frac2IntPri::Frac2IntPri(int dim_) : dim(dim_) {
    fractions = new Frac *[dim];
    priority = new int **[dim];
    frac_size = new int[dim];
    max_den = new int[dim];
}

Frac2IntPri::~Frac2IntPri() {
    if (fractions) {
        for (int i = 0; i < dim; i++) delete[] fractions[i];
        delete[] fractions;
    }

    if (priority) {
        for (int i = 0; i < dim; i++) {
            if (priority[i]) {
                for (int j = 0; j <= max_den[i]; j++) {
                    if (priority[i][j]) delete[] priority[i][j];
                }
                delete[] priority[i];
            }
        }
        delete[] priority;
    }

    delete[] frac_size;
    delete[] max_den;
}

// Build map for the i th dimension.
void Frac2IntPri::Build(int i, const int *den_arr, int den_size, int md) const {
    max_den[i] = md;
    frac_size[i] = GenFractions(i, den_arr, den_size, md);
    SetPriority(i);
}

// Convert a fraction vector to an integer priority vector.
void Frac2IntPri::Convert(vector<Frac> &frac_vec, vector<int> &pri_vec) const {
    for (int i = 0; i < dim; i++) {
        pri_vec[i] = BiSearch(fractions[i], frac_size[i], frac_vec[i]);
    }
}

// Convert a float vector to an integer priority vector.
void Frac2IntPri::Convert(vector<float> &f_vec, vector<int> &pri_vec) const {
    for (int i = 0; i < dim; i++) {
        pri_vec[i] = BiSearch(fractions[i], frac_size[i], f_vec[i]);
    }
}

// Convert an integer priority vector to a fraction vector.
bool Frac2IntPri::Convert(vector<int> &pri_vec, vector<Frac> &frac_vec) const {
    for (int i = 0; i < dim; i++) {
        if (pri_vec[i] >= frac_size[i]) return false;
        frac_vec[i] = fractions[i][pri_vec[i]];
    }
    return true;
}

void Frac2IntPri::Flush(const string &file){
    int mm_den = max_den[0];
    for (int i = 1; i < dim; i++) if (mm_den < max_den[i]) mm_den = max_den[i];
    int aux_arr[mm_den + 1], size;

    auto f = std::ofstream(file);
    f.write(reinterpret_cast<char *> (&dim), sizeof(int));

    for (int i = 0; i < dim; i++) {
        size = 0;
        for (int j = 0; j <= max_den[i]; j++) {
            if (priority[i][j]) aux_arr[size++] = j;
        }
        f.write(reinterpret_cast<char *>(&size), sizeof(int));
        f.write(reinterpret_cast<char *>(aux_arr), long(size * sizeof(int)));
    }
    f.close();
}

void Frac2IntPri::Load(std::basic_ifstream<char> &f) const {

    int size;

    for (int i = 0; i < dim; i++) {
        f.read(reinterpret_cast<char *>(&size), sizeof(int));

        int aux_arr[size];
        f.read(reinterpret_cast<char *>(aux_arr), long(size * sizeof(int)));
        Build(i, aux_arr, size, aux_arr[size-1]);
    }

}

Frac2IntPri* Frac2IntPri::Load(const string &file) {
    int dimension;
    auto f = std::ifstream(file);

    f.read(reinterpret_cast<char *>(&dimension), sizeof(int));

    Frac2IntPri * f2i = new Frac2IntPri(dimension);
    f2i->Load(f);
    f.close();

    return f2i;
}

bool Frac2IntPri::Equals(Frac2IntPri &f2i) const {
    // Equals basic elements
    if (dim != f2i.dim) return false;
    if (!ArrEq(max_den, f2i.max_den, dim)) return false;
    if (!ArrEq(frac_size, f2i.frac_size, dim)) return false;

    // Equals fractions
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < frac_size[i]; j++) {
            if (fractions[i][j].den != f2i.fractions[i][j].den || fractions[i][j].num != f2i.fractions[i][j].num) {
                return false;
            }
        }
    }
    return true;
}

// Generate and sort all fractions with denominators in den_arr.
int Frac2IntPri::GenFractions(int i, const int *den_arr, int den_size, int md) const {
    int num_of_frac, idx;
    bool bin[md + 1];
    Frac *frac_lst;

    // initialize
    memset(bin, false, (md + 1) * sizeof(bool));
    for (int j = 0; j < den_size; j++) bin[den_arr[j]] = true;

    num_of_frac = 1;  // 1 is for fraction with denominator 0
    for (int j = 1; j <= md; j++) if (bin[j]) num_of_frac += j + 1;

    // construct fractions
    frac_lst = new Frac[num_of_frac];
    frac_lst[0].num = 0;
    frac_lst[0].den = 0;

    idx = 1;
    for (int j = 1; j <= md; j++) {
        if (bin[j]) {
            for (int l = 0; l <= j; l++) {
                frac_lst[idx].den = j;
                frac_lst[idx++].num = l;
            }
        }
    }

    FracBinSort(frac_lst, md, num_of_frac);
    fractions[i] = frac_lst;
    return num_of_frac;
}

// Assign unique priority to fractions based on their values.
void Frac2IntPri::SetPriority(int i) const {
    int n, d, pn, pd, loc, idx, **pri;
    Frac *frac_lst = fractions[i];

    loc = 0;
    pri = new int *[max_den[i] + 1];
    for (int j = 0; j <= max_den[i]; j++) pri[j] = nullptr;

    pri[0] = new int[1];
    pri[0][0] = loc;

    idx = 1;
    while (frac_lst[idx].num == 0) {
        d = frac_lst[idx++].den;
        pri[d] = new int[d + 1];
        pri[d][0] = loc;
    }

    // init merge
    pd = frac_lst[idx].den;
    pn = frac_lst[idx++].num;
    frac_lst[++loc].den = pd;
    frac_lst[loc].num = pn;
    pri[pd][pn] = loc;

    // merge and assign
    for (; idx < frac_size[i]; idx++) {
        d = frac_lst[idx].den;
        n = frac_lst[idx].num;
        if (IsNotEqual(d, n, pd, pn)) {
            frac_lst[++loc].den = d;
            frac_lst[loc].num = n;
            pn = n;
            pd = d;
        }
        pri[d][n] = loc;
    }
    frac_size[i] = loc + 1;
    priority[i] = pri;
}

