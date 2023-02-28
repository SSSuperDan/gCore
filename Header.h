//
// Created by ldd on 2021/9/15.
//

#ifndef CSM4GMG_HEADER_H
#define CSM4GMG_HEADER_H

#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <functional>
#include <algorithm>
#include <dirent.h>
#include <unistd.h>
#include <iomanip>
#include <chrono>

#include <cmath>
#include <climits>
#include <random>
#include <string>
#include <cstring>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <exception>


#include "config.h"

using std::default_random_engine;
using std::mt19937;
using std::string;
using std::ifstream;
using std::ofstream;
using std::cerr;
using std::cout;

using std::setprecision;
using std::make_pair;
using std::to_string;
using std::floor;
using std::ceil;
using std::fabs;
using std::move;
using std::sort;
using std::endl;
using std::move;

using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::queue;
using std::pair;
using std::map;
using std::hash;
using std::random_device;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

#ifdef USE_MEMTRACK
#include "memtrack/memtrack.h"
using namespace MemTrack;
#endif


struct Frac {
    int den;
    int num;

    Frac() = default;
    Frac(int den_, int num_) : den(den_), num(num_) {}
};


#endif //CSM4GMG_HEADER_H



