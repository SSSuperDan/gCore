//
// Created by ldd on 2021/12/17.
//

#ifndef CSM4GMG_STRINGUTILS_H
#define CSM4GMG_STRINGUTILS_H

#include "../Header.h"

static void Split(const string &str, vector<string> &tokens, bool skip_first, char delim = ' ') {
    size_t last_pos, pos;

    tokens.clear();

    last_pos = str.find_first_not_of(delim, 0);
    if (skip_first) {
        last_pos = str.find(delim, last_pos);
        last_pos = str.find_first_not_of(delim, last_pos);
    }

    while (last_pos != string::npos) {
        pos = str.find(delim, last_pos);
        tokens.emplace_back(move(str.substr(last_pos, pos - last_pos)));
        last_pos = str.find_first_not_of(delim, pos);
    }
}

#endif //CSM4GMG_STRINGUTILS_H
