//
// Created by ldd on 2022/6/22.
//

#ifndef CSM4GMG_FILEUTILS_H
#define CSM4GMG_FILEUTILS_H

#include "../Header.h"

static void LoadColumn(const string &filename, int target_column, unordered_map<long long, string> &map, char sep = '\t') {
    string line, v_name;
    size_t pos1, pos2;
    long long v;

    auto fin = std::ifstream(filename);
    while (fin.peek() != EOF) {
        getline(fin, line);

        pos1 = line.find_first_of(sep);
        v = stoll(line.substr(0, pos1));

        for (int i = 1; i < target_column - 1; i++) {
            pos1 = line.find(sep, pos1 + 1);
        }
        pos2 = line.find(sep, pos1 + 1);
        if (pos2 == string::npos) pos2 = line.size();

        v_name = line.substr(pos1 + 1, pos2 - pos1);
        if (v_name[v_name.size() - 1] == '\r') v_name = v_name.substr(0, v_name.size() - 1);
        map[v] = v_name;
    }
    fin.close();

}

static string GetNextValidFile(const string &filename) {
    struct stat buffer{};
    string pref, suf;
    size_t pos1, pos2;
    int id;

    if (!stat(filename.c_str(), &buffer)) {
        pos1 = filename.find_last_of('_') + 1;
        pos2 = filename.find_last_of('.');
        try {
            id = std::stoi(filename.substr(pos1, pos2 - pos1));
            pref = filename.substr(0, pos1);
        } catch (std::invalid_argument &) {
            id = 0;
            pref = filename.substr(0, pos2) + '_';
        }

        suf = filename.substr(pos2);
        while (true) {
            string new_file = pref;
            new_file += std::to_string(++id);
            new_file += suf;

            if (stat(new_file.c_str(), &buffer)) {
                return new_file;
            }
        }
    } else return filename;
}

#endif //CSM4GMG_FILEUTILS_H
