//
// Created by ldd on 2021/10/6.
//

#ifndef CSM4GMG_VIZ_H
#define CSM4GMG_VIZ_H

#include "../Header.h"

static char default_viz_cmd[] = "dot -Tpng TmpFile/%s.dot -o TmpFile/%s.png";
static char default_eog_cmd[] = "eog TmpFile/%s.png";

static void RunViz(const string &dot_filename) {
    char viz_cmd[256], eog_cmd[256];

    sprintf(viz_cmd, default_viz_cmd, dot_filename.c_str(), dot_filename.c_str());
    sprintf(eog_cmd, default_eog_cmd, dot_filename.c_str());

    system(viz_cmd);
    system(eog_cmd);
}

#endif //CSM4GMG_VIZ_H
