//
// Created by ldd on 2021/9/22.
//

#ifndef CSM4GMG_TIMER_H
#define CSM4GMG_TIMER_H

#include "../Header.h"

class Timer {
public:
    Timer() {
        gettimeofday(&start, nullptr);
    }

    ~Timer() = default;

    void Stop() {
        gettimeofday(&end, nullptr);
    }

    [[nodiscard]] double GetTimeInMs() const {
        return 1000 * double(end.tv_sec - start.tv_sec) + double(end.tv_usec - start.tv_usec) / 1000;
    }

private:
    timeval start{};
    timeval end{};
};

#endif //CSM4GMG_TIMER_H
