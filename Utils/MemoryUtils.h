//
// Created by ldd on 2022/6/22.
//

#ifndef CSM4GMG_MEMORYUTILS_H
#define CSM4GMG_MEMORYUTILS_H
/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t GetPeakRSS() {
    struct rusage rusage{};
    getrusage(RUSAGE_SELF, &rusage);
    return (size_t) (rusage.ru_maxrss * 1024L);
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t GetCurrentRSS() {

    long rss = 0L;
    FILE *fp = nullptr;
    if ((fp = fopen("/proc/self/statm", "r")) == nullptr)
        return (size_t) 0L;      /* Can't open? */
    if (fscanf(fp, "%*s%ld", &rss) != 1) {
        fclose(fp);
        return (size_t) 0L;      /* Can't read? */
    }
    fclose(fp);
    return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);
}


#endif //CSM4GMG_MEMORYUTILS_H
