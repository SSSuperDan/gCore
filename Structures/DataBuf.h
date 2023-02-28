//
// Created by ldd on 2022/6/30.
//

#ifndef CSM4GMG_DATABUF_H
#define CSM4GMG_DATABUF_H

#include "../Header.h"

#define PAGE_SIZE (1<<12)

// Page containing massive small/medium-sized data lists.
template<class T>
struct Page {
    T data[PAGE_SIZE]{};
    Page *next_page{nullptr};
};

// Buffer for small/medium-sized data.
template<class T>
struct MDataBuf {

    int cur_pos{0};
    Page<T> *head;
    Page<T> *cur_page;

    inline T *Allocate(size_t len = 1) {
        T *space;

        if (cur_pos + len > PAGE_SIZE) {
            auto page = new Page<T>();
            cur_page->next_page = page;
            cur_page = page;
            cur_pos = 0;
        }

        space = cur_page->data + cur_pos;
        cur_pos += len;
        return space;
    }

    MDataBuf() {
        head = new Page<T>();
        cur_page = head;
    }

    ~MDataBuf() {
        cur_page = head;
        while (cur_page) {
            head = head->next_page;
            delete cur_page;
            cur_page = head;
        }
    }
};

// Page containing a single large-sized data list.
template<class T>
struct LPage {
    T *data{};
    LPage<T> *next_page{nullptr};

    LPage() = default;

    explicit LPage(int len) {
        data = new T[len];
    }

    ~LPage() {
        delete[] data;
    }
};

// Buffer for large-sized data.
template<class T>
struct LDataBuf {
    LPage<T> *head;
    LPage<T> *cur_page;

    inline T *Allocate(size_t len = 1) {
        auto page = new LPage<T>(len);
        cur_page->next_page = page;
        cur_page = page;
        return page->data;
    }

    LDataBuf() {
        head = new LPage<T>();
        cur_page = head;
    }

    ~LDataBuf() {
        cur_page = head;
        while (cur_page) {
            head = head->next_page;
            delete cur_page;
            cur_page = head;
        }
    }
};

template<class T>
struct DataBuf {

    MDataBuf<T> adj_buf;
    LDataBuf<T> ladj_buf;

    inline T *Allocate(size_t len) {
        return len > PAGE_SIZE ? ladj_buf.Allocate(len) : adj_buf.Allocate(len);
    }

    DataBuf() = default;

    ~DataBuf() = default;
};

#endif //CSM4GMG_DATABUF_H

