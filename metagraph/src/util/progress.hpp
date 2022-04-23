//
// Created by amir on 18/12/2020.
//
#pragma once
#include <istream>

namespace ts {


struct progress_bar {
    static size_t it;
    static size_t total;
    static size_t bar_len;
    static size_t bar_step;

    static void init(size_t total_iterations, size_t bar_len = 50);
    static void iter() ;
};


};
