#pragma once

#include <chrono>
#include <map>
#include <string>
#include <utility>
#include <cassert>
#include <vector>


namespace ts { // ts = Tensor Sketch

using namespace std::chrono;


class Timer {
  public:
    Timer(std::string name) :
            name(std::move(name)),
            birth(high_resolution_clock::now()){}

    Timer(const Timer &tt) :
            name(tt.name),
            birth(high_resolution_clock::now()){}
    ~Timer() {
        auto dur = high_resolution_clock::now() - birth;
        Timer::add_duration(name, dur);
    }

    static std::string summary();

  private:
    static void add_duration(const std::string &func_name, std::chrono::nanoseconds dur);


    std::string name;
    high_resolution_clock::time_point birth;

    static std::vector<std::map<std::string, std::chrono::nanoseconds>> durations_vec;
    static std::vector<std::map<std::string, size_t>> counts_vec;
};

} // namespace ts
