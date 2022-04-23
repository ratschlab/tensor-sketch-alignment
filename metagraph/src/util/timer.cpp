#include <omp.h>
#include <vector>
#include "timer.hpp"

namespace ts {

using namespace std::chrono;

std::vector<std::map<std::string, nanoseconds>> Timer::durations_vec
        = std::vector<std::map<std::string, nanoseconds>>(100);
std::vector<std::map<std::string, size_t>> Timer::counts_vec
        = std::vector<std::map<std::string, size_t>>(100);


void Timer::add_duration(const std::string &func_name, nanoseconds dur) {
    int tid = omp_get_thread_num();
    auto &durations = durations_vec[tid];

    if (durations.find(func_name) == durations.end()) { // doesn't contain `func_name`
        durations[func_name] = dur;
        counts_vec[tid][func_name] = 1;
    } else {
        durations[func_name] += dur;
        counts_vec[tid][func_name]++;
    }
}


std::string Timer::summary() {
    std::map<std::string, std::string> trans= {
            { "edit_distance", "ED" },
            { "minhash", "MH" },
            { "weighted_minhash", "WMH" },
            { "ordered_minhash", "OMH" },
            { "tensor_sketch", "TS" },
            { "tensor_slide_sketch", "TSS" },
            {"Int32Flattener", "I32FLAT"},
            {"DoubleFlattener", "FLAT"},
            {"seq2kmer", "S2K"}
    };
    std::map<std::string, size_t> total_counts;
    for (auto &counts : counts_vec) {
        for (auto const &[arg_name, arg_count] : counts) {
            if (total_counts.find(arg_name) != total_counts.end())
                total_counts[arg_name] += arg_count;
            else
                total_counts[arg_name] = arg_count;
        }
    }
    std::map<std::string, double> acc;
    for (auto &durations : Timer::durations_vec) {
        for (auto const &[arg_name, arg] : durations) {
            if (acc.find(arg_name) != acc.end()) {
                acc[arg_name] += arg.count();
            } else {
                acc[arg_name] = arg.count();
            }
        }
    }

    std::string str = "long name,short name, time, time sketch, time dist\n";

    for (auto const &[arg_name, arg] : acc) {
        double sk_time = (double)arg, dist_time;
        if (arg_name.find("hash") != std::string::npos && // contains *hash*
            arg_name.find("dist") == std::string::npos) { // doesn't contain *dist*
            sk_time += acc["seq2kmer"]; // add kmer computation time to MH* methods
        }
        if (arg_name == "edit_distance") {
            sk_time = sk_time /1e6/total_counts[arg_name];
            str += arg_name + "," + trans[arg_name] + "," + std::to_string(sk_time) + ",0,0\n";
        } else if (arg_name.find("dist") == std::string::npos && arg_name!="seq2kmer") {
            sk_time = sk_time /1e6/total_counts[arg_name] ;    // mean sketching time (ms)
            dist_time = acc[arg_name + "_dist"]/1e6/total_counts[arg_name + "_dist"]; // mean distance computation time (ms)
            str += arg_name + "," + trans[arg_name] +
                  "," + std::to_string(sk_time + dist_time) +
                  "," + std::to_string(sk_time) +
                  "," + std::to_string(dist_time) + '\n';
        }
    }
    return str;
}

} // namespace ts
