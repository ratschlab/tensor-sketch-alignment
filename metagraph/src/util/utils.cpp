#include "utils.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <gflags/gflags.h>
#include <numeric>

namespace ts {

std::string flag_values(char delimiter, bool skip_empty, bool include_flagfile) {
    const std::string short_letters = "fkmstw";
    std::vector<gflags::CommandLineFlagInfo> flags;
    gflags::GetAllFlags(&flags);
    std::string result;
    for (const auto &flag : flags) {
        if (skip_empty && flag.current_value.empty())
            continue;
        if (!include_flagfile && flag.name == "flagfile")
            continue;
        // Exclude short name flags.
        if (flag.name.size() == 1 && short_letters.find(flag.name[0]) != std::string::npos)
            continue;
        result += "--" + flag.name + "=" + flag.current_value + delimiter;
    }
    return result;
}

void write_output_meta() {
    std::string output_path;
    if (!gflags::GetCommandLineOption("o", &output_path))
        return;

    std::string meta_path = output_path + ".meta";
    std::ofstream meta(meta_path);
    meta << "#!/bin/sh\n";
    meta << "cd " << std::filesystem::current_path() << "\n";
    meta << gflags::GetArgv0() << " " << flag_values(' ', true, false) << "\n";

    std::filesystem::permissions(meta_path, std::filesystem::perms::owner_exec,
                                 std::filesystem::perm_options::add);
}

std::pair<double, double> avg_stddev(const std::vector<double> &v) {
    if (v.empty())
        return { 0, 0 };
    const double sum = std::accumulate(begin(v), end(v), 0.0);
    const double avg = sum / v.size();

    double var = 0;
    for (const auto &x : v)
        var += (x - avg) * (x - avg);

    return { avg, sqrt(var / v.size()) };
}

double median(const std::vector<double> &v) {
    assert(!v.empty());
    if (v.size() % 2)
        return v[v.size() / 2];
    return (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
}

} // namespace ts
