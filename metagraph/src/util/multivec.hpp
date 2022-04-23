#pragma once

#include "util/transformer.hpp"

#include <functional>
#include <type_traits>
#include <vector>

namespace ts { // ts = Tensor Sketch

template <class T>
using is_u_integral = typename std::enable_if<std::is_unsigned<T>::value>::type;

template <class T>
using Vec2D = std::vector<std::vector<T>>;

template <class T>
using Vec3D = std::vector<Vec2D<T>>;

template <class T>
using Vec4D = std::vector<Vec3D<T>>;


template <class T>
auto new2D(size_t d1, size_t d2, T val = 0) {
    return Vec2D<T>(d1, std::vector<T>(d2, val));
}
template <class T>
auto new3D(size_t d1, size_t d2, size_t d3, T val = 0) {
    return Vec3D<T>(d1, new2D(d2, d3, val));
}

template <class T>
void apply(std::vector<T> &vec, const transformer<T> &tr) {
    for (auto &v : vec) {
        v = tr.transform(v);
    }
}

template <class T>
void apply(Vec2D<T> &vec2D, const transformer<T> &tr) {
    for (auto &vec : vec2D) {
        apply(vec, tr);
    }
}

template <class T>
void apply(Vec3D<double> &vec3D, const transformer<T> &tr) {
    for (auto &vec2D : vec3D) {
        apply(vec2D, tr);
    }
}

} // namespace ts
