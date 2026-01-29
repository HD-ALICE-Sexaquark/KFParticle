#pragma once

#include <format>
#include <print>
#include <string_view>

#include <armadillo>

namespace KF::Utils {

template <size_t N>
arma::mat::fixed<N, N> StdVec2ArmaSymMat(const std::vector<double> &lower_tri) {

    std::size_t size{lower_tri.size()};
    if (N * (N + 1) / 2 != size) {
        throw std::invalid_argument("std::vector size must be a triangular number N*(N+1)/2");
    }

    arma::mat::fixed<N, N> lower;

    std::size_t idx{0};
    for (std::size_t i{0}; i < N; ++i) {
        for (std::size_t j{0}; j <= i; ++j) {
            lower(i, j) = static_cast<double>(lower_tri[idx++]);
        }
    }

    return arma::symmatl(lower);
}

[[maybe_unused]] static void PrintDouble(std::string_view fcn_name, std::string_view name, double val) {
    std::println(stdout, "({}) {} = {:13.6e}", fcn_name, name, val);
}

template <unsigned long long N>
[[maybe_unused]] static void Print(std::string_view fcn_name, std::string_view name, const arma::vec::fixed<N> &vec) {
    vec.t().print(std::format("({}) {} =", fcn_name, name));
}

template <unsigned long long N, unsigned long long M>
[[maybe_unused]] static void Print(std::string_view fcn_name, std::string_view name, const arma::mat::fixed<N, M> &mat) {
    mat.print(std::format("({}) {} =", fcn_name, name));
}

}  // namespace KF::Utils
