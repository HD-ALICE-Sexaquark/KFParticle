#pragma once

#include <format>
#include <print>
#include <string_view>

#include <armadillo>

namespace KF::Utils {

inline arma::mat vec_to_symmat(const std::vector<double> &lower_tri) {

    // Solve N*(N+1)/2 = size for N
    // N = (-1 + sqrt(1 + 8*size)) / 2
    std::size_t size{lower_tri.size()};
    double n_double{(-1. + std::sqrt(1. + 8. * double(size))) / 2.};
    std::size_t n{static_cast<std::size_t>(std::round(n_double))};

    // Validate that size is actually a triangular number
    if (n * (n + 1) / 2 != size) {
        throw std::invalid_argument("Vector size must be a triangular number N*(N+1)/2");
    }

    arma::mat lower(n, n, arma::fill::zeros);

    // Fill lower triangle row by row
    // Elements are assumed in order: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...
    std::size_t idx = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            lower(i, j) = static_cast<double>(lower_tri[idx++]);
        }
    }

    return arma::symmatl(lower);
}

[[maybe_unused]] static void PrintDouble(std::string_view fcn_name, std::string_view name, double val) {
    std::println(stdout, "({}) {} = {:13.6e}", fcn_name, name, val);
}

template <class S>
static void Print(std::string_view fcn_name, std::string_view name, const S &obj) {
    std::println(stdout, "({}) {} = {}", fcn_name, name, obj);
}

[[maybe_unused]] static void Print(std::string_view fcn_name, std::string_view name, const arma::vec &vec) {
    vec.t().print(std::format("({}) {} =", fcn_name, name));
}

[[maybe_unused]] static void Print(std::string_view fcn_name, std::string_view name, const arma::mat &mat) {
    mat.print(std::format("({}) {} =", fcn_name, name));
}

}  // namespace KF::Utils
