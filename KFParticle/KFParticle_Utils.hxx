#ifndef KFPARTICLE_UTILS_HXX
#define KFPARTICLE_UTILS_HXX

#include <cstddef>
#include <format>
#include <print>
#include <string_view>

#include "KFParticle_Math.hxx"

template <size_t N>
struct std::formatter<KF::Vector<N>> {
    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
    auto format(const KF::Vector<N> &vec, std::format_context &ctx) const {
        auto out = ctx.out();
        out = std::format_to(out, "(");
        for (size_t i{0}; i < N; ++i) {
            out = std::format_to(out, "{:13.6e}", vec[i]);
            if (i < N - 1) out = std::format_to(out, ", ");
        }
        out = std::format_to(out, ")");
        return out;
    }
};

template <size_t K>
struct std::formatter<KF::SymMatrix<K>> {
    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
    auto format(const KF::SymMatrix<K> &arr, std::format_context &ctx) const {
        auto out = ctx.out();
        out = std::format_to(out, "\n");
        size_t n_in_row{0};
        size_t max_n_row{1};
        for (size_t i{0}; i < K * (K + 1) / 2; ++i) {
            out = std::format_to(out, "{:13.6e}", arr[i]);
            ++n_in_row;
            if (n_in_row == max_n_row) {
                if (i < (K * (K + 1) / 2) - 1) out = std::format_to(out, "\n");
                n_in_row = 0;
                ++max_n_row;
            } else {
                out = std::format_to(out, "   ");
            }
        }
        return out;
    }
};

template <size_t N, size_t M>
struct std::formatter<KF::Matrix<N, M>> {
    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
    auto format(const KF::Matrix<N, M> &arr, std::format_context &ctx) const {
        auto out = ctx.out();
        out = std::format_to(out, "\n");
        for (size_t i{0}; i < N; ++i) {
            for (size_t j{0}; j < M; ++j) {
                out = std::format_to(out, "{:13.6e}", arr[i][j]);
                if (j < M - 1)
                    out = std::format_to(out, "   ");
                else if (i < N - 1)
                    out = std::format_to(out, "\n");
            }
        }
        return out;
    }
};

namespace KF::Utils {

[[maybe_unused]] static void PrintDouble(std::string_view fcn_name, std::string_view name, double val) {
    std::println(stdout, "({}) {} = {:13.6e}", fcn_name, name, val);
}

template <class S>
static void Print(std::string_view fcn_name, std::string_view name, const S &obj) {
    std::println(stdout, "({}) {} = {}", fcn_name, name, obj);
}

template <size_t N>
static void PrintSplitMatrix(std::string_view fcn_name, std::string_view name, const Vector<N> &arr1, const Vector<N> &arr2, const Vector<N> &arr3) {
    std::println(stdout, "({}) {} = ", fcn_name, name);
    for (size_t i{0}; i < N; ++i) {
        std::println(stdout, "{:13.6e}   {:13.6e}   {:13.6e}", arr1[i], arr2[i], arr3[i]);
    }
}

}  // namespace KF::Utils

#endif  // KFPARTICLE_UTILS_HXX
