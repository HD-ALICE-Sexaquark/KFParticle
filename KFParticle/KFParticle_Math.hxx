#ifndef KFPARTICLE_MATH_HXX
#define KFPARTICLE_MATH_HXX

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>

namespace KF {

namespace Const {
constexpr double Kappa{0.000299792458};  // (GeV/c) / (kG/cm)
constexpr double AbsAlmostZero{1.E-8};
constexpr double BigNumber{1.E8};
constexpr double Epsilon{1.E-6};
}  // namespace Const

template <size_t N>
using Vector = std::array<double, N>;

template <size_t N, size_t M>
using Matrix = std::array<std::array<double, N>, M>;

template <size_t N>
using SymMatrix = std::array<double, N *(N + 1) / 2>;

template <class S>
static void PrintValue(std::string_view fcn_name, std::string_view name, const S arr) {
    std::cout << "(" << fcn_name << ") " << name << " = " << arr << '\n';
}

template <size_t N>
static void PrintVector(std::string_view fcn_name, std::string_view name, const Vector<N> &arr) {
    std::cout << "(" << fcn_name << ") " << name << " = ";
    for (size_t i{0}; i < N; ++i) {
        std::cout << arr[i];
        if (i + 1 < N)
            std::cout << "    ";
        else
            std::cout << '\n';
    }
}

template <size_t N>
static void PrintSymMatrix(std::string_view fcn_name, std::string_view name, const SymMatrix<N> &arr) {
    std::cout << "(" << fcn_name << ") " << name << " =\n";
    size_t n_in_row{0};
    size_t max_n_row{1};
    for (size_t i{0}; i < N * (N + 1) / 2; ++i) {
        std::cout << arr[i];
        ++n_in_row;
        if (n_in_row == max_n_row) {
            std::cout << '\n';
            n_in_row = 0;
            ++max_n_row;
        } else {
            std::cout << "    ";
        }
    }
}

template <size_t N, size_t M>
static void PrintMatrix(std::string_view fcn_name, std::string_view name, const Matrix<N, M> &arr) {
    std::cout << "(" << fcn_name << ") " << name << " =\n";
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            std::cout << arr[i][j];
            if (j + 1 < M)
                std::cout << "    ";
            else
                std::cout << '\n';
        }
    }
}

template <size_t N>
static void PrintSplitMatrix(std::string_view fcn_name, std::string_view name, const Vector<N> &arr1, const Vector<N> &arr2, const Vector<N> &arr3) {
    std::cout << "(" << fcn_name << ") " << name << " =\n";
    for (size_t i{0}; i < N; ++i) {
        std::cout << "  " << arr1[i] << "    " << arr2[i] << "    " << arr3[i] << '\n';
    }
}

// Convert a pair of indices {i,j} of the covariance matrix to one index corresponding to the triangular form
template <typename D>
constexpr D IJ(D i, D j) {
    return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
}

template <size_t N, size_t M>
inline Vector<M> Slice(const Vector<N> &in, size_t begin_i = 0) {
    Vector<M> out{};
    for (size_t i{0}; i < M; ++i) {
        const size_t src{begin_i + i};
        if (src < N) out[i] = in[src];
    }
    return out;
}

template <size_t L, size_t K, size_t N, size_t M>
inline Matrix<N, M> Slice(const Matrix<L, K> &in, size_t begin_i = 0, size_t begin_j = 0) {
    Matrix<N, M> out{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            const size_t src_i{begin_i + i};
            const size_t src_j{begin_j + j};
            if (src_i < L && src_j < K) out[i][j] = in[src_i][src_j];
        }
    }
    return out;
}

namespace Math {

// Based on https://stackoverflow.com/a/64247207
template <class S>
constexpr std::pair<S, S> sincos(S arg) {
    return {std::sin(arg), std::cos(arg)};
}

// Return the dot product of vector `vec` with itself.
template <size_t N>
constexpr double SquaredNorm(const Vector<N> &vec) {
    double sum{0.};
    for (size_t i{0}; i < N; ++i) sum += vec[i] * vec[i];
    return sum;
}

// Return the norm of vector `vec`. It's equivalent to the square root of the dot product of vector `vec` with itself.
template <size_t N>
constexpr double Norm(const Vector<N> &vec) {
    return std::sqrt(SquaredNorm(vec));
}

// Symmetric 3x3 matrix a using modified Cholesky decomposition. The result is stored to the same matrix a.
// \param[in,out] a - 3x3 symmetric matrix
inline void InvertCholesky3(SymMatrix<3> &a) {

    Vector<3> d{};
    Matrix<3, 3> u{};

    for (size_t i{0}; i < 3; ++i) {
        double uud{0.};
        for (size_t j{0}; j < i; ++j) uud += u[j][i] * u[j][i] * d[j];
        uud = a[i * (i + 3) / 2] - uud;

        if (std::abs(uud) < Const::AbsAlmostZero) uud = Const::AbsAlmostZero;

        d[i] = uud / std::abs(uud);
        u[i][i] = std::sqrt(std::abs(uud));

        for (size_t j{i + 1}; j < 3; ++j) {
            uud = 0.;
            for (size_t k{0}; k < i; ++k) uud += u[k][i] * u[k][j] * d[k];
            uud = a[j * (j + 1) / 2 + i] - uud;
            u[i][j] = d[i] / u[i][i] * uud;
        }
    }

    Vector<3> u1{};

    for (size_t i{0}; i < 3; ++i) {
        u1[i] = u[i][i];
        u[i][i] = 1 / u[i][i];
    }
    for (size_t i{0}; i < 2; ++i) {
        u[i][i + 1] = -u[i][i + 1] * u[i][i] * u[i + 1][i + 1];
    }
    for (size_t i{0}; i < 1; ++i) {
        u[i][i + 2] = u[i][i + 1] * u1[i + 1] * u[i + 1][i + 2] - u[i][i + 2] * u[i][i] * u[i + 2][i + 2];
    }

    for (size_t i{0}; i < 3; ++i) a[i + 3] = u[i][2] * u[2][2] * d[2];
    for (size_t i{0}; i < 2; ++i) a[i + 1] = u[i][1] * u[1][1] * d[1] + u[i][2] * u[1][2] * d[2];
    a[0] = u[0][0] * u[0][0] * d[0] + u[0][1] * u[0][1] * d[1] + u[0][2] * u[0][2] * d[2];
}

// Return matrix multiplication Q * S * Q^T.
// Input arguments:
// - `Q` : square matrix
// - `S` : input symmetric matrix
template <size_t N>
static SymMatrix<N> MultQSQt(const Matrix<N, N> &Q, const SymMatrix<N> &S) {

    Matrix<N, N> SQT{};

    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < N; ++j) {
            for (size_t k{0}; k < N; ++k) SQT[i][j] += S[IJ(i, k)] * Q[j][k];
        }
    }

    SymMatrix<N> SOut{};

    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j <= i; ++j) {
            for (size_t k{0}; k < N; ++k) SOut[IJ(i, j)] += Q[i][k] * SQT[k][j];
        }
    }

    return SOut;
}

}  // namespace Math

}  // namespace KF

#endif  // KFPARTICLE_MATH_HXX
