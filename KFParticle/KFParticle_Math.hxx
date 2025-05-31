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
}  // namespace Const

template <size_t N>
using Vector = std::array<double, N>;

template <size_t N, size_t M>
using Matrix = std::array<std::array<double, N>, M>;

template <size_t N>
using SymMatrix = std::array<double, N *(N + 1) / 2>;

// Based on https://stackoverflow.com/a/64247207
template <class S>
inline std::pair<S, S> sincos(S arg) {
    return {std::sin(arg), std::cos(arg)};
}

// Return the dot product of vector `vec` with itself.
inline double squaredNorm(const Vector<3> &vec) { return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]; }

template <class S>
static void PrintValue(std::string_view fcn_name, std::string_view name, const S arr) {
    std::cout << "(" << fcn_name << ") " << name << " = " << arr << '\n';
}

template <int N>
static void PrintVector(std::string_view name, const Vector<N> &arr) {
    std::cout << name << " = ";
    for (int i{0}; i < N; ++i) {
        std::cout << arr[i];
        if (i + 1 < N)
            std::cout << "    ";
        else
            std::cout << '\n';
    }
}

template <int N>
static void PrintSymMatrix(std::string_view name, const SymMatrix<N> &arr) {
    std::cout << name << " =\n";
    int n_in_row{0};
    int max_n_row{1};
    for (int i{0}; i < N * (N + 1) / 2; ++i) {
        std::cout << arr[i];
        n_in_row++;
        if (n_in_row == max_n_row) {
            std::cout << '\n';
            n_in_row = 0;
            max_n_row++;
        } else {
            std::cout << "    ";
        }
    }
}

template <int N, int M>
static void PrintMatrix(std::string_view name, const Matrix<N, M> &arr) {
    std::cout << name << " =\n";
    for (int i{0}; i < N; ++i) {
        for (int j{0}; j < M; ++j) {
            std::cout << arr[i][j];
            if (j + 1 < M)
                std::cout << "    ";
            else
                std::cout << '\n';
        }
    }
}

template <int N>
static void PrintJoinedMatrix(std::string_view name, const Vector<N> &arr1, const Vector<N> &arr2, const Vector<N> &arr3) {
    std::cout << name << " =\n";
    for (int i{0}; i < N; ++i) {
        std::cout << "  " << arr1[i] << "    " << arr2[i] << "    " << arr3[i] << '\n';
    }
}

// Convert a pair of indices {i,j} of the covariance matrix to one index corresponding to the triangular form
inline int IJ(int i, int j) { return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }

// Symmetric 3x3 matrix a using modified Cholesky decomposition. The result is stored to the same matrix a.
// \param[in,out] a - 3x3 symmetric matrix
inline void InvertCholesky3(SymMatrix<3> &a) {

    Vector<3> d{};
    Matrix<3, 3> u{};

    for (int i{0}; i < 3; ++i) {
        double uud{0.};
        for (int j{0}; j < i; ++j) uud += u[j][i] * u[j][i] * d[j];
        uud = a[i * (i + 3) / 2] - uud;

        if (std::abs(uud) < Const::AbsAlmostZero) uud = Const::AbsAlmostZero;

        d[i] = uud / std::abs(uud);
        u[i][i] = std::sqrt(std::abs(uud));

        for (int j{i + 1}; j < 3; ++j) {
            uud = 0.;
            for (int k{0}; k < i; ++k) uud += u[k][i] * u[k][j] * d[k];
            uud = a[j * (j + 1) / 2 + i] - uud;
            u[i][j] = d[i] / u[i][i] * uud;
        }
    }

    Vector<3> u1{};

    for (int i{0}; i < 3; ++i) {
        u1[i] = u[i][i];
        u[i][i] = 1 / u[i][i];
    }
    for (int i{0}; i < 2; ++i) {
        u[i][i + 1] = -u[i][i + 1] * u[i][i] * u[i + 1][i + 1];
    }
    for (int i{0}; i < 1; ++i) {
        u[i][i + 2] = u[i][i + 1] * u1[i + 1] * u[i + 1][i + 2] - u[i][i + 2] * u[i][i] * u[i + 2][i + 2];
    }

    for (int i{0}; i < 3; ++i) a[i + 3] = u[i][2] * u[2][2] * d[2];
    for (int i{0}; i < 2; ++i) a[i + 1] = u[i][1] * u[1][1] * d[1] + u[i][2] * u[1][2] * d[2];
    a[0] = u[0][0] * u[0][0] * d[0] + u[0][1] * u[0][1] * d[1] + u[0][2] * u[0][2] * d[2];
}

// Matrix multiplication SOut = Q*S*Q^T, where Q - square matrix, S - symmetric matrix.
// \param[in] Q - square matrix
// \param[in] S - input symmetric matrix
template <int N>
static SymMatrix<N> MultQSQt(const Matrix<N, N> &Q, const SymMatrix<N> &S) {

    Matrix<N, N> SQT{};

    for (int i{0}; i < N; ++i) {
        for (int j{0}; j < N; ++j) {
            for (int k{0}; k < N; ++k) SQT[i][j] += S[IJ(i, k)] * Q[j][k];
        }
    }

    SymMatrix<N> SOut{};

    for (int i{0}; i < N; ++i) {
        for (int j{0}; j <= i; ++j) {
            for (int k{0}; k < N; ++k) SOut[IJ(i, j)] += Q[i][k] * SQT[k][j];
        }
    }

    return SOut;
}

}  // namespace KF

#endif  // KFPARTICLE_MATH_HXX
