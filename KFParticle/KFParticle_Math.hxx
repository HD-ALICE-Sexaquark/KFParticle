#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

#include "KFParticle_Const.hxx"

namespace KF {

template <size_t N>
struct alignas(KF_SIMD_ALIGN) Vector : std::array<double, N> {};

template <size_t Rows, size_t Cols>
struct alignas(KF_SIMD_ALIGN) Matrix : std::array<std::array<double, Cols>, Rows> {};

template <size_t K>
struct alignas(KF_SIMD_ALIGN) SymMatrix : std::array<double, (K * (K + 1)) / 2> {};

// Convert a pair of indices {i,j} of the covariance matrix to one index corresponding to the triangular form
template <typename D>
constexpr D IJ(D i, D j) {
    return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
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

// Return a symmetric 3x3 matrix using a modified Cholesky decomposition.
// Input:
// - `in` : a symmetric matrix (3x3)
constexpr SymMatrix<3> InvertCholesky3(const SymMatrix<3> &in) {

    SymMatrix<3> out{};
    Vector<3> d{};
    Matrix<3, 3> u{};

    for (size_t i{0}; i < 3; ++i) {
        double uud{0.};
        for (size_t j{0}; j < i; ++j) uud += u[j][i] * u[j][i] * d[j];
        uud = in[i * (i + 3) / 2] - uud;

        if (std::abs(uud) < Const::AbsAlmostZero) uud = Const::AbsAlmostZero;

        d[i] = uud / std::abs(uud);
        u[i][i] = std::sqrt(std::abs(uud));

        for (size_t j{i + 1}; j < 3; ++j) {
            uud = 0.;
            for (size_t k{0}; k < i; ++k) uud += u[k][i] * u[k][j] * d[k];
            uud = in[j * (j + 1) / 2 + i] - uud;
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

    for (size_t i{0}; i < 3; ++i) out[i + 3] = u[i][2] * u[2][2] * d[2];
    for (size_t i{0}; i < 2; ++i) out[i + 1] = u[i][1] * u[1][1] * d[1] + u[i][2] * u[1][2] * d[2];
    out[0] = u[0][0] * u[0][0] * d[0] + u[0][1] * u[0][1] * d[1] + u[0][2] * u[0][2] * d[2];

    return out;
}

// Extract a sub-vector (size M) from a vector `in` (size N) starting from index 0
template <size_t N, size_t M>
[[nodiscard]] inline Vector<M> Slice(const Vector<N> &in) {
    Vector<M> out{};
    for (size_t i{0}; i < M; ++i) {
        if (i < N) out[i] = in[i];
    }
    return out;
}

// Extract a sub-matrix (LxL symmetric) from a matrix `in` (KxK symmetric) starting from index 0
template <size_t K, size_t L>
[[nodiscard]] inline SymMatrix<L> Slice(const SymMatrix<K> &in) {
    SymMatrix<L> out{};
    constexpr size_t in_size{K * (K + 1) / 2};
    constexpr size_t out_size{L * (L + 1) / 2};
    for (size_t i{0}; i < out_size; ++i) {
        if (i < in_size) out[i] = in[i];
    }
    return out;
}

// Extract a sub-matrix (NxM) from a matrix `in` (LxK) starting from element (0,0)
template <size_t L, size_t K, size_t N, size_t M>
[[nodiscard]] inline Matrix<N, M> Slice(const Matrix<L, K> &in) {
    Matrix<N, M> out{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            if (i < L && j < K) out[i][j] = in[i][j];
        }
    }
    return out;
}

// Standard matrix multiplication
// C (NxK) = A (NxM) x B (MxK)
template <size_t N, size_t M, size_t K>
[[nodiscard]] constexpr Matrix<N, K> MultiplyMatrices(const Matrix<N, M> &A, const Matrix<M, K> &B) {
    Matrix<N, K> C{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t k{0}; k < M; ++k) {
            const double A_ik{A[i][k]};
            for (size_t j{0}; j < K; ++j) {
                C[i][j] += A_ik * B[k][j];
            }
        }
    }
    return C;
}

// Multiply symmetric matrices
// C (NxN) = S (N) x T (N)
// Note: result is non-symmetric
template <size_t N>
[[nodiscard]] constexpr Matrix<N, N> MultiplySymmetricMatrices(const SymMatrix<N> &S, const SymMatrix<N> &T) {
    Matrix<N, N> C{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t k{0}; k < N; ++k) {
            const double S_ik{S[IJ(i, k)]};
            for (size_t j{0}; j < N; ++j) {
                C[i][j] += S_ik * T[IJ(k, j)];
            }
        }
    }
    return C;
}

// Transpose a matrix
template <size_t N, size_t M>
[[nodiscard]] constexpr Matrix<M, N> Transpose(const Matrix<N, M> &A) {
    Matrix<M, N> B{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            B[j][i] = A[i][j];
        }
    }
    return B;
}

// Return square identity matrix
template <size_t N>
[[nodiscard]] constexpr Matrix<N, N> Identity() {
    Matrix<N, N> I{};
    for (size_t i{0}; i < N; ++i) I[i][i] = 1.;
    return I;
}

// Return addition of symmetric matrices:
// R (N) = S (N) + T (N)
template <size_t N>
[[nodiscard]] constexpr SymMatrix<N> AddMatrices(const SymMatrix<N> &S, const SymMatrix<N> &T) {
    SymMatrix<N> R{};
    constexpr size_t sym_size{(N * (N + 1)) / 2};
    for (size_t i{0}; i < sym_size; ++i) R[i] = S[i] + T[i];
    return R;
}

// Return addition of matrices:
// C = alpha * A + beta * B,
// where alpha and beta are scalars
template <size_t N, size_t M>
[[nodiscard]] constexpr Matrix<N, M> AddMatrices(const Matrix<N, M> &A, const Matrix<N, M> &B, double alpha = 1., double beta = 1.) {
    Matrix<N, M> C{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            C[i][j] = alpha * A[i][j] + beta * B[i][j];
        }
    }
    return C;
}

// Return a scaled matrix: B = scalar * A
template <size_t N, size_t M>
[[nodiscard]] constexpr Matrix<N, M> ScaleMatrix(const Matrix<N, M> &A, double scalar = 1.) {
    Matrix<N, M> B{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            B[i][j] = scalar * A[i][j];
        }
    }
    return B;
}

// Multiply symmetric matrix S with non-symmetric matrix A
// B (NxM) = S (N) x A (NxM)
template <size_t N, size_t M>
[[nodiscard]] constexpr Matrix<N, M> MultiplySymmWithNonSymm(const SymMatrix<N> &S, const Matrix<N, M> &A) {
    Matrix<N, M> B{};
    for (size_t k{0}; k < N; ++k) {
        for (size_t i{0}; i < N; ++i) {
            const double S_ik{S[IJ(i, k)]};
            for (size_t j{0}; j < M; ++j) {
                B[i][j] += S_ik * A[k][j];
            }
        }
    }
    return B;
}

// Multiply non-symmetric matrix A with symmetric matrix S on the right
// B (NxM) = A (NxM) × S (M)
template <size_t N, size_t M>
[[nodiscard]] constexpr Matrix<N, M> MultiplyNonSymmWithSymm(const Matrix<N, M> &A, const SymMatrix<M> &S) {
    Matrix<N, M> B{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            for (size_t k{0}; k < M; ++k) {
                B[i][j] += A[i][k] * S[IJ(k, j)];
            }
        }
    }
    return B;
}

// Return matrix multiplication Q x S x Q^T.
// Input arguments:
// - `Q` : square matrix
// - `S` : symmetric matrix
template <size_t N>
inline SymMatrix<N> MultiplyQSQT(const Matrix<N, N> &Q, const SymMatrix<N> &S) {
    Matrix<N, N> SQT{Math::MultiplySymmWithNonSymm(S, Math::Transpose(Q))};
    SymMatrix<N> QSQT{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j <= i; ++j) {  // only lower triangle
            double sum{0.};
            for (size_t k{0}; k < N; ++k) {
                sum += Q[i][k] * SQT[k][j];
            }
            QSQT[IJ(i, j)] = sum;
        }
    }
    return QSQT;
}

}  // namespace Math

}  // namespace KF
