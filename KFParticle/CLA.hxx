#pragma once

#include <array>
#include <cmath>
#include <cstdlib>
#include <format>

#if defined(__AVX512F__)
#define CLA_SIMD_ALIGN 64
#elif defined(__AVX2__)
#define CLA_SIMD_ALIGN 32
#else
#define CLA_SIMD_ALIGN alignof(double)
#endif

namespace CompactLinearAlgebra {

// ## Vector ## //

template <size_t N>
struct alignas(CLA_SIMD_ALIGN) Vector : std::array<double, N> {

    // Access operators.
    constexpr double &operator()(size_t i) { return (*this)[i]; }
    constexpr double operator()(size_t i) const { return (*this)[i]; }

    // Addition with another vector.
    [[nodiscard]] constexpr Vector<N> operator+(const Vector<N> &rhs) const {
        Vector<N> out{};
        for (size_t i{0}; i < N; ++i) out[i] = (*this)[i] + rhs[i];
        return out;
    }

    // Compound addition with another vector.
    constexpr Vector<N> &operator+=(const Vector<N> &rhs) {
        for (size_t i{0}; i < N; ++i) (*this)[i] += rhs[i];
        return *this;
    }

    // Subtraction.
    [[nodiscard]] constexpr Vector<N> operator-(const Vector<N> &rhs) const {
        Vector<N> out{};
        for (size_t i{0}; i < N; ++i) out[i] = (*this)[i] - rhs[i];
        return out;
    }

    // Compound subtraction.
    constexpr Vector<N> &operator-=(const Vector<N> &rhs) {
        for (size_t i{0}; i < N; ++i) (*this)[i] -= rhs[i];
        return *this;
    }

    // Unary negation.
    [[nodiscard]] constexpr Vector<N> operator-() const {
        Vector<N> out{};
        for (size_t i{0}; i < N; ++i) out[i] = -(*this)[i];
        return out;
    }

    // Multiplication with scalar (rhs).
    [[nodiscard]] constexpr Vector<N> operator*(double scalar) const {
        Vector<N> out{};
        for (size_t i{0}; i < N; ++i) out[i] = (*this)[i] * scalar;
        return out;
    }

    // Compound multiplication with scalar.
    constexpr Vector<N> &operator*=(double scalar) {
        for (size_t i{0}; i < N; ++i) (*this)[i] *= scalar;
        return *this;
    }

    // Division with scalar.
    [[nodiscard]] constexpr Vector<N> operator/(double scalar) const {
        Vector<N> out{};
        for (size_t i{0}; i < N; ++i) out[i] = (*this)[i] / scalar;
        return out;
    }

    // Compound division with scalar.
    constexpr Vector<N> &operator/=(double scalar) {
        for (size_t i{0}; i < N; ++i) (*this)[i] /= scalar;
        return *this;
    }

    // Extract a vector of size `M` from starting position `Pos`.
    template <size_t M, size_t Pos = 0>
    constexpr Vector<M> GetSlice() const {
        static_assert(Pos + M <= N, "GetSlice() exceeds vector bounds");
        Vector<M> out{};
        for (size_t i{0}; i < M; ++i) out[i] = (*this)[Pos + i];
        return out;
    }

    // In the vector, modify a sub-vector of size `M` starting from position `Pos`.
    template <size_t M, size_t Pos = 0>
    constexpr void SetSlice(const Vector<M> &src) {
        static_assert(Pos + M <= N, "SetSlice() exceeds vector bounds");
        for (size_t i{0}; i < M; ++i) (*this)[Pos + i] = src[i];
    }

    // In the vector, add to sub-vector of size `M` starting from position `Pos`.
    template <size_t M, size_t Pos = 0>
    constexpr void AddToSlice(const Vector<M> &src) {
        static_assert(Pos + M <= N, "AddToSlice() exceeds vector bounds");
        for (size_t i{0}; i < M; ++i) (*this)[Pos + i] += src[i];
    }

    // In the vector, subtract from sub-vector of size `M` starting from position `Pos`.
    template <size_t M, size_t Pos = 0>
    constexpr void SubtractFromSlice(const Vector<M> &src) {
        static_assert(Pos + M <= N, "SubtractFromSlice() exceeds vector bounds");
        for (size_t i{0}; i < M; ++i) (*this)[Pos + i] -= src[i];
    }

    // Return the dot product of the vector with itself.
    [[nodiscard]] constexpr double SquaredNorm() const {
        double sum{0.};
        for (size_t i{0}; i < N; ++i) sum += (*this)[i] * (*this)[i];
        return sum;
    }

    // Return the norm.
    // It's equivalent to the square root of the dot product of the vector with itself.
    [[nodiscard]] double Norm() const { return std::sqrt(SquaredNorm()); }
};

// Vector -- Free Functions //

// Vector multiplication with scalar (lhs).
template <size_t N>
[[nodiscard]] constexpr Vector<N> operator*(double scalar, const Vector<N> &vec) {
    return vec * scalar;
}

// Dot product between vectors.
template <size_t N>
[[nodiscard]] constexpr double Dot(const Vector<N> &a, const Vector<N> &b) {
    double sum{0.};
    for (size_t i{0}; i < N; ++i) sum += a[i] * b[i];
    return sum;
}

// ## Matrix ## //

template <size_t Rows, size_t Cols>
struct alignas(CLA_SIMD_ALIGN) Matrix : std::array<std::array<double, Cols>, Rows> {

    // Access operators.
    constexpr double operator()(size_t i, size_t j) const { return (*this)[i][j]; }
    constexpr double &operator()(size_t i, size_t j) { return (*this)[i][j]; }

    // Addition.
    [[nodiscard]] constexpr Matrix<Rows, Cols> operator+(const Matrix<Rows, Cols> &rhs) const {
        Matrix<Rows, Cols> out{};
        for (size_t i{0}; i < Rows; ++i) {
            for (size_t j{0}; j < Cols; ++j) {
                out[i][j] = (*this)[i][j] + rhs[i][j];
            }
        }
        return out;
    }

    // Compound addition.
    constexpr Matrix<Rows, Cols> &operator+=(const Matrix<Rows, Cols> &rhs) {
        for (size_t i{0}; i < Rows; ++i) {
            for (size_t j{0}; j < Cols; ++j) {
                (*this)[i][j] += rhs[i][j];
            }
        }
        return *this;
    }

    // Subtraction.
    [[nodiscard]] constexpr Matrix<Rows, Cols> operator-(const Matrix<Rows, Cols> &rhs) const {
        Matrix<Rows, Cols> out{};
        for (size_t i{0}; i < Rows; ++i) {
            for (size_t j{0}; j < Cols; ++j) {
                out[i][j] = (*this)[i][j] - rhs[i][j];
            }
        }
        return out;
    }

    // Compound subtraction.
    constexpr Matrix<Rows, Cols> &operator-=(const Matrix<Rows, Cols> &rhs) {
        for (size_t i{0}; i < Rows; ++i) {
            for (size_t j{0}; j < Cols; ++j) {
                (*this)[i][j] -= rhs[i][j];
            }
        }
        return *this;
    }

    // Unary negation.
    [[nodiscard]] constexpr Matrix<Rows, Cols> operator-() const {
        Matrix<Rows, Cols> out{};
        for (size_t i{0}; i < Rows; ++i) {
            for (size_t j{0}; j < Cols; ++j) {
                out[i][j] = -(*this)[i][j];
            }
        }
        return out;
    }

    // Matrix multiplication with scalar (rhs).
    [[nodiscard]] constexpr Matrix<Rows, Cols> operator*(double scalar) const {
        Matrix<Rows, Cols> out{};
        for (size_t i{0}; i < Rows; ++i) {
            for (size_t j{0}; j < Cols; ++j) {
                out[i][j] = (*this)[i][j] * scalar;
            }
        }
        return out;
    }

    // Compound scalar multiplication.
    constexpr Matrix<Rows, Cols> &operator*=(double scalar) {
        for (size_t i{0}; i < Rows; ++i) {
            for (size_t j{0}; j < Cols; ++j) {
                (*this)[i][j] *= scalar;
            }
        }
        return *this;
    }

    // Get a sub-matrix of size (`B_Rows` x `B_Cols`) starting at position `(RowPos, ColPos)`.
    template <size_t B_Rows, size_t B_Cols, size_t RowPos = 0, size_t ColPos = 0>
    constexpr Matrix<B_Rows, B_Cols> GetSlice() const {
        static_assert(RowPos + B_Rows <= Rows, "GetSlice() exceeds matrix row bounds");
        static_assert(ColPos + B_Cols <= Cols, "GetSlice() exceeds matrix column bounds");
        Matrix<B_Rows, B_Cols> out{};
        for (size_t i{0}; i < B_Rows; ++i) {
            for (size_t j{0}; j < B_Cols; ++j) {
                out[i][j] = (*this)[RowPos + i][ColPos + j];
            }
        }
        return out;
    }

    // Set a sub-matrix of size (`B_Rows` x `B_Cols`) starting at position `(RowPos, ColPos)`.
    template <size_t B_Rows, size_t B_Cols, size_t RowPos = 0, size_t ColPos = 0>
    constexpr void SetSlice(const Matrix<B_Rows, B_Cols> &src) {
        static_assert(RowPos + B_Rows <= Rows, "SetSlice() exceeds matrix row bounds");
        static_assert(ColPos + B_Cols <= Cols, "SetSlice() exceeds matrix column bounds");
        for (size_t i{0}; i < B_Rows; ++i) {
            for (size_t j{0}; j < B_Cols; ++j) {
                (*this)[RowPos + i][ColPos + j] = src[i][j];
            }
        }
    }

    // Add to a sub-matrix of size (`B_Rows` × `B_Cols`) starting at position `(RowPos, ColPos)`.
    template <size_t B_Rows, size_t B_Cols, size_t RowPos = 0, size_t ColPos = 0>
    constexpr void AddToSlice(const Matrix<B_Rows, B_Cols> &src) {
        static_assert(RowPos + B_Rows <= Rows, "AddToSlice() exceeds matrix row bounds");
        static_assert(ColPos + B_Cols <= Cols, "AddToSlice() exceeds matrix column bounds");
        for (size_t i{0}; i < B_Rows; ++i) {
            for (size_t j{0}; j < B_Cols; ++j) {
                (*this)[RowPos + i][ColPos + j] += src[i][j];
            }
        }
    }

    // Subtract from a sub-matrix of size (`B_Rows` × `B_Cols`) starting at position `(RowPos, ColPos)`.
    template <size_t B_Rows, size_t B_Cols, size_t RowPos = 0, size_t ColPos = 0>
    constexpr void SubtractFromSlice(const Matrix<B_Rows, B_Cols> &src) {
        static_assert(RowPos + B_Rows <= Rows, "SubtractFromSlice() exceeds matrix row bounds");
        static_assert(ColPos + B_Cols <= Cols, "SubtractFromSlice() exceeds matrix column bounds");
        for (size_t i{0}; i < B_Rows; ++i) {
            for (size_t j{0}; j < B_Cols; ++j) {
                (*this)[RowPos + i][ColPos + j] -= src[i][j];
            }
        }
    }

    // Set a row from a vector.
    template <size_t RowIdx>
    constexpr void SetRow(const Vector<Cols> &src) {
        static_assert(RowIdx < Rows, "SetRow() index out of bounds");
        for (size_t j{0}; j < Cols; ++j) (*this)[RowIdx][j] = src[j];
    }

    // Get a row as a vector.
    template <size_t RowIdx>
    [[nodiscard]] Vector<Cols> GetRow() const {
        static_assert(RowIdx < Rows, "GetRow() index out of bounds");
        Vector<Cols> out{};
        for (size_t j{0}; j < Cols; ++j) out[j] = (*this)[RowIdx][j];
        return out;
    }

    // Set a column from a vector.
    template <size_t ColIdx>
    constexpr void SetCol(const Vector<Rows> &src) {
        static_assert(ColIdx < Cols, "SetCol() index out of bounds");
        for (size_t i{0}; i < Rows; ++i) (*this)[i][ColIdx] = src[i];
    }

    // Get a column as a vector.
    template <size_t ColIdx>
    [[nodiscard]] Vector<Rows> GetCol() const {
        static_assert(ColIdx < Cols, "GetCol() index out of bounds");
        Vector<Rows> out{};
        for (size_t i{0}; i < Rows; ++i) out[i] = (*this)[i][ColIdx];
        return out;
    }

    // Transpose matrix.
    [[nodiscard]] constexpr Matrix<Cols, Rows> Transpose() const {
        Matrix<Cols, Rows> out{};
        for (size_t i{0}; i < Rows; ++i) {
            for (size_t j{0}; j < Cols; ++j) {
                out[j][i] = (*this)[i][j];
            }
        }
        return out;
    }
    [[nodiscard]] constexpr Matrix<Cols, Rows> T() const { return Transpose(); }
};

// Matrix -- Free Functions //

// Matrix multiplication with scalar (lhs).
template <size_t Rows, size_t Cols>
[[nodiscard]] constexpr Matrix<Rows, Cols> operator*(double scalar, const Matrix<Rows, Cols> &mat) {
    return mat * scalar;
}

// Matrix multiplication: (NxM) * (MxK) = (NxK).
template <size_t N, size_t M, size_t K>
[[nodiscard]] constexpr Matrix<N, K> operator*(const Matrix<N, M> &A, const Matrix<M, K> &B) {
    Matrix<N, K> out{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t k{0}; k < M; ++k) {
            const double a_ik = A[i][k];
            for (size_t j{0}; j < K; ++j) {
                out[i][j] += a_ik * B[k][j];
            }
        }
    }
    return out;
}

// Return an NxN identity matrix.
template <size_t N>
[[nodiscard]] constexpr Matrix<N, N> Identity() {
    Matrix<N, N> I{};
    for (size_t i{0}; i < N; ++i) I[i][i] = 1.;
    return I;
}

// ## Symmetric Matrix ## //

template <size_t K>
struct alignas(CLA_SIMD_ALIGN) SymMatrix : std::array<double, (K * (K + 1)) / 2> {

    // Single index access.
    constexpr double operator()(size_t idx) const { return (*this)[idx]; }
    constexpr double &operator()(size_t idx) { return (*this)[idx]; }

    // Two-index access (row, col).
    constexpr double operator()(size_t i, size_t j) const { return (j <= i) ? (*this)[i * (i + 1) / 2 + j] : (*this)[j * (j + 1) / 2 + i]; }
    constexpr double &operator()(size_t i, size_t j) { return (j <= i) ? (*this)[i * (i + 1) / 2 + j] : (*this)[j * (j + 1) / 2 + i]; }

    // Addition.
    [[nodiscard]] constexpr SymMatrix<K> operator+(const SymMatrix<K> &rhs) const {
        SymMatrix<K> out{};
        for (size_t i{0}; i < Size; ++i) out[i] = (*this)[i] + rhs[i];
        return out;
    }

    // Compound addition.
    constexpr SymMatrix<K> &operator+=(const SymMatrix<K> &rhs) {
        for (size_t i{0}; i < Size; ++i) (*this)[i] += rhs[i];
        return *this;
    }

    // Subtraction.
    [[nodiscard]] constexpr SymMatrix<K> operator-(const SymMatrix<K> &rhs) const {
        SymMatrix<K> out{};
        for (size_t i{0}; i < Size; ++i) out[i] = (*this)[i] - rhs[i];
        return out;
    }

    // Compound subtraction.
    constexpr SymMatrix<K> &operator-=(const SymMatrix<K> &rhs) {
        for (size_t i{0}; i < Size; ++i) (*this)[i] -= rhs[i];
        return *this;
    }

    // Unary negation.
    [[nodiscard]] constexpr SymMatrix<K> operator-() const {
        SymMatrix<K> out{};
        for (size_t i{0}; i < Size; ++i) out[i] = -(*this)[i];
        return out;
    }

    // Multiplication with scalar (rhs).
    [[nodiscard]] constexpr SymMatrix<K> operator*(double scalar) const {
        SymMatrix<K> out{};
        for (size_t i{0}; i < Size; ++i) out[i] = (*this)[i] * scalar;
        return out;
    }

    // Compound scalar multiplication.
    constexpr SymMatrix<K> &operator*=(double scalar) {
        for (size_t i{0}; i < Size; ++i) (*this)[i] *= scalar;
        return *this;
    }

    // Division with scalar (rhs).
    [[nodiscard]] constexpr SymMatrix<K> operator/(double scalar) const {
        SymMatrix<K> out{};
        for (size_t i{0}; i < Size; ++i) out[i] = (*this)[i] / scalar;
        return out;
    }

    // Compound scalar division.
    constexpr SymMatrix<K> &operator/=(double scalar) {
        for (size_t i{0}; i < Size; ++i) (*this)[i] /= scalar;
        return *this;
    }

    // Extract LxL symmetric sub-matrix starting at diagonal position `DiagPos`.
    template <size_t L, size_t DiagPos = 0>
    constexpr SymMatrix<L> GetSymSlice() const {
        static_assert(DiagPos + L <= K, "GetSlice() exceeds symmetric matrix bounds");
        SymMatrix<L> out{};
        for (size_t i{0}; i < L; ++i) {
            for (size_t j{0}; j <= i; ++j) {
                out(i, j) = (*this)(DiagPos + i, DiagPos + j);
            }
        }
        return out;
    }

    // Assign to LxL symmetric sub-matrix starting at diagonal position `DiagPos`.
    template <size_t L, size_t DiagPos = 0>
    constexpr void SetSymSlice(const SymMatrix<L> &src) {
        static_assert(DiagPos + L <= K, "SetSlice() exceeds symmetric matrix bounds");
        for (size_t i{0}; i < L; ++i) {
            for (size_t j{0}; j <= i; ++j) {
                (*this)(DiagPos + i, DiagPos + j) = src(i, j);
            }
        }
    }

    // Add LxL sym. sub-matrix to current sym. matrix starting at diagonal position `DiagPos`.
    template <size_t L, size_t DiagPos = 0>
    constexpr void AddToSymSlice(const SymMatrix<L> &src) {
        static_assert(DiagPos + L <= K, "AddToSlice() exceeds symmetric matrix bounds");
        for (size_t i{0}; i < L; ++i) {
            for (size_t j{0}; j <= i; ++j) {
                (*this)(DiagPos + i, DiagPos + j) += src(i, j);
            }
        }
    }

    // Subtract LxL sym. sub-matrix from sym. matrix starting at diagonal position `DiagPos`.
    template <size_t L, size_t DiagPos = 0>
    constexpr void SubtractFromSymSlice(const SymMatrix<L> &src) {
        static_assert(DiagPos + L <= K, "SubtractFromSlice() exceeds symmetric matrix bounds");
        for (size_t i{0}; i < L; ++i) {
            for (size_t j{0}; j <= i; ++j) {
                (*this)(DiagPos + i, DiagPos + j) -= src(i, j);
            }
        }
    }

    // Extract a rectangular (`R_Rows` x `R_Cols`) sub-matrix starting at `(RowPos, ColPos)`.
    template <size_t R_Rows, size_t R_Cols, size_t RowPos = 0, size_t ColPos = 0>
    constexpr Matrix<R_Rows, R_Cols> GetRectangularSlice() const {
        static_assert(RowPos + R_Rows <= K, "GetRectangularSlice() exceeds row bounds");
        static_assert(ColPos + R_Cols <= K, "GetRectangularSlice() exceeds column bounds");
        Matrix<R_Rows, R_Cols> out{};
        for (size_t i{0}; i < R_Rows; ++i) {
            for (size_t j{0}; j < R_Cols; ++j) {
                out[i][j] = (*this)(RowPos + i, ColPos + j);
            }
        }
        return out;
    }

    // Total number of elements of the general view of the matrix.
    static constexpr size_t Size = (K * (K + 1)) / 2;
};

// Symmetric Matrix -- Free Functions //

// Multiplication with scalar (lhs).
template <size_t K>
[[nodiscard]] constexpr SymMatrix<K> operator*(double scalar, const SymMatrix<K> &mat) {
    return mat * scalar;
}

// Multiplication between symmetric matrices. Result might not be symmetric.
template <size_t K>
[[nodiscard]] constexpr Matrix<K, K> operator*(const SymMatrix<K> &S, const SymMatrix<K> &T) {
    Matrix<K, K> out{};
    for (size_t i{0}; i < K; ++i) {
        for (size_t k{0}; k < K; ++k) {
            const double s_ik = S(i, k);
            for (size_t j{0}; j < K; ++j) {
                out[i][j] += s_ik * T(k, j);
            }
        }
    }
    return out;
}

// ## Mixed Operations ## //

// Matrix (NxM) * ColVector (Mx1) = ColVector (Nx1)
template <size_t N, size_t M>
[[nodiscard]] constexpr Vector<N> operator*(const Matrix<N, M> &A, const Vector<M> &v) {
    Vector<N> out{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            out[i] += A[i][j] * v[j];
        }
    }
    return out;
}

// RowVector (1xN) * Matrix (N×M) = RowVector (1xM)
template <size_t N, size_t M>
[[nodiscard]] constexpr Vector<M> operator*(const Vector<N> &v, const Matrix<N, M> &A) {
    Vector<M> out{};
    for (size_t j{0}; j < M; ++j) {
        for (size_t i{0}; i < N; ++i) {
            out[j] += v[i] * A[i][j];
        }
    }
    return out;
}

// SymMatrix (KxK) * ColVector (Kx1) = ColVector (Kx1).
template <size_t K>
[[nodiscard]] constexpr Vector<K> operator*(const SymMatrix<K> &S, const Vector<K> &v) {
    Vector<K> out{};
    for (size_t i{0}; i < K; ++i) {
        for (size_t j{0}; j < K; ++j) {
            out[i] += S(i, j) * v[j];
        }
    }
    return out;
}

// RowVector (1xK) * SymMatrix (KxK) = RowVector (1xK).
template <size_t K>
[[nodiscard]] constexpr Vector<K> operator*(const Vector<K> &v, const SymMatrix<K> &S) {
    return S * v;  // because matrix is symmetric, the order doesn't matter
}

// Outer product: ColVector (Nx1) * RowVector (1xM) = Matrix (NxM).
template <size_t N, size_t M>
[[nodiscard]] constexpr Matrix<N, M> Outer(const Vector<N> &a, const Vector<M> &b) {
    Matrix<N, M> out{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < M; ++j) {
            out[i][j] = a[i] * b[j];
        }
    }
    return out;
}

// SymMatrix (KxK) * Matrix (KxM) = Matrix (KxM).
template <size_t K, size_t M>
[[nodiscard]] constexpr Matrix<K, M> operator*(const SymMatrix<K> &S, const Matrix<K, M> &A) {
    Matrix<K, M> out{};
    for (size_t i{0}; i < K; ++i) {
        for (size_t k{0}; k < K; ++k) {
            const double s_ik = S(i, k);
            for (size_t j{0}; j < M; ++j) {
                out[i][j] += s_ik * A[k][j];
            }
        }
    }
    return out;
}

// Matrix (NxK) * SymMatrix (KxK) = Matrix (NxK).
template <size_t N, size_t K>
[[nodiscard]] constexpr Matrix<N, K> operator*(const Matrix<N, K> &A, const SymMatrix<K> &S) {
    Matrix<N, K> out{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t k{0}; k < K; ++k) {
            const double a_ik = A[i][k];
            for (size_t j{0}; j < K; ++j) {
                out[i][j] += a_ik * S(k, j);
            }
        }
    }
    return out;
}

// Quadratic form.
// RowVector (1xK) * SymMatrix (KxK) * ColVector (Kx1) = scalar
template <size_t K>
[[nodiscard]] constexpr double VTSV(const Vector<K> &v, const SymMatrix<K> &S) {
    double sum{0.};
    // diagonal terms
    for (size_t i{0}; i < K; ++i) sum += v[i] * S(i, i) * v[i];
    // off-diagonal terms (counted twice due to symmetry)
    for (size_t i{1}; i < K; ++i) {
        for (size_t j{0}; j < i; ++j) {
            sum += 2.0 * v[i] * S(i, j) * v[j];
        }
    }
    return sum;
}

// Similarity transform.
// A (NxK) * S (KxK, symm.) * A^T (KxN) = (NxN, symm.)
template <size_t N, size_t K>
[[nodiscard]] constexpr SymMatrix<N> ASAT(const Matrix<N, K> &A, const SymMatrix<K> &S) {
    SymMatrix<N> out{};
    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j <= i; ++j) {
            double sum{0.};
            for (size_t k{0}; k < K; ++k) {
                for (size_t l{0}; l < K; ++l) {
                    sum += A[i][k] * S(k, l) * A[j][l];
                }
            }
            out(i, j) = sum;
        }
    }
    return out;
}

// Sum a square matrix with its transpose. The result is a symmetric matrix.
// A (KxK) + A^T (KxK) = S (KxK, symm.)
template <size_t K>
[[nodiscard]] constexpr SymMatrix<K> ItselfPlusItsTranspose(const Matrix<K, K> &A) {
    SymMatrix<K> out{};
    for (size_t i{0}; i < K; ++i) {
        for (size_t j{0}; j <= i; ++j) {
            out(i, j) = A[i][j] + A[j][i];
        }
    }
    return out;
}

// Do the outer product of a vector with itself. The result is a symmetric matrix.
// v^T (Kx1) * v (1xK) = S (KxK, symm.)
template <size_t K>
[[nodiscard]] constexpr SymMatrix<K> OuterProductWithItself(const Vector<K> &v) {
    SymMatrix<K> out{};
    for (size_t i{0}; i < K; ++i) {
        for (size_t j{0}; j <= i; ++j) {
            out(i, j) = v[i] * v[j];
        }
    }
    return out;
}

// ## Complex Operations ## //

// Invert 3x3 symmetric matrix via a modified Cholesky decomposition.
// - `src` : [input] a symmetric 3x3 matrix.
// - `inv` : [output by reference] invert of `S`.
inline void Invert_3x3SymMatrix_ModCholesky(const SymMatrix<3> &src, SymMatrix<3> &inv) {

    constexpr double kRegularization = 1E-8;

    Vector<3> d{};
    Matrix<3, 3> u{};

    // Decomposition: S = U^T * D * U (modified Cholesky)
    for (size_t i{0}; i < 3; ++i) {
        double diag{0.};
        for (size_t j{0}; j < i; ++j) diag += u[j][i] * u[j][i] * d[j];
        diag = src[i * (i + 3) / 2] - diag;

        if (std::abs(diag) < kRegularization) diag = kRegularization;

        d[i] = diag / std::abs(diag);  // sign: +1 or -1
        u[i][i] = std::sqrt(std::abs(diag));

        for (size_t j{i + 1}; j < 3; ++j) {
            double off{0.};
            for (size_t k{0}; k < i; ++k) off += u[k][i] * u[k][j] * d[k];
            off = src[j * (j + 1) / 2 + i] - off;
            u[i][j] = d[i] / u[i][i] * off;
        }
    }

    // Inversion of U (in-place)
    Vector<3> u_diag{u[0][0], u[1][1], u[2][2]};

    u[0][0] = 1.0 / u[0][0];
    u[1][1] = 1.0 / u[1][1];
    u[2][2] = 1.0 / u[2][2];

    u[0][1] = -u[0][1] * u[0][0] * u[1][1];
    u[1][2] = -u[1][2] * u[1][1] * u[2][2];

    u[0][2] = u[0][1] * u_diag[1] * u[1][2] - u[0][2] * u[0][0] * u[2][2];

    // Assemble S^{-1} = U^{-1} * D^{-1} * (U^{-1})^T
    inv[3] = u[0][2] * u[2][2] * d[2];
    inv[4] = u[1][2] * u[2][2] * d[2];
    inv[5] = u[2][2] * u[2][2] * d[2];

    inv[1] = u[0][1] * u[1][1] * d[1] + u[0][2] * u[1][2] * d[2];
    inv[2] = u[1][1] * u[1][1] * d[1] + u[1][2] * u[1][2] * d[2];

    inv[0] = u[0][0] * u[0][0] * d[0] + u[0][1] * u[0][1] * d[1] + u[0][2] * u[0][2] * d[2];
}

}  // namespace CompactLinearAlgebra

// ## Print Utilities ## //

template <size_t N>
struct std::formatter<CompactLinearAlgebra::Vector<N>> {

    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

    auto format(const CompactLinearAlgebra::Vector<N> &vec, std::format_context &ctx) const {
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

template <size_t N, size_t M>
struct std::formatter<CompactLinearAlgebra::Matrix<N, M>> {

    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

    auto format(const CompactLinearAlgebra::Matrix<N, M> &arr, std::format_context &ctx) const {
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

template <size_t K>
struct std::formatter<CompactLinearAlgebra::SymMatrix<K>> {

    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

    auto format(const CompactLinearAlgebra::SymMatrix<K> &arr, std::format_context &ctx) const {
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
