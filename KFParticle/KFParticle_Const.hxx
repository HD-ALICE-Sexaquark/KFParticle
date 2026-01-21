#pragma once

#if defined(__AVX512F__)
#define KF_SIMD_ALIGN 64
#elif defined(__AVX2__)
#define KF_SIMD_ALIGN 32
#else
#define KF_SIMD_ALIGN alignof(double)
#endif

namespace KF::Const {
constexpr double Kappa{0.000299792458};  // (GeV/c) / (kG/cm)
constexpr double AbsAlmostZero{1.E-8};
constexpr double Epsilon{1.E-6};
constexpr double Initial_C_xx{100.};
constexpr double Initial_C_yy{100.};
constexpr double Initial_C_zz{100.};
constexpr double Initial_C_SS{1.};
}  // namespace KF::Const
