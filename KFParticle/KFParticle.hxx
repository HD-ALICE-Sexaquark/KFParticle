/*
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *               2007-2019 Sergey Gorbunov
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <print>
#include <utility>
#include <vector>

#include "KFParticle_Math.hxx"
#include "KFParticle_Utils.hxx"

namespace KF {

struct alignas(KF_SIMD_ALIGN) PCA {
    PCA() = default;
    PCA(double x, double y, double z, double px, double py, double pz) : xyz{x, y, z}, dir{px, py, pz} {};
    Vector<3> xyz{};
    Vector<3> dir{};
};
struct alignas(KF_SIMD_ALIGN) Cache {
    PCA pca;
    double theta{0.};
    double sin{0.};
    double cos{0.};
    double sB{0.};
    double cB{0.};
    double ds{0.};
};

namespace Result {
struct alignas(KF_SIMD_ALIGN) Minimization : Cache {
    Vector<6> ds_dr{};
    Vector<6> ds_dr1{};
};
struct alignas(KF_SIMD_ALIGN) Transport {
    Matrix<6, 6> jacob{};
    Matrix<6, 6> corr{};
    SymMatrix<8> C{};
    Vector<8> P{};
};
struct alignas(KF_SIMD_ALIGN) Measurement {
    Measurement(const Vector<8> &p1, const Vector<8> &p2, const SymMatrix<8> &c1, const SymMatrix<8> &c2)  //
        : C1{c1}, C2{c2}, P1{p1}, P2{p2} {}

    SymMatrix<8> C1{};
    SymMatrix<8> C2{};
    Matrix<3, 3> D{};
    Vector<8> P1{};
    Vector<8> P2{};
};
struct alignas(KF_SIMD_ALIGN) MassConstraint {
    SymMatrix<8> C{};
    Matrix<7, 7> jacob{};
    Vector<8> P{};
};
}  // namespace Result

// @class KFParticleBase
// @brief The base of KFParticle class, describes particle objects.
// @author S.Gorbunov, I.Kisel, M.Zyzak
// @date 05.02.2019
// @version 1.0
//
// Contains the main mathematics of the KFParticle.
class alignas(KF_SIMD_ALIGN) Particle {
   public:
    Particle(const Particle &) = default;
    Particle(Particle &&) = default;
    Particle &operator=(const Particle &) = default;
    Particle &operator=(Particle &&) = default;

    Particle() { Initialize(); }
    Particle(const Vector<6> &p, const SymMatrix<6> &cov, int charge, double mass) { Initialize(p, cov, charge, mass); }
    Particle(const Vector<7> &p, const SymMatrix<7> &cov, int charge) { Initialize(p, cov, charge); }
    ~Particle() = default;

    [[nodiscard]] double X() const noexcept { return fP[0]; }     // return X coordinate of the particle
    [[nodiscard]] double Y() const noexcept { return fP[1]; }     // return Y coordinate of the particle
    [[nodiscard]] double Z() const noexcept { return fP[2]; }     // return Z coordinate of the particle
    [[nodiscard]] double Px() const noexcept { return fP[3]; }    // return X component of the momentum
    [[nodiscard]] double Py() const noexcept { return fP[4]; }    // return Y component of the momentum
    [[nodiscard]] double Pz() const noexcept { return fP[5]; }    // return Z component of the momentum
    [[nodiscard]] double E() const noexcept { return fP[6]; }     // return energy of the particle
    [[nodiscard]] double S() const noexcept { return fP[7]; }     // return dS=l/p, l - decay length, defined if production vertex is set
    [[nodiscard]] int Charge() const noexcept { return fQ; }      // return charge of the particle
    [[nodiscard]] double Chi2() const noexcept { return fChi2; }  // return Chi2 of the fit
    [[nodiscard]] int NDF() const noexcept { return fNDF; }       // return number of degrees of freedom
    [[nodiscard]] double Chi2NDF() const { return fChi2 / static_cast<double>(fNDF); }  // return Chi2/ndf

    [[nodiscard]] double P2() const { return Math::SquaredNorm<3>({Px(), Py(), Pz()}); }
    [[nodiscard]] double P() const { return std::sqrt(P2()); }
    [[nodiscard]] double Pt() const { return Math::Norm<2>({Px(), Py()}); }
    [[nodiscard]] double Mass() const {
        double mass2{E() * E() - P2()};
        if (mass2 < 0.) return -1.;  // protection
        return std::sqrt(mass2);
    }

    // Pseudorapidity.
    [[nodiscard]] double Eta() const { return std::atanh(Pz() / P() + Const::Epsilon); }

    // Rapidity.
    [[nodiscard]] double Rapidity() const { return std::log((E() + Pz()) / (E() - Pz() + Const::Epsilon)) / 2.; }

    // Radius (cm) in cylindrical coordinates.
    [[nodiscard]] double Radius2D() const { return Math::Norm<2>({X(), Y()}); }

    // Radius (cm) in spherical coordinates.
    [[nodiscard]] double Radius3D() const { return Math::Norm<3>({X(), Y(), Z()}); }

    // Return point of closest approach (PCA) of a certain daughter after minimization.
    [[nodiscard]] PCA GetPCA(size_t index_daughter) const {
        if (fPCAs.size() <= index_daughter) return {0., 0., 0., 0., 0., 0.};  // protection
        return fPCAs[index_daughter];
    }

    // Return distance of closest approach (DCA) (cm) between added daughter and fitted vertex.
    [[nodiscard]] std::optional<double> GetDCA(size_t index_daughter) const {
        if (fPCAs.size() <= index_daughter) return std::nullopt;  // protection
        Vector<3> diff{};
        for (size_t i{0}; i < 3; ++i) {
            diff[i] = fP[i] - fPCAs[index_daughter].xyz[i];
        }
        return Math::Norm<3>(diff);
    }

    // Return distance of closest approach (DCA) (cm) between added daughter1 and added daughter2.
    [[nodiscard]] std::optional<double> GetDCA(size_t index_daughter1, size_t index_daughter2) const {
        if (fPCAs.size() <= index_daughter1 || fPCAs.size() <= index_daughter2) return std::nullopt;  // protection
        Vector<3> diff{};
        for (size_t i{0}; i < 3; ++i) {
            diff[i] = fPCAs[index_daughter1].xyz[i] - fPCAs[index_daughter2].xyz[i];
        }
        return Math::Norm<3>(diff);
    }

    // Return distance of closest approach (DCA) (cm) in XY plane between added daughter and fitted vertex.
    [[nodiscard]] std::optional<double> GetDCAxy(size_t index_daughter) const {
        if (fPCAs.size() <= index_daughter) return std::nullopt;  // protection
        Vector<2> diff{};
        for (size_t i{0}; i < 2; ++i) {
            diff[i] = fP[i] - fPCAs[index_daughter].xyz[i];
        }
        return Math::Norm<2>(diff);
    }

    // Return distance of closest approach (DCA) (cm) in XY plane between added daughter1 and added daughter2.
    [[nodiscard]] std::optional<double> GetDCAxy(size_t index_daughter1, size_t index_daughter2) const {
        if (fPCAs.size() <= index_daughter1 || fPCAs.size() <= index_daughter2) return std::nullopt;  // protection
        Vector<2> diff{};
        for (size_t i{0}; i < 2; ++i) {
            diff[i] = fPCAs[index_daughter1].xyz[i] - fPCAs[index_daughter2].xyz[i];
        }
        return Math::Norm<2>(diff);
    }

    [[nodiscard]] double GetParameter(size_t i) const { return fP[i]; }   // return P[i] parameter
    [[nodiscard]] double GetCovariance(size_t i) const { return fC[i]; }  // return C[i] element of the covariance matrix in the lower triangular form
    [[nodiscard]] double GetCovariance(size_t i, size_t j) const { return fC[IJ(i, j)]; }  // return C[i,j] element of the covariance matrix

    void AddDaughterWithEnergyFit(const Particle &daughter, double bz, double chi2_threshold = 1E4);
    void AddDaughter(const Particle &daughter, double bz) {
        if (fNDF < -1) {  // first daughter -> just copy
            fNDF += 2;
            fQ = daughter.Charge();
            for (size_t i{0}; i < 7; ++i) fP[i] = daughter.fP[i];
            for (size_t i{0}; i < 28; ++i) fC[i] = daughter.fC[i];
            return;
        }
        AddDaughterWithEnergyFit(daughter, bz);
    }
    void AddProductionVertex(const Vector<3> &prod_vtx, const SymMatrix<3> &cov, double bz, double chi2_threshold = 1E4);
    void AddMassConstraint(double target_mass);

    void Print() const {
        std::println(stdout, "(X,Y,Z,S)    = ({:13.6e}, {:13.6e}, {:13.6e}, {:13.6e})", fP[0], fP[1], fP[2], fP[7]);
        std::println(stdout, "(Px,Py,Pz,E) = ({:13.6e}, {:13.6e}, {:13.6e}, {:13.6e})", fP[3], fP[4], fP[5], fP[6]);
        std::println(stdout, "Mass         = {:13.6e}", Mass());
        std::println(stdout, "Radius2D     = {:13.6e}", Radius2D());
        std::println(stdout, "Chi2/NDF     = {:13.6e} / {} = {:13.6e}", fChi2, fNDF, Chi2NDF());
        std::println(stdout, "CovMatrix    = {}", fC);
        std::println(stdout, "DCAxy_Dau    = {:13.6e}", GetDCAxy(0, 1));
        std::println(stdout, "DCAxy_Neg    = {:13.6e}", GetDCAxy(0));
        std::println(stdout, "DCAxy_Pos    = {:13.6e}", GetDCAxy(1));
    }

   protected:
    // Set Cxx=Cyy=Czz=100 and Css=1
    // Note: it will modify the state of the current `KF::Particle`
    void Initialize() {
        fC[0] = Const::Initial_C_xx;
        fC[2] = Const::Initial_C_yy;
        fC[5] = Const::Initial_C_zz;
        fC[35] = Const::Initial_C_SS;
    }

    // Set the parameters of the particle:
    // Input arguments:
    // - `param`  : position and momentum { X, Y, Z, Px, Py, Pz }
    // - `cov`    : lower-triangular part of the symmetric 6x6 covariance matrix
    // - `charge` : charge of the particle in elementary charge units
    // - `mass`   : the mass hypothesis
    // Note: it will modify the state of the current `KF::Particle`
    void Initialize(const Vector<6> &param, const SymMatrix<6> &cov, int charge, double mass) {
#if KF_DEBUG
        std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

        for (size_t i{0}; i < 6; ++i) fP[i] = param[i];
        double energy{std::sqrt(mass * mass + fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5])};
        fP[6] = energy;
        fP[7] = 0.;

        double h0{fP[3] / energy};
        double h1{fP[4] / energy};
        double h2{fP[5] / energy};

        for (size_t i{0}; i < 21; ++i) fC[i] = cov[i];
        fC[21] = h0 * fC[6] + h1 * fC[10] + h2 * fC[15];
        fC[22] = h0 * fC[7] + h1 * fC[11] + h2 * fC[16];
        fC[23] = h0 * fC[8] + h1 * fC[12] + h2 * fC[17];
        fC[24] = h0 * fC[9] + h1 * fC[13] + h2 * fC[18];
        fC[25] = h0 * fC[13] + h1 * fC[14] + h2 * fC[19];
        fC[26] = h0 * fC[18] + h1 * fC[19] + h2 * fC[20];
        fC[27] = (h0 * h0 * fC[9] + h1 * h1 * fC[14] + h2 * h2 * fC[20] + 2 * (h0 * h1 * fC[13] + h0 * h2 * fC[18] + h1 * h2 * fC[19]));
        fC[35] = 1.;

        fQ = charge;
#if KF_DEBUG
        Utils::Print(__FUNCTION__, "CovMatrix", fC);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
    }

    // Set the parameters of the particle:
    // Input arguments:
    // - `param`  : position, momentum and energy { X, Y, Z, Px, Py, Pz, E }
    // - `cov`    : lower-triangular part of the symmetric 7x7 covariance matrix
    // - `charge` : charge of the particle in elementary charge units
    // Note: it will modify the state of the current `KF::Particle`
    void Initialize(const Vector<7> &param, const SymMatrix<7> &cov, int charge) {
#if KF_DEBUG
        std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif
        for (size_t i{0}; i < 7; ++i) fP[i] = param[i];
        fP[7] = 0.;
        for (size_t i{0}; i < 28; ++i) fC[i] = cov[i];
        fC[35] = 1.;
        fQ = charge;
#if KF_DEBUG
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
    }

    [[nodiscard]] Result::Measurement GetMeasurement(const Particle &daughter, double bz) const;

    [[nodiscard]] Result::Minimization MinimizeLinePoint(const Vector<3> &v) const;
    [[nodiscard]] Result::Minimization MinimizeHelixPoint(const Vector<3> &v, double bz) const;
    [[nodiscard]] Result::Minimization Minimize(const Vector<3> &v, double bz = 0) const {
        if (std::abs(fQ) < Const::AbsAlmostZero || std::abs(bz) < Const::AbsAlmostZero) {
            return MinimizeLinePoint(v);
        }
        return MinimizeHelixPoint(v, bz);
    }

    [[nodiscard]] std::pair<Result::Minimization, Result::Minimization> MinimizeLineLine(const Particle &p) const;
    [[nodiscard]] std::pair<Result::Minimization, Result::Minimization> MinimizeHelixHelix(const Particle &p, double bz) const;
    [[nodiscard]] std::pair<Result::Minimization, Result::Minimization> Minimize(const Particle &p, double bz = 0.) const {
        if ((std::abs(fQ) < Const::AbsAlmostZero && std::abs(p.fQ) < Const::AbsAlmostZero) || std::abs(bz) < Const::AbsAlmostZero) {
            return MinimizeLineLine(p);
        }
        return MinimizeHelixHelix(p, bz);
    }

    [[nodiscard]] Result::Transport TransportBz(const Result::Minimization &min, double bz) const;
    [[nodiscard]] Result::Transport TransportLine(const Result::Minimization &min) const;
    [[nodiscard]] Result::Transport Transport(const Result::Minimization &min, double bz = 0.) const {
        if (std::abs(fQ) < Const::AbsAlmostZero || std::abs(bz) < Const::AbsAlmostZero) {
            return TransportLine(min);
        }
        return TransportBz(min, bz);
    }

    SymMatrix<8> fC{};  // lower-triangular part of the symmetric 8x8 covariance matrix

    // Registered Points of Closest Approach (PCAs).
    // 0) If there's no or a single daughter has been added <-> no fitted vertex => size = 0
    // 1) After a second daughter has been added <-> there is a fitted vertex => size = 2
    // 2) After that, for any additional daughter or production vertex => size += 1
    std::vector<PCA> fPCAs;

    Vector<8> fP{};    // particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}
    double fChi2{0.};  // chi2
    int fNDF{-3};      // number of degrees of freedom
    int fQ{0};         // charge of the particle in units of elementary charge
};

}  // namespace KF
