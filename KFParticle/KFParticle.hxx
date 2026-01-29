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

#include <armadillo>

#include "KFParticle_Const.hxx"
#include "KFParticle_Utils.hxx"

namespace KF {

struct PCA {
    PCA() = default;
    PCA(double x, double y, double z, double px, double py, double pz) : xyz{x, y, z}, dir{px, py, pz} {};
    arma::vec::fixed<3> xyz;
    arma::vec::fixed<3> dir;
};
struct Cache {
    PCA pca;
    double theta{0.};
    double sin{0.};
    double cos{0.};
    double sB{0.};
    double cB{0.};
    double ds{0.};
};

namespace Result {
struct Minimization : Cache {
    arma::vec::fixed<6> ds_dr;
    arma::vec::fixed<6> ds_dr1;
};
struct Transport {
    arma::mat::fixed<8, 8> C;
    arma::mat::fixed<6, 6> jacob;
    arma::mat::fixed<6, 6> corr;
    arma::vec::fixed<8> P;
};
struct Measurement {
    Measurement(const arma::vec::fixed<8> &p1, const arma::vec::fixed<8> &p2, const arma::mat::fixed<8, 8> &c1, const arma::mat::fixed<8, 8> &c2)
        : C1{c1}, C2{c2}, P1{p1}, P2{p2} {}

    arma::mat::fixed<8, 8> C1;
    arma::mat::fixed<8, 8> C2;
    arma::mat::fixed<3, 3> D;
    arma::vec::fixed<8> P1;
    arma::vec::fixed<8> P2;
};
struct MassConstraint {
    arma::mat::fixed<8, 8> jacob;
    arma::mat::fixed<7, 7> C;
    arma::vec::fixed<8> P;
};
}  // namespace Result

// @class KFParticleBase
// @brief The base of KFParticle class, describes particle objects.
// @author S.Gorbunov, I.Kisel, M.Zyzak
// @date 05.02.2019
// @version 1.0
//
// Contains the main mathematics of the KFParticle.
class Particle {
   public:
    Particle(const Particle &) = default;
    Particle(Particle &&) noexcept = default;
    Particle &operator=(const Particle &) = default;
    Particle &operator=(Particle &&) noexcept = default;

    Particle() { Initialize(); }
    Particle(const arma::vec::fixed<6> &p, const arma::mat::fixed<6, 6> &cov, int charge, double mass) { Initialize(p, cov, charge, mass); }
    Particle(const arma::vec::fixed<7> &p, const arma::mat::fixed<7, 7> &cov, int charge) { Initialize(p, cov, charge); }
    ~Particle() = default;

    [[nodiscard]] double X() const noexcept { return fP(0); }     // return X coordinate of the particle
    [[nodiscard]] double Y() const noexcept { return fP(1); }     // return Y coordinate of the particle
    [[nodiscard]] double Z() const noexcept { return fP(2); }     // return Z coordinate of the particle
    [[nodiscard]] double Px() const noexcept { return fP(3); }    // return X component of the momentum
    [[nodiscard]] double Py() const noexcept { return fP(4); }    // return Y component of the momentum
    [[nodiscard]] double Pz() const noexcept { return fP(5); }    // return Z component of the momentum
    [[nodiscard]] double E() const noexcept { return fP(6); }     // return energy of the particle
    [[nodiscard]] double S() const noexcept { return fP(7); }     // return dS=l/p, l - decay length, defined if production vertex is set
    [[nodiscard]] double Chi2() const noexcept { return fChi2; }  // return Chi2 of the fit
    [[nodiscard]] int NDF() const noexcept { return fNDF; }       // return number of degrees of freedom
    [[nodiscard]] int Charge() const noexcept { return fQ; }      // return charge of the particle

    [[nodiscard]] double Chi2NDF() const { return fChi2 / static_cast<double>(fNDF); }  // return Chi2/ndf
    [[nodiscard]] double P() const { return std::hypot(Px(), Py(), Pz()); }
    [[nodiscard]] double Pt() const { return std::hypot(Px(), Py()); }
    [[nodiscard]] std::optional<double> Mass() const {
        double m_squared{(E() - P()) * (E() + P())};
        if (m_squared < 0.) return std::nullopt;  // protection
        return std::sqrt(m_squared);
    }

    // Pseudorapidity.
    [[nodiscard]] double Eta() const { return std::atanh(Pz() / P() + Const::Epsilon); }

    // Rapidity.
    [[nodiscard]] double Rapidity() const { return std::log((E() + Pz()) / (E() - Pz() + Const::Epsilon)) / 2.; }

    // Radius (cm) in cylindrical coordinates.
    [[nodiscard]] double Radius2D() const { return std::hypot(X(), Y()); }

    // Radius (cm) in spherical coordinates.
    [[nodiscard]] double Radius3D() const { return std::hypot(X(), Y(), Z()); }

    // Return point of closest approach (PCA) of a certain daughter after minimization.
    [[nodiscard]] PCA GetPCA(size_t index_daughter) const {
        if (fPCAs.size() <= index_daughter) return {0., 0., 0., 0., 0., 0.};  // protection
        return fPCAs[index_daughter];
    }

    // Return distance of closest approach (DCA) (cm) between added daughter and fitted vertex.
    [[nodiscard]] std::optional<double> GetDCA(size_t index_daughter) const {
        if (fPCAs.size() <= index_daughter) return std::nullopt;  // protection
        arma::vec::fixed<3> diff = fP.head(3) - fPCAs[index_daughter].xyz;
        return arma::norm(diff);
    }

    // Return distance of closest approach (DCA) (cm) between added daughter1 and added daughter2.
    [[nodiscard]] std::optional<double> GetDCA(size_t index_daughter1, size_t index_daughter2) const {
        if (fPCAs.size() <= index_daughter1 || fPCAs.size() <= index_daughter2) return std::nullopt;  // protection
        arma::vec::fixed<3> diff = fPCAs[index_daughter1].xyz - fPCAs[index_daughter2].xyz;
        return norm(diff);
    }

    // Return distance of closest approach (DCA) (cm) in XY plane between added daughter and fitted vertex.
    [[nodiscard]] std::optional<double> GetDCAxy(size_t index_daughter) const {
        if (fPCAs.size() <= index_daughter) return std::nullopt;  // protection
        arma::vec::fixed<2> diff = fP.head(2) - fPCAs[index_daughter].xyz.head(2);
        return arma::norm(diff);
    }

    // Return distance of closest approach (DCA) (cm) in XY plane between added daughter1 and added daughter2.
    [[nodiscard]] std::optional<double> GetDCAxy(size_t index_daughter1, size_t index_daughter2) const {
        if (fPCAs.size() <= index_daughter1 || fPCAs.size() <= index_daughter2) return std::nullopt;  // protection
        arma::vec::fixed<2> diff = fPCAs[index_daughter1].xyz.head(2) - fPCAs[index_daughter2].xyz.head(2);
        return norm(diff);
    }

    [[nodiscard]] double GetParameter(size_t i) const { return fP(i); }
    [[nodiscard]] double GetCovariance(size_t i, size_t j) const { return fC(i, j); }

    [[nodiscard]] double SigmaX2() const { return fC(0, 0); }
    [[nodiscard]] double SigmaXY() const { return fC(1, 0); }
    [[nodiscard]] double SigmaY2() const { return fC(1, 1); }
    [[nodiscard]] double SigmaXZ() const { return fC(2, 0); }
    [[nodiscard]] double SigmaYZ() const { return fC(2, 1); }
    [[nodiscard]] double SigmaZ2() const { return fC(2, 2); }
    [[nodiscard]] double SigmaXPx() const { return fC(3, 0); }
    [[nodiscard]] double SigmaYPx() const { return fC(3, 1); }
    [[nodiscard]] double SigmaZPx() const { return fC(3, 2); }
    [[nodiscard]] double SigmaPx2() const { return fC(3, 3); }
    [[nodiscard]] double SigmaXPy() const { return fC(4, 0); }
    [[nodiscard]] double SigmaYPy() const { return fC(4, 1); }
    [[nodiscard]] double SigmaZPy() const { return fC(4, 2); }
    [[nodiscard]] double SigmaPxPy() const { return fC(4, 3); }
    [[nodiscard]] double SigmaPy2() const { return fC(4, 4); }
    [[nodiscard]] double SigmaXPz() const { return fC(5, 0); }
    [[nodiscard]] double SigmaYPz() const { return fC(5, 1); }
    [[nodiscard]] double SigmaZPz() const { return fC(5, 2); }
    [[nodiscard]] double SigmaPxPz() const { return fC(5, 3); }
    [[nodiscard]] double SigmaPyPz() const { return fC(5, 4); }
    [[nodiscard]] double SigmaPz2() const { return fC(5, 5); }
    [[nodiscard]] double SigmaXE() const { return fC(6, 0); }
    [[nodiscard]] double SigmaYE() const { return fC(6, 1); }
    [[nodiscard]] double SigmaZE() const { return fC(6, 2); }
    [[nodiscard]] double SigmaPxE() const { return fC(6, 3); }
    [[nodiscard]] double SigmaPyE() const { return fC(6, 4); }
    [[nodiscard]] double SigmaPzE() const { return fC(6, 5); }
    [[nodiscard]] double SigmaE2() const { return fC(6, 6); }
    [[nodiscard]] double SigmaXS() const { return fC(7, 0); }
    [[nodiscard]] double SigmaYS() const { return fC(7, 1); }
    [[nodiscard]] double SigmaZS() const { return fC(7, 2); }
    [[nodiscard]] double SigmaPxS() const { return fC(7, 3); }
    [[nodiscard]] double SigmaPyS() const { return fC(7, 4); }
    [[nodiscard]] double SigmaPzS() const { return fC(7, 5); }
    [[nodiscard]] double SigmaES() const { return fC(7, 6); }
    [[nodiscard]] double SigmaS2() const { return fC(7, 7); }

    bool AddDaughterWithEnergyFit(const Particle &daughter, double bz, double chi2_threshold = 1E4);
    bool AddDaughter(const Particle &daughter, double bz) {
        if (fNDF < -1) {  // first daughter -> just copy
            fP = daughter.fP;
            fC = daughter.fC;
            fQ = daughter.Charge();
            fNDF += 2;
            return true;
        }
        return AddDaughterWithEnergyFit(daughter, bz);
    }
    bool AddProductionVertex(const arma::vec::fixed<3> &prod_vtx, const arma::mat::fixed<3, 3> &cov, double bz, double chi2_threshold = 1E4);
    bool AddMassConstraint(double target_mass);

    void Print() const {
        std::println(stdout, "(X,Y,Z,S)    = ({:13.6e}, {:13.6e}, {:13.6e}, {:13.6e})", fP[0], fP[1], fP[2], fP[7]);
        std::println(stdout, "(Px,Py,Pz,E) = ({:13.6e}, {:13.6e}, {:13.6e}, {:13.6e})", fP[3], fP[4], fP[5], fP[6]);
        std::println(stdout, "Mass         = {:13.6e}", (*Mass()));
        std::println(stdout, "Radius2D     = {:13.6e}", Radius2D());
        std::println(stdout, "Chi2/NDF     = {:13.6e} / {} = {:13.6e}", fChi2, fNDF, Chi2NDF());
        fC.print("CovMatrix    =");
        std::println(stdout, "DCAxy_Dau    = {:13.6e}", (*GetDCAxy(0, 1)));
        std::println(stdout, "DCAxy_Neg    = {:13.6e}", (*GetDCAxy(0)));
        std::println(stdout, "DCAxy_Pos    = {:13.6e}", (*GetDCAxy(1)));
    }

   protected:
    // Set Cxx=Cyy=Czz=100 and Css=1
    // Note: it will modify the state of the current `KF::Particle`
    void Initialize() {
        fC(0, 0) = Const::Initial_C_xx;
        fC(1, 1) = Const::Initial_C_yy;
        fC(2, 2) = Const::Initial_C_zz;
        fC(7, 7) = Const::Initial_C_SS;
    }

    // Set the parameters of the particle:
    // Input arguments:
    // - `param`  : position and momentum { X, Y, Z, Px, Py, Pz }
    // - `cov`    : covariance matrix -- 6x6 symm.
    // - `charge` : charge of the particle in elementary charge units
    // - `mass`   : the mass hypothesis
    // Note: it will modify the state of the current `KF::Particle`
    void Initialize(const arma::vec::fixed<6> &param, const arma::mat::fixed<6, 6> &cov, int charge, double mass) {
#if KF_DEBUG
        std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

        fP.head(6) = param;
        double momentum{std::hypot(fP(3), fP(4), fP(5))};
        double energy{std::hypot(mass, momentum)};
        fP(6) = energy;
        fP(7) = 0.;

        double h0{fP(3) / energy};
        double h1{fP(4) / energy};
        double h2{fP(5) / energy};

        fC.submat(0, 0, arma::size(6, 6)) = cov;
        fC(6, 0) = h0 * fC(3, 0) + h1 * fC(4, 0) + h2 * fC(5, 0);
        fC(6, 1) = h0 * fC(3, 1) + h1 * fC(4, 1) + h2 * fC(5, 1);
        fC(6, 2) = h0 * fC(3, 2) + h1 * fC(4, 2) + h2 * fC(5, 2);
        fC(6, 3) = h0 * fC(3, 3) + h1 * fC(4, 3) + h2 * fC(5, 3);
        fC(6, 4) = h0 * fC(4, 3) + h1 * fC(4, 4) + h2 * fC(5, 4);
        fC(6, 5) = h0 * fC(5, 3) + h1 * fC(5, 4) + h2 * fC(5, 5);
        fC(6, 6) = (h0 * h0 * fC(3, 3) + h1 * h1 * fC(4, 4) + h2 * h2 * fC(5, 5) +  //
                    2 * (h0 * h1 * fC(4, 3) + h0 * h2 * fC(5, 3) + h1 * h2 * fC(5, 4)));
        fC(7, 7) = 1.;
        fC = arma::symmatl(fC);  // force symmetry

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
    void Initialize(const arma::vec::fixed<7> &param, const arma::mat::fixed<7, 7> &cov, int charge) {
#if KF_DEBUG
        std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif
        fP.head(7) = param;
        fP(7) = 0.;

        fC.submat(0, 0, arma::size(7, 7)) = cov;
        fC(7, 7) = 1.;

        fQ = charge;
#if KF_DEBUG
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
    }

    [[nodiscard]] Result::Measurement GetMeasurement(const Particle &daughter, double bz) const;

    [[nodiscard]] Result::Minimization MinimizeLinePoint(const arma::vec::fixed<3> &v) const;
    [[nodiscard]] Result::Minimization MinimizeHelixPoint(const arma::vec::fixed<3> &v, double bz) const;
    [[nodiscard]] Result::Minimization Minimize(const arma::vec::fixed<3> &v, double bz = 0) const {
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

    arma::mat::fixed<8, 8> fC;  // symmetric 8x8 covariance matrix

    // Registered Points of Closest Approach (PCAs).
    // 0) If there's no or a single daughter has been added <-> no fitted vertex => size = 0
    // 1) After a second daughter has been added <-> there is a fitted vertex => size = 2
    // 2) After that, for any additional daughter or production vertex => size += 1
    std::vector<PCA> fPCAs;

    arma::vec::fixed<8> fP;  // particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}
    double fChi2{0.};        // chi2
    int fNDF{-3};            // number of degrees of freedom
    int fQ{0};               // charge of the particle in units of elementary charge
};

}  // namespace KF
