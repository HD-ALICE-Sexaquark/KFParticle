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

#ifndef KFPARTICLE_HXX
#define KFPARTICLE_HXX

#include <cmath>
#include <iostream>
#include <utility>

#include "KFParticle_Math.hxx"

namespace KF {

struct alignas(32) PCA {
    PCA() = default;
    PCA(double x, double y, double z, double px, double py, double pz) : pos{x, y, z}, dir{px, py, pz} {};
    Vector<3> pos{};
    Vector<3> dir{};
};
struct alignas(32) Cache {
    PCA pca;
    double theta{0.};
    double sin{0.};
    double cos{0.};
    double sB{0.};
    double cB{0.};
    double ds{0.};
};

namespace Result {
struct alignas(32) Minimization : Cache {
    Vector<6> ds_dr{};
    Vector<6> ds_dr1{};
};
struct alignas(32) Transport {
    Matrix<6, 6> jacob{};
    Matrix<6, 6> corr{};
    SymMatrix<8> C{};
    Vector<8> P{};
};
struct alignas(32) Measurement {
    SymMatrix<8> C1{};
    SymMatrix<8> C2{};
    Matrix<3, 3> D{};
    Vector<8> P1{};
    Vector<8> P2{};
};
struct alignas(32) MassConstraint {
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
class alignas(32) Particle {
   public:
    Particle(const Particle &) = default;
    Particle(Particle &&) = delete;
    Particle &operator=(const Particle &) = default;
    Particle &operator=(Particle &&) = delete;

    Particle() { Initialize(); }
    Particle(const Vector<6> &p, const SymMatrix<6> &cov, int charge, double mass) { Initialize(p, cov, charge, mass); }
    ~Particle() = default;

    double X() const { return fP[0]; }        // return X coordinate of the particle
    double Y() const { return fP[1]; }        // return Y coordinate of the particle
    double Z() const { return fP[2]; }        // return Z coordinate of the particle
    double Px() const { return fP[3]; }       // return X component of the momentum
    double Py() const { return fP[4]; }       // return Y component of the momentum
    double Pz() const { return fP[5]; }       // return Z component of the momentum
    double E() const { return fP[6]; }        // return energy of the particle
    double S() const { return fP[7]; }        // return dS=l/p, l - decay length, defined if production vertex is set
    int GetQ() const { return fQ; }           // return charge of the particle
    double GetChi2() const { return fChi2; }  // return Chi2 of the fit
    int GetNDF() const { return fNDF; }       // return number of degrees of freedom

    double P2() const { return squaredNorm({fP[3], fP[4], fP[5]}); };
    double Mass() const {
        double mass2{E() * E() - P2()};
        if (mass2 < 0.) return -1.;
        return std::sqrt(mass2);
    }

    double GetParameter(int i) const { return fP[i]; }                 // return P[i] parameter
    double GetCovariance(int i) const { return fC[i]; }                // return C[i] element of the covariance matrix in the lower triangular form
    double GetCovariance(int i, int j) const { return fC[IJ(i, j)]; }  // return C[i,j] element of the covariance matrix
    SymMatrix<6> Cov_6x6() const { return Slice<36, 21>(fC); }

    std::pair<PCA, PCA> AddDaughterWithEnergyFit(const Particle &daughter, double bz, double chi2_threshold = 1E4);
    std::pair<PCA, PCA> AddDaughter(const Particle &daughter, double bz) {
        if (fNDF < -1) {  // first daughter -> just copy
            fNDF += 2;
            fQ = daughter.GetQ();
            for (int i{0}; i < 7; ++i) fP[i] = daughter.fP[i];
            for (int i{0}; i < 28; ++i) fC[i] = daughter.fC[i];
            return {{X(), Y(), Z(), Px(), Py(), Pz()},  //
                    {X(), Y(), Z(), Px(), Py(), Pz()}};
        }
        return AddDaughterWithEnergyFit(daughter, bz);
    }
    PCA AddProductionVertex(const Vector<3> &prod_vtx, const SymMatrix<3> &cov, double bz, double chi2_threshold = 1E4);
    void AddMassConstraint(double target_mass);

    void Print() {
        std::cout << "(X,Y,Z,S)    = " << fP[0] << "    " << fP[1] << "    " << fP[2] << "    " << fP[7] << '\n';
        std::cout << "(Px,Py,Pz,E) = " << fP[3] << "    " << fP[4] << "    " << fP[5] << "    " << fP[6] << '\n';
        std::cout << "Mass         = " << Mass() << '\n';
        std::cout << "Radius       = " << std::sqrt(fP[0] * fP[0] + fP[1] * fP[1]) << '\n';
        std::cout << "Chi2/NDF     = " << fChi2 << "/" << fNDF << '\n';
        std::cout << "CovMatrix    = ";
        size_t n_in_row{0};
        size_t max_n_row{1};
        for (size_t i{0}; i < 8 * 9 / 2; ++i) {
            std::cout << fC[i];
            ++n_in_row;
            if (n_in_row == max_n_row) {
                std::cout << '\n';
                if (i + 1 < 8 * 9 / 2) std::cout << "               ";
                n_in_row = 0;
                ++max_n_row;
            } else {
                std::cout << "    ";
            }
        }
    }

   protected:
    // Set Cxx=Cyy=Czz=100 and Css=1
    // Note: it will modify the state of the current `KF::Particle`
    void Initialize() {
        fC[0] = 100.;
        fC[2] = 100.;
        fC[5] = 100.;
        fC[35] = 1.;
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
        std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

        for (int i{0}; i < 6; ++i) fP[i] = param[i];
        double energy{std::sqrt(mass * mass + fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5])};
        fP[6] = energy;
        fP[7] = 0.;

        double h0{fP[3] / energy};
        double h1{fP[4] / energy};
        double h2{fP[5] / energy};

        for (int i{0}; i < 21; ++i) fC[i] = cov[i];
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
        PrintSymMatrix<8>(__FUNCTION__, "CovMatrix", fC);
        std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    }

    Result::Measurement GetMeasurement(const Particle &daughter, double bz) const;

    Result::Minimization MinimizeLinePoint(const Vector<3> &v) const;
    Result::Minimization MinimizeHelixPoint(const Vector<3> &v, double bz) const;
    Result::Minimization Minimize(const Vector<3> &v, double bz = 0) const {
        if (std::abs(fQ) < Const::AbsAlmostZero || std::abs(bz) < Const::AbsAlmostZero) {
            return MinimizeLinePoint(v);
        }
        return MinimizeHelixPoint(v, bz);
    }

    std::pair<Result::Minimization, Result::Minimization> MinimizeLineLine(const Particle &p) const;
    std::pair<Result::Minimization, Result::Minimization> MinimizeHelixHelix(const Particle &p, double bz) const;
    std::pair<Result::Minimization, Result::Minimization> Minimize(const Particle &p, double bz = 0.) const {
        if ((std::abs(fQ) < Const::AbsAlmostZero && std::abs(p.fQ) < Const::AbsAlmostZero) || std::abs(bz) < Const::AbsAlmostZero) {
            return MinimizeLineLine(p);
        }
        return MinimizeHelixHelix(p, bz);
    }

    Result::Transport TransportBz(const Result::Minimization &min, double bz) const;
    Result::Transport TransportLine(const Result::Minimization &min) const;
    Result::Transport Transport(const Result::Minimization &min, double bz = 0.) const {
        if (std::abs(fQ) < Const::AbsAlmostZero || std::abs(bz) < Const::AbsAlmostZero) {
            return TransportLine(min);
        }
        return TransportBz(min, bz);
    }

    SymMatrix<8> fC{};  // lower-triangular part of the symmetric 8x8 covariance matrix
    Vector<8> fP{};     // particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}
    double fChi2{0.};   // chi2
    int fNDF{-3};       // number of degrees of freedom
    int fQ{0};          // the charge of the particle in units of elementary charge
};

}  // namespace KF

#endif  // KFPARTICLE_HXX
