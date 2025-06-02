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

#ifndef KFPARTICLE_H
#define KFPARTICLE_H

#include <cmath>
#include <iostream>
#include <utility>

#include "KFParticle_Math.hxx"

namespace KF {

struct alignas(32) Cache {
    Vector<3> dir{};
    Vector<3> pca{};
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

    void Initialize();
    void Initialize(const Vector<6> &param, const SymMatrix<6> &cov, int charge, double mass);

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

    // define the construction method for the current particle (see description of fConstructMethod)
    void SetConstructMethod(int m) { fConstructMethod = m; }
    void SetMassHypo(double m) { fMassHypo = m; }            // set the mass hypothesis to the particle, is used when fConstructMethod = 2
    const double &GetMassHypo() const { return fMassHypo; }  // return the mass hypothesis
    const double &GetSumDaughterMass() const { return fSumDaughterMass; }  // return the sum of masses of the daughters

    // Accessors
    Vector<3> XYZ() const { return {fP[0], fP[1], fP[2]}; }                             // PENDING: copies
    Vector<6> XYZPxPyPz() const { return {fP[0], fP[1], fP[2], fP[3], fP[4], fP[5]}; }  // PENDING: copies
    double GetX() const { return fP[0]; }                                               // return X coordinate of the particle, fP[0]
    double GetY() const { return fP[1]; }                                               // return Y coordinate of the particle, fP[1]
    double GetZ() const { return fP[2]; }                                               // return Z coordinate of the particle, fP[2]
    double GetPx() const { return fP[3]; }                                              // return X component of the momentum, fP[3]
    double GetPy() const { return fP[4]; }                                              // return Y component of the momentum, fP[4]
    double GetPz() const { return fP[5]; }                                              // return Z component of the momentum, fP[5]
    double GetE() const { return fP[6]; }                                               // return energy of the particle, fP[6]
    double GetS() const { return fP[7]; }     // return dS=l/p, l - decay length, fP[7], defined if production vertex is set
    int GetQ() const { return fQ; }           // return charge of the particle
    double GetChi2() const { return fChi2; }  // return Chi2 of the fit
    int GetNDF() const { return fNDF; }       // return number of decrease of freedom

    const double &X() const { return fP[0]; }     // return X coordinate of the particle, fP[0]
    const double &Y() const { return fP[1]; }     // return Y coordinate of the particle, fP[1]
    const double &Z() const { return fP[2]; }     // return Z coordinate of the particle, fP[2]
    const double &Px() const { return fP[3]; }    // return X component of the momentum, fP[3]
    const double &Py() const { return fP[4]; }    // return Y component of the momentum, fP[4]
    const double &Pz() const { return fP[5]; }    // return Z component of the momentum, fP[5]
    const double &E() const { return fP[6]; }     // return energy of the particle, fP[6]
    const double &S() const { return fP[7]; }     // return dS=l/p, l - decay length, fP[7], defined if production vertex is set
    const int &Q() const { return fQ; }           // return charge of the particle
    const double &Chi2() const { return fChi2; }  // return Chi2 of the fit
    const int &NDF() const { return fNDF; }       // return number of decrease of freedom

    double GetParameter(int i) const { return fP[i]; }                 // return P[i] parameter
    double GetCovariance(int i) const { return fC[i]; }                // return C[i] element of the covariance matrix in the lower triangular form
    double GetCovariance(int i, int j) const { return fC[IJ(i, j)]; }  // return C[i,j] element of the covariance matrix
    SymMatrix<6> Cov_6x6() const {
        SymMatrix<6> cov;
        for (int i{0}; i < (6 * (6 + 1) / 2); ++i) cov[i] = fC[i];
        return cov;
    }

    //  MODIFIERS
    double &X() { return fP[0]; }     // modifier of X coordinate of the particle, fP[0]
    double &Y() { return fP[1]; }     // modifier of Y coordinate of the particle, fP[1]
    double &Z() { return fP[2]; }     // modifier of Z coordinate of the particle, fP[2]
    double &Px() { return fP[3]; }    // modifier of X component of the momentum, fP[3]
    double &Py() { return fP[4]; }    // modifier of Y component of the momentum, fP[4]
    double &Pz() { return fP[5]; }    // modifier of Z component of the momentum, fP[5]
    double &E() { return fP[6]; }     // modifier of energy of the particle, fP[6]
    double &S() { return fP[7]; }     // modifier of dS=l/p, l - decay length, fP[7], defined if production vertex is set
    int &Q() { return fQ; }           // modifier of charge of the particle
    double &Chi2() { return fChi2; }  // modifier of Chi2 of the fit
    int &NDF() { return fNDF; }       // modifier of number of decrease of freedom

    double &Parameter(int i) { return fP[i]; }                 // modifier of P[i] parameter
    double &Covariance(int i) { return fC[i]; }                // modifier of C[i] element of the covariance matrix in the lower triangular form
    double &Covariance(int i, int j) { return fC[IJ(i, j)]; }  // modifier of C[i,j] element of the covariance matrix

    void AddDaughterWithEnergyFit(const Particle &daughter, double bz, double chi2_threshold = 1E4);
    void AddDaughterWithEnergyFitMC(const Particle &daughter, double bz, double chi2_threshold = 1E4);
    void AddDaughter(const Particle &daughter, double bz) {
        if (fNDF < -1) {  // first daughter -> just copy
            fNDF += 2;
            fQ = daughter.GetQ();
            for (int i{0}; i < 7; ++i) fP[i] = daughter.fP[i];
            for (int i{0}; i < 28; ++i) fC[i] = daughter.fC[i];
            fMassHypo = daughter.fMassHypo;
            fSumDaughterMass = daughter.fSumDaughterMass;
            return;
        }

        if (fConstructMethod == 0)
            AddDaughterWithEnergyFit(daughter, bz);
        else if (fConstructMethod == 2)
            AddDaughterWithEnergyFitMC(daughter, bz);

        fSumDaughterMass += daughter.fSumDaughterMass;
        fMassHypo = -1.;
    }

    void AddProductionVertex(const Vector<3> &prod_vtx, const SymMatrix<3> &cov, double bz, double chi2_threshold = 1E4);

    Result::MassConstraint SetMassConstraint(double mass) const;
    void SetLinearMassConstraint(double mass, double sigma_mass = 0);

    void Print() {
        std::cout << "(X,Y,Z)      = " << fP[0] << "    " << fP[1] << "    " << fP[2] << '\n';
        std::cout << "Radius       = " << std::sqrt(fP[0] * fP[0] + fP[1] * fP[1]) << '\n';
        std::cout << "(Px,Py,Pz,E) = " << fP[3] << "    " << fP[4] << "    " << fP[5] << "    " << fP[6] << '\n';
        std::cout << "Mass         = " << std::sqrt(fP[6] * fP[6] - fP[3] * fP[3] - fP[4] * fP[4] - fP[5] * fP[5]) << '\n';
        std::cout << "Chi2/NDF     = " << fChi2 << "/" << fNDF << '\n';
    }

   protected:
    double &Cij(int i, int j) { return fC[IJ(i, j)]; }
    SymMatrix<8> fC{};            // low-triangle covariance matrix of fP
    Vector<8> fP{};               // particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}
    double fChi2{0.};             // chi2
    double fSumDaughterMass{0.};  // sum of the daughter particles masses Needed to set the constraint on the minimum mass during particle
    double fMassHypo{-1.};        // the mass hypothesis, used for the constraints during particle construction
    int fNDF{-3};                 // number of degrees of freedom
    int fQ{0};                    // the charge of the particle in units of elementary charge

    // Determine particle construction method.
    // 0 - Energy considered as an independent variable, fitted independently from momentum, without any constraints on mass
    // 2 - Energy considered as an independent variable, fitted independently from momentum, with constraints on mass of daughter particle
    int fConstructMethod{0};
};

}  // namespace KF

#endif
