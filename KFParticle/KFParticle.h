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

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

namespace KF {

template <size_t N>
using Vector = std::array<double, N>;

template <size_t N, size_t M>
using Matrix = std::array<std::array<double, N>, M>;

template <size_t N>
using SymMatrix = std::array<double, N *(N + 1) / 2>;

namespace Const {
constexpr double Kappa{0.000299792458};  // (GeV/c) / (kG/cm)
constexpr double AbsAlmostZero{1.E-8};
}  // namespace Const

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

    bool GetMeasurement(double bz, const Particle &daughter, Vector<8> &m, SymMatrix<8> &V, Matrix<3, 3> &D);

    void Construct(double bz, const Particle *v_daughters[], int n_daughters, const Particle *parent = nullptr, double mass = -1.);

    double GetDStoPointLine(const Vector<3> &xyz, Vector<6> &ds_dr) const;
    double GetDStoPointBz(double bz, const Vector<3> &xyz, Vector<6> &ds_dr) const;

    void GetDStoParticleLine(const Particle &p, double ds[2], Vector<6> dsdr[4]) const;
    void GetDStoParticleBz(double bz, const Particle &p, double dS[2], Vector<6> ds_dr[4]) const;

    void TransportToDS(double bz, double ds, const Vector<6> &dsdr);
    void TransportBz(double bz, double ds, const Vector<6> &dsdr, Vector<8> &P, SymMatrix<8> &C, const Vector<6> &dsdr1, Matrix<6, 6> &jacob,
                     Matrix<6, 6> &corr) const;
    void TransportToDecayVertex(double bz);
    void TransportToProductionVertex(double bz);
    void TransportLine(double ds, const Vector<6> &ds_dr, Vector<8> &P, SymMatrix<8> &C, const Vector<6> &ds_dr1, Matrix<6, 6> &jacob,
                       Matrix<6, 6> &corr) const;

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
        for (int i{0}; i < (6 * (6 + 1) / 2); i++) cov[i] = fC[i];
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

    void AddDaughter(double bz, const Particle &daughter);
    void AddDaughterWithEnergyFit(double bz, const Particle &daughter);
    void AddDaughterWithEnergyFitMC(double bz, const Particle &daughter);

    void SetNonlinearMassConstraint(double mass);
    void SetMassConstraint(double mass, double sigma_mass = 0);
    void SetMassConstraint(Vector<8> &mP, SymMatrix<8> &mC, Matrix<7, 7> &mJ, double mass) const;

    void SetProductionVertex(const Particle &vtx, double bz);
    void SetNoDecayLength(double bz);  // set no decay length for resonances

    static void InvertCholesky3(SymMatrix<3> &a);

    template <int N>
    SymMatrix<N> MultQSQt(const Matrix<N, N> &Q, const SymMatrix<N> &S) const;

    template <int N>
    void PrintVector(std::string_view name, const Vector<N> arr) const {
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
    void PrintSymMatrix(std::string_view name, const SymMatrix<N> arr) const {
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
    void PrintMatrix(std::string_view name, const Matrix<N, M> &arr) const {
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
    void PrintJoinedMatrix(std::string_view name, const Vector<N> &arr1, const Vector<N> &arr2, const Vector<N> &arr3) const {
        std::cout << name << " =\n";
        for (int i{0}; i < N; ++i) {
            std::cout << "  " << arr1[i] << "    " << arr2[i] << "    " << arr3[i] << '\n';
        }
    }

    void Print() {
        std::cout << "# V0" << '\n';
        std::cout << "(X,Y,Z)      = " << fP[0] << "    " << fP[1] << "    " << fP[2] << '\n';
        std::cout << "Radius       = " << std::sqrt(fP[0] * fP[0] + fP[1] * fP[1]) << '\n';
        std::cout << "(Px,Py,Pz,E) = " << fP[3] << "    " << fP[4] << "    " << fP[5] << "    " << fP[6] << '\n';
        std::cout << "Mass         = " << std::sqrt(fP[6] * fP[6] - fP[3] * fP[3] - fP[4] * fP[4] - fP[5] * fP[5]) << '\n';
        std::cout << "Chi2/NDF     = " << fChi2 << '\n';
    }

    template <size_t N>
    static std::array<double, N> Zero() {
        std::array<double, N> vec;
        vec.fill(0);
        return vec;
    }

    template <size_t N, size_t M>
    static Matrix<N, M> Zero() {
        Matrix<N, M> mat;
        for (auto &row : mat) row.fill(0.);
        return mat;
    }

   protected:
    // Convert a pair of indices {i,j} of the covariance matrix to one index corresponding to the triangular form
    static int IJ(int i, int j) { return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }
    // Return an element of the covariance matrix with {i,j} indices
    double &Cij(int i, int j) { return fC[IJ(i, j)]; }

    Vector<8> fP{0., 0., 0., 0., 0., 0., 0., 0.};  // particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}
    SymMatrix<8> fC{Zero<8 * 9 / 2>()};            // low-triangle covariance matrix of fP
    double fChi2{0.};                              // chi^2
    int fNDF{-3};                                  // number of degrees of freedom
    double fSFromDecay{0.};                        // distance from the decay vertex to the current position
    double fSumDaughterMass{0.};      // sum of the daughter particles masses Needed to set the constraint on the minimum mass during particle
                                      // construction
    double fMassHypo{-1.};            // the mass hypothesis, used for the constraints during particle construction
    bool fAtProductionVertex{false};  // flag shows if particle is at the production point
    int fQ{0};                        // the charge of the particle in units of elementary charge

    // Determine particle construction method.
    // 0 - Energy considered as an independent variable, fitted independently from momentum, without any constraints on mass
    // 2 - Energy considered as an independent variable, fitted independently from momentum, with constraints on mass of daughter particle
    int fConstructMethod{0};
};

}  // namespace KF

#endif
