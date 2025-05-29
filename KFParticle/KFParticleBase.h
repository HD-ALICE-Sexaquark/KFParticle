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

#ifndef KFPARTICLEBASE_H
#define KFPARTICLEBASE_H

#include <cmath>
#include <iostream>

#include "KFPTrack.h"

constexpr float kCLight{0.000299792458};
constexpr float LocalSmall{1.E-8};

// @class KFParticleBase
// @brief The base of KFParticle class, describes particle objects.
// @author S.Gorbunov, I.Kisel, M.Zyzak
// @date 05.02.2019
// @version 1.0
//
// Contains the main mathematics of the KFParticle.
class alignas(32) KFParticleBase {
   public:
    KFParticleBase() : fChi2{0}, fNDF{-3}, fSFromDecay{0}, fSumDaughterMass{0}, fMassHypo{-1}, fAtProductionVertex{0}, fQ{0}, fConstructMethod{0} {
        Initialize();
    }
    KFParticleBase(const KFPTrack &track, const float mass) {
        track.XvYvZv(fP);
        track.PxPyPz(fP + 3);
        float energy = std::sqrt(mass * mass + fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]);
        fP[6] = energy;
        fP[7] = 0.;

        track.GetCovarianceXYZPxPyPz(fC);
        float h0 = fP[3] / energy;
        float h1 = fP[4] / energy;
        float h2 = fP[5] / energy;
        fC[21] = h0 * fC[6] + h1 * fC[10] + h2 * fC[15];
        fC[22] = h0 * fC[7] + h1 * fC[11] + h2 * fC[16];
        fC[23] = h0 * fC[8] + h1 * fC[12] + h2 * fC[17];
        fC[24] = h0 * fC[9] + h1 * fC[13] + h2 * fC[18];
        fC[25] = h0 * fC[13] + h1 * fC[14] + h2 * fC[19];
        fC[26] = h0 * fC[18] + h1 * fC[19] + h2 * fC[20];
        fC[27] = (h0 * h0 * fC[9] + h1 * h1 * fC[14] + h2 * h2 * fC[20] + 2 * (h0 * h1 * fC[13] + h0 * h2 * fC[18] + h1 * h2 * fC[19]));
        for (int i = 28; i < 36; ++i) fC[i] = 0.;
        fC[35] = 1.;

        fQ = track.Charge();
        fChi2 = track.GetChi2();
        fNDF = track.GetNDF();
        fSFromDecay = 0;
        fSumDaughterMass = 0;
        fMassHypo = 0;
        fAtProductionVertex = false;
        fConstructMethod = 0;
    }

    ~KFParticleBase() = default;

    void Initialize(const float param[], const float cov[], int charge, float mass);
    void Initialize();

    void Construct(float bz, const KFParticleBase *v_daughters[], int n_daughters, const KFParticleBase *prod_vtx = 0, float mass = -1.);

    float GetDStoPointLine(const float xyz[3], float dsdr[6]) const;
    float GetDStoPointBz(float bz, const float xyz[3], float dsdr[6]) const;

    void GetDStoParticleLine(const KFParticleBase &p, float ds[2], float dsdr[4][6]) const;
    void GetDStoParticleBz(float bz, const KFParticleBase &p, float ds[2], float dsdr[4][6]) const;

    void TransportToDS(float bz, float ds, const float *dsdr);
    void TransportBz(float bz, float ds, const float *dsdr, float p[], float c[], float *dsdr1 = 0, float *f = 0, float *f1 = 0) const;
    void TransportToDecayVertex(float bz);
    void TransportToProductionVertex(float bz);

    // define the construction method for the current particle (see description of fConstructMethod)
    void SetConstructMethod(int m) { fConstructMethod = m; }
    void SetMassHypo(float m) { fMassHypo = m; }            // set the mass hypothesis to the particle, is used when fConstructMethod = 2
    const float &GetMassHypo() const { return fMassHypo; }  // return the mass hypothesis
    const float &GetSumDaughterMass() const { return fSumDaughterMass; }  // return the sum of masses of the daughters

    // Accessors
    float GetX() const { return fP[0]; }     // return X coordinate of the particle, fP[0]
    float GetY() const { return fP[1]; }     // return Y coordinate of the particle, fP[1]
    float GetZ() const { return fP[2]; }     // return Z coordinate of the particle, fP[2]
    float GetPx() const { return fP[3]; }    // return X component of the momentum, fP[3]
    float GetPy() const { return fP[4]; }    // return Y component of the momentum, fP[4]
    float GetPz() const { return fP[5]; }    // return Z component of the momentum, fP[5]
    float GetE() const { return fP[6]; }     // return energy of the particle, fP[6]
    float GetS() const { return fP[7]; }     // return dS=l/p, l - decay length, fP[7], defined if production vertex is set
    char GetQ() const { return fQ; }         // return charge of the particle
    float GetChi2() const { return fChi2; }  // return Chi2 of the fit
    int GetNDF() const { return fNDF; }      // return number of decrease of freedom

    const float &X() const { return fP[0]; }     // return X coordinate of the particle, fP[0]
    const float &Y() const { return fP[1]; }     // return Y coordinate of the particle, fP[1]
    const float &Z() const { return fP[2]; }     // return Z coordinate of the particle, fP[2]
    const float &Px() const { return fP[3]; }    // return X component of the momentum, fP[3]
    const float &Py() const { return fP[4]; }    // return Y component of the momentum, fP[4]
    const float &Pz() const { return fP[5]; }    // return Z component of the momentum, fP[5]
    const float &E() const { return fP[6]; }     // return energy of the particle, fP[6]
    const float &S() const { return fP[7]; }     // return dS=l/p, l - decay length, fP[7], defined if production vertex is set
    const char &Q() const { return fQ; }         // return charge of the particle
    const float &Chi2() const { return fChi2; }  // return Chi2 of the fit
    const int &NDF() const { return fNDF; }      // return number of decrease of freedom

    float GetParameter(int i) const { return fP[i]; }                 // return P[i] parameter
    float GetCovariance(int i) const { return fC[i]; }                // return C[i] element of the covariance matrix in the lower triangular form
    float GetCovariance(int i, int j) const { return fC[IJ(i, j)]; }  // return C[i,j] element of the covariance matrix

    //  MODIFIERS
    float &X() { return fP[0]; }     // modifier of X coordinate of the particle, fP[0]
    float &Y() { return fP[1]; }     // modifier of Y coordinate of the particle, fP[1]
    float &Z() { return fP[2]; }     // modifier of Z coordinate of the particle, fP[2]
    float &Px() { return fP[3]; }    // modifier of X component of the momentum, fP[3]
    float &Py() { return fP[4]; }    // modifier of Y component of the momentum, fP[4]
    float &Pz() { return fP[5]; }    // modifier of Z component of the momentum, fP[5]
    float &E() { return fP[6]; }     // modifier of energy of the particle, fP[6]
    float &S() { return fP[7]; }     // modifier of dS=l/p, l - decay length, fP[7], defined if production vertex is set
    char &Q() { return fQ; }         // modifier of charge of the particle
    float &Chi2() { return fChi2; }  // modifier of Chi2 of the fit
    int &NDF() { return fNDF; }      // modifier of number of decrease of freedom

    float &Parameter(int i) { return fP[i]; }                 // modifier of P[i] parameter
    float &Covariance(int i) { return fC[i]; }                // modifier of C[i] element of the covariance matrix in the lower triangular form
    float &Covariance(int i, int j) { return fC[IJ(i, j)]; }  // modifier of C[i,j] element of the covariance matrix

    // CONSTRUCTION OF THE PARTICLE BY ITS DAUGHTERS AND MOTHER USING THE KALMAN FILTER METHOD

    void AddDaughter(float bz, const KFParticleBase &daughter);

    void AddDaughterWithEnergyFit(float bz, const KFParticleBase &daughter);
    void AddDaughterWithEnergyFitMC(float bz, const KFParticleBase &daughter);

    void SetProductionVertex(const KFParticleBase &vtx, float bz);

    void SetNonlinearMassConstraint(float mass);
    void SetMassConstraint(float mass, float sigma_mass = 0);

    void SetNoDecayLength(float bz);  // set no decay length for resonances

    // OTHER UTILITIES

    static void InvertCholetsky3(float a[6]);
    static void MultQSQt(const float q[], const float s[], float s_out[], const int k_n);

    template <typename D, int N>
    void PrintVector(std::string_view name, D arr[N]) const {
        std::cout << "  " << name << " = ";
        for (int i{0}; i < N; ++i) {
            std::cout << arr[i];
            if (i + 1 < N)
                std::cout << "    ";
            else
                std::cout << '\n';
        }
    }

    template <typename D, int N>
    void PrintSymMatrix(std::string_view name, D arr[N]) const {
        std::cout << "  " << name << " =\n";
        int n_in_row{0};
        int max_n_row{1};
        for (int i{0}; i < N; ++i) {
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

    template <typename D, int N>
    void PrintNonSymMatrix(std::string_view name, D arr[N]) const {
        std::cout << "  " << name << " =\n";
        int counter{0};
        for (int i{0}; i < N; ++i) {
            std::cout << arr[i];
            counter++;
            if (counter < std::sqrt(N))
                std::cout << "    ";
            else {
                std::cout << '\n';
                counter = 0;
            }
        }
    }
    template <typename D, int N, int M>
    void PrintMatrix(std::string_view name, D arr[N][M]) const {
        std::cout << "  " << name << " =\n";
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

    template <typename D, int N>
    void PrintJoinedMatrix(std::string_view name, D arr1[N], D arr2[N], D arr3[N]) const {
        std::cout << "  " << name << " =\n";
        for (int i = 0; i < N; ++i) {
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

   protected:
    // Convert a pair of indices {i,j} of the covariance matrix to one index corresponding to the triangular form
    static int IJ(int i, int j) { return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }
    // Return an element of the covariance matrix with {i,j} indices
    float &Cij(int i, int j) { return fC[IJ(i, j)]; }
    void TransportLine(float S, const float *dsdr, float P[], float C[], float *dsdr1, float *F, float *F1) const;
    bool GetMeasurement(float bz, const KFParticleBase &daughter, float m[], float V[], float D[3][3]);
    void SetMassConstraint(float *mP, float *mC, float mJ[7][7], float mass);

    float fP[8];               // particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}
    float fC[36];              // low-triangle covariance matrix of fP
    float fChi2;               // chi^2
    int fNDF;                  // number of degrees of freedom
    float fSFromDecay;         // distance from the decay vertex to the current position
    float fSumDaughterMass;    // sum of the daughter particles masses Needed to set the constraint on the minimum mass during particle construction
    float fMassHypo;           // the mass hypothesis, used for the constraints during particle construction
    bool fAtProductionVertex;  // flag shows if particle is at the production point
    char fQ;                   // the charge of the particle in the units of the elementary charge

    // Determine particle construction method.
    // 0 - Energy considered as an independent variable, fitted independently from momentum, without any constraints on mass
    // 2 - Energy considered as an independent variable, fitted independently from momentum, with constraints on mass of daughter particle
    char fConstructMethod;
};

#endif
