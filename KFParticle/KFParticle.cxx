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

#include "KFParticle.h"
#include "KFParticle_Math.hxx"

#include <algorithm>
#include <tuple>

namespace KF {

// Set Cxx=Cyy=Czz=100. and Css=1.
void Particle::Initialize() {
    fC[0] = 100.;
    fC[2] = 100.;
    fC[5] = 100.;
    fC[35] = 1.;
}

// Set the parameters of the particle:
//  \param[in] param[6] = { X, Y, Z, Px, Py, Pz } - position and momentum
//  \param[in] cov[21]  - lower-triangular part of the 6x6 covariance matrix:
//           (  0  .  .  .  .  . )
//           (  1  2  .  .  .  . )
// Cov[21] = (  3  4  5  .  .  . )
//           (  6  7  8  9  .  . )
//           ( 10 11 12 13 14  . )
//           ( 15 16 17 18 19 20 )
//  \param[in] charge - charge of the particle in elementary charge units
//  \param[in] mass - the mass hypothesis
void Particle::Initialize(const Vector<6>& param, const SymMatrix<6>& cov, int charge, double mass) {
#if KF_DEBUG
    std::cout << "-- starting (Initialize) --" << '\n';
#endif

    for (int i{0}; i < 6; ++i) fP[i] = param[i];
    double energy = std::sqrt(mass * mass + fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]);
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
    fSumDaughterMass = mass;
    fMassHypo = mass;
#if KF_DEBUG
    PrintSymMatrix<8>("(Initialize) CovMatrix", fC);
    std::cout << "-- finished (Initialize) --" << '\n';
#endif
}

// Obtain the measurements from the current particle and the daughter to be added for the Kalman filter
// mathematics. If these are two first daughters they are transported to the point of the closest approach,
// if the third or higher daughter is added it is transported to the DCA point of the already constructed
// vertex. The correlations are taken into account in the covariance matrices of both measurements,
// the correlation matrix of two measurements is also calculated.
// Input:
// - daughter - the daughter particle to be added, stays unchanged
// - bz - z-component of magnetic field
// Return:
// - P[8] - the output parameters of the daughter particle at the DCA point
// - V[36] - the output covariance matrix of the daughter parameters, takes into account the correlation
// - D[3][3] - the correlation matrix between the current and daughter particles
Result::Measurement Particle::GetMeasurement(const Particle& daughter, double bz) const {
#if KF_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    if (fNDF == -1) {
        auto [min1, min2] = MinimizeHelixHelix(daughter, bz);

        auto tpr1 = TransportBz(min1, bz);
        auto tpr2 = daughter.TransportBz(min2, bz);

#if KF_DEBUG
        PrintMatrix<6, 6>("(GetMeasurement) F1", tpr1.jacob);
        PrintMatrix<6, 6>("(GetMeasurement) F2", tpr1.corr);
        PrintMatrix<6, 6>("(GetMeasurement) F3", tpr2.corr);
        PrintMatrix<6, 6>("(GetMeasurement) F4", tpr2.jacob);
#endif
        SymMatrix<6> V0Tmp{MultQSQt<6>(tpr1.corr, daughter.Cov_6x6())};
        SymMatrix<6> V1Tmp{MultQSQt<6>(tpr2.corr, Cov_6x6())};
#if KF_DEBUG
        PrintSymMatrix<6>("(GetMeasurement) V0Tmp", V0Tmp);
        PrintSymMatrix<6>("(GetMeasurement) V1Tmp", V1Tmp);
#endif

        Result::Measurement meas{tpr1.C, tpr2.C, Matrix<3, 3>{}, tpr1.P, tpr2.P};

        for (int iC{0}; iC < 21; ++iC) {
            meas.C1[iC] += V0Tmp[iC];
            meas.C2[iC] += V1Tmp[iC];
        }

        Matrix<6, 6> C1F1T;
        for (int i{0}; i < 6; ++i) {
            for (int j{0}; j < 6; ++j) {
                C1F1T[i][j] = 0.;
                for (int k{0}; k < 6; ++k) {
                    C1F1T[i][j] += fC[IJ(i, k)] * tpr1.jacob[j][k];
                }
            }
        }
        Matrix<6, 6> F3C1F1T;
        for (int i{0}; i < 6; ++i) {
            for (int j{0}; j < 6; ++j) {
                F3C1F1T[i][j] = 0.;
                for (int k{0}; k < 6; ++k) {
                    F3C1F1T[i][j] += tpr2.corr[i][k] * C1F1T[k][j];
                }
            }
        }
        Matrix<6, 6> C2F2T;
        for (int i{0}; i < 6; ++i) {
            for (int j{0}; j < 6; ++j) {
                C2F2T[i][j] = 0.;
                for (int k{0}; k < 6; ++k) {
                    C2F2T[i][j] += daughter.fC[IJ(i, k)] * tpr1.corr[j][k];
                }
            }
        }
        for (int i{0}; i < 3; ++i) {
            for (int j{0}; j < 3; ++j) {
                meas.D[i][j] = F3C1F1T[i][j];
                for (int k{0}; k < 6; ++k) {
                    meas.D[i][j] += tpr2.jacob[i][k] * C2F2T[k][j];
                }
            }
        }
#if KF_DEBUG
        PrintVector<8>("(GetMeasurement) meas.P1", meas.P1);
        PrintVector<8>("(GetMeasurement) meas.P2", meas.P2);
        PrintSymMatrix<8>("(GetMeasurement) meas.C1", meas.C1);
        PrintSymMatrix<8>("(GetMeasurement) meas.C2", meas.C2);
        PrintMatrix<6, 6>("(GetMeasurement) C1F1T", C1F1T);
        PrintMatrix<6, 6>("(GetMeasurement) F3C1F1T", F3C1F1T);
        PrintMatrix<6, 6>("(GetMeasurement) C2F2T", C2F2T);
        PrintMatrix<3, 3>("(GetMeasurement) D", meas.D);
        std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif

        return meas;
    }
    Result::Measurement meas;  // ULTRA PENDING

    Vector<6> ds_dr{0., 0., 0., 0., 0., 0.};
    // double ds{daughter.GetDStoPointBz(bz, XYZ(), ds_dr)}; // PENDING

    Vector<6> dsdp{-ds_dr[0], -ds_dr[1], -ds_dr[2], 0., 0., 0.};

    Matrix<6, 6> F{};
    Matrix<6, 6> F1{};
    // daughter.TransportBz(bz, ds, ds_dr, m, V, dsdp, F, F1); // PENDING

    Matrix<3, 6> VFT;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 6; ++j) {
            VFT[i][j] = 0.;
            for (int k{0}; k < 3; ++k) {
                VFT[i][j] += fC[IJ(i, k)] * F1[j][k];
            }
        }
    }

    Matrix<6, 6> FVFT;
    for (int i{0}; i < 6; ++i) {
        for (int j{0}; j < 6; ++j) {
            FVFT[i][j] = 0.;
            for (int k{0}; k < 3; ++k) {
                FVFT[i][j] += F1[i][k] * VFT[k][j];
            }
        }
    }

    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            for (int k{0}; k < 3; ++k) {
                meas.D[i][j] += fC[IJ(j, k)] * F1[i][k];
            }
        }
    }
    meas.C2[0] += FVFT[0][0];
    meas.C2[1] += FVFT[1][0];
    meas.C2[2] += FVFT[1][1];
    meas.C2[3] += FVFT[2][0];
    meas.C2[4] += FVFT[2][1];
    meas.C2[5] += FVFT[2][2];
    return meas;
}

// Add daughter to the current particle. Depending on the selected construction method uses:
// 1) Either simplifyed fast mathematics which consideres momentum and energy as
// independent variables and thus ignores constraint on the fixed mass (fConstructMethod = 0).
// In this case the mass of the daughter particle can be corrupted when the constructed vertex
// is added as the measurement and the mass of the output short-lived particle can become
// unphysical - smaller then the threshold. Implemented in the
// AddDaughterWithEnergyFit() function
// 2) Or slower but correct mathematics which requires that the masses of daughter particles
// stays fixed in the construction process (fConstructMethod = 2). Implemented in the
// AddDaughterWithEnergyFitMC() function.
// \param[in] daughter - the daughter particle
void Particle::AddDaughter(const Particle& daughter, double bz) {

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

// Add daughter to the current particle. Uses simplifyed fast mathematics which consideres momentum
// and energy as independent variables and thus ignores constraint on the fixed mass.
// In this case the mass of the daughter particle can be corrupted when the constructed vertex
// is added as the measurement and the mass of the output short-lived particle can become
// unphysical - smaller then the threshold.
// Input:
// - daughter       : the daughter particle
// - bz             : z-component of magnetic field
// - chi2_threshold : do an early cut of chi2
// CONSIDERATION: it will modify the state of the current KF::Particle
void Particle::AddDaughterWithEnergyFit(const Particle& daughter, double bz, double chi2_threshold) {
#if KF_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    auto meas = GetMeasurement(daughter, bz);

    SymMatrix<3> mS{meas.C1[0] + meas.C2[0],                           //
                    meas.C1[1] + meas.C2[1], meas.C1[2] + meas.C2[2],  //
                    meas.C1[3] + meas.C2[3], meas.C1[4] + meas.C2[4], meas.C1[5] + meas.C2[5]};

    InvertCholesky3(mS);
#if KF_DEBUG
    PrintVector<8>("(AddDaughterWithEnergyFit) meas.P1", meas.P1);
    PrintSymMatrix<8>("(AddDaughterWithEnergyFit) meas.C1", meas.C1);
    PrintVector<8>("(AddDaughterWithEnergyFit) meas.P2", meas.P2);
    PrintSymMatrix<8>("(AddDaughterWithEnergyFit) meas.C2", meas.C2);
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) meas.D", meas.D);
    PrintSymMatrix<3>("(AddDaughterWithEnergyFit) mS", mS);
#endif

    Vector<3> zeta{meas.P2[0] - meas.P1[0], meas.P2[1] - meas.P1[1], meas.P2[2] - meas.P1[2]};

    double dChi2{(mS[0] * zeta[0] + mS[1] * zeta[1] + mS[3] * zeta[2]) * zeta[0] + (mS[1] * zeta[0] + mS[2] * zeta[1] + mS[4] * zeta[2]) * zeta[1] +
                 (mS[3] * zeta[0] + mS[4] * zeta[1] + mS[5] * zeta[2]) * zeta[2]};
#if KF_DEBUG
    PrintVector<3>("(AddDaughterWithEnergyFit) zeta", zeta);
    PrintValue(__FUNCTION__, "dChi2", dChi2);
#endif
    if (dChi2 > chi2_threshold) return;

    // update current particle state //

    for (int i{0}; i < 8; ++i) fP[i] = meas.P1[i];
    for (int i{0}; i < 28; ++i) fC[i] = meas.C1[i];

    // add daughter's momentum to particle momentum //

    fP[3] += meas.P2[3];
    fP[4] += meas.P2[4];
    fP[5] += meas.P2[5];
    fP[6] += meas.P2[6];

    fC[9] += meas.C2[9];
    fC[13] += meas.C2[13];
    fC[14] += meas.C2[14];
    fC[18] += meas.C2[18];
    fC[19] += meas.C2[19];
    fC[20] += meas.C2[20];
    fC[24] += meas.C2[24];
    fC[25] += meas.C2[25];
    fC[26] += meas.C2[26];
    fC[27] += meas.C2[27];

#if KF_DEBUG
    PrintVector<8>("(AddDaughterWithEnergyFit) fP (before Kalman gain)", fP);
    PrintSymMatrix<8>("(AddDaughterWithEnergyFit) fC (before Kalman gain)", fC);
#endif

    // CHt = CH' - D'
    // Kalman gain K = mCH'*S
    // New estimation of the vertex position r += K*zeta
    // New covariance matrix C -= K*(mCH')'

    Vector<7> mCHt0{
        meas.C1[0], meas.C1[1], meas.C1[3], meas.C1[6] - meas.C2[6], meas.C1[10] - meas.C2[10], meas.C1[15] - meas.C2[15], meas.C1[21] - meas.C2[21]};
    Vector<7> mCHt1{
        meas.C1[1], meas.C1[2], meas.C1[4], meas.C1[7] - meas.C2[7], meas.C1[11] - meas.C2[11], meas.C1[16] - meas.C2[16], meas.C1[22] - meas.C2[22]};
    Vector<7> mCHt2{
        meas.C1[3], meas.C1[4], meas.C1[5], meas.C1[8] - meas.C2[8], meas.C1[12] - meas.C2[12], meas.C1[17] - meas.C2[17], meas.C1[23] - meas.C2[23]};

    Vector<7> k0;
    Vector<7> k1;
    Vector<7> k2;
    for (int i{0}; i < 7; ++i) {
        k0[i] = mCHt0[i] * mS[0] + mCHt1[i] * mS[1] + mCHt2[i] * mS[3];
        k1[i] = mCHt0[i] * mS[1] + mCHt1[i] * mS[2] + mCHt2[i] * mS[4];
        k2[i] = mCHt0[i] * mS[3] + mCHt1[i] * mS[4] + mCHt2[i] * mS[5];
    }

    for (int i{0}; i < 7; ++i) fP[i] = fP[i] + k0[i] * zeta[0] + k1[i] * zeta[1] + k2[i] * zeta[2];

    for (int i{0}, k{0}; i < 7; ++i) {
        for (int j{0}; j <= i; ++j, ++k) {
            fC[k] = fC[k] - (k0[i] * mCHt0[j] + k1[i] * mCHt1[j] + k2[i] * mCHt2[j]);
        }
    }

#if KF_DEBUG
    PrintJoinedMatrix<7>("(AddDaughterWithEnergyFit) mCH", mCHt0, mCHt1, mCHt2);
    PrintJoinedMatrix<7>("(AddDaughterWithEnergyFit) KGain", k0, k1, k2);
    PrintVector<8>("(AddDaughterWithEnergyFit) fP (after Kalman gain)", fP);
    PrintSymMatrix<8>("(AddDaughterWithEnergyFit) fC (after Kalman gain)", fC);
#endif

    Matrix<3, 3> K;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            K[i][j] = 0.;
            for (int k{0}; k < 3; ++k) K[i][j] += meas.C1[IJ(i, k)] * mS[IJ(k, j)];
        }
    }

    Matrix<3, 3> K2;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) K2[i][j] = -K[j][i];
        K2[i][i] += 1.;
    }

    Matrix<3, 3> A;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            A[i][j] = 0.;
            for (int k{0}; k < 3; ++k) A[i][j] += meas.D[i][k] * K2[k][j];
        }
    }

    Matrix<3, 3> M;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            M[i][j] = 0.;
            for (int k{0}; k < 3; ++k) M[i][j] += K[i][k] * A[k][j];
        }
    }

    fC[0] += 2. * M[0][0];
    fC[1] += M[0][1] + M[1][0];
    fC[2] += 2. * M[1][1];
    fC[3] += M[0][2] + M[2][0];
    fC[4] += M[1][2] + M[2][1];
    fC[5] += 2. * M[2][2];
#if KF_DEBUG
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) K", K);
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) K2", K2);
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) A", A);
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) M", M);
    PrintSymMatrix<8>("(AddDaughterWithEnergyFit) fC (the end)", fC);
#endif

    // Calculate Chi^2

    fNDF += 2;
    fQ += daughter.GetQ();
    fChi2 += dChi2;
#if KF_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// Add daughter to the current particle. Uses slower but correct mathematics,
// which requires that the masses of daughter particles stays fixed in the construction process.
// Input:
// - daughter       : the daughter particle
// - bz             : z-component of magnetic field
// - chi2_threshold : do an early cut of chi2
// CONSIDERATION: it will modify the state of the current KF::Particle
void Particle::AddDaughterWithEnergyFitMC(const Particle& daughter, double bz, double chi2_threshold) {
#if KF_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    auto meas = GetMeasurement(daughter, bz);

    SymMatrix<3> mS{meas.C1[0] + meas.C2[0],                           //
                    meas.C1[1] + meas.C2[1], meas.C1[2] + meas.C2[2],  //
                    meas.C1[3] + meas.C2[3], meas.C1[4] + meas.C2[4], meas.C1[5] + meas.C2[5]};
    InvertCholesky3(mS);

    Vector<3> zeta{meas.P2[0] - meas.P1[0], meas.P2[1] - meas.P1[1], meas.P2[2] - meas.P1[2]};
    double dChi2{(mS[0] * zeta[0] + mS[1] * zeta[1] + mS[3] * zeta[2]) * zeta[0] + (mS[1] * zeta[0] + mS[2] * zeta[1] + mS[4] * zeta[2]) * zeta[1] +
                 (mS[3] * zeta[0] + mS[4] * zeta[1] + mS[5] * zeta[2]) * zeta[2]};
    if (dChi2 > chi2_threshold) return;

    Matrix<3, 3> K;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            K[i][j] = 0.;
            for (int k{0}; k < 3; ++k) K[i][j] += meas.C1[IJ(i, k)] * mS[IJ(k, j)];
        }
    }

    // CHt = CH'

    Vector<7> mCHt0{meas.C1[0], meas.C1[1], meas.C1[3], meas.C1[6], meas.C1[10], meas.C1[15], meas.C1[21]};
    Vector<7> mCHt1{meas.C1[1], meas.C1[2], meas.C1[4], meas.C1[7], meas.C1[11], meas.C1[16], meas.C1[22]};
    Vector<7> mCHt2{meas.C1[3], meas.C1[4], meas.C1[5], meas.C1[8], meas.C1[12], meas.C1[17], meas.C1[23]};

    // Kalman gain K = mCH'*S

    Vector<7> k0;
    Vector<7> k1;
    Vector<7> k2;
    for (int i{0}; i < 7; ++i) {
        k0[i] = mCHt0[i] * mS[0] + mCHt1[i] * mS[1] + mCHt2[i] * mS[3];
        k1[i] = mCHt0[i] * mS[1] + mCHt1[i] * mS[2] + mCHt2[i] * mS[4];
        k2[i] = mCHt0[i] * mS[3] + mCHt1[i] * mS[4] + mCHt2[i] * mS[5];
    }

    // last iteration -> update the particle

    // VHt = VH'

    Vector<7> mVHt0{meas.C2[0], meas.C2[1], meas.C2[3], meas.C2[6], meas.C2[10], meas.C2[15], meas.C2[21]};
    Vector<7> mVHt1{meas.C2[1], meas.C2[2], meas.C2[4], meas.C2[7], meas.C2[11], meas.C2[16], meas.C2[22]};
    Vector<7> mVHt2{meas.C2[3], meas.C2[4], meas.C2[5], meas.C2[8], meas.C2[12], meas.C2[17], meas.C2[23]};

    // Kalman gain Km = mCH'*S

    Vector<7> km0;
    Vector<7> km1;
    Vector<7> km2;
    for (int i{0}; i < 7; ++i) {
        km0[i] = mVHt0[i] * mS[0] + mVHt1[i] * mS[1] + mVHt2[i] * mS[3];
        km1[i] = mVHt0[i] * mS[1] + mVHt1[i] * mS[2] + mVHt2[i] * mS[4];
        km2[i] = mVHt0[i] * mS[3] + mVHt1[i] * mS[4] + mVHt2[i] * mS[5];
    }

    for (int i{0}; i < 7; ++i) fP[i] = meas.P1[i] + k0[i] * zeta[0] + k1[i] * zeta[1] + k2[i] * zeta[2];

    for (int i{0}; i < 7; ++i) meas.P2[i] = meas.P2[i] - km0[i] * zeta[0] - km1[i] * zeta[1] - km2[i] * zeta[2];

    for (int i{0}, k{0}; i < 7; ++i) {
        for (int j{0}; j <= i; ++j, ++k) {
            fC[k] = meas.C1[k] - (k0[i] * mCHt0[j] + k1[i] * mCHt1[j] + k2[i] * mCHt2[j]);
        }
    }

    for (int i{0}, k{0}; i < 7; ++i) {
        for (int j{0}; j <= i; ++j, ++k) meas.C2[k] = meas.C2[k] - (km0[i] * mVHt0[j] + km1[i] * mVHt1[j] + km2[i] * mVHt2[j]);
    }

    Matrix<7, 7> mDf;

    for (int i{0}; i < 7; ++i) {
        for (int j{0}; j < 7; ++j) mDf[i][j] = km0[i] * mCHt0[j] + km1[i] * mCHt1[j] + km2[i] * mCHt2[j];
    }

    Matrix<7, 7> mJ1{};
    Matrix<7, 7> mJ2{};

    double mMassParticle{fP[6] * fP[6] - fP[3] * fP[3] - fP[4] * fP[4] - fP[5] * fP[5]};
    double mMassDaughter{meas.P2[6] * meas.P2[6] - meas.P2[3] * meas.P2[3] - meas.P2[4] * meas.P2[4] - meas.P2[5] * meas.P2[5]};
    if (mMassParticle > 0.) mMassParticle = std::sqrt(mMassParticle);
    if (mMassDaughter > 0.) mMassDaughter = std::sqrt(mMassDaughter);

    if (fMassHypo > -0.5)
        SetMassConstraint(fP, fC, mJ1, fMassHypo);
    else if ((mMassParticle < fSumDaughterMass) || (fP[6] < 0.))
        SetMassConstraint(fP, fC, mJ1, fSumDaughterMass);

    if (daughter.fMassHypo > -0.5)
        SetMassConstraint(meas.P2, meas.C2, mJ2, daughter.fMassHypo);
    else if ((mMassDaughter < daughter.fSumDaughterMass) || (meas.P2[6] < 0.))
        SetMassConstraint(meas.P2, meas.C2, mJ2, daughter.fSumDaughterMass);

    Matrix<7, 7> mDJ;
    for (int i{0}; i < 7; ++i) {
        for (int j{0}; j < 7; ++j) {
            mDJ[i][j] = 0.;
            for (int k{0}; k < 7; ++k) mDJ[i][j] += mDf[i][k] * mJ1[j][k];
        }
    }

    for (int i{0}; i < 7; ++i) {
        for (int j{0}; j < 7; ++j) {
            mDf[i][j] = 0.;
            for (int l{0}; l < 7; ++l) mDf[i][j] += mJ2[i][l] * mDJ[l][j];
        }
    }

    // Add the daughter momentum to the particle momentum

    fP[3] += meas.P2[3];
    fP[4] += meas.P2[4];
    fP[5] += meas.P2[5];
    fP[6] += meas.P2[6];

    fC[9] += meas.C2[9];
    fC[13] += meas.C2[13];
    fC[14] += meas.C2[14];
    fC[18] += meas.C2[18];
    fC[19] += meas.C2[19];
    fC[20] += meas.C2[20];
    fC[24] += meas.C2[24];
    fC[25] += meas.C2[25];
    fC[26] += meas.C2[26];
    fC[27] += meas.C2[27];

    fC[6] += mDf[3][0];
    fC[7] += mDf[3][1];
    fC[8] += mDf[3][2];
    fC[10] += mDf[4][0];
    fC[11] += mDf[4][1];
    fC[12] += mDf[4][2];
    fC[15] += mDf[5][0];
    fC[16] += mDf[5][1];
    fC[17] += mDf[5][2];
    fC[21] += mDf[6][0];
    fC[22] += mDf[6][1];
    fC[23] += mDf[6][2];

    fC[9] += mDf[3][3] + mDf[3][3];
    fC[13] += mDf[4][3] + mDf[3][4];
    fC[14] += mDf[4][4] + mDf[4][4];
    fC[18] += mDf[5][3] + mDf[3][5];
    fC[19] += mDf[5][4] + mDf[4][5];
    fC[20] += mDf[5][5] + mDf[5][5];
    fC[24] += mDf[6][3] + mDf[3][6];
    fC[25] += mDf[6][4] + mDf[4][6];
    fC[26] += mDf[6][5] + mDf[5][6];
    fC[27] += mDf[6][6] + mDf[6][6];

    Matrix<3, 3> K2;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) K2[i][j] = -K[j][i];
        K2[i][i] += 1.;
    }

    Matrix<3, 3> A;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            A[i][j] = 0.;
            for (int k{0}; k < 3; ++k) A[i][j] += meas.D[i][k] * K2[k][j];
        }
    }

    Matrix<3, 3> M;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            M[i][j] = 0.;
            for (int k{0}; k < 3; ++k) M[i][j] += K[i][k] * A[k][j];
        }
    }

    fC[0] += 2. * M[0][0];
    fC[1] += M[0][1] + M[1][0];
    fC[2] += 2. * M[1][1];
    fC[3] += M[0][2] + M[2][0];
    fC[4] += M[1][2] + M[2][1];
    fC[5] += 2. * M[2][2];

    // calculate Chi^2

    fNDF += 2;
    fQ += daughter.GetQ();
    fChi2 += dChi2;
#if KF_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// Add a vertex as a point-like measurement to the current particle.
// The eights parameter of the state vector is filled with the decay
// length to the momentum ratio (s = l/p).
// The corresponding covariances are calculated as well.
// The parameters of the particle are stored at the position of the production vertex.
// \param[in] vtx - the assumed production vertex
void Particle::SetProductionVertex(const Particle& vtx, double bz) {
    /*
        Vector<6> m{vtx.XYZPxPyPz()};
        SymMatrix<6> mV{vtx.Cov_6x6()};

        Vector<3> decayPoint{fP[0], fP[1], fP[2]};
        SymMatrix<3> decayPointCov{fC[0], fC[1], fC[2], fC[3], fC[4], fC[5]};

        auto D = Zero<6, 6>();

        bool noS{fC[35] <= 0};  // no decay length allowed

        if (noS) {
            TransportToDecayVertex(bz);
            fP[7] = 0.;
            for (int i{28}; i < 36; ++i) fC[i] = 0.;
        } else {
            Vector<6> dsdr{0., 0., 0., 0., 0., 0.};
            double ds{GetDStoPointBz(bz, vtx.XYZ(), dsdr)};

            Vector<6> dsdp{-dsdr[0], -dsdr[1], -dsdr[2], 0., 0., 0.};

            auto F = Zero<6, 6>();
            auto F1 = Zero<6, 6>();
            TransportBz(bz, ds, dsdr, fP, fC, dsdp, F, F1);

            SymMatrix<6> CTmp{MultQSQt<6>(F1, mV)};

            for (int iC{0}; iC < 6; ++iC) fC[iC] += CTmp[iC];

            for (int i{0}; i < 6; ++i) {
                for (int j{0}; j < 3; ++j) {
                    for (int k{0}; k < 3; ++k) {
                        D[i][j] += mV[IJ(j, k)] * F1[i][k];
                    }
                }
            }
        }

        SymMatrix<3> mS{fC[0] + mV[0],                 //
                        fC[1] + mV[1], fC[2] + mV[2],  //
                        fC[3] + mV[3], fC[4] + mV[4], fC[5] + mV[5]};
        InvertCholesky3(mS);

        Vector<3> res{m[0] - X(), m[1] - Y(), m[2] - Z()};

        Matrix<3, 3> K;
        for (int i{0}; i < 3; ++i) {
            for (int j{0}; j < 3; ++j) {
                K[i][j] = 0.;
                for (int k{0}; k < 3; ++k) K[i][j] += fC[IJ(i, k)] * mS[IJ(k, j)];
            }
        }

        Vector<7> mCHt0{fC[0], fC[1], fC[3], fC[6], fC[10], fC[15], fC[21]};
        Vector<7> mCHt1{fC[1], fC[2], fC[4], fC[7], fC[11], fC[16], fC[22]};
        Vector<7> mCHt2{fC[3], fC[4], fC[5], fC[8], fC[12], fC[17], fC[23]};

        Vector<7> k0;
        Vector<7> k1;
        Vector<7> k2;
        for (int i{0}; i < 7; ++i) {
            k0[i] = mCHt0[i] * mS[0] + mCHt1[i] * mS[1] + mCHt2[i] * mS[3];
            k1[i] = mCHt0[i] * mS[1] + mCHt1[i] * mS[2] + mCHt2[i] * mS[4];
            k2[i] = mCHt0[i] * mS[3] + mCHt1[i] * mS[4] + mCHt2[i] * mS[5];
        }

        for (int i{0}; i < 7; ++i) fP[i] = fP[i] + k0[i] * res[0] + k1[i] * res[1] + k2[i] * res[2];

        for (int i{0}, k{0}; i < 7; ++i) {
            for (int j{0}; j <= i; ++j, ++k) {
                fC[k] = fC[k] - (k0[i] * mCHt0[j] + k1[i] * mCHt1[j] + k2[i] * mCHt2[j]);
            }
        }

        Matrix<3, 3> K2;
        for (int i{0}; i < 3; ++i) {
            for (int j{0}; j < 3; ++j) K2[i][j] = -K[j][i];
            K2[i][i] += 1.;
        }

        Matrix<3, 3> A;
        for (int i{0}; i < 3; ++i) {
            for (int j{0}; j < 3; ++j) {
                A[i][j] = 0.;
                for (int k{0}; k < 3; ++k) A[i][j] += D[k][i] * K2[k][j];
            }
        }

        Matrix<3, 3> M;
        for (int i{0}; i < 3; ++i) {
            for (int j{0}; j < 3; ++j) {
                M[i][j] = 0.;
                for (int k{0}; k < 3; ++k) M[i][j] += K[i][k] * A[k][j];
            }
        }

        fC[0] += 2. * M[0][0];
        fC[1] += M[0][1] + M[1][0];
        fC[2] += 2. * M[1][1];
        fC[3] += M[0][2] + M[2][0];
        fC[4] += M[1][2] + M[2][1];
        fC[5] += 2. * M[2][2];

        fChi2 += (mS[0] * res[0] + mS[1] * res[1] + mS[3] * res[2]) * res[0] + (mS[1] * res[0] + mS[2] * res[1] + mS[4] * res[2]) * res[1] +
                 (mS[3] * res[0] + mS[4] * res[1] + mS[5] * res[2]) * res[2];
        fNDF += 2;

        if (noS) {
            fP[7] = 0.;
            for (int i{28}; i < 36; ++i) fC[i] = 0.;
            fSFromDecay = 0.;
        } else {
            Vector<6> dsdr{0., 0., 0., 0., 0., 0.};
            fP[7] = GetDStoPointBz(bz, decayPoint, dsdr);

            Vector<6> dsdp{-dsdr[0], -dsdr[1], -dsdr[2], 0., 0., 0.};

            auto F = Zero<6, 6>();
            auto F1 = Zero<6, 6>();
            Vector<8> tmpP;
            SymMatrix<8> tmpC;
            TransportBz(bz, fP[7], dsdr, tmpP, tmpC, dsdp, F, F1);

            fC[35] = 0.;
            for (int iDsDr{0}; iDsDr < 6; ++iDsDr) {
                double dsdrC{0.};
                double dsdpV{0.};

                for (int k{0}; k < 6; ++k) dsdrC += dsdr[k] * fC[IJ(k, iDsDr)];  // (-dsdr[k])*fC[k,j]

                fC[iDsDr + 28] = dsdrC;
                fC[35] += dsdrC * dsdr[iDsDr];
                if (iDsDr < 3) {
                    for (int k{0}; k < 3; ++k) dsdpV -= dsdr[k] * decayPointCov[IJ(k, iDsDr)];
                    fC[35] -= dsdpV * dsdr[iDsDr];
                }
            }
            fSFromDecay = -fP[7];
        }

        fAtProductionVertex = true;
        */
}

// Set the exact nonlinear mass constraint on the state vector mP with the covariance matrix mC.
// \param[in,out] mP - the state vector to be modified
// \param[in,out] mC - the corresponding covariance matrix
// \param[in,out] mJ - the Jacobian between initial and modified parameters
// \param[in] mass - the mass to be set on the state vector mP
void Particle::SetMassConstraint(Vector<8>& mP, SymMatrix<8>& mC, Matrix<7, 7>& mJ, double mass) const {

    // Set nonlinear mass constraint (Mass) on the state vector mP with a covariance matrix mC.

    double energy2{mP[6] * mP[6]};
    double p2{mP[3] * mP[3] + mP[4] * mP[4] + mP[5] * mP[5]};
    double mass2{mass * mass};

    double a{energy2 - p2 + 2. * mass2};
    double b{-2. * (energy2 + p2)};
    double c{energy2 - p2 - mass2};

    double lambda{0.};
    if (std::abs(b) > Const::AbsAlmostZero) lambda = -c / b;

    double d{4. * energy2 * p2 - mass2 * (energy2 - p2 - 2. * mass2)};
    if (d >= 0 && std::abs(a) > Const::AbsAlmostZero) lambda = (energy2 + p2 - std::sqrt(d)) / a;

    if (mP[6] < 0)           // If energy < 0 we need a lambda < 0
        lambda = -1000000.;  // Empirical, a better solution should be found

    for (int iIter{0}; iIter < 100; ++iIter) {
        double lambda2{lambda * lambda};
        double lambda4{lambda2 * lambda2};

        double lambda0{lambda};

        double f{-mass2 * lambda4 + a * lambda2 + b * lambda + c};
        double df{-4. * mass2 * lambda2 * lambda + 2. * a * lambda + b};
        if (std::abs(df) < Const::AbsAlmostZero) break;
        lambda -= f / df;
        if (std::abs(lambda0 - lambda) < Const::AbsAlmostZero) break;
    }

    double lpi{1. / (1. + lambda)};
    double lmi{1. / (1. - lambda)};
    double lp2i{lpi * lpi};
    double lm2i{lmi * lmi};

    double lambda2{lambda * lambda};

    double dfl{-4. * mass2 * lambda2 * lambda + 2. * a * lambda + b};
    Vector<4> dfx{-2. * (1. + lambda) * (1. + lambda) * mP[3], -2. * (1. + lambda) * (1. + lambda) * mP[4],
                  -2. * (1. + lambda) * (1. + lambda) * mP[5], 2. * (1. - lambda) * (1. - lambda) * mP[6]};
    Vector<4> dlx{1., 1., 1., 1.};
    if (std::abs(dfl) > 1.e-10) {
        for (int i{0}; i < 4; ++i) dlx[i] = -dfx[i] / dfl;
    }

    Vector<4> dxx{mP[3] * lm2i, mP[4] * lm2i, mP[5] * lm2i, -mP[6] * lp2i};

    // mJ = Zero<7, 7>(); // PENDING
    mJ[0][0] = 1.;
    mJ[1][1] = 1.;
    mJ[2][2] = 1.;

    for (int i{3}; i < 7; ++i) {
        for (int j{3}; j < 7; ++j) mJ[i][j] = dlx[j - 3] * dxx[i - 3];
    }

    for (int i{3}; i < 6; ++i) mJ[i][i] += lmi;
    mJ[6][6] += lpi;

    Matrix<7, 7> mCJ;
    for (int i{0}; i < 7; ++i) {
        for (int j{0}; j < 7; ++j) {
            mCJ[i][j] = 0.;
            for (int k{0}; k < 7; ++k) {
                mCJ[i][j] += mC[IJ(i, k)] * mJ[j][k];
            }
        }
    }

    for (int i{0}; i < 7; ++i) {
        for (int j{0}; j <= i; ++j) {
            mC[IJ(i, j)] = 0.;
            for (int l{0}; l < 7; ++l) {
                mC[IJ(i, j)] += mJ[i][l] * mCJ[l][j];
            }
        }
    }

    mP[3] *= lmi;
    mP[4] *= lmi;
    mP[5] *= lmi;
    mP[6] *= lpi;
}

// Set the exact nonlinear mass constraint on the current particle.
// \param[in] mass - the mass to be set on the particle
void Particle::SetNonlinearMassConstraint(double mass) {

    double px{fP[3]};
    double py{fP[4]};
    double pz{fP[5]};
    double energy{fP[6]};

    double residual{energy * energy - px * px - py * py - pz * pz - mass * mass};
    double dm2{4. * (px * px * fC[9] + py * py * fC[14] + pz * pz * fC[20] + energy * energy * fC[27] +
                     2. * (px * py * fC[13] + pz * (px * fC[18] + py * fC[19]) - energy * (px * fC[24] + py * fC[25] + pz * fC[26])))};
    fChi2 += residual * residual / dm2;
    fNDF += 1;

    Matrix<7, 7> mJ{};
    SetMassConstraint(fP, fC, mJ, mass);
    fMassHypo = mass;
    fSumDaughterMass = mass;
}

// Set linearised mass constraint on the current particle. The constraint can be set with an uncertainty.
// \param[in] Mass - the mass to be set on the state vector mP
// \param[in] SigmaMass - uncertainty of the constraint
void Particle::SetMassConstraint(double mass, double sigma_mass) {

    fMassHypo = mass;
    fSumDaughterMass = mass;

    double m2{mass * mass};                   // measurement, weighted by mass
    double s2{m2 * sigma_mass * sigma_mass};  // sigma^2

    double p2{fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]};
    double e0{std::sqrt(m2 + p2)};

    Vector<8> mH{0., 0., 0., -2 * fP[3], -2 * fP[4], -2 * fP[5], 2 * fP[6], 0.};

    double zeta{e0 * e0 - e0 * fP[6]};
    zeta = m2 - (fP[6] * fP[6] - p2);

    Vector<8> mCHt;
    double s2_est{0.};
    for (int i{0}; i < 8; ++i) {
        for (int j{0}; j < 8; ++j) mCHt[i] += Cij(i, j) * mH[j];
        s2_est += mH[i] * mCHt[i];
    }

    if (s2_est < 1.e-20)
        return;  // calculated mass error is already 0,
                 // the particle can not be constrained on mass

    double w2{1. / (s2 + s2_est)};
    fChi2 += zeta * zeta * w2;
    fNDF += 1;
    for (int i{0}, ii{0}; i < 8; ++i) {
        double ki{mCHt[i] * w2};
        fP[i] += ki * zeta;
        for (int j{0}; j <= i; ++j) fC[++ii] -= ki * mCHt[j];
    }
}

// Return dS = l/p parameter, where
// 1) l - signed distance to the DCA point with the input xyz point;
// 2) p - momentum of the particle;
// assuming the straigth line trajectory. Is used for particles with charge 0 or in case of zero magnetic field.
// Also calculate partial derivatives dsdr of the parameter dS over the state vector of the current particle.
// Input:
// - xyz[3] - point where particle should be transported
// Return:
// - dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
// TO CONSIDER: p2 and (p2*p2) cannot be zero
Result::MinPart2Vtx Particle::MinimizeLinePoint(const Vector<3>& xyz) const {

    double x0{fP[0]};
    double y0{fP[1]};
    double z0{fP[2]};
    double px0{fP[3]};
    double py0{fP[4]};
    double pz0{fP[5]};
    double dx{xyz[0] - x0};
    double dy{xyz[1] - y0};
    double dz{xyz[2] - z0};

    double p2{px0 * px0 + py0 * py0 + pz0 * pz0};
    double a{px0 * dx + py0 * dy + pz0 * dz};

    Result::MinPart2Vtx min;
    min.ds = a / p2;
    min.ds_dr[0] = -px0 / p2;
    min.ds_dr[1] = -py0 / p2;
    min.ds_dr[2] = -pz0 / p2;
    min.ds_dr[3] = (dx * p2 - 2. * px0 * a) / (p2 * p2);
    min.ds_dr[4] = (dy * p2 - 2. * py0 * a) / (p2 * p2);
    min.ds_dr[5] = (dz * p2 - 2. * pz0 * a) / (p2 * p2);
    min.pca[0] = x0 + px0 * min.ds;
    min.pca[1] = y0 + py0 * min.ds;
    min.pca[2] = z0 + pz0 * min.ds;

    return min;
}

// return dS = l/p parameter, where
// 1) l - signed distance to the DCA point with the input xyz point;
// 2) p - momentum of the particle;
// under the assumption of the constant homogeneous field Bz.
// Also calculate partial derivatives dsdr of the parameter dS over the state vector of the current particle.
// \param[in] bz - magnetic field Bz
// \param[in] xyz[3] - point, to which particle should be transported
// \param[out] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
Result::MinPart2Vtx Particle::MinimizeHelixPoint(const Vector<3>& xyz, double bz) const {

    Result::MinPart2Vtx min;

    // 1. find min. distance in XY plane //

    double x{fP[0]};
    double y{fP[1]};
    double z{fP[2]};
    double px{fP[3]};
    double py{fP[4]};
    double pz{fP[5]};

    double bq{bz * fQ * Const::Kappa};
    double pt2{px * px + py * py};

    double dx{xyz[0] - x};
    double dy{xyz[1] - y};
    double dz{xyz[2] - z};
    double a{dx * px + dy * py};

    double abq{bq * a};

    min.ds = std::atan2(abq, pt2 + bq * (dy * px - dx * py)) / bq;

    // 2. add z-component as small correction //

    double bs{bq * min.ds};

    double s{std::sin(bs)};
    double c{std::cos(bs)};

    if (std::abs(bq) < Const::AbsAlmostZero) bq = Const::AbsAlmostZero;
    double bbq{bq * (dx * py - dy * px) - pt2};

    double den{abq * abq + bbq * bbq};
    den = den < Const::AbsAlmostZero ? Const::AbsAlmostZero : den;

    min.ds_dr[0] = (px * bbq - py * abq) / den;
    min.ds_dr[1] = (px * abq + py * bbq) / den;
    min.ds_dr[2] = 0.;
    min.ds_dr[3] = -(dx * bbq + dy * abq + 2. * px * a) / den;
    min.ds_dr[4] = (dx * abq - dy * bbq - 2. * py * a) / den;
    min.ds_dr[5] = 0.;

    double sz{0.};
    double cCoeff{(bbq * c - abq * s) - pz * pz};
    if (std::abs(cCoeff) > Const::AbsAlmostZero) sz = (min.ds * pz - dz) * pz / cCoeff;

    Vector<6> dcdr{-bq * py * c - bbq * s * bq * min.ds_dr[0] + px * bq * s - abq * c * bq * min.ds_dr[0],
                   bq * px * c - bbq * s * bq * min.ds_dr[1] + py * bq * s - abq * c * bq * min.ds_dr[1],
                   0.,
                   (-bq * dy - 2. * px) * c - bbq * s * bq * min.ds_dr[3] - dx * bq * s - abq * c * bq * min.ds_dr[3],
                   (bq * dx - 2. * py) * c - bbq * s * bq * min.ds_dr[4] - dy * bq * s - abq * c * bq * min.ds_dr[4],
                   -2. * pz};

    for (int iP{0}; iP < 6; ++iP) min.ds_dr[iP] += pz * pz / cCoeff * min.ds_dr[iP] - sz / cCoeff * dcdr[iP];
    min.ds_dr[2] += pz / cCoeff;
    min.ds_dr[5] += (2. * pz * min.ds - dz) / cCoeff;

    min.ds += sz;

    /*
    bs = bq * dS;
    s = std::sin(bs);
    c = std::cos(bs);

    double sB = s / bq;
    double cB = (1. - c) / bq;

    double p[5];
    p[0] = x + sB * px + cB * py;
    p[1] = y - cB * px + sB * py;
    p[2] = z + dS * pz;
    p[3] = c * px + s * py;
    p[4] = -s * px + c * py;

    dx = xyz[0] - p[0];
    dy = xyz[1] - p[1];
    dz = xyz[2] - p[2];
    a = dx * p[3] + dy * p[4] + dz * pz;
    abq = bq * a;

    dS += std::atan2(abq, p2 + bq * (dy * p[3] - dx * p[4])) / bq;
    */

    return min;
}

// Calculate dS = l/p parameters for two particles, where
// 1) l - signed distance to the DCA point with the other particle;
// 2) p - momentum of the particle;
// under the assumption of the constant homogeneous field Bz. dS[0] is the transport parameter for the current particle,
// dS[1] - for the particle "p".
// Also calculate partial derivatives dsdr of the parameters dS[0] and dS[1] over the state vectors of the particles:
// 1) dsdr[0][6] = d(dS[0])/d(param1);
// 2) dsdr[1][6] = d(dS[0])/d(param2);
// 3) dsdr[2][6] = d(dS[1])/d(param1);
// 4) dsdr[3][6] = d(dS[1])/d(param2);
// where param1 are parameters of the current particle (if the pointer is not provided it is initialised with fP) and
// param2 are parameters of the second particle "p" (if the pointer is not provided it is initialised with p.fP). Parameters
// param1 and param2 should be either provided both or both set to null pointers.
// Input:
// \param[in] Bz - magnetic field Bz
// \param[in] p - second particle
// Return:
// \param[out] dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
// \param[out] dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
// \param[in] param1 - optional parameter, is used in case if the parameters of the current particles are rotated
// to other coordinate system (see GetDStoParticleBy() function), otherwise fP are used
// \param[in] param2 - optional parameter, is used in case if the parameters of the second particles are rotated
// to other coordinate system (see GetDStoParticleBy() function), otherwise p.fP are used
std::pair<Result::MinPart2Part, Result::MinPart2Part> Particle::MinimizeHelixHelix(const Particle& p, double bz) const {
#if KF_DEBUG
    std::cout << "-- starting (MinimizeHelixHelix) --" << '\n';
#endif

    Result::MinPart2Part min1;
    Result::MinPart2Part min2;

    // 1 -- find points of closest approach (PCAs) in XY plane //

    double bq1{bz * fQ * Const::Kappa};
    double bq2{bz * p.fQ * Const::Kappa};

    bool isStraight1{std::abs(bq1) < Const::AbsAlmostZero};
    bool isStraight2{std::abs(bq2) < Const::AbsAlmostZero};
    if (isStraight1 && isStraight2) return MinimizeLineLine(p);

    double px01{Px()};
    double py01{Py()};
    double pz01{Pz()};
    double px02{p.Px()};
    double py02{p.Py()};
    double pz02{p.Pz()};
    double pt12{px01 * px01 + py01 * py01};
    double pt22{px02 * px02 + py02 * py02};
    double x01{X()};
    double y01{Y()};
    double z01{Z()};
    double x02{p.X()};
    double y02{p.Y()};
    double z02{p.Z()};

    double dx0{x01 - x02};
    double dy0{y01 - y02};
    double dr02{dx0 * dx0 + dy0 * dy0};
    double drp1{dx0 * px01 + dy0 * py01};
    double dxyp1{dx0 * py01 - dy0 * px01};
    double drp2{dx0 * px02 + dy0 * py02};
    double dxyp2{dx0 * py02 - dy0 * px02};
    double p1p2{px01 * px02 + py01 * py02};
    double dp1p2{px01 * py02 - px02 * py01};

    double k11{bq2 * drp1 - dp1p2};
    double k21{bq1 * (bq2 * dxyp1 - p1p2) + bq2 * pt12};
    double k12{bq1 * drp2 - dp1p2};
    double k22{bq2 * (bq1 * dxyp2 + p1p2) - bq1 * pt22};

    double kp{dxyp1 * bq2 - dxyp2 * bq1 - p1p2};
    double kd{dr02 * bq1 * bq2 / 2. + kp};
    double c1{-bq1 * kd - pt12 * bq2};
    double c2{bq2 * kd + pt22 * bq1};

    double d1{std::max(pt12 * pt22 - kd * kd, 0.)};
    d1 = std::sqrt(d1);

    // 1.a -- select solution with minimum distance in 3D  //

    double dx{0.};
    double dy{0.};
    double dz{0.};

    int w_sign{+1};  // winner sign
    double dca_sq{std::numeric_limits<double>::max()};

    for (auto sign : {+1, -1}) {
        // particle 1 //
        Cache tmp1;
        if (!isStraight1) {
            tmp1.theta = std::atan2(bq1 * (k11 * c1 + sign * k21 * d1), sign * bq1 * k11 * d1 * bq1 - k21 * c1);
            std::tie(tmp1.sin, tmp1.cos) = sincos(tmp1.theta);
            tmp1.sB = tmp1.sin / bq1;
            tmp1.cB = (1. - tmp1.cos) / bq1;
            tmp1.ds = tmp1.theta / bq1;

            tmp1.pca[0] = x01 + tmp1.sB * px01 + tmp1.cB * py01;
            tmp1.pca[1] = y01 - tmp1.cB * px01 + tmp1.sB * py01;
            tmp1.pca[2] = z01 + tmp1.ds * pz01;
            tmp1.dir[0] = tmp1.cos * px01 + tmp1.sin * py01;
            tmp1.dir[1] = -tmp1.sin * px01 + tmp1.cos * py01;
            tmp1.dir[2] = pz01;
        } else {
            tmp1.ds = (k11 * c1 + sign * k21 * d1) / (-k21 * c1);
            tmp1.pca[0] = x01 + px01 * tmp1.ds;
            tmp1.pca[1] = y01 + py01 * tmp1.ds;
            tmp1.pca[2] = z01 + pz01 * tmp1.ds;
        }

        // particle 2 //
        Cache tmp2;
        if (!isStraight2) {
            tmp2.theta = std::atan2(bq2 * (k12 * c2 + sign * k22 * d1), sign * bq2 * k12 * d1 * bq2 - k22 * c2);
            std::tie(tmp2.sin, tmp2.cos) = sincos(tmp2.theta);
            tmp2.sB = tmp2.sin / bq2;
            tmp2.cB = (1. - tmp2.cos) / bq2;
            tmp2.ds = tmp2.theta / bq2;

            tmp2.pca[0] = x02 + tmp2.sB * px02 + tmp2.cB * py02;
            tmp2.pca[1] = y02 - tmp2.cB * px02 + tmp2.sB * py02;
            tmp2.pca[2] = z02 + tmp2.ds * pz02;
            tmp2.dir[0] = tmp2.cos * px02 + tmp2.sin * py02;
            tmp2.dir[1] = -tmp2.sin * px02 + tmp2.cos * py02;
            tmp2.dir[2] = pz02;
        } else {
            tmp2.ds = (k12 * c2 + sign * k22 * d1) / (-k22 * c2);
            tmp2.pca[0] = x02 + px02 * tmp2.ds;
            tmp2.pca[1] = y02 + py02 * tmp2.ds;
            tmp2.pca[2] = z02 + pz02 * tmp2.ds;
        }

        Vector<3> tmp_diff{tmp2.pca};
        for (int i{0}; i < 3; ++i) tmp_diff[i] -= tmp1.pca[i];
        double tmp_dca_sq{squaredNorm(tmp_diff)};

        // store //
        if (tmp_dca_sq < dca_sq) {
            static_cast<Cache&>(min1) = tmp1;
            static_cast<Cache&>(min2) = tmp2;

            dx = tmp_diff[0];
            dy = tmp_diff[1];
            dz = tmp_diff[2];

            w_sign = sign;
            dca_sq = tmp_dca_sq;
        }
    }
#if KF_DEBUG
    std::cout << "(MinimizeHelixHelix) min1.ds (no z-correction) = " << min1.ds << '\n';
    PrintVector<3>("(MinimizeHelixHelix) min1.(x,y,z)", min1.pca);
    std::cout << "(MinimizeHelixHelix) min2.ds (no z-correction) = " << min2.ds << '\n';
    PrintVector<3>("(MinimizeHelixHelix) min2.(x,y,z)", min2.pca);
    std::cout << "(MinimizeHelixHelix) dx = " << dx << '\n';
    std::cout << "(MinimizeHelixHelix) dy = " << dy << '\n';
    std::cout << "(MinimizeHelixHelix) dz = " << dz << '\n';
    std::cout << "(MinimizeHelixHelix) w_sign = " << w_sign << '\n';
    PrintValue(__FUNCTION__, "dca", std::sqrt(dca_sq));
#endif

    // 1.b -- handle derivatives //

    Vector<6> dk11dr1{bq2 * px01, bq2 * py01, 0., bq2 * dx0 - py02, bq2 * dy0 + px02, 0.};
    Vector<6> dk11dr2{-bq2 * px01, -bq2 * py01, 0., py01, -px01, 0.};
    Vector<6> dk12dr1{bq1 * px02, bq1 * py02, 0., -py02, px02, 0.};
    Vector<6> dk12dr2{-bq1 * px02, -bq1 * py02, 0., bq1 * dx0 + py01, bq1 * dy0 - px01, 0.};
    Vector<6> dk21dr1{
        bq1 * bq2 * py01, -bq1 * bq2 * px01, 0., 2. * bq2 * px01 + bq1 * (-(bq2 * dy0) - px02), 2. * bq2 * py01 + bq1 * (bq2 * dx0 - py02), 0.};
    Vector<6> dk21dr2{-(bq1 * bq2 * py01), bq1 * bq2 * px01, 0., -(bq1 * px01), -(bq1 * py01), 0.};
    Vector<6> dk22dr1{bq1 * bq2 * py02, -(bq1 * bq2 * px02), 0., bq2 * px02, bq2 * py02, 0.};
    Vector<6> dk22dr2{
        -(bq1 * bq2 * py02), bq1 * bq2 * px02, 0., bq2 * (-(bq1 * dy0) + px01) - 2. * bq1 * px02, bq2 * (bq1 * dx0 + py01) - 2. * bq1 * py02, 0.};

    Vector<6> dkddr1{
        bq1 * bq2 * dx0 + bq2 * py01 - bq1 * py02, bq1 * bq2 * dy0 - bq2 * px01 + bq1 * px02, 0., -bq2 * dy0 - px02, bq2 * dx0 - py02, 0.};
    Vector<6> dkddr2{
        -bq1 * bq2 * dx0 - bq2 * py01 + bq1 * py02, -bq1 * bq2 * dy0 + bq2 * px01 - bq1 * px02, 0., bq1 * dy0 - px01, -bq1 * dx0 - py01, 0.};

    Vector<6> dc1dr1{-(bq1 * (bq1 * bq2 * dx0 + bq2 * py01 - bq1 * py02)), -(bq1 * (bq1 * bq2 * dy0 - bq2 * px01 + bq1 * px02)), 0.,
                     -2. * bq2 * px01 - bq1 * (-bq2 * dy0 - px02),         -2. * bq2 * py01 - bq1 * (bq2 * dx0 - py02),          0.};
    Vector<6> dc1dr2{-bq1 * (-bq1 * bq2 * dx0 - bq2 * py01 + bq1 * py02),
                     -bq1 * (-bq1 * bq2 * dy0 + bq2 * px01 - bq1 * px02),
                     0.,
                     -bq1 * (bq1 * dy0 - px01),
                     -bq1 * (-bq1 * dx0 - py01),
                     0.};
    Vector<6> dc2dr1{bq2 * (bq1 * bq2 * dx0 + bq2 * py01 - bq1 * py02),
                     bq2 * (bq1 * bq2 * dy0 - bq2 * px01 + bq1 * px02),
                     0.,
                     bq2 * (-bq2 * dy0 - px02),
                     bq2 * (bq2 * dx0 - py02),
                     0.};
    Vector<6> dc2dr2{bq2 * (-bq1 * bq2 * dx0 - bq2 * py01 + bq1 * py02), bq2 * (-bq1 * bq2 * dy0 + bq2 * px01 - bq1 * px02), 0.,
                     bq2 * (bq1 * dy0 - px01) + 2 * bq1 * px02,          bq2 * (-(bq1 * dx0) - py01) + 2 * bq1 * py02,       0.};
    Vector<6> dd1dr1{};
    Vector<6> dd1dr2{};
    if (d1 > 0.) {
        for (int i{0}; i < 6; ++i) {
            dd1dr1[i] = -kd * dkddr1[i] / d1;
            dd1dr2[i] = -kd * dkddr2[i] / d1;
        }
        dd1dr1[3] += px01 * pt22 / d1;
        dd1dr1[4] += py01 * pt22 / d1;
        dd1dr2[3] += px02 * pt12 / d1;
        dd1dr2[4] += py02 * pt12 / d1;
    }
#if KF_DEBUG
    PrintVector<6>("(MinimizeHelixHelix) dk11dr1", dk11dr1);
    PrintVector<6>("(MinimizeHelixHelix) dk11dr2", dk11dr2);
    PrintVector<6>("(MinimizeHelixHelix) dk12dr1", dk12dr1);
    PrintVector<6>("(MinimizeHelixHelix) dk12dr2", dk12dr2);
    PrintVector<6>("(MinimizeHelixHelix) dk21dr1", dk21dr1);
    PrintVector<6>("(MinimizeHelixHelix) dk21dr2", dk21dr2);
    PrintVector<6>("(MinimizeHelixHelix) dk22dr1", dk22dr1);
    PrintVector<6>("(MinimizeHelixHelix) dk22dr2", dk22dr2);
    PrintVector<6>("(MinimizeHelixHelix) dkddr1", dkddr1);
    PrintVector<6>("(MinimizeHelixHelix) dkddr2", dkddr2);
    PrintVector<6>("(MinimizeHelixHelix) dc1dr1", dc1dr1);
    PrintVector<6>("(MinimizeHelixHelix) dc1dr2", dc1dr2);
    PrintVector<6>("(MinimizeHelixHelix) dc2dr1", dc2dr1);
    PrintVector<6>("(MinimizeHelixHelix) dc2dr2", dc2dr2);
    PrintVector<6>("(MinimizeHelixHelix) dd1dr1", dd1dr1);
    PrintVector<6>("(MinimizeHelixHelix) dd1dr2", dd1dr2);
#endif

    if (!isStraight1) {
        double a{bq1 * (k11 * c1 + w_sign * k21 * d1)};
        double b{w_sign * bq1 * k11 * d1 * bq1 - k21 * c1};
        double c{b * b + a * a};
        double d{c > 0. ? (1. / bq1 * 1. / c) : 0.};

        for (int iP{0}; iP < 6; ++iP) {
            double dadr1{bq1 * (dk11dr1[iP] * c1 + k11 * dc1dr1[iP] + w_sign * dk21dr1[iP] * d1 + w_sign * k21 * dd1dr1[iP])};
            double dadr2{bq1 * (dk11dr2[iP] * c1 + k11 * dc1dr2[iP] + w_sign * dk21dr2[iP] * d1 + w_sign * k21 * dd1dr2[iP])};
            double dbdr1{w_sign * bq1 * bq1 * (dk11dr1[iP] * d1 + k11 * dd1dr1[iP]) - (dk21dr1[iP] * c1 + k21 * dc1dr1[iP])};
            double dbdr2{w_sign * bq1 * bq1 * (dk11dr2[iP] * d1 + k11 * dd1dr2[iP]) - (dk21dr2[iP] * c1 + k21 * dc1dr2[iP])};

            min1.ds_dr[iP] = d * (dadr1 * b - dbdr1 * a);
            min1.ds_dr1[iP] = d * (dadr2 * b - dbdr2 * a);
        }
    } else {
        double a{k11 * c1 + w_sign * k21 * d1};
        double b{-k21 * c1};

        for (int iP{0}; iP < 6; ++iP) {
            double dadr1{dk11dr1[iP] * c1 + k11 * dc1dr1[iP] + w_sign * dk21dr1[iP] * d1 + w_sign * k21 * dd1dr1[iP]};
            double dadr2{dk11dr2[iP] * c1 + k11 * dc1dr2[iP] + w_sign * dk21dr2[iP] * d1 + w_sign * k21 * dd1dr2[iP]};
            double dbdr1{-dk21dr1[iP] * c1 - k21 * dc1dr1[iP]};
            double dbdr2{-dk21dr2[iP] * c1 - k21 * dc1dr2[iP]};

            min1.ds_dr[iP] = dadr1 / b - dbdr1 * a / (b * b);
            min1.ds_dr1[iP] = dadr2 / b - dbdr2 * a / (b * b);
        }
    }

    if (!isStraight2) {
        double a{bq2 * (k12 * c2 + w_sign * k22 * d1)};
        double b{w_sign * bq2 * k12 * d1 * bq2 - k22 * c2};
        double c{b * b + a * a};
        double d{c > 0. ? (1. / bq2 * 1. / c) : 0.};

        for (int iP{0}; iP < 6; ++iP) {
            double dadr1{bq2 * (dk12dr1[iP] * c2 + k12 * dc2dr1[iP] + w_sign * dk22dr1[iP] * d1 + w_sign * k22 * dd1dr1[iP])};
            double dadr2{bq2 * (dk12dr2[iP] * c2 + k12 * dc2dr2[iP] + w_sign * dk22dr2[iP] * d1 + w_sign * k22 * dd1dr2[iP])};
            double dbdr1{w_sign * bq2 * bq2 * (dk12dr1[iP] * d1 + k12 * dd1dr1[iP]) - (dk22dr1[iP] * c2 + k22 * dc2dr1[iP])};
            double dbdr2{w_sign * bq2 * bq2 * (dk12dr2[iP] * d1 + k12 * dd1dr2[iP]) - (dk22dr2[iP] * c2 + k22 * dc2dr2[iP])};

            min2.ds_dr1[iP] = d * (dadr1 * b - dbdr1 * a);
            min2.ds_dr[iP] = d * (dadr2 * b - dbdr2 * a);
        }
    } else {
        double a{k12 * c2 + w_sign * k22 * d1};
        double b{-k22 * c2};

        for (int iP{0}; iP < 6; ++iP) {
            double dadr1{dk12dr1[iP] * c2 + k12 * dc2dr1[iP] + w_sign * dk22dr1[iP] * d1 + w_sign * k22 * dd1dr1[iP]};
            double dadr2{dk12dr2[iP] * c2 + k12 * dc2dr2[iP] + w_sign * dk22dr2[iP] * d1 + w_sign * k22 * dd1dr2[iP]};
            double dbdr1{-dk22dr1[iP] * c2 - k22 * dc2dr1[iP]};
            double dbdr2{-dk22dr2[iP] * c2 - k22 * dc2dr2[iP]};

            min2.ds_dr1[iP] = (std::abs(b) > Const::AbsAlmostZero ? dadr1 / b : 0.) - (b * b > Const::AbsAlmostZero ? dbdr1 * a / (b * b) : 0.);
            min2.ds_dr[iP] = (std::abs(b) > Const::AbsAlmostZero ? dadr2 / b : 0.) - (b * b > Const::AbsAlmostZero ? dbdr2 * a / (b * b) : 0.);
        }
    }

    double px1{min1.dir[0]};
    double py1{min1.dir[1]};
    double px2{min2.dir[0]};
    double py2{min2.dir[1]};
#if KF_DEBUG
    std::cout << "(MinimizeHelixHelix) px1 = " << px1 << '\n';
    std::cout << "(MinimizeHelixHelix) py1 = " << py1 << '\n';
    std::cout << "(MinimizeHelixHelix) px2 = " << px2 << '\n';
    std::cout << "(MinimizeHelixHelix) py2 = " << py2 << '\n';
#endif

    // 2 -- add z-component as small correction //

    double p12{px1 * px1 + py1 * py1 + pz01 * pz01};
    double p22{px2 * px2 + py2 * py2 + pz02 * pz02};
    double lp1p2{px1 * px2 + py1 * py2 + pz01 * pz02};

    double detp{lp1p2 * lp1p2 - p12 * p22};
    if (std::abs(detp) < Const::AbsAlmostZero || detp * detp < Const::AbsAlmostZero) return {min1, min2};  // protection

    // 2.a -- update derivatives //

    double ldrp1{px1 * dx + py1 * dy + pz01 * dz};
    double ldrp2{px2 * dx + py2 * dy + pz02 * dz};
    double a1{ldrp2 * lp1p2 - ldrp1 * p22};
    double a2{ldrp2 * p12 - ldrp1 * lp1p2};
    double lp1p2_ds0{bq1 * (px2 * py1 - py2 * px1)};
    double lp1p2_ds1{bq2 * (px1 * py2 - py1 * px2)};
    double ldrp1_ds0{-p12 + bq1 * (py1 * dx - px1 * dy)};
    double ldrp1_ds1{lp1p2};
    double ldrp2_ds0{-lp1p2};
    double ldrp2_ds1{p22 + bq2 * (py2 * dx - px2 * dy)};
    double detp_ds0{2. * lp1p2 * lp1p2_ds0};
    double detp_ds1{2. * lp1p2 * lp1p2_ds1};
    double a1_ds0{ldrp2_ds0 * lp1p2 + ldrp2 * lp1p2_ds0 - ldrp1_ds0 * p22};
    double a1_ds1{ldrp2_ds1 * lp1p2 + ldrp2 * lp1p2_ds1 - ldrp1_ds1 * p22};
    double a2_ds0{ldrp2_ds0 * p12 - ldrp1_ds0 * lp1p2 - ldrp1 * lp1p2_ds0};
    double a2_ds1{ldrp2_ds1 * p12 - ldrp1_ds1 * lp1p2 - ldrp1 * lp1p2_ds1};

    double dsl1ds0{a1_ds0 / detp - a1 * detp_ds0 / (detp * detp)};
    double dsl1ds1{a1_ds1 / detp - a1 * detp_ds1 / (detp * detp)};
    double dsl2ds0{a2_ds0 / detp - a2 * detp_ds0 / (detp * detp)};
    double dsl2ds1{a2_ds1 / detp - a2 * detp_ds1 / (detp * detp)};
#if KF_DEBUG
    std::cout << "(MinimizeHelixHelix) a1 = " << a1 << '\n';
    std::cout << "(MinimizeHelixHelix) a2 = " << a2 << '\n';
    std::cout << "(MinimizeHelixHelix) lp1p2_ds0 = " << lp1p2_ds0 << '\n';
    std::cout << "(MinimizeHelixHelix) lp1p2_ds1 = " << lp1p2_ds1 << '\n';
    std::cout << "(MinimizeHelixHelix) ldrp1_ds0 = " << ldrp1_ds0 << '\n';
    std::cout << "(MinimizeHelixHelix) ldrp1_ds1 = " << ldrp1_ds1 << '\n';
    std::cout << "(MinimizeHelixHelix) ldrp2_ds0 = " << ldrp2_ds0 << '\n';
    std::cout << "(MinimizeHelixHelix) ldrp2_ds1 = " << ldrp2_ds1 << '\n';
    std::cout << "(MinimizeHelixHelix) detp_ds0 = " << detp_ds0 << '\n';
    std::cout << "(MinimizeHelixHelix) detp_ds1 = " << detp_ds1 << '\n';
    std::cout << "(MinimizeHelixHelix) a1_ds0 = " << a1_ds0 << '\n';
    std::cout << "(MinimizeHelixHelix) a1_ds1 = " << a1_ds1 << '\n';
    std::cout << "(MinimizeHelixHelix) a2_ds0 = " << a2_ds0 << '\n';
    std::cout << "(MinimizeHelixHelix) a2_ds1 = " << a2_ds1 << '\n';
    std::cout << "(MinimizeHelixHelix) dsl1ds0 = " << dsl1ds0 << '\n';
    std::cout << "(MinimizeHelixHelix) dsl1ds1 = " << dsl1ds1 << '\n';
    std::cout << "(MinimizeHelixHelix) dsl2ds0 = " << dsl2ds0 << '\n';
    std::cout << "(MinimizeHelixHelix) dsl2ds1 = " << dsl2ds1 << '\n';
#endif

    Vector<6> dsldr0{};
    Vector<6> dsldr1{};
    Vector<6> dsldr2{};
    Vector<6> dsldr3{};
    for (int iP{0}; iP < 6; ++iP) {
        dsldr0[iP] = dsl1ds0 * min1.ds_dr[iP] + dsl1ds1 * min2.ds_dr1[iP];
        dsldr1[iP] = dsl1ds0 * min1.ds_dr1[iP] + dsl1ds1 * min2.ds_dr[iP];
        dsldr2[iP] = dsl2ds0 * min1.ds_dr[iP] + dsl2ds1 * min2.ds_dr1[iP];
        dsldr3[iP] = dsl2ds0 * min1.ds_dr1[iP] + dsl2ds1 * min2.ds_dr[iP];
    }

    for (int iP{0}; iP < 6; ++iP) {
        min1.ds_dr[iP] += dsldr0[iP];
        min1.ds_dr1[iP] += dsldr1[iP];
        min2.ds_dr1[iP] += dsldr2[iP];
        min2.ds_dr[iP] += dsldr3[iP];
    }
#if KF_DEBUG
    PrintVector<6>("(MinimizeHelixHelix) min1.ds_dr (after z-correction 1)", min1.ds_dr);
    PrintVector<6>("(MinimizeHelixHelix) min1.ds_dr1 (after z-correction 1)", min1.ds_dr1);
    PrintVector<6>("(MinimizeHelixHelix) min2.ds_dr1 (after z-correction 1)", min2.ds_dr1);
    PrintVector<6>("(MinimizeHelixHelix) min2.ds_dr (after z-correction 1)", min2.ds_dr);
#endif

    Vector<6> lp1p2_dr0{0., 0., 0., min1.cos * px2 - py2 * min1.sin, min1.cos * py2 + px2 * min1.sin, pz02};
    Vector<6> lp1p2_dr1{0., 0., 0., min2.cos * px1 - py1 * min2.sin, min2.cos * py1 + px1 * min2.sin, pz01};
    Vector<6> ldrp1_dr0{-px1,
                        -py1,
                        -pz01,
                        min1.cB * py1 - px1 * min1.sB + min1.cos * dx - min1.sin * dy,
                        -min1.cB * px1 - py1 * min1.sB + min1.sin * dx + min1.cos * dy,
                        -min1.ds * pz01 + dz};
    Vector<6> ldrp1_dr1{px1, py1, pz01, -min2.cB * py1 + px1 * min2.sB, min2.cB * px1 + py1 * min2.sB, min2.ds * pz01};
    Vector<6> ldrp2_dr0{-px2, -py2, -pz02, min1.cB * py2 - px2 * min1.sB, -min1.cB * px2 - py2 * min1.sB, -min1.ds * pz02};
    Vector<6> ldrp2_dr1{px2,
                        py2,
                        pz02,
                        -min2.cB * py2 + px2 * min2.sB + min2.cos * dx - min2.sin * dy,
                        min2.cB * px2 + py2 * min2.sB + min2.sin * dx + min2.cos * dy,
                        dz + min2.ds * pz02};
    Vector<6> p12_dr0{0., 0., 0., 2. * px01, 2. * py01, 2. * pz01};
    Vector<6> p22_dr1{0., 0., 0., 2. * px02, 2. * py02, 2. * pz02};
#if KF_DEBUG
    PrintVector<6>("(MinimizeHelixHelix) lp1p2_dr0", lp1p2_dr0);
    PrintVector<6>("(MinimizeHelixHelix) lp1p2_dr1", lp1p2_dr1);
    PrintVector<6>("(MinimizeHelixHelix) ldrp1_dr0", ldrp1_dr0);
    PrintVector<6>("(MinimizeHelixHelix) ldrp1_dr1", ldrp1_dr1);
    PrintVector<6>("(MinimizeHelixHelix) ldrp2_dr0", ldrp2_dr0);
    PrintVector<6>("(MinimizeHelixHelix) ldrp2_dr1", ldrp2_dr1);
    PrintVector<6>("(MinimizeHelixHelix) p12_dr0", p12_dr0);
    PrintVector<6>("(MinimizeHelixHelix) p22_dr1", p22_dr1);
#endif

    for (int iP{0}; iP < 6; ++iP) {
        double a1_dr0{ldrp2_dr0[iP] * lp1p2 + ldrp2 * lp1p2_dr0[iP] - ldrp1_dr0[iP] * p22};
        double a1_dr1{ldrp2_dr1[iP] * lp1p2 + ldrp2 * lp1p2_dr1[iP] - ldrp1_dr1[iP] * p22 - ldrp1 * p22_dr1[iP]};
        double a2_dr0{ldrp2_dr0[iP] * p12 + ldrp2 * p12_dr0[iP] - ldrp1_dr0[iP] * lp1p2 - ldrp1 * lp1p2_dr0[iP]};
        double a2_dr1{ldrp2_dr1[iP] * p12 - ldrp1_dr1[iP] * lp1p2 - ldrp1 * lp1p2_dr1[iP]};
        double detp_dr0{2. * lp1p2 * lp1p2_dr0[iP] - p12_dr0[iP] * p22};
        double detp_dr1{2. * lp1p2 * lp1p2_dr1[iP] - p12 * p22_dr1[iP]};

        min1.ds_dr[iP] += a1_dr0 / detp - a1 * detp_dr0 / (detp * detp);
        min1.ds_dr1[iP] += a1_dr1 / detp - a1 * detp_dr1 / (detp * detp);
        min2.ds_dr1[iP] += a2_dr0 / detp - a2 * detp_dr0 / (detp * detp);
        min2.ds_dr[iP] += a2_dr1 / detp - a2 * detp_dr1 / (detp * detp);
    }
#if KF_DEBUG
    PrintVector<6>("(MinimizeHelixHelix) min1.ds_dr (after z-correction 2)", min1.ds_dr);
    PrintVector<6>("(MinimizeHelixHelix) min1.ds_dr1 (after z-correction 2)", min1.ds_dr1);
    PrintVector<6>("(MinimizeHelixHelix) min2.ds_dr1 (after z-correction 2)", min2.ds_dr1);
    PrintVector<6>("(MinimizeHelixHelix) min2.ds_dr (after z-correction 2)", min2.ds_dr);
#endif

    // 2.b -- update ds //

    min1.ds += a1 / detp;
    min2.ds += a2 / detp;
#if KF_DEBUG
    std::cout << "(MinimizeHelixHelix) min1.ds (after z-correction) = " << min1.ds << '\n';
    std::cout << "(MinimizeHelixHelix) min2.ds (after z-correction) = " << min2.ds << '\n';
#endif

    // 2.c -- update rest of cache properties //

    min1.theta = bq1 * min1.ds;
    std::tie(min1.sin, min1.cos) = sincos(min1.theta);
    min1.sB = min1.sin / bq1;
    min1.cB = (1. - min1.cos) / bq1;

    min1.pca[0] = x01 + min1.sB * px01 + min1.cB * py01;
    min1.pca[1] = y01 - min1.cB * px01 + min1.sB * py01;
    min1.pca[2] = z01 + min1.ds * pz01;
    min1.dir[0] = min1.cos * px01 + min1.sin * py01;
    min1.dir[1] = -min1.sin * px01 + min1.cos * py01;
    min1.dir[2] = pz01;

    min2.theta = bq2 * min2.ds;
    std::tie(min2.sin, min2.cos) = sincos(min2.theta);
    min2.sB = min2.sin / bq2;
    min2.cB = (1. - min2.cos) / bq2;

    min2.pca[0] = x02 + min2.sB * px02 + min2.cB * py02;
    min2.pca[1] = y02 - min2.cB * px02 + min2.sB * py02;
    min2.pca[2] = z02 + min2.ds * pz02;
    min2.dir[0] = min2.cos * px02 + min2.sin * py02;
    min2.dir[1] = -min2.sin * px02 + min2.cos * py02;
    min2.dir[2] = pz02;
#if KF_DEBUG
    PrintVector<3>("(MinimizeHelixHelix) min1.(x,y,z) (after z-correction)", min1.pca);
    PrintVector<3>("(MinimizeHelixHelix) min2.(x,y,z) (after z-correction)", min2.pca);
    std::cout << "-- finished (MinimizeHelixHelix) --" << '\n';
#endif

    return {min1, min2};
}

// Calculate dS = l/p parameters for two particles, where
// 1) l - signed distance to the DCA point with the other particle;
// 2) p - momentum of the particle;
// under the assumption of the straight line trajectory. Is used for particles with charge 0 or in case of zero magnetic field.
// dS[0] is the transport parameter for the current particle, dS[1] - for the particle "p".
// Also calculate partial derivatives dsdr of the parameters dS[0] and dS[1] over the state vectors of the particles:
// 1) dsdr[0][6] = d(dS[0])/d(param1);
// 2) dsdr[1][6] = d(dS[0])/d(param2);
// 3) dsdr[2][6] = d(dS[1])/d(param1);
// 4) dsdr[3][6] = d(dS[1])/d(param2);
// where param1 are parameters of the current particle fP and
// param2 are parameters of the second particle p.fP.
// Input:
// - p - second particle
// Return:
// - dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
// - dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
std::pair<Result::MinPart2Part, Result::MinPart2Part> Particle::MinimizeLineLine(const Particle& p) const {
#if KF_DEBUG
    std::cout << "-- started (" << __FUNCTION__ << ") --" << '\n';
#endif
    Result::MinPart2Part min1;
    Result::MinPart2Part min2;

    double x01{fP[0]};
    double y01{fP[1]};
    double z01{fP[2]};
    double px01{fP[3]};
    double py01{fP[4]};
    double pz01{fP[5]};

    double x02{p.fP[0]};
    double y02{p.fP[1]};
    double z02{p.fP[2]};
    double px02{p.fP[3]};
    double py02{p.fP[4]};
    double pz02{p.fP[5]};

    double p12{px01 * px01 + py01 * py01 + pz01 * pz01};
    double p22{px02 * px02 + py02 * py02 + pz02 * pz02};
    double p1p2{px01 * px02 + py01 * py02 + pz01 * pz02};

    double drp1{px01 * (x02 - x01) + py01 * (y02 - y01) + pz01 * (z02 - z01)};
    double drp2{px02 * (x02 - x01) + py02 * (y02 - y01) + pz02 * (z02 - z01)};

    double detp{p1p2 * p1p2 - p12 * p22};
    if (std::abs(detp) < Const::AbsAlmostZero) {
        return {min1, min2};  // PENDING
    }

    min1.ds = (drp2 * p1p2 - drp1 * p22) / detp;
    min1.pca[0] = x01 + px01 * min1.ds;
    min1.pca[1] = y01 + py01 * min1.ds;
    min1.pca[2] = z01 + pz01 * min1.ds;
    min2.ds = (drp2 * p12 - drp1 * p1p2) / detp;
    min2.pca[0] = x02 + px02 * min2.ds;
    min2.pca[1] = y02 + py02 * min2.ds;
    min2.pca[2] = z02 + pz02 * min2.ds;

    Vector<6> drp1_dr1{-px01, -py01, -pz01, -x01 + x02, -y01 + y02, -z01 + z02};
    Vector<6> drp1_dr2{px01, py01, pz01, 0., 0., 0.};
    Vector<6> drp2_dr1{-px02, -py02, -pz02, 0., 0., 0.};
    Vector<6> drp2_dr2{px02, py02, pz02, -x01 + x02, -y01 + y02, -z01 + z02};
    Vector<6> dp1p2_dr1{0., 0., 0., px02, py02, pz02};
    Vector<6> dp1p2_dr2{0., 0., 0., px01, py01, pz01};
    Vector<6> dp12_dr1{0., 0., 0., 2. * px01, 2. * py01, 2. * pz01};
    Vector<6> dp12_dr2{0., 0., 0., 0., 0., 0.};
    Vector<6> dp22_dr1{0., 0., 0., 0., 0., 0.};
    Vector<6> dp22_dr2{0., 0., 0., 2. * px02, 2. * py02, 2. * pz02};
    Vector<6> ddetp_dr1{0., 0., 0., -2 * p22 * px01 + 2. * p1p2 * px02, -2 * p22 * py01 + 2. * p1p2 * py02, -2 * p22 * pz01 + 2. * p1p2 * pz02};
    Vector<6> ddetp_dr2{0., 0., 0., 2. * p1p2 * px01 - 2. * p12 * px02, 2. * p1p2 * py01 - 2. * p12 * py02, 2. * p1p2 * pz01 - 2. * p12 * pz02};

    double a1{drp2 * p1p2 - drp1 * p22};
    double a2{drp2 * p12 - drp1 * p1p2};

    for (int i{0}; i < 6; ++i) {
        double da1_dr1{drp2_dr1[i] * p1p2 + drp2 * dp1p2_dr1[i] - drp1_dr1[i] * p22 - drp1 * dp22_dr1[i]};
        double da1_dr2{drp2_dr2[i] * p1p2 + drp2 * dp1p2_dr2[i] - drp1_dr2[i] * p22 - drp1 * dp22_dr2[i]};
        double da2_dr1{drp2_dr1[i] * p12 + drp2 * dp12_dr1[i] - drp1_dr1[i] * p1p2 - drp1 * dp1p2_dr1[i]};
        double da2_dr2{drp2_dr2[i] * p12 + drp2 * dp12_dr2[i] - drp1_dr2[i] * p1p2 - drp1 * dp1p2_dr2[i]};

        min1.ds_dr[i] = da1_dr1 / detp - a1 * ddetp_dr1[i] / (detp * detp);
        min1.ds_dr1[i] = da1_dr2 / detp - a1 * ddetp_dr2[i] / (detp * detp);
        min2.ds_dr1[i] = da2_dr1 / detp - a2 * ddetp_dr1[i] / (detp * detp);
        min2.ds_dr[i] = da2_dr2 / detp - a2 * ddetp_dr2[i] / (detp * detp);
    }
#if KF_DEBUG
    std::cout << "(MinimizeLineLine) min1.ds = " << min1.ds << '\n';
    PrintVector<3>("(MinimizeLineLine) min1.pca", min1.pca);
    PrintVector<6>("(MinimizeLineLine) min1.ds_dr", min1.ds_dr);
    PrintVector<6>("(MinimizeLineLine) min1.ds_dr1", min1.ds_dr1);
    std::cout << "(MinimizeLineLine) min2.ds = " << min2.ds << '\n';
    PrintVector<3>("(MinimizeLineLine) min2.pca", min2.pca);
    PrintVector<6>("(MinimizeLineLine) min2.ds_dr", min2.ds_dr);
    PrintVector<6>("(MinimizeLineLine) min2.ds_dr1", min2.ds_dr1);
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif

    return {min1, min2};
}

// Transport the parameters and their covariance matrix of the current particle assuming constant homogeneous
// magnetic field Bz on the length defined by the transport parameter dS = l/p, where l is the signed distance and p is
// the momentum of the current particle.
// The obtained parameters and covariance matrix are stored to the arrays P and
// C respectively. P and C can be set to the parameters fP and covariance matrix fC of the current particle. In this
// case the particle parameters will be modified. Dependence of the transport parameter dS on the state vector of the
// current particle is taken into account in the covariance matrix using partial derivatives dsdr = d(dS)/d(fP). If
// a pointer to F is initialised the transport jacobian F = d(fP new)/d(fP old) is stored.
// Since dS can depend on the state vector r1 of other particle or vertex, the corelation matrix
// F1 = d(fP new)/d(r1) can be optionally calculated if a pointer F1 is provided.
// Parameters F and F1 should be either both initialised or both set to null pointer.
// \param[in] Bz - z-component of the constant homogeneous magnetic field Bz
// \param[in] dS - transport parameter which defines the distance to which particle should be transported
// \param[in] dsdr[6] = ds/dr - partial derivatives of the parameter dS over the state vector of the current particle
// \param[out] P[8] - array, where transported parameters should be stored
// \param[out] C[36] - array, where transported covariance matrix (8x8) should be stored in the lower triangular form
// \param[in] dsdr1[6] = ds/dr - partial derivatives of the parameter dS over the state vector of another particle
// or vertex
// \param[out] F[36] - transport jacobian, 6x6 matrix F = d(fP new)/d(fP old)
// \param[out] F1[36] - correlation 6x6 matrix between the current particle and particle or vertex
// with the state vector r1, to which the current particle is being transported, F1 = d(fP new)/d(r1)
Result::Transport Particle::TransportBz(const Result::MinPart2Part& min, double bz) const {
#if KF_DEBUG
    std::cout << "-- starting (TransportBz) --" << '\n';
#endif
    Result::Transport tpr;

    double bq{bz * fQ * Const::Kappa};

    double px0{fP[3]};
    double py0{fP[4]};
    double pz0{fP[5]};

    tpr.P[0] = fP[0] + min.sB * px0 + min.cB * py0;
    tpr.P[1] = fP[1] - min.cB * px0 + min.sB * py0;
    tpr.P[2] = fP[2] + min.ds * pz0;
    tpr.P[3] = min.cos * px0 + min.sin * py0;
    tpr.P[4] = -min.sin * px0 + min.cos * py0;
    tpr.P[5] = fP[5];
    tpr.P[6] = fP[6];
    tpr.P[7] = fP[7];

    Matrix<8, 8> mJ{};
    for (int i{0}; i < 8; ++i) mJ[i][i] = 1.;
    mJ[0][3] = min.sB;
    mJ[0][4] = min.cB;
    mJ[1][3] = -min.cB;
    mJ[1][4] = min.sB;
    mJ[2][5] = min.ds;
    mJ[3][3] = min.cos;
    mJ[3][4] = min.sin;
    mJ[4][3] = -min.sin;
    mJ[4][4] = min.cos;

    Matrix<6, 6> mJds{};
    mJds[0][3] = min.cos;
    mJds[0][4] = min.sin;
    mJds[1][3] = -min.sin;
    mJds[1][4] = min.cos;
    mJds[2][5] = 1.;
    mJds[3][3] = -bq * min.sin;
    mJds[3][4] = bq * min.cos;
    mJds[4][3] = -bq * min.cos;
    mJds[4][4] = -bq * min.sin;

    for (int i1{0}; i1 < 6; ++i1) {
        for (int i2{0}; i2 < 6; ++i2) {
            mJ[i1][i2] += mJds[i1][3] * px0 * min.ds_dr[i2] + mJds[i1][4] * py0 * min.ds_dr[i2] + mJds[i1][5] * pz0 * min.ds_dr[i2];
        }
    }

    tpr.C = MultQSQt<8>(mJ, fC);

    for (int i{0}; i < 6; ++i) {
        for (int j{0}; j < 6; ++j) tpr.jacob[i][j] = mJ[i][j];
    }

    for (int i1{0}; i1 < 6; ++i1) {
        for (int i2{0}; i2 < 6; ++i2) {
            tpr.corr[i1][i2] = mJds[i1][3] * px0 * min.ds_dr1[i2] + mJds[i1][4] * py0 * min.ds_dr1[i2] + mJds[i1][5] * pz0 * min.ds_dr1[i2];
        }
    }
#if KF_DEBUG
    PrintVector<8>("(TransportBz) State", tpr.P);
    PrintSymMatrix<8>("(TransportBz) Cov", tpr.C);
    PrintMatrix<6, 6>("(TransportBz) Jacob", tpr.jacob);
    PrintMatrix<6, 6>("(TransportBz) Corr", tpr.corr);
    std::cout << "-- finished (TransportBz) --" << '\n';
#endif

    return tpr;
}

// Transports the parameters and their covariance matrix of the current particle assuming the straight line trajectory
// on the length defined by the transport parameter dS = l/p, where l is the signed distance and p is
// the momentum of the current particle. The obtained parameters and covariance matrix are stored to the arrays P and
// C respectively. P and C can be set to the parameters fP and covariance matrix fC of the current particle. In this
// case the particle parameters will be modified. Dependence of the transport parameter dS on the state vector of the
// current particle is taken into account in the covariance matrix using partial derivatives dsdr = d(dS)/d(fP). If
// a pointer to F is initialised the transport jacobian F = d(fP new)/d(fP old) is stored.
// Since dS can depend on the state vector r1 of other particle or vertex, the corelation matrix
// F1 = d(fP new)/d(r1) can be optionally calculated if a pointer F1 is provided.
// *Parameters F and F1 should be either both initialised or both set to null pointer.
// \param[in] dS - transport parameter which defines the distance to which particle should be transported
// \param[in] dsdr[6] = ds/dr - partial derivatives of the parameter dS over the state vector of the current particle
// \param[out] P[8] - array, where transported parameters should be stored
// \param[out] C[36] - array, where transported covariance matrix (8x8) should be stored in the lower triangular form
// \param[in] dsdr1[6] = ds/dr - partial derivatives of the parameter dS over the state vector of another particle
// or vertex
// \param[out] F[36] - optional parameter, transport jacobian, 6x6 matrix F = d(fP new)/d(fP old)
// \param[out] F1[36] - optional parameter, corelation 6x6 matrix betweeen the current particle and particle or vertex
// with the state vector r1, to which the current particle is being transported, F1 = d(fP new)/d(r1)
Result::Transport Particle::TransportLine(const Result::MinPart2Part& min) const {

    Result::Transport tpr;

    Matrix<8, 8> mJ{};
    mJ[0][0] = 1.;
    mJ[0][3] = min.ds;
    mJ[1][1] = 1.;
    mJ[1][4] = min.ds;
    mJ[2][2] = 1.;
    mJ[2][5] = min.ds;
    mJ[3][3] = 1.;
    mJ[4][4] = 1.;
    mJ[5][5] = 1.;
    mJ[6][6] = 1.;
    mJ[7][7] = 1.;

    double px{fP[3]};
    double py{fP[4]};
    double pz{fP[5]};

    tpr.P[0] = fP[0] + min.ds * fP[3];
    tpr.P[1] = fP[1] + min.ds * fP[4];
    tpr.P[2] = fP[2] + min.ds * fP[5];
    tpr.P[3] = fP[3];
    tpr.P[4] = fP[4];
    tpr.P[5] = fP[5];
    tpr.P[6] = fP[6];
    tpr.P[7] = fP[7];

    Matrix<6, 6> mJds{};
    mJds[0][3] = 1.;
    mJds[1][4] = 1.;
    mJds[2][5] = 1.;

    for (int i1{0}; i1 < 6; ++i1) {
        for (int i2{0}; i2 < 6; ++i2) {
            mJ[i1][i2] += mJds[i1][3] * px * min.ds_dr[i2] + mJds[i1][4] * py * min.ds_dr[i2] + mJds[i1][5] * pz * min.ds_dr[i2];
        }
    }
    tpr.C = MultQSQt<8>(mJ, fC);

    for (int i{0}; i < 6; ++i) {
        for (int j{0}; j < 6; ++j) tpr.jacob[i][j] = mJ[i][j];
    }
    for (int i1{0}; i1 < 6; ++i1) {
        for (int i2{0}; i2 < 6; ++i2) {
            tpr.corr[i1][i2] = mJds[i1][3] * px * min.ds_dr1[i2] + mJds[i1][4] * py * min.ds_dr1[i2] + mJds[i1][5] * pz * min.ds_dr1[i2];
        }
    }

    return tpr;
}

}  // namespace KF
