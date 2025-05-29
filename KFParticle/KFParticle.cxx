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

#include <algorithm>

namespace KF {

// Set Cxx=Cyy=Czz=100. and Css=1.
void Particle::Initialize() {
    fC[0] = 100.;
    fC[2] = 100.;
    fC[5] = 100.;
    fC[35] = 1.;
}

// set the parameters of the particle:
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
// the correlation matrix of two measurements is also calculated. Parameters of the current particle are
// modified by this function, the daughter is not changed, its parameters are stored to the output arrays
// after modifications.
// \param[in] daughter - the daughter particle to be added, stays unchanged
// \param[out] m[8] - the output parameters of the daughter particle at the DCA point
// \param[out] V[36] - the output covariance matrix of the daughter parameters, takes into account the correlation
// \param[out] D[3][3] - the correlation matrix between the current and daughter particles
bool Particle::GetMeasurement(double bz, const Particle& daughter, Vector<8>& m, SymMatrix<8>& V, Matrix<3, 3>& D) {

    if (fNDF == -1) {
        double ds[2]{0., 0.};
        Vector<6> ds_dr[4];
        auto F1 = Zero<6, 6>();
        auto F2 = Zero<6, 6>();
        auto F3 = Zero<6, 6>();
        auto F4 = Zero<6, 6>();
        GetDStoParticleBz(bz, daughter, ds, ds_dr);

        SymMatrix<6> C{Cov_6x6()};
#if KF_DEBUG
        PrintSymMatrix<6>("(GetMeasurement) C", C);
#endif

        TransportBz(bz, ds[0], ds_dr[0], fP, fC, ds_dr[1], F1, F2);
        daughter.TransportBz(bz, ds[1], ds_dr[3], m, V, ds_dr[2], F4, F3);
#if KF_DEBUG
        PrintMatrix<6, 6>("(GetMeasurement) F1", F1);
        PrintMatrix<6, 6>("(GetMeasurement) F2", F2);
        PrintMatrix<6, 6>("(GetMeasurement) F3", F3);
        PrintMatrix<6, 6>("(GetMeasurement) F4", F4);
#endif
        SymMatrix<6> V0Tmp{MultQSQt<6>(F2, daughter.Cov_6x6())};
        SymMatrix<6> V1Tmp{MultQSQt<6>(F3, C)};
#if KF_DEBUG
        PrintSymMatrix<6>("(GetMeasurement) V0Tmp", V0Tmp);
        PrintSymMatrix<6>("(GetMeasurement) V1Tmp", V1Tmp);
#endif

        for (int iC{0}; iC < 21; ++iC) {
            fC[iC] += V0Tmp[iC];
            V[iC] += V1Tmp[iC];
        }

        Matrix<6, 6> C1F1T;
        for (int i{0}; i < 6; ++i) {
            for (int j{0}; j < 6; ++j) {
                C1F1T[i][j] = 0.;
                for (int k{0}; k < 6; ++k) {
                    C1F1T[i][j] += C[IJ(i, k)] * F1[j][k];
                }
            }
        }
        Matrix<6, 6> F3C1F1T;
        for (int i{0}; i < 6; ++i) {
            for (int j{0}; j < 6; ++j) {
                F3C1F1T[i][j] = 0.;
                for (int k{0}; k < 6; ++k) {
                    F3C1F1T[i][j] += F3[i][k] * C1F1T[k][j];
                }
            }
        }
        Matrix<6, 6> C2F2T;
        for (int i{0}; i < 6; ++i) {
            for (int j{0}; j < 6; ++j) {
                C2F2T[i][j] = 0.;
                for (int k{0}; k < 6; ++k) {
                    C2F2T[i][j] += daughter.fC[IJ(i, k)] * F2[j][k];
                }
            }
        }
        for (int i{0}; i < 3; ++i) {
            for (int j{0}; j < 3; ++j) {
                D[i][j] = F3C1F1T[i][j];
                for (int k{0}; k < 6; ++k) {
                    D[i][j] += F4[i][k] * C2F2T[k][j];
                }
            }
        }
#if KF_DEBUG
        PrintSymMatrix<8>("(GetMeasurement) fC", fC);
        PrintSymMatrix<8>("(GetMeasurement) V", V);
        PrintMatrix<6, 6>("(GetMeasurement) C1F1T", C1F1T);
        PrintMatrix<6, 6>("(GetMeasurement) F3C1F1T", F3C1F1T);
        PrintMatrix<6, 6>("(GetMeasurement) C2F2T", C2F2T);
        PrintMatrix<3, 3>("(GetMeasurement) D", D);
#endif
    } else {
        Vector<6> ds_dr{0., 0., 0., 0., 0., 0};
        double ds{daughter.GetDStoPointBz(bz, XYZ(), ds_dr)};

        Vector<6> dsdp{-ds_dr[0], -ds_dr[1], -ds_dr[2], 0., 0., 0.};

        auto F = Zero<6, 6>();
        auto F1 = Zero<6, 6>();
        daughter.TransportBz(bz, ds, ds_dr, m, V, dsdp, F, F1);

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
                D[i][j] = 0.;
                for (int k{0}; k < 3; ++k) {
                    D[i][j] += fC[IJ(j, k)] * F1[i][k];
                }
            }
        }
        V[0] += FVFT[0][0];
        V[1] += FVFT[1][0];
        V[2] += FVFT[1][1];
        V[3] += FVFT[2][0];
        V[4] += FVFT[2][1];
        V[5] += FVFT[2][2];
    }

    return true;
}

// Adds daughter to the current particle. Depending on the selected construction method uses:
// 1) Either simplifyed fast mathematics which consideres momentum and energy as
// independent variables and thus ignores constraint on the fixed mass (fConstructMethod = 0).
// In this case the mass of the daughter particle can be corrupted when the constructed vertex
// is added as the measurement and the mass of the output short-lived particle can become
// unphysical - smaller then the threshold. Implemented in the
// AddDaughterWithEnergyFit() function
// 2) Or slower but correct mathematics which requires that the masses of daughter particles
// stays fixed in the construction process (fConstructMethod = 2). Implemented in the
// AddDaughterWithEnergyFitMC() function.
// \param[in] Daughter - the daughter particle
void Particle::AddDaughter(double bz, const Particle& daughter) {

    if (fNDF < -1) {  // first daughter -> just copy
        fNDF = -1;
        fQ = daughter.GetQ();
        for (int i{0}; i < 7; ++i) fP[i] = daughter.fP[i];
        for (int i{0}; i < 28; ++i) fC[i] = daughter.fC[i];
        fSFromDecay = 0.;
        fMassHypo = daughter.fMassHypo;
        fSumDaughterMass = daughter.fSumDaughterMass;
        return;
    }

    if (fConstructMethod == 0)
        AddDaughterWithEnergyFit(bz, daughter);
    else if (fConstructMethod == 2)
        AddDaughterWithEnergyFitMC(bz, daughter);

    fSumDaughterMass += daughter.fSumDaughterMass;
    fMassHypo = -1.;
}

// Adds daughter to the current particle. Uses simplifyed fast mathematics which consideres momentum
// and energy as independent variables and thus ignores constraint on the fixed mass.
// In this case the mass of the daughter particle can be corrupted when the constructed vertex
// is added as the measurement and the mass of the output short-lived particle can become
// unphysical - smaller then the threshold.
// \param[in] Daughter - the daughter particle
void Particle::AddDaughterWithEnergyFit(double bz, const Particle& daughter) {

    Vector<8> m{Zero<8>()};
    SymMatrix<8> mV{Zero<36>()};

    auto D = Zero<3, 3>();
    if (!GetMeasurement(bz, daughter, m, mV, D)) return;

    SymMatrix<3> mS{fC[0] + mV[0],                 //
                    fC[1] + mV[1], fC[2] + mV[2],  //
                    fC[3] + mV[3], fC[4] + mV[4], fC[5] + mV[5]};

    InvertCholesky3(mS);
#if KF_DEBUG
    PrintVector<8>("(AddDaughterWithEnergyFit) m", m);
    PrintSymMatrix<8>("(AddDaughterWithEnergyFit) mV", mV);
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) D", D);
    PrintSymMatrix<3>("(AddDaughterWithEnergyFit) mS", mS);
#endif
    // Residual (measured - estimated)

    Vector<3> zeta{m[0] - fP[0], m[1] - fP[1], m[2] - fP[2]};

    double dChi2{(mS[0] * zeta[0] + mS[1] * zeta[1] + mS[3] * zeta[2]) * zeta[0] + (mS[1] * zeta[0] + mS[2] * zeta[1] + mS[4] * zeta[2]) * zeta[1] +
                 (mS[3] * zeta[0] + mS[4] * zeta[1] + mS[5] * zeta[2]) * zeta[2]};
    // if (dChi2 > 1e9) return;  // PENDING: consider chi2 threshold

    Matrix<3, 3> K;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            K[i][j] = 0.;
            for (int k{0}; k < 3; ++k) K[i][j] += fC[IJ(i, k)] * mS[IJ(k, j)];
        }
    }
#if KF_DEBUG
    PrintVector<3>("(AddDaughterWithEnergyFit) zeta", zeta);
    std::cout << "(AddDaughterWithEnergyFit) dChi2 = " << dChi2 << '\n';
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) K", K);
#endif

    // CHt = CH' - D'

    Vector<7> mCHt0{fC[0], fC[1], fC[3], fC[6] - mV[6], fC[10] - mV[10], fC[15] - mV[15], fC[21] - mV[21]};
    Vector<7> mCHt1{fC[1], fC[2], fC[4], fC[7] - mV[7], fC[11] - mV[11], fC[16] - mV[16], fC[22] - mV[22]};
    Vector<7> mCHt2{fC[3], fC[4], fC[5], fC[8] - mV[8], fC[12] - mV[12], fC[17] - mV[17], fC[23] - mV[23]};

    // Kalman gain K = mCH'*S

    Vector<7> k0;
    Vector<7> k1;
    Vector<7> k2;
    for (int i{0}; i < 7; ++i) {
        k0[i] = mCHt0[i] * mS[0] + mCHt1[i] * mS[1] + mCHt2[i] * mS[3];
        k1[i] = mCHt0[i] * mS[1] + mCHt1[i] * mS[2] + mCHt2[i] * mS[4];
        k2[i] = mCHt0[i] * mS[3] + mCHt1[i] * mS[4] + mCHt2[i] * mS[5];
    }
#if KF_DEBUG
    PrintJoinedMatrix<7>("(AddDaughterWithEnergyFit) mCH", mCHt0, mCHt1, mCHt2);
    PrintJoinedMatrix<7>("(AddDaughterWithEnergyFit) KGain", k0, k1, k2);
#endif

    // Add the daughter momentum to the particle momentum

    fP[3] += m[3];
    fP[4] += m[4];
    fP[5] += m[5];
    fP[6] += m[6];

    fC[9] += mV[9];
    fC[13] += mV[13];
    fC[14] += mV[14];
    fC[18] += mV[18];
    fC[19] += mV[19];
    fC[20] += mV[20];
    fC[24] += mV[24];
    fC[25] += mV[25];
    fC[26] += mV[26];
    fC[27] += mV[27];

#if KF_DEBUG
    PrintVector<8>("(AddDaughterWithEnergyFit) fP", fP);
    PrintSymMatrix<8>("(AddDaughterWithEnergyFit) fC", fC);
#endif

    // New estimation of the vertex position r += K*zeta

    for (int i{0}; i < 7; ++i) fP[i] = fP[i] + k0[i] * zeta[0] + k1[i] * zeta[1] + k2[i] * zeta[2];

    // New covariance matrix C -= K*(mCH')'

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
            for (int k{0}; k < 3; ++k) A[i][j] += D[i][k] * K2[k][j];
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
    PrintVector<8>("(AddDaughterWithEnergyFit) fP", fP);
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) K2", K2);
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) A", A);
    PrintMatrix<3, 3>("(AddDaughterWithEnergyFit) M", M);
    PrintSymMatrix<8>("(AddDaughterWithEnergyFit) fC", fC);
    std::cout << "-- finished (AddDaughterWithEnergyFit) --" << '\n';
#endif

    // Calculate Chi^2

    fNDF += 2;
    fQ += daughter.GetQ();
    fSFromDecay = 0.;
    fChi2 += dChi2;
}

// Add daughter to the current particle. Uses slower but correct mathematics
// which requires that the masses of daughter particles
// stays fixed in the construction process.
// \param[in] daughter - the daughter particle
void Particle::AddDaughterWithEnergyFitMC(double bz, const Particle& daughter) {

    Vector<8> m;
    SymMatrix<8> mV;

    Matrix<3, 3> D;
    GetMeasurement(bz, daughter, m, mV, D);

    SymMatrix<3> mS{fC[0] + mV[0],                 //
                    fC[1] + mV[1], fC[2] + mV[2],  //
                    fC[3] + mV[3], fC[4] + mV[4], fC[5] + mV[5]};
    InvertCholesky3(mS);

    // Residual (measured - estimated)
    Vector<3> zeta{m[0] - fP[0], m[1] - fP[1], m[2] - fP[2]};

    Matrix<3, 3> K;
    for (int i{0}; i < 3; ++i) {
        for (int j{0}; j < 3; ++j) {
            K[i][j] = 0.;
            for (int k{0}; k < 3; ++k) K[i][j] += fC[IJ(i, k)] * mS[IJ(k, j)];
        }
    }

    // CHt = CH'

    Vector<7> mCHt0{fC[0], fC[1], fC[3], fC[6], fC[10], fC[15], fC[21]};
    Vector<7> mCHt1{fC[1], fC[2], fC[4], fC[7], fC[11], fC[16], fC[22]};
    Vector<7> mCHt2{fC[3], fC[4], fC[5], fC[8], fC[12], fC[17], fC[23]};

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

    Vector<7> mVHt0{mV[0], mV[1], mV[3], mV[6], mV[10], mV[15], mV[21]};
    Vector<7> mVHt1{mV[1], mV[2], mV[4], mV[7], mV[11], mV[16], mV[22]};
    Vector<7> mVHt2{mV[3], mV[4], mV[5], mV[8], mV[12], mV[17], mV[23]};

    // Kalman gain Km = mCH'*S

    Vector<7> km0;
    Vector<7> km1;
    Vector<7> km2;
    for (int i{0}; i < 7; ++i) {
        km0[i] = mVHt0[i] * mS[0] + mVHt1[i] * mS[1] + mVHt2[i] * mS[3];
        km1[i] = mVHt0[i] * mS[1] + mVHt1[i] * mS[2] + mVHt2[i] * mS[4];
        km2[i] = mVHt0[i] * mS[3] + mVHt1[i] * mS[4] + mVHt2[i] * mS[5];
    }

    for (int i{0}; i < 7; ++i) fP[i] = fP[i] + k0[i] * zeta[0] + k1[i] * zeta[1] + k2[i] * zeta[2];

    for (int i{0}; i < 7; ++i) m[i] = m[i] - km0[i] * zeta[0] - km1[i] * zeta[1] - km2[i] * zeta[2];

    for (int i{0}, k{0}; i < 7; ++i) {
        for (int j{0}; j <= i; ++j, ++k) {
            fC[k] = fC[k] - (k0[i] * mCHt0[j] + k1[i] * mCHt1[j] + k2[i] * mCHt2[j]);
        }
    }

    for (int i{0}, k{0}; i < 7; ++i) {
        for (int j{0}; j <= i; ++j, ++k) mV[k] = mV[k] - (km0[i] * mVHt0[j] + km1[i] * mVHt1[j] + km2[i] * mVHt2[j]);
    }

    Matrix<7, 7> mDf;

    for (int i{0}; i < 7; ++i) {
        for (int j{0}; j < 7; ++j) mDf[i][j] = km0[i] * mCHt0[j] + km1[i] * mCHt1[j] + km2[i] * mCHt2[j];
    }

    auto mJ1 = Zero<7, 7>();
    auto mJ2 = Zero<7, 7>();

    double mMassParticle{fP[6] * fP[6] - (fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5])};
    double mMassDaughter{m[6] * m[6] - (m[3] * m[3] + m[4] * m[4] + m[5] * m[5])};
    if (mMassParticle > 0.) mMassParticle = std::sqrt(mMassParticle);
    if (mMassDaughter > 0.) mMassDaughter = std::sqrt(mMassDaughter);

    if (fMassHypo > -0.5)
        SetMassConstraint(fP, fC, mJ1, fMassHypo);
    else if ((mMassParticle < fSumDaughterMass) || (fP[6] < 0.))
        SetMassConstraint(fP, fC, mJ1, fSumDaughterMass);

    if (daughter.fMassHypo > -0.5)
        SetMassConstraint(m, mV, mJ2, daughter.fMassHypo);
    else if ((mMassDaughter < daughter.fSumDaughterMass) || (m[6] < 0.))
        SetMassConstraint(m, mV, mJ2, daughter.fSumDaughterMass);

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

    fP[3] += m[3];
    fP[4] += m[4];
    fP[5] += m[5];
    fP[6] += m[6];

    fC[9] += mV[9];
    fC[13] += mV[13];
    fC[14] += mV[14];
    fC[18] += mV[18];
    fC[19] += mV[19];
    fC[20] += mV[20];
    fC[24] += mV[24];
    fC[25] += mV[25];
    fC[26] += mV[26];
    fC[27] += mV[27];

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
            for (int k{0}; k < 3; ++k) A[i][j] += D[i][k] * K2[k][j];
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
    fSFromDecay = 0.;
    fChi2 += (mS[0] * zeta[0] + mS[1] * zeta[1] + mS[3] * zeta[2]) * zeta[0] + (mS[1] * zeta[0] + mS[2] * zeta[1] + mS[4] * zeta[2]) * zeta[1] +
             (mS[3] * zeta[0] + mS[4] * zeta[1] + mS[5] * zeta[2]) * zeta[2];
}

// Adds a vertex as a point-like measurement to the current particle.
// The eights parameter of the state vector is filled with the decay
// length to the momentum ratio (s = l/p).
// The corresponding covariances are calculated as well.
// The parameters of the particle are stored at the position of the production vertex.
// \param[in] vtx - the assumed production vertex
void Particle::SetProductionVertex(const Particle& vtx, double bz) {

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

    mJ = Zero<7, 7>();
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

    Matrix<7, 7> mJ;
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

    auto mCHt = Zero<8>();
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

// set constraint on the zero decay length. When the production point is set
// the measurement from this particle is created at the decay point.
void Particle::SetNoDecayLength(double bz) {

    TransportToDecayVertex(bz);

    Vector<8> h{0., 0., 0., 0., 0., 0., 0., 1.};

    double zeta{-fP[7]};
    for (int i{0}; i < 8; ++i) zeta -= h[i] * (fP[i] - fP[i]);

    double s{fC[35]};
    if (s > 1.e-20) {
        s = 1. / s;
        fChi2 += zeta * zeta * s;
        fNDF += 1;
        for (int i{0}, ii{0}; i < 7; ++i) {
            double ki{fC[28 + i] * s};
            fP[i] += ki * zeta;
            for (int j{0}; j <= i; ++j) fC[++ii] -= ki * fC[28 + j];
        }
    }
    fP[7] = 0.;
    for (int i{28}; i < 36; ++i) fC[i] = 0.;
}

// Constructs a short-lived particle from a set of daughter particles:
// 1) all parameters of the "this" objects are initialised;
// 2) daughters are added one after another;
// 3) if Parent pointer is not null, the production vertex is set to it;
// 4) if Mass hypothesis >=0 the mass constraint is set.
// \param[in] v_daughters - array of daughter particles
// \param[in] n_daughters - number of daughter particles in the input array
// \param[in] parent - optional parrent particle
// \param[in] mass - optional mass hypothesis
void Particle::Construct(double bz, const Particle* v_daughters[], int n_daughters, const Particle* parent, double mass) {

    fAtProductionVertex = false;
    fSFromDecay = 0;
    fSumDaughterMass = 0;

    for (int i{0}; i < 36; ++i) fC[i] = 0.;
    fC[35] = 1.;

    fNDF = -3;
    fChi2 = 0.;
    fQ = 0;

    for (int itr{0}; itr < n_daughters; ++itr) {
        AddDaughter(bz, *v_daughters[itr]);
    }

    if (mass >= 0) SetMassConstraint(mass);
    if (parent) SetProductionVertex(*parent, bz);
}

// Transports the particle to its decay vertex
void Particle::TransportToDecayVertex(double bz) {
    Vector<6> dsdr{0., 0., 0., 0., 0., 0.};
    if (fSFromDecay != 0) TransportToDS(bz, -fSFromDecay, dsdr);
    fAtProductionVertex = false;
}

// Transports the particle to its production vertex
void Particle::TransportToProductionVertex(double bz) {
    Vector<6> dsdr{0., 0., 0., 0., 0., 0.};
    if (fSFromDecay != -fP[7]) TransportToDS(bz, -fSFromDecay - fP[7], dsdr);
    fAtProductionVertex = true;
}

// Transport the particle on a certain distance. The distance is defined by the dS=l/p parameter, where
// 1) l - signed distance
// 2) p - momentum of the particle
// \param[in] dS = l/p - distance normalised to the momentum of the particle to be transported on
// \param[in] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
void Particle::TransportToDS(double bz, double ds, const Vector<6>& dsdr) {
    Vector<6> dummy_dsdr;
    Matrix<6, 6> dummy_jacob;
    Matrix<6, 6> dummy_corr;
    TransportBz(bz, ds, dsdr, fP, fC, dummy_dsdr, dummy_jacob, dummy_corr);
    fSFromDecay += ds;
}

// Return dS = l/p parameter, where
// 1) l - signed distance to the DCA point with the input xyz point;
// 2) p - momentum of the particle;
// assuming the straigth line trajectory. Is used for particles with charge 0 or in case of zero magnetic field.
// Also calculate partial derivatives dsdr of the parameter dS over the state vector of the current particle.
// \param[in] xyz[3] - point where particle should be transported
// \param[out] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
double Particle::GetDStoPointLine(const Vector<3>& xyz, Vector<6>& ds_dr) const {

    double p2{fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]};

    double a{fP[3] * (xyz[0] - fP[0]) + fP[4] * (xyz[1] - fP[1]) + fP[5] * (xyz[2] - fP[2])};
    ds_dr[0] = -fP[3] / p2;
    ds_dr[1] = -fP[4] / p2;
    ds_dr[2] = -fP[5] / p2;
    ds_dr[3] = ((xyz[0] - fP[0]) * p2 - 2. * fP[3] * a) / (p2 * p2);
    ds_dr[4] = ((xyz[1] - fP[1]) * p2 - 2. * fP[4] * a) / (p2 * p2);
    ds_dr[5] = ((xyz[2] - fP[2]) * p2 - 2. * fP[5] * a) / (p2 * p2);

    return a / p2;
}

// return dS = l/p parameter, where
// 1) l - signed distance to the DCA point with the input xyz point;
// 2) p - momentum of the particle;
// under the assumption of the constant homogeneous field Bz.
// Also calculate partial derivatives dsdr of the parameter dS over the state vector of the current particle.
// \param[in] bz - magnetic field Bz
// \param[in] xyz[3] - point, to which particle should be transported
// \param[out] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
double Particle::GetDStoPointBz(double bz, const Vector<3>& xyz, Vector<6>& ds_dr) const {

    double ds{0.};

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

    ds = std::atan2(abq, pt2 + bq * (dy * px - dx * py)) / bq;

    // 2. add z-component as small correction //

    double bs{bq * ds};

    double s{std::sin(bs)};
    double c{std::cos(bs)};

    if (std::abs(bq) < Const::AbsAlmostZero) bq = Const::AbsAlmostZero;
    double bbq{bq * (dx * py - dy * px) - pt2};

    double den{abq * abq + bbq * bbq};
    den = den < Const::AbsAlmostZero ? Const::AbsAlmostZero : den;

    ds_dr[0] = (px * bbq - py * abq) / den;
    ds_dr[1] = (px * abq + py * bbq) / den;
    ds_dr[2] = 0.;
    ds_dr[3] = -(dx * bbq + dy * abq + 2. * px * a) / den;
    ds_dr[4] = (dx * abq - dy * bbq - 2. * py * a) / den;
    ds_dr[5] = 0.;

    double sz{0.};
    double cCoeff{(bbq * c - abq * s) - pz * pz};
    if (std::abs(cCoeff) > Const::AbsAlmostZero) sz = (ds * pz - dz) * pz / cCoeff;

    Vector<6> dcdr{-bq * py * c - bbq * s * bq * ds_dr[0] + px * bq * s - abq * c * bq * ds_dr[0],
                   bq * px * c - bbq * s * bq * ds_dr[1] + py * bq * s - abq * c * bq * ds_dr[1],
                   0.,
                   (-bq * dy - 2. * px) * c - bbq * s * bq * ds_dr[3] - dx * bq * s - abq * c * bq * ds_dr[3],
                   (bq * dx - 2. * py) * c - bbq * s * bq * ds_dr[4] - dy * bq * s - abq * c * bq * ds_dr[4],
                   -2. * pz};

    for (int iP{0}; iP < 6; ++iP) ds_dr[iP] += pz * pz / cCoeff * ds_dr[iP] - sz / cCoeff * dcdr[iP];
    ds_dr[2] += pz / cCoeff;
    ds_dr[5] += (2. * pz * ds - dz) / cCoeff;

    ds += sz;

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

    return ds;
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
// \param[in] Bz - magnetic field Bz
// \param[in] p - second particle
// \param[out] dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
// \param[out] dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
// \param[in] param1 - optional parameter, is used in case if the parameters of the current particles are rotated
// to other coordinate system (see GetDStoParticleBy() function), otherwise fP are used
// \param[in] param2 - optional parameter, is used in case if the parameters of the second particles are rotated
// to other coordinate system (see GetDStoParticleBy() function), otherwise p.fP are used
void Particle::GetDStoParticleBz(double bz, const Particle& p, double ds[2], Vector<6> ds_dr[4]) const {

    // 1. find minimum distance in XY plane //

    double bq1{bz * fQ * Const::Kappa};
    double bq2{bz * p.fQ * Const::Kappa};

    bool isStraight1{std::abs(bq1) < Const::AbsAlmostZero};
    bool isStraight2{std::abs(bq2) < Const::AbsAlmostZero};

    if (isStraight1 && isStraight2) {
        GetDStoParticleLine(p, ds, ds_dr);
        return;
    }

    double px1{fP[3]};
    double py1{fP[4]};
    double pz1{fP[5]};

    double px2{p.fP[3]};
    double py2{p.fP[4]};
    double pz2{p.fP[5]};

    double pt12{px1 * px1 + py1 * py1};
    double pt22{px2 * px2 + py2 * py2};

    double x01{fP[0]};
    double y01{fP[1]};
    double z01{fP[2]};

    double x02{p.fP[0]};
    double y02{p.fP[1]};
    double z02{p.fP[2]};

    double ds1[2]{0., 0.};
    double ds2[2]{0., 0.};

    double dx0{x01 - x02};
    double dy0{y01 - y02};
    double dr02{dx0 * dx0 + dy0 * dy0};
    double drp1{dx0 * px1 + dy0 * py1};
    double dxyp1{dx0 * py1 - dy0 * px1};
    double drp2{dx0 * px2 + dy0 * py2};
    double dxyp2{dx0 * py2 - dy0 * px2};
    double p1p2{px1 * px2 + py1 * py2};
    double dp1p2{px1 * py2 - px2 * py1};

    double k11{bq2 * drp1 - dp1p2};
    double k21{bq1 * (bq2 * dxyp1 - p1p2) + bq2 * pt12};
    double k12{bq1 * drp2 - dp1p2};
    double k22{bq2 * (bq1 * dxyp2 + p1p2) - bq1 * pt22};

    double kp{dxyp1 * bq2 - dxyp2 * bq1 - p1p2};
    double kd{dr02 * bq1 * bq2 / 2. + kp};
    double c1{-(bq1 * kd + pt12 * bq2)};
    double c2{bq2 * kd + pt22 * bq1};

    double d1{std::max(pt12 * pt22 - kd * kd, 0.)};
    d1 = std::sqrt(d1);

    double ds1dR1[2][6];
    double ds2dR2[2][6];

    double ds1dR2[2][6];
    double ds2dR1[2][6];

    Vector<6> dk11dr1{bq2 * px1, bq2 * py1, 0., bq2 * dx0 - py2, bq2 * dy0 + px2, 0.};
    Vector<6> dk11dr2{-bq2 * px1, -bq2 * py1, 0., py1, -px1, 0.};
    Vector<6> dk12dr1{bq1 * px2, bq1 * py2, 0., -py2, px2, 0.};
    Vector<6> dk12dr2{-bq1 * px2, -bq1 * py2, 0., bq1 * dx0 + py1, bq1 * dy0 - px1, 0.};
    Vector<6> dk21dr1{
        bq1 * bq2 * py1, -bq1 * bq2 * px1, 0., 2. * bq2 * px1 + bq1 * (-(bq2 * dy0) - px2), 2. * bq2 * py1 + bq1 * (bq2 * dx0 - py2), 0.};
    Vector<6> dk21dr2{-(bq1 * bq2 * py1), bq1 * bq2 * px1, 0., -(bq1 * px1), -(bq1 * py1), 0.};
    Vector<6> dk22dr1{bq1 * bq2 * py2, -(bq1 * bq2 * px2), 0., bq2 * px2, bq2 * py2, 0.};
    Vector<6> dk22dr2{
        -(bq1 * bq2 * py2), bq1 * bq2 * px2, 0., bq2 * (-(bq1 * dy0) + px1) - 2. * bq1 * px2, bq2 * (bq1 * dx0 + py1) - 2. * bq1 * py2, 0.};

    Vector<6> dkddr1{bq1 * bq2 * dx0 + bq2 * py1 - bq1 * py2, bq1 * bq2 * dy0 - bq2 * px1 + bq1 * px2, 0., -bq2 * dy0 - px2, bq2 * dx0 - py2, 0.};
    Vector<6> dkddr2{-bq1 * bq2 * dx0 - bq2 * py1 + bq1 * py2, -bq1 * bq2 * dy0 + bq2 * px1 - bq1 * px2, 0., bq1 * dy0 - px1, -bq1 * dx0 - py1, 0.};

    Vector<6> dc1dr1{-(bq1 * (bq1 * bq2 * dx0 + bq2 * py1 - bq1 * py2)), -(bq1 * (bq1 * bq2 * dy0 - bq2 * px1 + bq1 * px2)), 0.,
                     -2. * bq2 * px1 - bq1 * (-(bq2 * dy0) - px2),       -2. * bq2 * py1 - bq1 * (bq2 * dx0 - py2),          0.};
    Vector<6> dc1dr2{-(bq1 * (-(bq1 * bq2 * dx0) - bq2 * py1 + bq1 * py2)),
                     -(bq1 * (-(bq1 * bq2 * dy0) + bq2 * px1 - bq1 * px2)),
                     0.,
                     -(bq1 * (bq1 * dy0 - px1)),
                     -(bq1 * (-(bq1 * dx0) - py1)),
                     0.};

    Vector<6> dc2dr1{bq2 * (bq1 * bq2 * dx0 + bq2 * py1 - bq1 * py2),
                     bq2 * (bq1 * bq2 * dy0 - bq2 * px1 + bq1 * px2),
                     0.,
                     bq2 * (-(bq2 * dy0) - px2),
                     bq2 * (bq2 * dx0 - py2),
                     0.};
    Vector<6> dc2dr2{bq2 * (-(bq1 * bq2 * dx0) - bq2 * py1 + bq1 * py2), bq2 * (-(bq1 * bq2 * dy0) + bq2 * px1 - bq1 * px2), 0.,
                     bq2 * (bq1 * dy0 - px1) + 2. * bq1 * px2,           bq2 * (-(bq1 * dx0) - py1) + 2. * bq1 * py2,        0.};
    Vector<6> dd1dr1{0., 0., 0., 0., 0., 0.};
    Vector<6> dd1dr2{0., 0., 0., 0., 0., 0.};
    if (d1 > 0) {
        for (int i{0}; i < 6; ++i) {
            dd1dr1[i] = -kd / d1 * dkddr1[i];
            dd1dr2[i] = -kd / d1 * dkddr2[i];
        }
        dd1dr1[3] += px1 / d1 * pt22;
        dd1dr1[4] += py1 / d1 * pt22;
        dd1dr2[3] += px2 / d1 * pt12;
        dd1dr2[4] += py2 / d1 * pt12;
    }
#if KF_DEBUG
    PrintVector<6>("(GetDStoParticleBz) dk11dr1", dk11dr1);
    PrintVector<6>("(GetDStoParticleBz) dk11dr2", dk11dr2);
    PrintVector<6>("(GetDStoParticleBz) dk12dr1", dk12dr1);
    PrintVector<6>("(GetDStoParticleBz) dk12dr2", dk12dr2);
    PrintVector<6>("(GetDStoParticleBz) dk21dr1", dk21dr1);
    PrintVector<6>("(GetDStoParticleBz) dk21dr2", dk21dr2);
    PrintVector<6>("(GetDStoParticleBz) dk22dr1", dk22dr1);
    PrintVector<6>("(GetDStoParticleBz) dk22dr2", dk22dr2);
    PrintVector<6>("(GetDStoParticleBz) dkddr1", dkddr1);
    PrintVector<6>("(GetDStoParticleBz) dkddr2", dkddr2);
    PrintVector<6>("(GetDStoParticleBz) dc1dr1", dc1dr1);
    PrintVector<6>("(GetDStoParticleBz) dc1dr2", dc1dr2);
    PrintVector<6>("(GetDStoParticleBz) dc2dr1", dc2dr1);
    PrintVector<6>("(GetDStoParticleBz) dc2dr2", dc2dr2);
    PrintVector<6>("(GetDStoParticleBz) dd1dr1", dd1dr1);
    PrintVector<6>("(GetDStoParticleBz) dd1dr2", dd1dr2);
#endif

    if (!isStraight1) {
        ds1[0] = std::atan2(bq1 * (k11 * c1 + k21 * d1), bq1 * k11 * d1 * bq1 - k21 * c1) / bq1;
        ds1[1] = std::atan2(bq1 * (k11 * c1 - k21 * d1), -bq1 * k11 * d1 * bq1 - k21 * c1) / bq1;

        double a{bq1 * (k11 * c1 + k21 * d1)};
        double b{bq1 * k11 * d1 * bq1 - k21 * c1};
        for (int iP{0}; iP < 6; ++iP) {
            if ((b * b + a * a) > 0) {
                double dadr1{bq1 * (dk11dr1[iP] * c1 + k11 * dc1dr1[iP] + dk21dr1[iP] * d1 + k21 * dd1dr1[iP])};
                double dadr2{bq1 * (dk11dr2[iP] * c1 + k11 * dc1dr2[iP] + dk21dr2[iP] * d1 + k21 * dd1dr2[iP])};
                double dbdr1{bq1 * bq1 * (dk11dr1[iP] * d1 + k11 * dd1dr1[iP]) - (dk21dr1[iP] * c1 + k21 * dc1dr1[iP])};
                double dbdr2{bq1 * bq1 * (dk11dr2[iP] * d1 + k11 * dd1dr2[iP]) - (dk21dr2[iP] * c1 + k21 * dc1dr2[iP])};

                ds1dR1[0][iP] = 1. / bq1 * 1. / (b * b + a * a) * (dadr1 * b - dbdr1 * a);
                ds1dR2[0][iP] = 1. / bq1 * 1. / (b * b + a * a) * (dadr2 * b - dbdr2 * a);
            } else {
                ds1dR1[0][iP] = 0.;
                ds1dR2[0][iP] = 0.;
            }
        }

        a = bq1 * (k11 * c1 - k21 * d1);
        b = -bq1 * k11 * d1 * bq1 - k21 * c1;
        for (int iP{0}; iP < 6; ++iP) {
            if ((b * b + a * a) > 0) {
                double dadr1{bq1 * (dk11dr1[iP] * c1 + k11 * dc1dr1[iP] - (dk21dr1[iP] * d1 + k21 * dd1dr1[iP]))};
                double dadr2{bq1 * (dk11dr2[iP] * c1 + k11 * dc1dr2[iP] - (dk21dr2[iP] * d1 + k21 * dd1dr2[iP]))};
                double dbdr1{-bq1 * bq1 * (dk11dr1[iP] * d1 + k11 * dd1dr1[iP]) - (dk21dr1[iP] * c1 + k21 * dc1dr1[iP])};
                double dbdr2{-bq1 * bq1 * (dk11dr2[iP] * d1 + k11 * dd1dr2[iP]) - (dk21dr2[iP] * c1 + k21 * dc1dr2[iP])};

                ds1dR1[1][iP] = 1 / bq1 * 1 / (b * b + a * a) * (dadr1 * b - dbdr1 * a);
                ds1dR2[1][iP] = 1 / bq1 * 1 / (b * b + a * a) * (dadr2 * b - dbdr2 * a);
            } else {
                ds1dR1[1][iP] = 0;
                ds1dR2[1][iP] = 0;
            }
        }
    }
    if (!isStraight2) {
        ds2[0] = std::atan2(bq2 * k12 * c2 + k22 * d1 * bq2, (bq2 * k12 * d1 * bq2 - k22 * c2)) / bq2;
        ds2[1] = std::atan2(bq2 * k12 * c2 - k22 * d1 * bq2, (-bq2 * k12 * d1 * bq2 - k22 * c2)) / bq2;

        double a{bq2 * (k12 * c2 + k22 * d1)};
        double b{bq2 * k12 * d1 * bq2 - k22 * c2};
        for (int iP{0}; iP < 6; ++iP) {
            if ((b * b + a * a) > 0) {
                double dadr1{bq2 * (dk12dr1[iP] * c2 + k12 * dc2dr1[iP] + dk22dr1[iP] * d1 + k22 * dd1dr1[iP])};
                double dadr2{bq2 * (dk12dr2[iP] * c2 + k12 * dc2dr2[iP] + dk22dr2[iP] * d1 + k22 * dd1dr2[iP])};
                double dbdr1{bq2 * bq2 * (dk12dr1[iP] * d1 + k12 * dd1dr1[iP]) - (dk22dr1[iP] * c2 + k22 * dc2dr1[iP])};
                double dbdr2{bq2 * bq2 * (dk12dr2[iP] * d1 + k12 * dd1dr2[iP]) - (dk22dr2[iP] * c2 + k22 * dc2dr2[iP])};

                ds2dR1[0][iP] = 1 / bq2 * 1 / (b * b + a * a) * (dadr1 * b - dbdr1 * a);
                ds2dR2[0][iP] = 1 / bq2 * 1 / (b * b + a * a) * (dadr2 * b - dbdr2 * a);
            } else {
                ds2dR1[0][iP] = 0;
                ds2dR2[0][iP] = 0;
            }
        }

        a = bq2 * (k12 * c2 - k22 * d1);
        b = -bq2 * k12 * d1 * bq2 - k22 * c2;
        for (int iP{0}; iP < 6; ++iP) {
            if ((b * b + a * a) > 0) {
                double dadr1{bq2 * (dk12dr1[iP] * c2 + k12 * dc2dr1[iP] - (dk22dr1[iP] * d1 + k22 * dd1dr1[iP]))};
                double dadr2{bq2 * (dk12dr2[iP] * c2 + k12 * dc2dr2[iP] - (dk22dr2[iP] * d1 + k22 * dd1dr2[iP]))};
                double dbdr1{-bq2 * bq2 * (dk12dr1[iP] * d1 + k12 * dd1dr1[iP]) - (dk22dr1[iP] * c2 + k22 * dc2dr1[iP])};
                double dbdr2{-bq2 * bq2 * (dk12dr2[iP] * d1 + k12 * dd1dr2[iP]) - (dk22dr2[iP] * c2 + k22 * dc2dr2[iP])};

                ds2dR1[1][iP] = 1 / bq2 * 1 / (b * b + a * a) * (dadr1 * b - dbdr1 * a);
                ds2dR2[1][iP] = 1 / bq2 * 1 / (b * b + a * a) * (dadr2 * b - dbdr2 * a);
            } else {
                ds2dR1[1][iP] = 0;
                ds2dR2[1][iP] = 0;
            }
        }
    }
    if (isStraight1 && pt12 > 0.) {
        ds1[0] = (k11 * c1 + k21 * d1) / (-k21 * c1);
        ds1[1] = (k11 * c1 - k21 * d1) / (-k21 * c1);

        double a{k11 * c1 + k21 * d1};
        double b{-k21 * c1};

        for (int iP{0}; iP < 6; ++iP) {
            if (b * b > 0) {
                double dadr1{(dk11dr1[iP] * c1 + k11 * dc1dr1[iP] + dk21dr1[iP] * d1 + k21 * dd1dr1[iP])};
                double dadr2{(dk11dr2[iP] * c1 + k11 * dc1dr2[iP] + dk21dr2[iP] * d1 + k21 * dd1dr2[iP])};
                double dbdr1{-(dk21dr1[iP] * c1 + k21 * dc1dr1[iP])};
                double dbdr2{-(dk21dr2[iP] * c1 + k21 * dc1dr2[iP])};

                ds1dR1[0][iP] = dadr1 / b - dbdr1 * a / (b * b);
                ds1dR2[0][iP] = dadr2 / b - dbdr2 * a / (b * b);
            } else {
                ds1dR1[0][iP] = 0;
                ds1dR2[0][iP] = 0;
            }
        }

        a = k11 * c1 - k21 * d1;
        for (int iP{0}; iP < 6; ++iP) {
            if (b * b > 0) {
                double dadr1{(dk11dr1[iP] * c1 + k11 * dc1dr1[iP] - dk21dr1[iP] * d1 - k21 * dd1dr1[iP])};
                double dadr2{(dk11dr2[iP] * c1 + k11 * dc1dr2[iP] - dk21dr2[iP] * d1 - k21 * dd1dr2[iP])};
                double dbdr1{-(dk21dr1[iP] * c1 + k21 * dc1dr1[iP])};
                double dbdr2{-(dk21dr2[iP] * c1 + k21 * dc1dr2[iP])};

                ds1dR1[1][iP] = dadr1 / b - dbdr1 * a / (b * b);
                ds1dR2[1][iP] = dadr2 / b - dbdr2 * a / (b * b);
            } else {
                ds1dR1[1][iP] = 0;
                ds1dR2[1][iP] = 0;
            }
        }
    }
    if (isStraight2 && pt22 > 0.) {
        ds2[0] = (k12 * c2 + k22 * d1) / (-k22 * c2);
        ds2[1] = (k12 * c2 - k22 * d1) / (-k22 * c2);

        double a{k12 * c2 + k22 * d1};
        double b{-k22 * c2};

        for (int iP{0}; iP < 6; ++iP) {
            if (b * b > 0) {
                double dadr1{(dk12dr1[iP] * c2 + k12 * dc2dr1[iP] + dk22dr1[iP] * d1 + k22 * dd1dr1[iP])};
                double dadr2{(dk12dr2[iP] * c2 + k12 * dc2dr2[iP] + dk22dr2[iP] * d1 + k22 * dd1dr2[iP])};
                double dbdr1{-(dk22dr1[iP] * c2 + k22 * dc2dr1[iP])};
                double dbdr2{-(dk22dr2[iP] * c2 + k22 * dc2dr2[iP])};

                ds2dR1[0][iP] = dadr1 / b - dbdr1 * a / (b * b);
                ds2dR2[0][iP] = dadr2 / b - dbdr2 * a / (b * b);
            } else {
                ds2dR1[0][iP] = 0;
                ds2dR2[0][iP] = 0;
            }
        }

        a = k12 * c2 - k22 * d1;
        for (int iP{0}; iP < 6; ++iP) {
            if (b * b > 0) {
                double dadr1{dk12dr1[iP] * c2 + k12 * dc2dr1[iP] - dk22dr1[iP] * d1 - k22 * dd1dr1[iP]};
                double dadr2{dk12dr2[iP] * c2 + k12 * dc2dr2[iP] - dk22dr2[iP] * d1 - k22 * dd1dr2[iP]};
                double dbdr1{-dk22dr1[iP] * c2 - k22 * dc2dr1[iP]};
                double dbdr2{-dk22dr2[iP] * c2 - k22 * dc2dr2[iP]};

                ds2dR1[1][iP] = dadr1 / b - dbdr1 * a / (b * b);
                ds2dR2[1][iP] = dadr2 / b - dbdr2 * a / (b * b);
            } else {
                ds2dR1[1][iP] = 0;
                ds2dR2[1][iP] = 0;
            }
        }
    }

    // select a point which is close to the primary vertex (with the smallest r)

    double dr2[2];
    double tmp_x1[2];
    double tmp_y1[2];
    double tmp_z1[2];
    double tmp_x2[2];
    double tmp_y2[2];
    double tmp_z2[2];
    double tmp_dx[2];
    double tmp_dy[2];
    double tmp_dz[2];
    int winner{0};
    for (int iP{0}; iP < 2; ++iP) {
        double bs1{bq1 * ds1[iP]};
        double bs2{bq2 * ds2[iP]};
        double sss{std::sin(bs1)};
        double ccc{std::cos(bs1)};

        double sB{sss / bq1};
        double cB{(1. - ccc) / bq1};

        tmp_x1[iP] = fP[0] + sB * px1 + cB * py1;
        tmp_y1[iP] = fP[1] - cB * px1 + sB * py1;
        tmp_z1[iP] = fP[2] + ds1[iP] * fP[5];

        sss = std::sin(bs2);
        ccc = std::cos(bs2);

        sB = sss / bq2;
        cB = (1. - ccc) / bq2;

        tmp_x2[iP] = p.fP[0] + sB * px2 + cB * py2;
        tmp_y2[iP] = p.fP[1] - cB * px2 + sB * py2;
        tmp_z2[iP] = p.fP[2] + ds2[iP] * p.fP[5];

        tmp_dx[iP] = tmp_x1[iP] - tmp_x2[iP];
        tmp_dy[iP] = tmp_y1[iP] - tmp_y2[iP];
        tmp_dz[iP] = tmp_z1[iP] - tmp_z2[iP];

        dr2[iP] = tmp_dx[iP] * tmp_dx[iP] + tmp_dy[iP] * tmp_dy[iP] + tmp_dz[iP] * tmp_dz[iP];
    }

    bool isFirstRoot{dr2[0] < dr2[1]};
    if (isFirstRoot) {
        ds[0] = ds1[0];
        ds[1] = ds2[0];
        for (int iP{0}; iP < 6; ++iP) {
            ds_dr[0][iP] = ds1dR1[0][iP];
            ds_dr[1][iP] = ds1dR2[0][iP];
            ds_dr[2][iP] = ds2dR1[0][iP];
            ds_dr[3][iP] = ds2dR2[0][iP];
        }
    } else {
        ds[0] = ds1[1];
        ds[1] = ds2[1];
        for (int iP{0}; iP < 6; ++iP) {
            ds_dr[0][iP] = ds1dR1[1][iP];
            ds_dr[1][iP] = ds1dR2[1][iP];
            ds_dr[2][iP] = ds2dR1[1][iP];
            ds_dr[3][iP] = ds2dR2[1][iP];
        }
        winner = 1;
    }
#if KF_DEBUG
    std::cout << "(GetDStoParticleBz) min1.ds = " << ds[0] << '\n';
    std::cout << "(GetDStoParticleBz) min1.(x,y,z) = " << tmp_x1[winner] << ", " << tmp_y1[winner] << ", " << tmp_z1[winner] << '\n';
    std::cout << "(GetDStoParticleBz) min2.ds = " << ds[1] << '\n';
    std::cout << "(GetDStoParticleBz) min2.(x,y,z) = " << tmp_x2[winner] << ", " << tmp_y2[winner] << ", " << tmp_z2[winner] << '\n';
    std::cout << "(GetDStoParticleBz) dx = " << tmp_dx[winner] << '\n';
    std::cout << "(GetDStoParticleBz) dy = " << tmp_dy[winner] << '\n';
    std::cout << "(GetDStoParticleBz) dz = " << tmp_dz[winner] << '\n';
    std::cout << "(GetDStoParticleBz) dca_3d = " << std::sqrt(dr2[winner]) << '\n';
#endif

    // 2. add z-component as small correction //

    double bs1{bq1 * ds[0]};
    double bs2{bq2 * ds[1]};
    double sss{std::sin(bs1)};
    double ccc{std::cos(bs1)};

    double sB{sss / bq1};
    double cB{(1. - ccc) / bq1};

    double x1{x01 + sB * px1 + cB * py1};
    double y1{y01 - cB * px1 + sB * py1};
    double z1{z01 + ds[0] * pz1};
    double ppx1{ccc * px1 + sss * py1};
    double ppy1{-sss * px1 + ccc * py1};
    double ppz1{pz1};

    double sss1{std::sin(bs2)};
    double ccc1{std::cos(bs2)};

    double sB1{sss1 / bq2};
    double cB1{(1. - ccc1) / bq2};

    double x2{x02 + sB1 * px2 + cB1 * py2};
    double y2{y02 - cB1 * px2 + sB1 * py2};
    double z2{z02 + ds[1] * pz2};
    double ppx2{ccc1 * px2 + sss1 * py2};
    double ppy2{-sss1 * px2 + ccc1 * py2};
    double ppz2{pz2};

    double p12{ppx1 * ppx1 + ppy1 * ppy1 + ppz1 * ppz1};
    double p22{ppx2 * ppx2 + ppy2 * ppy2 + ppz2 * ppz2};
    double lp1p2{ppx1 * ppx2 + ppy1 * ppy2 + ppz1 * ppz2};

    double dx{x2 - x1};
    double dy{y2 - y1};
    double dz{z2 - z1};

    double ldrp1{ppx1 * dx + ppy1 * dy + ppz1 * dz};
    double ldrp2{ppx2 * dx + ppy2 * dy + ppz2 * dz};

    double detp{lp1p2 * lp1p2 - p12 * p22};
    if (std::abs(detp) < 1.e-4) detp = 1;  // PENDING

    // ds_dr calculation
    double a1{ldrp2 * lp1p2 - ldrp1 * p22};
    double a2{ldrp2 * p12 - ldrp1 * lp1p2};
    double lp1p2_ds0{bq1 * (ppx2 * ppy1 - ppy2 * ppx1)};
    double lp1p2_ds1{bq2 * (ppx1 * ppy2 - ppy1 * ppx2)};
    double ldrp1_ds0{-p12 + bq1 * (ppy1 * dx - ppx1 * dy)};
    double ldrp1_ds1{lp1p2};
    double ldrp2_ds0{-lp1p2};
    double ldrp2_ds1{p22 + bq2 * (ppy2 * dx - ppx2 * dy)};
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
    std::cout << "(GetDStoParticleBz) a1 = " << a1 << '\n';
    std::cout << "(GetDStoParticleBz) a2 = " << a2 << '\n';
    std::cout << "(GetDStoParticleBz) lp1p2_ds0 = " << lp1p2_ds0 << '\n';
    std::cout << "(GetDStoParticleBz) lp1p2_ds1 = " << lp1p2_ds1 << '\n';
    std::cout << "(GetDStoParticleBz) ldrp1_ds0 = " << ldrp1_ds0 << '\n';
    std::cout << "(GetDStoParticleBz) ldrp1_ds1 = " << ldrp1_ds1 << '\n';
    std::cout << "(GetDStoParticleBz) ldrp2_ds0 = " << ldrp2_ds0 << '\n';
    std::cout << "(GetDStoParticleBz) ldrp2_ds1 = " << ldrp2_ds1 << '\n';
    std::cout << "(GetDStoParticleBz) detp_ds0 = " << detp_ds0 << '\n';
    std::cout << "(GetDStoParticleBz) detp_ds1 = " << detp_ds1 << '\n';
    std::cout << "(GetDStoParticleBz) a1_ds0 = " << a1_ds0 << '\n';
    std::cout << "(GetDStoParticleBz) a1_ds1 = " << a1_ds1 << '\n';
    std::cout << "(GetDStoParticleBz) a2_ds0 = " << a2_ds0 << '\n';
    std::cout << "(GetDStoParticleBz) a2_ds1 = " << a2_ds1 << '\n';
    std::cout << "(GetDStoParticleBz) dsl1ds0 = " << dsl1ds0 << '\n';
    std::cout << "(GetDStoParticleBz) dsl1ds1 = " << dsl1ds1 << '\n';
    std::cout << "(GetDStoParticleBz) dsl2ds0 = " << dsl2ds0 << '\n';
    std::cout << "(GetDStoParticleBz) dsl2ds1 = " << dsl2ds1 << '\n';
#endif
    double dsldr[4][6];
    for (int iP{0}; iP < 6; ++iP) {
        dsldr[0][iP] = dsl1ds0 * ds_dr[0][iP] + dsl1ds1 * ds_dr[2][iP];
        dsldr[1][iP] = dsl1ds0 * ds_dr[1][iP] + dsl1ds1 * ds_dr[3][iP];
        dsldr[2][iP] = dsl2ds0 * ds_dr[0][iP] + dsl2ds1 * ds_dr[2][iP];
        dsldr[3][iP] = dsl2ds0 * ds_dr[1][iP] + dsl2ds1 * ds_dr[3][iP];
    }

    for (int iDS{0}; iDS < 4; ++iDS) {
        for (int iP{0}; iP < 6; ++iP) ds_dr[iDS][iP] += dsldr[iDS][iP];
    }

    Vector<6> lp1p2_dr0{0., 0., 0., ccc * ppx2 - ppy2 * sss, ccc * ppy2 + ppx2 * sss, pz2};
    Vector<6> lp1p2_dr1{0., 0., 0., ccc1 * ppx1 - ppy1 * sss1, ccc1 * ppy1 + ppx1 * sss1, pz1};
    Vector<6> ldrp1_dr0{
        -ppx1, -ppy1, -pz1, cB * ppy1 - ppx1 * sB + ccc * dx - sss * dy, -cB * ppx1 - ppy1 * sB + sss * dx + ccc * dy, -ds[0] * pz1 + dz};
    Vector<6> ldrp1_dr1{ppx1, ppy1, pz1, -cB1 * ppy1 + ppx1 * sB1, cB1 * ppx1 + ppy1 * sB1, ds[1] * pz1};
    Vector<6> ldrp2_dr0{-ppx2, -ppy2, -pz2, cB * ppy2 - ppx2 * sB, -cB * ppx2 - ppy2 * sB, -ds[0] * pz2};
    Vector<6> ldrp2_dr1{
        ppx2, ppy2, pz2, -cB1 * ppy2 + ppx2 * sB1 + ccc1 * dx - sss1 * dy, cB1 * ppx2 + ppy2 * sB1 + sss1 * dx + ccc1 * dy, dz + ds[1] * pz2};
    Vector<6> p12_dr0{0., 0., 0., 2. * px1, 2. * py1, 2. * pz1};
    Vector<6> p22_dr1{0., 0., 0., 2. * px2, 2. * py2, 2. * pz2};
    for (int iP{0}; iP < 6; ++iP) {
        double a1_dr0{ldrp2_dr0[iP] * lp1p2 + ldrp2 * lp1p2_dr0[iP] - ldrp1_dr0[iP] * p22};
        double a1_dr1{ldrp2_dr1[iP] * lp1p2 + ldrp2 * lp1p2_dr1[iP] - ldrp1_dr1[iP] * p22 - ldrp1 * p22_dr1[iP]};
        double a2_dr0{ldrp2_dr0[iP] * p12 + ldrp2 * p12_dr0[iP] - ldrp1_dr0[iP] * lp1p2 - ldrp1 * lp1p2_dr0[iP]};
        double a2_dr1{ldrp2_dr1[iP] * p12 - ldrp1_dr1[iP] * lp1p2 - ldrp1 * lp1p2_dr1[iP]};
        double detp_dr0{2. * lp1p2 * lp1p2_dr0[iP] - p12_dr0[iP] * p22};
        double detp_dr1{2. * lp1p2 * lp1p2_dr1[iP] - p12 * p22_dr1[iP]};

        ds_dr[0][iP] += a1_dr0 / detp - a1 * detp_dr0 / (detp * detp);
        ds_dr[1][iP] += a1_dr1 / detp - a1 * detp_dr1 / (detp * detp);
        ds_dr[2][iP] += a2_dr0 / detp - a2 * detp_dr0 / (detp * detp);
        ds_dr[3][iP] += a2_dr1 / detp - a2 * detp_dr1 / (detp * detp);
    }

    ds[0] += a1 / detp;
    ds[1] += a2 / detp;
#if KF_DEBUG
    PrintVector<6>("(GetDStoParticleBz) lp1p2_dr0", lp1p2_dr0);
    PrintVector<6>("(GetDStoParticleBz) lp1p2_dr1", lp1p2_dr1);
    PrintVector<6>("(GetDStoParticleBz) ldrp1_dr0", ldrp1_dr0);
    PrintVector<6>("(GetDStoParticleBz) ldrp1_dr1", ldrp1_dr1);
    PrintVector<6>("(GetDStoParticleBz) ldrp2_dr0", ldrp2_dr0);
    PrintVector<6>("(GetDStoParticleBz) ldrp2_dr1", ldrp2_dr1);
    PrintVector<6>("(GetDStoParticleBz) p12_dr0", p12_dr0);
    PrintVector<6>("(GetDStoParticleBz) p22_dr1", p22_dr1);
    PrintVector<6>("(GetDStoParticleBz) min1.dsdr", ds_dr[0]);
    PrintVector<6>("(GetDStoParticleBz) min1.dsdr1", ds_dr[1]);
    PrintVector<6>("(GetDStoParticleBz) min2.dsdr1", ds_dr[2]);
    PrintVector<6>("(GetDStoParticleBz) min2.dsdr", ds_dr[3]);
    std::cout << "(GetDStoParticleBz) min1.ds (after z-correction) = " << ds[0] << '\n';
    std::cout << "(GetDStoParticleBz) min2.ds (after z-correction) = " << ds[1] << '\n';
#endif
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
// \param[in] p - second particle
// \param[out] dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
// \param[out] dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
void Particle::GetDStoParticleLine(const Particle& p, double ds[2], Vector<6> dsdr[4]) const {

    double p12{fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]};
    double p22{p.fP[3] * p.fP[3] + p.fP[4] * p.fP[4] + p.fP[5] * p.fP[5]};
    double p1p2{fP[3] * p.fP[3] + fP[4] * p.fP[4] + fP[5] * p.fP[5]};

    double drp1{fP[3] * (p.fP[0] - fP[0]) + fP[4] * (p.fP[1] - fP[1]) + fP[5] * (p.fP[2] - fP[2])};
    double drp2{p.fP[3] * (p.fP[0] - fP[0]) + p.fP[4] * (p.fP[1] - fP[1]) + p.fP[5] * (p.fP[2] - fP[2])};

    double detp{p1p2 * p1p2 - p12 * p22};
    if (std::abs(detp) < 1.e-4) detp = 1;  // PENDING

    ds[0] = (drp2 * p1p2 - drp1 * p22) / detp;
    ds[1] = (drp2 * p12 - drp1 * p1p2) / detp;

    double x01{fP[0]};
    double y01{fP[1]};
    double z01{fP[2]};
    double px1{fP[3]};
    double py1{fP[4]};
    double pz1{fP[5]};

    double x02{p.fP[0]};
    double y02{p.fP[1]};
    double z02{p.fP[2]};
    double px2{p.fP[3]};
    double py2{p.fP[4]};
    double pz2{p.fP[5]};

    Vector<6> drp1_dr1{-px1, -py1, -pz1, -x01 + x02, -y01 + y02, -z01 + z02};
    Vector<6> drp1_dr2{px1, py1, pz1, 0., 0., 0.};
    Vector<6> drp2_dr1{-px2, -py2, -pz2, 0., 0., 0.};
    Vector<6> drp2_dr2{px2, py2, pz2, -x01 + x02, -y01 + y02, -z01 + z02};
    Vector<6> dp1p2_dr1{0., 0., 0., px2, py2, pz2};
    Vector<6> dp1p2_dr2{0., 0., 0., px1, py1, pz1};
    Vector<6> dp12_dr1{0., 0., 0., 2. * px1, 2. * py1, 2. * pz1};
    Vector<6> dp12_dr2{0., 0., 0., 0., 0., 0.};
    Vector<6> dp22_dr1{0., 0., 0., 0., 0., 0.};
    Vector<6> dp22_dr2{0., 0., 0., 2. * px2, 2. * py2, 2. * pz2};
    Vector<6> ddetp_dr1{0., 0., 0., -2 * p22 * px1 + 2. * p1p2 * px2, -2 * p22 * py1 + 2. * p1p2 * py2, -2 * p22 * pz1 + 2. * p1p2 * pz2};
    Vector<6> ddetp_dr2{0., 0., 0., 2. * p1p2 * px1 - 2. * p12 * px2, 2. * p1p2 * py1 - 2. * p12 * py2, 2. * p1p2 * pz1 - 2. * p12 * pz2};

    double a1{drp2 * p1p2 - drp1 * p22};
    double a2{drp2 * p12 - drp1 * p1p2};
    for (int i{0}; i < 6; ++i) {
        double da1_dr1{drp2_dr1[i] * p1p2 + drp2 * dp1p2_dr1[i] - drp1_dr1[i] * p22 - drp1 * dp22_dr1[i]};
        double da1_dr2{drp2_dr2[i] * p1p2 + drp2 * dp1p2_dr2[i] - drp1_dr2[i] * p22 - drp1 * dp22_dr2[i]};
        double da2_dr1{drp2_dr1[i] * p12 + drp2 * dp12_dr1[i] - drp1_dr1[i] * p1p2 - drp1 * dp1p2_dr1[i]};
        double da2_dr2{drp2_dr2[i] * p12 + drp2 * dp12_dr2[i] - drp1_dr2[i] * p1p2 - drp1 * dp1p2_dr2[i]};

        dsdr[0][i] = da1_dr1 / detp - a1 * ddetp_dr1[i] / (detp * detp);
        dsdr[1][i] = da1_dr2 / detp - a1 * ddetp_dr2[i] / (detp * detp);
        dsdr[2][i] = da2_dr1 / detp - a2 * ddetp_dr1[i] / (detp * detp);
        dsdr[3][i] = da2_dr2 / detp - a2 * ddetp_dr2[i] / (detp * detp);
    }
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
void Particle::TransportBz(double bz, double ds, const Vector<6>& dsdr, Vector<8>& P, SymMatrix<8>& C, const Vector<6>& dsdr1, Matrix<6, 6>& jacob,
                           Matrix<6, 6>& corr) const {
#if KF_DEBUG
    std::cout << "-- starting (TransportBz) --" << '\n';
#endif
    double bq{bz * fQ * Const::Kappa};
    double bs{bq * ds};
    double s{std::sin(bs)};
    double c{std::cos(bs)};
    double sB{s / bq};
    double cB{(1 - c) / bq};

    double px{fP[3]};
    double py{fP[4]};
    double pz{fP[5]};

    P[0] = fP[0] + sB * px + cB * py;
    P[1] = fP[1] - cB * px + sB * py;
    P[2] = fP[2] + ds * pz;
    P[3] = c * px + s * py;
    P[4] = -s * px + c * py;
    P[5] = fP[5];
    P[6] = fP[6];
    P[7] = fP[7];

    auto mJ = Zero<8, 8>();
    for (int i{0}; i < 8; ++i) mJ[i][i] = 1.;
    mJ[0][3] = sB;
    mJ[0][4] = cB;
    mJ[1][3] = -cB;
    mJ[1][4] = sB;
    mJ[2][5] = ds;
    mJ[3][3] = c;
    mJ[3][4] = s;
    mJ[4][3] = -s;
    mJ[4][4] = c;

    auto mJds = Zero<6, 6>();
    mJds[0][3] = c;
    mJds[0][4] = s;
    mJds[1][3] = -s;
    mJds[1][4] = c;
    mJds[2][5] = 1.;
    mJds[3][3] = -bq * s;
    mJds[3][4] = bq * c;
    mJds[4][3] = -bq * c;
    mJds[4][4] = -bq * s;

    for (int i1{0}; i1 < 6; ++i1) {
        for (int i2{0}; i2 < 6; ++i2) mJ[i1][i2] += mJds[i1][3] * px * dsdr[i2] + mJds[i1][4] * py * dsdr[i2] + mJds[i1][5] * pz * dsdr[i2];
    }

    C = MultQSQt<8>(mJ, fC);

    for (int i{0}; i < 6; ++i) {
        for (int j{0}; j < 6; ++j) jacob[i][j] = mJ[i][j];
    }

    for (int i1{0}; i1 < 6; ++i1) {
        for (int i2{0}; i2 < 6; ++i2) corr[i1][i2] = mJds[i1][3] * px * dsdr1[i2] + mJds[i1][4] * py * dsdr1[i2] + mJds[i1][5] * pz * dsdr1[i2];
    }
#if KF_DEBUG
    PrintVector<8>("(TransportBz) State", P);
    PrintSymMatrix<8>("(TransportBz) Cov", C);
    PrintMatrix<6, 6>("(TransportBz) Jacob", jacob);
    PrintMatrix<6, 6>("(TransportBz) Corr", corr);
    std::cout << "-- finished (TransportBz) --" << '\n';
#endif
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
// *Parameters F and F1 should be either both initialised or
//     both set to null pointer.
// \param[in] dS - transport parameter which defines the distance to which particle should be transported
// \param[in] dsdr[6] = ds/dr - partial derivatives of the parameter dS over the state vector of the current particle
// \param[out] P[8] - array, where transported parameters should be stored
// \param[out] C[36] - array, where transported covariance matrix (8x8) should be stored in the lower triangular form
// \param[in] dsdr1[6] = ds/dr - partial derivatives of the parameter dS over the state vector of another particle
// or vertex
// \param[out] F[36] - optional parameter, transport jacobian, 6x6 matrix F = d(fP new)/d(fP old)
// \param[out] F1[36] - optional parameter, corelation 6x6 matrix betweeen the current particle and particle or vertex
// with the state vector r1, to which the current particle is being transported, F1 = d(fP new)/d(r1)
void Particle::TransportLine(double ds, const Vector<6>& ds_dr, Vector<8>& P, SymMatrix<8>& C, const Vector<6>& ds_dr1, Matrix<6, 6>& jacob,
                             Matrix<6, 6>& corr) const {

    auto mJ = Zero<8, 8>();
    mJ[0][0] = 1.;
    mJ[0][3] = ds;
    mJ[1][1] = 1.;
    mJ[1][4] = ds;
    mJ[2][2] = 1.;
    mJ[2][5] = ds;
    mJ[3][3] = 1.;
    mJ[4][4] = 1.;
    mJ[5][5] = 1.;
    mJ[6][6] = 1.;
    mJ[7][7] = 1.;

    double px{fP[3]};
    double py{fP[4]};
    double pz{fP[5]};

    P[0] = fP[0] + ds * fP[3];
    P[1] = fP[1] + ds * fP[4];
    P[2] = fP[2] + ds * fP[5];
    P[3] = fP[3];
    P[4] = fP[4];
    P[5] = fP[5];
    P[6] = fP[6];
    P[7] = fP[7];

    auto mJds = Zero<6, 6>();
    mJds[0][3] = 1.;
    mJds[1][4] = 1.;
    mJds[2][5] = 1.;

    for (int i1{0}; i1 < 6; ++i1) {
        for (int i2{0}; i2 < 6; ++i2) {
            mJ[i1][i2] += mJds[i1][3] * px * ds_dr[i2] + mJds[i1][4] * py * ds_dr[i2] + mJds[i1][5] * pz * ds_dr[i2];
        }
    }
    C = MultQSQt<8>(mJ, fC);

    for (int i{0}; i < 6; ++i) {
        for (int j{0}; j < 6; ++j) jacob[i][j] = mJ[i][j];
    }
    for (int i1{0}; i1 < 6; ++i1) {
        for (int i2{0}; i2 < 6; ++i2) {
            corr[i1][i2] = mJds[i1][3] * px * ds_dr1[i2] + mJds[i1][4] * py * ds_dr1[i2] + mJds[i1][5] * pz * ds_dr1[i2];
        }
    }
}

// Symmetric 3x3 matrix a using modified Cholesky decomposition. The result is stored to the same matrix a.
// \param[in,out] a - 3x3 symmetric matrix
void Particle::InvertCholesky3(SymMatrix<3>& a) {

    auto d = Zero<3>();
    auto u = Zero<3, 3>();

    for (int i{0}; i < 3; ++i) {
        double uud{0.};
        for (int j{0}; j < i; ++j) uud += u[j][i] * u[j][i] * d[j];
        uud = a[i * (i + 3) / 2] - uud;

        if (std::abs(uud) < Const::AbsAlmostZero) uud = Const::AbsAlmostZero;

        d[i] = uud / std::abs(uud);
        u[i][i] = std::sqrt(std::abs(uud));

        for (int j{i + 1}; j < 3; ++j) {
            uud = 0.;
            for (int k{0}; k < i; ++k) uud += u[k][i] * u[k][j] * d[k];
            uud = a[j * (j + 1) / 2 + i] - uud;
            u[i][j] = d[i] / u[i][i] * uud;
        }
    }

    Vector<3> u1;

    for (int i{0}; i < 3; ++i) {
        u1[i] = u[i][i];
        u[i][i] = 1 / u[i][i];
    }
    for (int i{0}; i < 2; ++i) {
        u[i][i + 1] = -u[i][i + 1] * u[i][i] * u[i + 1][i + 1];
    }
    for (int i{0}; i < 1; ++i) {
        u[i][i + 2] = u[i][i + 1] * u1[i + 1] * u[i + 1][i + 2] - u[i][i + 2] * u[i][i] * u[i + 2][i + 2];
    }

    for (int i{0}; i < 3; ++i) a[i + 3] = u[i][2] * u[2][2] * d[2];
    for (int i{0}; i < 2; ++i) a[i + 1] = u[i][1] * u[1][1] * d[1] + u[i][2] * u[1][2] * d[2];
    a[0] = u[0][0] * u[0][0] * d[0] + u[0][1] * u[0][1] * d[1] + u[0][2] * u[0][2] * d[2];
}

// Matrix multiplication SOut = Q*S*Q^T, where Q - square matrix, S - symmetric matrix.
// \param[in] Q - square matrix
// \param[in] S - input symmetric matrix
template <int N>
SymMatrix<N> Particle::MultQSQt(const Matrix<N, N>& Q, const SymMatrix<N>& S) const {

    Matrix<N, N> SQT;

    for (int i{0}; i < N; ++i) {
        for (int j{0}; j < N; ++j) {
            SQT[i][j] = 0.;
            for (int k{0}; k < N; ++k) SQT[i][j] += S[IJ(i, k)] * Q[j][k];
        }
    }

    SymMatrix<N> SOut;

    for (int i{0}; i < N; ++i) {
        for (int j{0}; j <= i; ++j) {
            SOut[IJ(i, j)] = 0.;
            for (int k{0}; k < N; ++k) SOut[IJ(i, j)] += Q[i][k] * SQT[k][j];
        }
    }

    return SOut;
}

}  // namespace KF
