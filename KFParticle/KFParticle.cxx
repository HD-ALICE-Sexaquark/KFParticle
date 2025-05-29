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

// set the parameters of the particle:
//  \param[in] param[6] = { X, Y, Z, Px, Py, Pz } - position and momentum
//  \param[in] cov[21]  - lower-triangular part of the covariance matrix:@n
//           (  0  .  .  .  .  . )
//           (  1  2  .  .  .  . )
// Cov[21] = (  3  4  5  .  .  . )
//           (  6  7  8  9  .  . )
//           ( 10 11 12 13 14  . )
//           ( 15 16 17 18 19 20 )
//  \param[in] charge - charge of the particle in elementary charge units
//  \param[in] mass - the mass hypothesis
void KFParticle::Initialize(const float param[], const float cov[], int charge, float mass) {

    for (int i = 0; i < 6; ++i) fP[i] = param[i];
    for (int i = 0; i < 21; ++i) fC[i] = cov[i];

    float energy = std::sqrt(mass * mass + fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]);
    fP[6] = energy;
    fP[7] = 0.;
    fQ = charge;
    fNDF = 0;
    fChi2 = 0.;
    fAtProductionVertex = false;
    fSFromDecay = 0.;

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

    fSumDaughterMass = mass;
    fMassHypo = mass;
}

// Initialises the parameters by default:
// 1) all parameters are set to 0;
// 2) all elements of the covariance matrix are set to 0 except Cxx=Cyy=Czz=100;
// 3) Q = 0;
// 4) chi2 is set to 0;
// 5) NDF = -3, since 3 parameters should be fitted: X, Y, Z.
void KFParticle::Initialize() {
    for (int i = 0; i < 8; ++i) fP[i] = 0.;
    for (int i = 0; i < 36; ++i) fC[i] = 0.;
    fC[0] = fC[2] = fC[5] = 100.;
    fC[35] = 1.;
    fNDF = -3;
    fChi2 = 0.;
    fQ = 0;
    fSFromDecay = 0.;
    fAtProductionVertex = false;
    fSumDaughterMass = 0.;
    fMassHypo = -1.;
}

// Obtains the measurements from the current particle and the daughter to be added for the Kalman filter
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
bool KFParticle::GetMeasurement(float bz, const KFParticle& daughter, float m[], float V[], float D[3][3]) {

    if (fNDF == -1) {
        float ds[2]{0., 0.};
        float dsdr[4][6];
        float F1[36];
        float F2[36];
        float F3[36];
        float F4[36];
        for (int i1 = 0; i1 < 36; ++i1) {
            F1[i1] = 0.;
            F2[i1] = 0.;
            F3[i1] = 0.;
            F4[i1] = 0.;
        }
        GetDStoParticleBz(bz, daughter, ds, dsdr);

        if (std::abs(ds[0] * fP[5]) > 1000. || std::abs(ds[1] * daughter.fP[5]) > 1000.) return 0;

        float V0Tmp[36] = {0.};
        float V1Tmp[36] = {0.};

        float C[36];
        for (int iC = 0; iC < 36; ++iC) C[iC] = fC[iC];

        TransportBz(bz, ds[0], dsdr[0], fP, fC, dsdr[1], F1, F2);
        daughter.TransportBz(bz, ds[1], dsdr[3], m, V, dsdr[2], F4, F3);
#if KF_DEBUG
        PrintNonSymMatrix<float, 36>("(AddDaughterWithEnergyFit) F1", F1);
        PrintNonSymMatrix<float, 36>("(AddDaughterWithEnergyFit) F2", F2);
        PrintNonSymMatrix<float, 36>("(AddDaughterWithEnergyFit) F3", F3);
        PrintNonSymMatrix<float, 36>("(AddDaughterWithEnergyFit) F4", F4);
#endif
        MultQSQt(F2, daughter.fC, V0Tmp, 6);
        MultQSQt(F3, C, V1Tmp, 6);

        for (int iC = 0; iC < 21; ++iC) {
            fC[iC] += V0Tmp[iC];
            V[iC] += V1Tmp[iC];
        }

        float C1F1T[6][6];
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                C1F1T[i][j] = 0.;
                for (int k = 0; k < 6; ++k) {
                    C1F1T[i][j] += C[IJ(i, k)] * F1[j * 6 + k];
                }
            }
        float F3C1F1T[6][6];
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                F3C1F1T[i][j] = 0.;
                for (int k = 0; k < 6; ++k) {
                    F3C1F1T[i][j] += F3[i * 6 + k] * C1F1T[k][j];
                }
            }
        float C2F2T[6][6];
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                C2F2T[i][j] = 0.;
                for (int k = 0; k < 6; ++k) {
                    C2F2T[i][j] += daughter.fC[IJ(i, k)] * F2[j * 6 + k];
                }
            }
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                D[i][j] = F3C1F1T[i][j];
                for (int k = 0; k < 6; ++k) {
                    D[i][j] += F4[i * 6 + k] * C2F2T[k][j];
                }
            }
    } else {
        float dsdr[6];
        float dS = daughter.GetDStoPointBz(bz, fP, dsdr);

        float dsdp[6] = {-dsdr[0], -dsdr[1], -dsdr[2], 0, 0, 0};

        float F[36];
        float F1[36];
        for (int i2 = 0; i2 < 36; ++i2) {
            F[i2] = 0.;
            F1[i2] = 0.;
        }
        daughter.TransportBz(bz, dS, dsdr, m, V, dsdp, F, F1);

        //     float V1Tmp[36] = {0.};
        //     MultQSQt(F1, fC, V1Tmp, 6);

        //     for(int iC=0; iC<21; ++iC)
        //       V[iC] += V1Tmp[iC];

        float VFT[3][6];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 6; ++j) {
                VFT[i][j] = 0.;
                for (int k = 0; k < 3; ++k) {
                    VFT[i][j] += fC[IJ(i, k)] * F1[j * 6 + k];
                }
            }

        float FVFT[6][6];
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                FVFT[i][j] = 0.;
                for (int k = 0; k < 3; ++k) {
                    FVFT[i][j] += F1[i * 6 + k] * VFT[k][j];
                }
            }

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                D[i][j] = 0.;
                for (int k = 0; k < 3; ++k) {
                    D[i][j] += fC[IJ(j, k)] * F1[i * 6 + k];
                }
            }

        V[0] += FVFT[0][0];
        V[1] += FVFT[1][0];
        V[2] += FVFT[1][1];
        V[3] += FVFT[2][0];
        V[4] += FVFT[2][1];
        V[5] += FVFT[2][2];

        //     if(fNDF > 100)
        //     {
        //       float dx = fP[0] - m[0];
        //       float dy = fP[1] - m[1];
        //       float dz = fP[2] - m[2];
        //       float sigmaS = 3.f*std::sqrt( (dx*dx + dy*dy + dz*dz) / (m[3]*m[3] + m[4]*m[4] + m[5]*m[5]) );
        //
        //       float h[3] = { m[3]*sigmaS, m[4]*sigmaS, m[5]*sigmaS };
        //       V[0]+= h[0]*h[0];
        //       V[1]+= h[1]*h[0];
        //       V[2]+= h[1]*h[1];
        //       V[3]+= h[2]*h[0];
        //       V[4]+= h[2]*h[1];
        //       V[5]+= h[2]*h[2];
        //     }
    }

    return 1;
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
void KFParticle::AddDaughter(float bz, const KFParticle& daughter) {

    if (fNDF < -1) {  // first daughter -> just copy
        fNDF = -1;
        fQ = daughter.GetQ();
        for (int i = 0; i < 7; ++i) fP[i] = daughter.fP[i];
        for (int i = 0; i < 28; ++i) fC[i] = daughter.fC[i];
        fSFromDecay = 0;
        fMassHypo = daughter.fMassHypo;
        fSumDaughterMass = daughter.fSumDaughterMass;
        return;
    }

    if (static_cast<int>(fConstructMethod) == 0)
        AddDaughterWithEnergyFit(bz, daughter);
    else if (static_cast<int>(fConstructMethod) == 2)
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
void KFParticle::AddDaughterWithEnergyFit(float bz, const KFParticle& daughter) {

    float m[8];
    float mV[36];

    float D[3][3];
    if (!GetMeasurement(bz, daughter, m, mV, D)) return;

    float mS[6]{fC[0] + mV[0], fC[1] + mV[1], fC[2] + mV[2], fC[3] + mV[3], fC[4] + mV[4], fC[5] + mV[5]};

    InvertCholetsky3(mS);
#if KF_DEBUG
    PrintVector<float, 8>("(AddDaughterWithEnergyFit) m", m);
    PrintSymMatrix<float, 36>("(AddDaughterWithEnergyFit) mV", mV);
    PrintMatrix<float, 3, 3>("(AddDaughterWithEnergyFit) D", D);
    PrintSymMatrix<float, 6>("(AddDaughterWithEnergyFit) mS", mS);
#endif
    // Residual (measured - estimated)

    float zeta[3] = {m[0] - fP[0], m[1] - fP[1], m[2] - fP[2]};

    float dChi2 = (mS[0] * zeta[0] + mS[1] * zeta[1] + mS[3] * zeta[2]) * zeta[0] + (mS[1] * zeta[0] + mS[2] * zeta[1] + mS[4] * zeta[2]) * zeta[1] +
                  (mS[3] * zeta[0] + mS[4] * zeta[1] + mS[5] * zeta[2]) * zeta[2];
    if (dChi2 > 1e9) return;
    //     if(fNDF > 100 && dChi2 > 9) return;

    float K[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            K[i][j] = 0;
            for (int k = 0; k < 3; ++k) K[i][j] += fC[IJ(i, k)] * mS[IJ(k, j)];
        }
#if KF_DEBUG
    PrintVector<float, 3>("(AddDaughterWithEnergyFit) zeta", zeta);
    std::cout << "(AddDaughterWithEnergyFit) dChi2 = " << dChi2 << '\n';
    PrintMatrix<float, 3, 3>("(AddDaughterWithEnergyFit) K", K);
#endif

    // CHt = CH' - D'
    float mCHt0[7];
    float mCHt1[7];
    float mCHt2[7];

    mCHt0[0] = fC[0];
    mCHt1[0] = fC[1];
    mCHt2[0] = fC[3];
    mCHt0[1] = fC[1];
    mCHt1[1] = fC[2];
    mCHt2[1] = fC[4];
    mCHt0[2] = fC[3];
    mCHt1[2] = fC[4];
    mCHt2[2] = fC[5];
    mCHt0[3] = fC[6] - mV[6];
    mCHt1[3] = fC[7] - mV[7];
    mCHt2[3] = fC[8] - mV[8];
    mCHt0[4] = fC[10] - mV[10];
    mCHt1[4] = fC[11] - mV[11];
    mCHt2[4] = fC[12] - mV[12];
    mCHt0[5] = fC[15] - mV[15];
    mCHt1[5] = fC[16] - mV[16];
    mCHt2[5] = fC[17] - mV[17];
    mCHt0[6] = fC[21] - mV[21];
    mCHt1[6] = fC[22] - mV[22];
    mCHt2[6] = fC[23] - mV[23];

    // Kalman gain K = mCH'*S

    float k0[7];
    float k1[7];
    float k2[7];

    for (int i = 0; i < 7; ++i) {
        k0[i] = mCHt0[i] * mS[0] + mCHt1[i] * mS[1] + mCHt2[i] * mS[3];
        k1[i] = mCHt0[i] * mS[1] + mCHt1[i] * mS[2] + mCHt2[i] * mS[4];
        k2[i] = mCHt0[i] * mS[3] + mCHt1[i] * mS[4] + mCHt2[i] * mS[5];
    }
#if KF_DEBUG
    PrintJoinedMatrix<float, 7>("(AddDaughterWithEnergyFit) mCH", mCHt0, mCHt1, mCHt2);
    PrintJoinedMatrix<float, 7>("(AddDaughterWithEnergyFit) KGain", k0, k1, k2);
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
    PrintVector<float, 8>("(AddDaughterWithEnergyFit) fP", fP);
    PrintSymMatrix<float, 36>("(AddDaughterWithEnergyFit) fC", fC);
#endif

    // New estimation of the vertex position r += K*zeta

    for (int i = 0; i < 7; ++i) fP[i] = fP[i] + k0[i] * zeta[0] + k1[i] * zeta[1] + k2[i] * zeta[2];

    // New covariance matrix C -= K*(mCH')'

    for (int i = 0, k = 0; i < 7; ++i) {
        for (int j = 0; j <= i; ++j, ++k) {
            fC[k] = fC[k] - (k0[i] * mCHt0[j] + k1[i] * mCHt1[j] + k2[i] * mCHt2[j]);
        }
    }

    float K2[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) K2[i][j] = -K[j][i];
        K2[i][i] += 1.;
    }

    float A[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[i][j] = 0.;
            for (int k = 0; k < 3; ++k) {
                A[i][j] += D[i][k] * K2[k][j];
            }
        }
    }

    double M[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            M[i][j] = 0.;
            for (int k = 0; k < 3; ++k) {
                M[i][j] += K[i][k] * A[k][j];
            }
        }
    }

    fC[0] += 2. * M[0][0];
    fC[1] += M[0][1] + M[1][0];
    fC[2] += 2. * M[1][1];
    fC[3] += M[0][2] + M[2][0];
    fC[4] += M[1][2] + M[2][1];
    fC[5] += 2. * M[2][2];
#if KF_DEBUG
    PrintVector<float, 8>("(AddDaughterWithEnergyFit) fP", fP);
    PrintMatrix<float, 3, 3>("(AddDaughterWithEnergyFit) K2", K2);
    PrintMatrix<float, 3, 3>("(AddDaughterWithEnergyFit) A", A);
    PrintMatrix<double, 3, 3>("(AddDaughterWithEnergyFit) M", M);
    PrintSymMatrix<float, 36>("(AddDaughterWithEnergyFit) fC", fC);
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
// \param[in] Daughter - the daughter particle
void KFParticle::AddDaughterWithEnergyFitMC(float bz, const KFParticle& Daughter) {

    float m[8];
    float mV[36];

    float D[3][3];
    GetMeasurement(bz, Daughter, m, mV, D);

    float mS[6]{fC[0] + mV[0], fC[1] + mV[1], fC[2] + mV[2], fC[3] + mV[3], fC[4] + mV[4], fC[5] + mV[5]};
    InvertCholetsky3(mS);

    // Residual (measured - estimated)
    float zeta[3]{m[0] - fP[0], m[1] - fP[1], m[2] - fP[2]};

    float K[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            K[i][j] = 0.;
            for (int k = 0; k < 3; ++k) K[i][j] += fC[IJ(i, k)] * mS[IJ(k, j)];
        }
    }

    // CHt = CH'

    float mCHt0[7];
    float mCHt1[7];
    float mCHt2[7];

    mCHt0[0] = fC[0];
    mCHt1[0] = fC[1];
    mCHt2[0] = fC[3];
    mCHt0[1] = fC[1];
    mCHt1[1] = fC[2];
    mCHt2[1] = fC[4];
    mCHt0[2] = fC[3];
    mCHt1[2] = fC[4];
    mCHt2[2] = fC[5];
    mCHt0[3] = fC[6];
    mCHt1[3] = fC[7];
    mCHt2[3] = fC[8];
    mCHt0[4] = fC[10];
    mCHt1[4] = fC[11];
    mCHt2[4] = fC[12];
    mCHt0[5] = fC[15];
    mCHt1[5] = fC[16];
    mCHt2[5] = fC[17];
    mCHt0[6] = fC[21];
    mCHt1[6] = fC[22];
    mCHt2[6] = fC[23];

    // Kalman gain K = mCH'*S

    float k0[7];
    float k1[7];
    float k2[7];

    for (int i = 0; i < 7; ++i) {
        k0[i] = mCHt0[i] * mS[0] + mCHt1[i] * mS[1] + mCHt2[i] * mS[3];
        k1[i] = mCHt0[i] * mS[1] + mCHt1[i] * mS[2] + mCHt2[i] * mS[4];
        k2[i] = mCHt0[i] * mS[3] + mCHt1[i] * mS[4] + mCHt2[i] * mS[5];
    }

    // last itearation -> update the particle

    // VHt = VH'

    float mVHt0[7];
    float mVHt1[7];
    float mVHt2[7];

    mVHt0[0] = mV[0];
    mVHt1[0] = mV[1];
    mVHt2[0] = mV[3];
    mVHt0[1] = mV[1];
    mVHt1[1] = mV[2];
    mVHt2[1] = mV[4];
    mVHt0[2] = mV[3];
    mVHt1[2] = mV[4];
    mVHt2[2] = mV[5];
    mVHt0[3] = mV[6];
    mVHt1[3] = mV[7];
    mVHt2[3] = mV[8];
    mVHt0[4] = mV[10];
    mVHt1[4] = mV[11];
    mVHt2[4] = mV[12];
    mVHt0[5] = mV[15];
    mVHt1[5] = mV[16];
    mVHt2[5] = mV[17];
    mVHt0[6] = mV[21];
    mVHt1[6] = mV[22];
    mVHt2[6] = mV[23];

    // Kalman gain Km = mCH'*S

    float km0[7];
    float km1[7];
    float km2[7];

    for (int i = 0; i < 7; ++i) {
        km0[i] = mVHt0[i] * mS[0] + mVHt1[i] * mS[1] + mVHt2[i] * mS[3];
        km1[i] = mVHt0[i] * mS[1] + mVHt1[i] * mS[2] + mVHt2[i] * mS[4];
        km2[i] = mVHt0[i] * mS[3] + mVHt1[i] * mS[4] + mVHt2[i] * mS[5];
    }

    for (int i = 0; i < 7; ++i) fP[i] = fP[i] + k0[i] * zeta[0] + k1[i] * zeta[1] + k2[i] * zeta[2];

    for (int i = 0; i < 7; ++i) m[i] = m[i] - km0[i] * zeta[0] - km1[i] * zeta[1] - km2[i] * zeta[2];

    for (int i = 0, k = 0; i < 7; ++i) {
        for (int j = 0; j <= i; ++j, ++k) {
            fC[k] = fC[k] - (k0[i] * mCHt0[j] + k1[i] * mCHt1[j] + k2[i] * mCHt2[j]);
        }
    }

    for (int i = 0, k = 0; i < 7; ++i) {
        for (int j = 0; j <= i; ++j, ++k) {
            mV[k] = mV[k] - (km0[i] * mVHt0[j] + km1[i] * mVHt1[j] + km2[i] * mVHt2[j]);
        }
    }

    float mDf[7][7];

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            mDf[i][j] = km0[i] * mCHt0[j] + km1[i] * mCHt1[j] + km2[i] * mCHt2[j];
        }
    }

    float mJ1[7][7];
    float mJ2[7][7];
    for (int iPar1 = 0; iPar1 < 7; ++iPar1) {
        for (int iPar2 = 0; iPar2 < 7; ++iPar2) {
            mJ1[iPar1][iPar2] = 0.;
            mJ2[iPar1][iPar2] = 0.;
        }
    }

    float mMassParticle = fP[6] * fP[6] - (fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]);
    float mMassDaughter = m[6] * m[6] - (m[3] * m[3] + m[4] * m[4] + m[5] * m[5]);
    if (mMassParticle > 0) mMassParticle = std::sqrt(mMassParticle);
    if (mMassDaughter > 0) mMassDaughter = std::sqrt(mMassDaughter);

    if (fMassHypo > -0.5)
        SetMassConstraint(fP, fC, mJ1, fMassHypo);
    else if ((mMassParticle < fSumDaughterMass) || (fP[6] < 0))
        SetMassConstraint(fP, fC, mJ1, fSumDaughterMass);

    if (Daughter.fMassHypo > -0.5)
        SetMassConstraint(m, mV, mJ2, Daughter.fMassHypo);
    else if ((mMassDaughter < Daughter.fSumDaughterMass) || (m[6] < 0))
        SetMassConstraint(m, mV, mJ2, Daughter.fSumDaughterMass);

    float mDJ[7][7];
    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            mDJ[i][j] = 0.;
            for (int k = 0; k < 7; ++k) {
                mDJ[i][j] += mDf[i][k] * mJ1[j][k];
            }
        }
    }

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            mDf[i][j] = 0.;
            for (int l = 0; l < 7; ++l) {
                mDf[i][j] += mJ2[i][l] * mDJ[l][j];
            }
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

    float K2[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) K2[i][j] = -K[j][i];
        K2[i][i] += 1.;
    }

    float A[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            A[i][j] = 0.;
            for (int k = 0; k < 3; ++k) {
                A[i][j] += D[i][k] * K2[k][j];
            }
        }

    double M[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            M[i][j] = 0.;
            for (int k = 0; k < 3; ++k) {
                M[i][j] += K[i][k] * A[k][j];
            }
        }

    fC[0] += 2 * M[0][0];
    fC[1] += M[0][1] + M[1][0];
    fC[2] += 2 * M[1][1];
    fC[3] += M[0][2] + M[2][0];
    fC[4] += M[1][2] + M[2][1];
    fC[5] += 2 * M[2][2];

    // Calculate Chi^2

    fNDF += 2;
    fQ += Daughter.GetQ();
    fSFromDecay = 0;
    fChi2 += (mS[0] * zeta[0] + mS[1] * zeta[1] + mS[3] * zeta[2]) * zeta[0] + (mS[1] * zeta[0] + mS[2] * zeta[1] + mS[4] * zeta[2]) * zeta[1] +
             (mS[3] * zeta[0] + mS[4] * zeta[1] + mS[5] * zeta[2]) * zeta[2];
}

// Adds a vertex as a point-like measurement to the current particle.
// The eights parameter of the state vector is filled with the decay
// length to the momentum ratio (s = l/p). The corresponding covariances
// are calculated as well. The parameters of the particle are stored
// at the position of the production vertex.
// \param[in] Vtx - the assumed producation vertex
void KFParticle::SetProductionVertex(const KFParticle& vtx, float bz) {

    const float* m = vtx.fP;
    const float* mV = vtx.fC;

    float decayPoint[3]{fP[0], fP[1], fP[2]};
    float decayPointCov[6]{fC[0], fC[1], fC[2], fC[3], fC[4], fC[5]};

    float D[6][6];
    for (int iD1 = 0; iD1 < 6; ++iD1) {
        for (int iD2 = 0; iD2 < 6; ++iD2) D[iD1][iD2] = 0.;
    }

    bool noS = (fC[35] <= 0);  // no decay length allowed

    if (noS) {
        TransportToDecayVertex(bz);
        fP[7] = 0;
        fC[28] = fC[29] = fC[30] = fC[31] = fC[32] = fC[33] = fC[34] = fC[35] = 0;
    } else {
        float dsdr[6] = {0., 0., 0., 0., 0., 0.};
        float dS = GetDStoPointBz(bz, vtx.fP, dsdr);

        float dsdp[6] = {-dsdr[0], -dsdr[1], -dsdr[2], 0, 0, 0};

        float F[36], F1[36];
        for (int i2 = 0; i2 < 36; ++i2) {
            F[i2] = 0;
            F1[i2] = 0;
        }
        TransportBz(bz, dS, dsdr, fP, fC, dsdp, F, F1);

        float CTmp[36] = {0.};
        MultQSQt(F1, mV, CTmp, 6);

        for (int iC = 0; iC < 6; ++iC) fC[iC] += CTmp[iC];

        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 3; ++j) {
                D[i][j] = 0;
                for (int k = 0; k < 3; ++k) {
                    D[i][j] += mV[IJ(j, k)] * F1[i * 6 + k];
                }
            }
        }
    }

    float mS[6]{fC[0] + mV[0], fC[1] + mV[1], fC[2] + mV[2], fC[3] + mV[3], fC[4] + mV[4], fC[5] + mV[5]};
    InvertCholetsky3(mS);

    float res[3]{m[0] - X(), m[1] - Y(), m[2] - Z()};

    float K[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            K[i][j] = 0;
            for (int k = 0; k < 3; ++k) K[i][j] += fC[IJ(i, k)] * mS[IJ(k, j)];
        }

    float mCHt0[7];
    float mCHt1[7];
    float mCHt2[7];
    mCHt0[0] = fC[0];
    mCHt1[0] = fC[1];
    mCHt2[0] = fC[3];
    mCHt0[1] = fC[1];
    mCHt1[1] = fC[2];
    mCHt2[1] = fC[4];
    mCHt0[2] = fC[3];
    mCHt1[2] = fC[4];
    mCHt2[2] = fC[5];
    mCHt0[3] = fC[6];
    mCHt1[3] = fC[7];
    mCHt2[3] = fC[8];
    mCHt0[4] = fC[10];
    mCHt1[4] = fC[11];
    mCHt2[4] = fC[12];
    mCHt0[5] = fC[15];
    mCHt1[5] = fC[16];
    mCHt2[5] = fC[17];
    mCHt0[6] = fC[21];
    mCHt1[6] = fC[22];
    mCHt2[6] = fC[23];

    float k0[7];
    float k1[7];
    float k2[7];
    for (int i = 0; i < 7; ++i) {
        k0[i] = mCHt0[i] * mS[0] + mCHt1[i] * mS[1] + mCHt2[i] * mS[3];
        k1[i] = mCHt0[i] * mS[1] + mCHt1[i] * mS[2] + mCHt2[i] * mS[4];
        k2[i] = mCHt0[i] * mS[3] + mCHt1[i] * mS[4] + mCHt2[i] * mS[5];
    }

    for (int i = 0; i < 7; ++i) fP[i] = fP[i] + k0[i] * res[0] + k1[i] * res[1] + k2[i] * res[2];

    for (int i = 0, k = 0; i < 7; ++i) {
        for (int j = 0; j <= i; ++j, ++k) {
            fC[k] = fC[k] - (k0[i] * mCHt0[j] + k1[i] * mCHt1[j] + k2[i] * mCHt2[j]);
        }
    }

    float K2[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) K2[i][j] = -K[j][i];
        K2[i][i] += 1.;
    }

    float A[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[i][j] = 0;
            for (int k = 0; k < 3; ++k) {
                A[i][j] += D[k][i] * K2[k][j];
            }
        }
    }

    double M[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            M[i][j] = 0;
            for (int k = 0; k < 3; ++k) {
                M[i][j] += K[i][k] * A[k][j];
            }
        }
    }

    fC[0] += 2 * M[0][0];
    fC[1] += M[0][1] + M[1][0];
    fC[2] += 2 * M[1][1];
    fC[3] += M[0][2] + M[2][0];
    fC[4] += M[1][2] + M[2][1];
    fC[5] += 2 * M[2][2];

    fChi2 += (mS[0] * res[0] + mS[1] * res[1] + mS[3] * res[2]) * res[0] + (mS[1] * res[0] + mS[2] * res[1] + mS[4] * res[2]) * res[1] +
             (mS[3] * res[0] + mS[4] * res[1] + mS[5] * res[2]) * res[2];
    fNDF += 2;

    if (noS) {
        fP[7] = 0.;
        fC[28] = fC[29] = fC[30] = fC[31] = fC[32] = fC[33] = fC[34] = fC[35] = 0.;
        fSFromDecay = 0.;
    } else {
        float dsdr[6]{0., 0., 0., 0., 0., 0.};
        fP[7] = GetDStoPointBz(bz, decayPoint, dsdr);

        float dsdp[6]{-dsdr[0], -dsdr[1], -dsdr[2], 0., 0., 0.};

        float F[36];
        float F1[36];
        for (int i2 = 0; i2 < 36; ++i2) {
            F[i2] = 0.;
            F1[i2] = 0.;
        }
        float tmpP[8];
        float tmpC[36];
        TransportBz(bz, fP[7], dsdr, tmpP, tmpC, dsdp, F, F1);

        fC[35] = 0.;
        for (int iDsDr = 0; iDsDr < 6; ++iDsDr) {
            float dsdrC = 0.;
            float dsdpV = 0.;

            for (int k = 0; k < 6; ++k) dsdrC += dsdr[k] * fC[IJ(k, iDsDr)];  // (-dsdr[k])*fC[k,j]

            fC[iDsDr + 28] = dsdrC;
            fC[35] += dsdrC * dsdr[iDsDr];
            if (iDsDr < 3) {
                for (int k = 0; k < 3; ++k) dsdpV -= dsdr[k] * decayPointCov[IJ(k, iDsDr)];
                fC[35] -= dsdpV * dsdr[iDsDr];
            }
        }
        fSFromDecay = -fP[7];
    }

    fAtProductionVertex = true;
}

// set the exact nonlinear mass constraint on the state vector mP with the covariance matrix mC.
// \param[in,out] mP - the state vector to be modified
// \param[in,out] mC - the corresponding covariance matrix
// \param[in,out] mJ - the Jacobian between initial and modified parameters
// \param[in] mass - the mass to be set on the state vector mP
void KFParticle::SetMassConstraint(float* mP, float* mC, float mJ[7][7], float mass) {

    // Set nonlinear mass constraint (Mass) on the state vector mP with a covariance matrix mC.

    float energy2 = mP[6] * mP[6];
    float p2 = mP[3] * mP[3] + mP[4] * mP[4] + mP[5] * mP[5];
    float mass2 = mass * mass;

    float a = energy2 - p2 + 2. * mass2;
    float b = -2. * (energy2 + p2);
    float c = energy2 - p2 - mass2;

    float lambda = 0.;
    if (std::abs(b) > 1.e-10) lambda = -c / b;

    float d = 4. * energy2 * p2 - mass2 * (energy2 - p2 - 2. * mass2);
    if (d >= 0 && std::abs(a) > 1.e-10) lambda = (energy2 + p2 - std::sqrt(d)) / a;

    if (mP[6] < 0)           // If energy < 0 we need a lambda < 0
        lambda = -1000000.;  // Empirical, a better solution should be found

    int iIter = 0;
    for (iIter = 0; iIter < 100; ++iIter) {
        float lambda2 = lambda * lambda;
        float lambda4 = lambda2 * lambda2;

        float lambda0 = lambda;

        float f = -mass2 * lambda4 + a * lambda2 + b * lambda + c;
        float df = -4. * mass2 * lambda2 * lambda + 2. * a * lambda + b;
        if (std::abs(df) < 1.e-10) break;
        lambda -= f / df;
        if (std::abs(lambda0 - lambda) < 1.e-8) break;
    }

    float lpi = 1. / (1. + lambda);
    float lmi = 1. / (1. - lambda);
    float lp2i = lpi * lpi;
    float lm2i = lmi * lmi;

    float lambda2 = lambda * lambda;

    float dfl = -4. * mass2 * lambda2 * lambda + 2. * a * lambda + b;
    float dfx[7] = {0};  //,0,0,0};
    dfx[0] = -2. * (1. + lambda) * (1. + lambda) * mP[3];
    dfx[1] = -2. * (1. + lambda) * (1. + lambda) * mP[4];
    dfx[2] = -2. * (1. + lambda) * (1. + lambda) * mP[5];
    dfx[3] = 2. * (1. - lambda) * (1. - lambda) * mP[6];
    float dlx[4] = {1, 1, 1, 1};
    if (std::abs(dfl) > 1.e-10) {
        for (int i = 0; i < 4; ++i) dlx[i] = -dfx[i] / dfl;
    }

    float dxx[4]{mP[3] * lm2i, mP[4] * lm2i, mP[5] * lm2i, -mP[6] * lp2i};

    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 7; ++j) mJ[i][j] = 0.;
    mJ[0][0] = 1.;
    mJ[1][1] = 1.;
    mJ[2][2] = 1.;

    for (int i = 3; i < 7; ++i) {
        for (int j = 3; j < 7; ++j) mJ[i][j] = dlx[j - 3] * dxx[i - 3];
    }

    for (int i = 3; i < 6; ++i) mJ[i][i] += lmi;
    mJ[6][6] += lpi;

    float mCJ[7][7];

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            mCJ[i][j] = 0.;
            for (int k = 0; k < 7; ++k) {
                mCJ[i][j] += mC[IJ(i, k)] * mJ[j][k];
            }
        }
    }

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j <= i; ++j) {
            mC[IJ(i, j)] = 0.;
            for (int l = 0; l < 7; ++l) {
                mC[IJ(i, j)] += mJ[i][l] * mCJ[l][j];
            }
        }
    }

    mP[3] *= lmi;
    mP[4] *= lmi;
    mP[5] *= lmi;
    mP[6] *= lpi;
}

// set the exact nonlinear mass constraint on the current particle.
// \param[in] mass - the mass to be set on the particle
void KFParticle::SetNonlinearMassConstraint(float mass) {

    float px = fP[3];
    float py = fP[4];
    float pz = fP[5];
    float energy = fP[6];

    float residual = (energy * energy - px * px - py * py - pz * pz) - mass * mass;
    float dm2 =
        float(4.f) * (px * px * fC[9] + py * py * fC[14] + pz * pz * fC[20] + energy * energy * fC[27] +
                      float(2.f) * (px * py * fC[13] + pz * (px * fC[18] + py * fC[19]) - energy * (px * fC[24] + py * fC[25] + pz * fC[26])));
    float dChi2 = residual * residual / dm2;
    fChi2 += dChi2;
    fNDF += 1;

    float mJ[7][7];
    SetMassConstraint(fP, fC, mJ, mass);
    fMassHypo = mass;
    fSumDaughterMass = mass;
}

// Set linearised mass constraint on the current particle. The constraint can be set with an uncertainty.
// \param[in] Mass - the mass to be set on the state vector mP
// \param[in] SigmaMass - uncertainty of the constraint
void KFParticle::SetMassConstraint(float mass, float sigma_mass) {

    fMassHypo = mass;
    fSumDaughterMass = mass;

    float m2 = mass * mass;                   // measurement, weighted by Mass
    float s2 = m2 * sigma_mass * sigma_mass;  // sigma^2

    float p2 = fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5];
    float e0 = std::sqrt(m2 + p2);

    float mH[8]{0., 0., 0., 0., 0., 0., 0., 0.};
    mH[3] = -2 * fP[3];
    mH[4] = -2 * fP[4];
    mH[5] = -2 * fP[5];
    mH[6] = 2 * fP[6];  // e0;

    float zeta = e0 * e0 - e0 * fP[6];
    zeta = m2 - (fP[6] * fP[6] - p2);

    float mCHt[8];
    float s2_est = 0.;
    for (int i = 0; i < 8; ++i) {
        mCHt[i] = 0.;
        for (int j = 0; j < 8; ++j) mCHt[i] += Cij(i, j) * mH[j];
        s2_est += mH[i] * mCHt[i];
    }

    if (s2_est < 1.e-20)
        return;  // calculated mass error is already 0,
                 // the particle can not be constrained on mass

    float w2 = 1. / (s2 + s2_est);
    fChi2 += zeta * zeta * w2;
    fNDF += 1;
    for (int i = 0, ii = 0; i < 8; ++i) {
        float ki = mCHt[i] * w2;
        fP[i] += ki * zeta;
        for (int j = 0; j <= i; ++j) fC[++ii] -= ki * mCHt[j];
    }
}

// set constraint on the zero decay length. When the production point is set
// the measurement from this particle is created at the decay point.
void KFParticle::SetNoDecayLength(float bz) {

    TransportToDecayVertex(bz);

    float h[8];
    h[0] = h[1] = h[2] = h[3] = h[4] = h[5] = h[6] = 0;
    h[7] = 1;

    float zeta = 0. - fP[7];
    for (int i = 0; i < 8; ++i) zeta -= h[i] * (fP[i] - fP[i]);

    float s = fC[35];
    if (s > 1.e-20) {
        s = 1. / s;
        fChi2 += zeta * zeta * s;
        fNDF += 1;
        for (int i = 0, ii = 0; i < 7; ++i) {
            float ki = fC[28 + i] * s;
            fP[i] += ki * zeta;
            for (int j = 0; j <= i; ++j) fC[++ii] -= ki * fC[28 + j];
        }
    }
    fP[7] = 0;
    fC[28] = fC[29] = fC[30] = fC[31] = fC[32] = fC[33] = fC[34] = fC[35] = 0;
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
void KFParticle::Construct(float bz, const KFParticle* v_daughters[], int n_daughters, const KFParticle* parent, float mass) {

    fAtProductionVertex = false;
    fSFromDecay = 0;
    fSumDaughterMass = 0;

    for (int i = 0; i < 36; ++i) fC[i] = 0.;
    fC[35] = 1.;

    fNDF = -3;
    fChi2 = 0.;
    fQ = 0;

    for (int itr = 0; itr < n_daughters; ++itr) {
        AddDaughter(bz, *v_daughters[itr]);
    }

    if (mass >= 0) SetMassConstraint(mass);
    if (parent) SetProductionVertex(*parent, bz);
}

// Transports the particle to its decay vertex
void KFParticle::TransportToDecayVertex(float bz) {
    float dsdr[6] = {0.};
    if (fSFromDecay != 0) TransportToDS(bz, -fSFromDecay, dsdr);
    fAtProductionVertex = false;
}

// Transports the particle to its production vertex
void KFParticle::TransportToProductionVertex(float bz) {
    float dsdr[6] = {0.};
    if (fSFromDecay != -fP[7]) TransportToDS(bz, -fSFromDecay - fP[7], dsdr);
    fAtProductionVertex = true;
}

// Transport the particle on a certain distane. The distance is defined by the dS=l/p parameter, where
// 1) l - signed distance;
// 2) p - momentum of the particle.
// \param[in] dS = l/p - distance normalised to the momentum of the particle to be transported on
// \param[in] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
void KFParticle::TransportToDS(float bz, float ds, const float* dsdr) {
    TransportBz(bz, ds, dsdr, fP, fC);
    fSFromDecay += ds;
}

// return dS = l/p parameter, where
// 1) l - signed distance to the DCA point with the input xyz point;
// 2) p - momentum of the particle;
// assuming the straigth line trajectory. Is used for particles with charge 0 or in case of zero magnetic field.
// Also calculate partial derivatives dsdr of the parameter dS over the state vector of the current particle.
// \param[in] xyz[3] - point where particle should be transported
// \param[out] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
float KFParticle::GetDStoPointLine(const float xyz[3], float dsdr[6]) const {

    float p2 = fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5];
    if (p2 < 1.e-4) p2 = 1;

    float a = fP[3] * (xyz[0] - fP[0]) + fP[4] * (xyz[1] - fP[1]) + fP[5] * (xyz[2] - fP[2]);
    dsdr[0] = -fP[3] / p2;
    dsdr[1] = -fP[4] / p2;
    dsdr[2] = -fP[5] / p2;
    dsdr[3] = ((xyz[0] - fP[0]) * p2 - 2. * fP[3] * a) / (p2 * p2);
    dsdr[4] = ((xyz[1] - fP[1]) * p2 - 2. * fP[4] * a) / (p2 * p2);
    dsdr[5] = ((xyz[2] - fP[2]) * p2 - 2. * fP[5] * a) / (p2 * p2);

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
float KFParticle::GetDStoPointBz(float bz, const float xyz[3], float dsdr[6]) const {

    float dS(0.);

    float x = fP[0];
    float y = fP[1];
    float z = fP[2];
    float px = fP[3];
    float py = fP[4];
    float pz = fP[5];

    float bq = bz * fQ * kCLight;
    float pt2 = px * px + py * py;
    float p2 = pt2 + pz * pz;

    float dx = xyz[0] - x;
    float dy = xyz[1] - y;
    float dz = xyz[2] - z;
    float a = dx * px + dy * py;

    float abq = bq * a;

    bool mask = std::abs(bq) < LocalSmall;
    if (mask && p2 > 1.e-4f) {
        dS = (a + dz * pz) / p2;

        dsdr[0] = -px / p2;
        dsdr[1] = -py / p2;
        dsdr[2] = -pz / p2;
        dsdr[3] = (dx * p2 - 2. * px * (a + dz * pz)) / (p2 * p2);
        dsdr[4] = (dy * p2 - 2. * py * (a + dz * pz)) / (p2 * p2);
        dsdr[5] = (dz * p2 - 2. * pz * (a + dz * pz)) / (p2 * p2);
    }
    if (mask) return dS;

    dS = std::atan2(abq, pt2 + bq * (dy * px - dx * py)) / bq;

    float bs = bq * dS;

    float s = std::sin(bs);
    float c = std::cos(bs);

    if (std::abs(bq) < LocalSmall) bq = LocalSmall;
    float bbq = bq * (dx * py - dy * px) - pt2;

    float den = abq * abq + bbq * bbq;
    den = den < LocalSmall ? LocalSmall : den;

    dsdr[0] = (px * bbq - py * abq) / den;
    dsdr[1] = (px * abq + py * bbq) / den;
    dsdr[2] = 0.;
    dsdr[3] = -(dx * bbq + dy * abq + 2. * px * a) / den;
    dsdr[4] = (dx * abq - dy * bbq - 2. * py * a) / den;
    dsdr[5] = 0.;

    float sz{0.};
    float cCoeff = (bbq * c - abq * s) - pz * pz;
    if (std::abs(cCoeff) > LocalSmall) sz = (dS * pz - dz) * pz / cCoeff;

    float dcdr[6] = {0.};
    dcdr[0] = -bq * py * c - bbq * s * bq * dsdr[0] + px * bq * s - abq * c * bq * dsdr[0];
    dcdr[1] = bq * px * c - bbq * s * bq * dsdr[1] + py * bq * s - abq * c * bq * dsdr[1];
    dcdr[3] = (-bq * dy - 2 * px) * c - bbq * s * bq * dsdr[3] - dx * bq * s - abq * c * bq * dsdr[3];
    dcdr[4] = (bq * dx - 2 * py) * c - bbq * s * bq * dsdr[4] - dy * bq * s - abq * c * bq * dsdr[4];
    dcdr[5] = -2 * pz;

    for (int iP = 0; iP < 6; ++iP) dsdr[iP] += pz * pz / cCoeff * dsdr[iP] - sz / cCoeff * dcdr[iP];
    dsdr[2] += pz / cCoeff;
    dsdr[5] += (2. * pz * dS - dz) / cCoeff;

    dS += sz;

    /*
    bs = bq * dS;
    s = std::sin(bs);
    c = std::cos(bs);

    float sB = s / bq;
    float cB = (1. - c) / bq;

    float p[5];
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

    return dS;
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
void KFParticle::GetDStoParticleBz(float bz, const KFParticle& p, float dS[2], float dsdr[4][6]) const {

    // in XY plane first root
    float bq1 = bz * fQ * kCLight;
    float bq2 = bz * p.fQ * kCLight;

    bool isStraight1 = std::abs(bq1) < 1.E-8;
    bool isStraight2 = std::abs(bq2) < 1.E-8;

    if (isStraight1 && isStraight2) {
        GetDStoParticleLine(p, dS, dsdr);
        return;
    }

    float px1 = fP[3];
    float py1 = fP[4];
    float pz1 = fP[5];

    float px2 = p.fP[3];
    float py2 = p.fP[4];
    float pz2 = p.fP[5];

    float pt12 = px1 * px1 + py1 * py1;
    float pt22 = px2 * px2 + py2 * py2;

    float x01 = fP[0];
    float y01 = fP[1];
    float z01 = fP[2];

    float x02 = p.fP[0];
    float y02 = p.fP[1];
    float z02 = p.fP[2];

    float dS1[2]{0., 0.};
    float dS2[2]{0., 0.};

    float dx0 = x01 - x02;
    float dy0 = y01 - y02;
    float dr02 = dx0 * dx0 + dy0 * dy0;
    float drp1 = dx0 * px1 + dy0 * py1;
    float dxyp1 = dx0 * py1 - dy0 * px1;
    float drp2 = dx0 * px2 + dy0 * py2;
    float dxyp2 = dx0 * py2 - dy0 * px2;
    float p1p2 = px1 * px2 + py1 * py2;
    float dp1p2 = px1 * py2 - px2 * py1;

    float k11 = bq2 * drp1 - dp1p2;
    float k21 = bq1 * (bq2 * dxyp1 - p1p2) + bq2 * pt12;
    float k12 = bq1 * drp2 - dp1p2;
    float k22 = bq2 * (bq1 * dxyp2 + p1p2) - bq1 * pt22;

    float kp = dxyp1 * bq2 - dxyp2 * bq1 - p1p2;
    float kd = dr02 * bq1 * bq2 / 2. + kp;
    float c1 = -(bq1 * kd + pt12 * bq2);
    float c2 = bq2 * kd + pt22 * bq1;

    float d1 = std::max<float>(pt12 * pt22 - kd * kd, 0.);
    d1 = std::sqrt(d1);

    // find two points of closest approach in XY plane

    float dS1dR1[2][6];
    float dS2dR2[2][6];

    float dS1dR2[2][6];
    float dS2dR1[2][6];

    float dk11dr1[6]{bq2 * px1, bq2 * py1, 0, bq2 * dx0 - py2, bq2 * dy0 + px2, 0};
    float dk11dr2[6]{-bq2 * px1, -bq2 * py1, 0, py1, -px1, 0};
    float dk12dr1[6]{bq1 * px2, bq1 * py2, 0, -py2, px2, 0};
    float dk12dr2[6]{-bq1 * px2, -bq1 * py2, 0, bq1 * dx0 + py1, bq1 * dy0 - px1, 0};
    float dk21dr1[6]{bq1 * bq2 * py1, -bq1 * bq2 * px1, 0, 2 * bq2 * px1 + bq1 * (-(bq2 * dy0) - px2), 2 * bq2 * py1 + bq1 * (bq2 * dx0 - py2), 0};
    float dk21dr2[6]{-(bq1 * bq2 * py1), bq1 * bq2 * px1, 0, -(bq1 * px1), -(bq1 * py1), 0};
    float dk22dr1[6]{bq1 * bq2 * py2, -(bq1 * bq2 * px2), 0, bq2 * px2, bq2 * py2, 0};
    float dk22dr2[6]{-(bq1 * bq2 * py2), bq1 * bq2 * px2, 0, bq2 * (-(bq1 * dy0) + px1) - 2 * bq1 * px2, bq2 * (bq1 * dx0 + py1) - 2 * bq1 * py2, 0};

    float dkddr1[6]{bq1 * bq2 * dx0 + bq2 * py1 - bq1 * py2, bq1 * bq2 * dy0 - bq2 * px1 + bq1 * px2, 0, -bq2 * dy0 - px2, bq2 * dx0 - py2, 0};
    float dkddr2[6]{-bq1 * bq2 * dx0 - bq2 * py1 + bq1 * py2, -bq1 * bq2 * dy0 + bq2 * px1 - bq1 * px2, 0, bq1 * dy0 - px1, -bq1 * dx0 - py1, 0};

    float dc1dr1[6]{-(bq1 * (bq1 * bq2 * dx0 + bq2 * py1 - bq1 * py2)), -(bq1 * (bq1 * bq2 * dy0 - bq2 * px1 + bq1 * px2)), 0,
                    -2 * bq2 * px1 - bq1 * (-(bq2 * dy0) - px2),        -2 * bq2 * py1 - bq1 * (bq2 * dx0 - py2),           0};
    float dc1dr2[6]{-(bq1 * (-(bq1 * bq2 * dx0) - bq2 * py1 + bq1 * py2)),
                    -(bq1 * (-(bq1 * bq2 * dy0) + bq2 * px1 - bq1 * px2)),
                    0,
                    -(bq1 * (bq1 * dy0 - px1)),
                    -(bq1 * (-(bq1 * dx0) - py1)),
                    0};

    float dc2dr1[6]{bq2 * (bq1 * bq2 * dx0 + bq2 * py1 - bq1 * py2),
                    bq2 * (bq1 * bq2 * dy0 - bq2 * px1 + bq1 * px2),
                    0,
                    bq2 * (-(bq2 * dy0) - px2),
                    bq2 * (bq2 * dx0 - py2),
                    0};
    float dc2dr2[6]{bq2 * (-(bq1 * bq2 * dx0) - bq2 * py1 + bq1 * py2),
                    bq2 * (-(bq1 * bq2 * dy0) + bq2 * px1 - bq1 * px2),
                    0,
                    bq2 * (bq1 * dy0 - px1) + 2 * bq1 * px2,
                    bq2 * (-(bq1 * dx0) - py1) + 2 * bq1 * py2,
                    0};

    float dd1dr1[6]{0, 0, 0, 0, 0, 0};
    float dd1dr2[6]{0, 0, 0, 0, 0, 0};
    if (d1 > 0) {
        for (int i = 0; i < 6; ++i) {
            dd1dr1[i] = -kd / d1 * dkddr1[i];
            dd1dr2[i] = -kd / d1 * dkddr2[i];
        }
        dd1dr1[3] += px1 / d1 * pt22;
        dd1dr1[4] += py1 / d1 * pt22;
        dd1dr2[3] += px2 / d1 * pt12;
        dd1dr2[4] += py2 / d1 * pt12;
    }

    if (!isStraight1) {
        dS1[0] = std::atan2(bq1 * (k11 * c1 + k21 * d1), bq1 * k11 * d1 * bq1 - k21 * c1) / bq1;
        dS1[1] = std::atan2(bq1 * (k11 * c1 - k21 * d1), -bq1 * k11 * d1 * bq1 - k21 * c1) / bq1;

        float a = bq1 * (k11 * c1 + k21 * d1);
        float b = bq1 * k11 * d1 * bq1 - k21 * c1;
        for (int iP = 0; iP < 6; ++iP) {
            if ((b * b + a * a) > 0) {
                float dadr1 = bq1 * (dk11dr1[iP] * c1 + k11 * dc1dr1[iP] + dk21dr1[iP] * d1 + k21 * dd1dr1[iP]);
                float dadr2 = bq1 * (dk11dr2[iP] * c1 + k11 * dc1dr2[iP] + dk21dr2[iP] * d1 + k21 * dd1dr2[iP]);
                float dbdr1 = bq1 * bq1 * (dk11dr1[iP] * d1 + k11 * dd1dr1[iP]) - (dk21dr1[iP] * c1 + k21 * dc1dr1[iP]);
                float dbdr2 = bq1 * bq1 * (dk11dr2[iP] * d1 + k11 * dd1dr2[iP]) - (dk21dr2[iP] * c1 + k21 * dc1dr2[iP]);

                dS1dR1[0][iP] = 1. / bq1 * 1. / (b * b + a * a) * (dadr1 * b - dbdr1 * a);
                dS1dR2[0][iP] = 1. / bq1 * 1. / (b * b + a * a) * (dadr2 * b - dbdr2 * a);
            } else {
                dS1dR1[0][iP] = 0.;
                dS1dR2[0][iP] = 0.;
            }
        }

        a = bq1 * (k11 * c1 - k21 * d1);
        b = -bq1 * k11 * d1 * bq1 - k21 * c1;
        for (int iP = 0; iP < 6; ++iP) {
            if ((b * b + a * a) > 0) {
                float dadr1 = bq1 * (dk11dr1[iP] * c1 + k11 * dc1dr1[iP] - (dk21dr1[iP] * d1 + k21 * dd1dr1[iP]));
                float dadr2 = bq1 * (dk11dr2[iP] * c1 + k11 * dc1dr2[iP] - (dk21dr2[iP] * d1 + k21 * dd1dr2[iP]));
                float dbdr1 = -bq1 * bq1 * (dk11dr1[iP] * d1 + k11 * dd1dr1[iP]) - (dk21dr1[iP] * c1 + k21 * dc1dr1[iP]);
                float dbdr2 = -bq1 * bq1 * (dk11dr2[iP] * d1 + k11 * dd1dr2[iP]) - (dk21dr2[iP] * c1 + k21 * dc1dr2[iP]);

                dS1dR1[1][iP] = 1 / bq1 * 1 / (b * b + a * a) * (dadr1 * b - dbdr1 * a);
                dS1dR2[1][iP] = 1 / bq1 * 1 / (b * b + a * a) * (dadr2 * b - dbdr2 * a);
            } else {
                dS1dR1[1][iP] = 0;
                dS1dR2[1][iP] = 0;
            }
        }
    }
    if (!isStraight2) {
        dS2[0] = std::atan2(bq2 * k12 * c2 + k22 * d1 * bq2, (bq2 * k12 * d1 * bq2 - k22 * c2)) / bq2;
        dS2[1] = std::atan2(bq2 * k12 * c2 - k22 * d1 * bq2, (-bq2 * k12 * d1 * bq2 - k22 * c2)) / bq2;

        float a = bq2 * (k12 * c2 + k22 * d1);
        float b = bq2 * k12 * d1 * bq2 - k22 * c2;
        for (int iP = 0; iP < 6; ++iP) {
            if ((b * b + a * a) > 0) {
                float dadr1 = bq2 * (dk12dr1[iP] * c2 + k12 * dc2dr1[iP] + dk22dr1[iP] * d1 + k22 * dd1dr1[iP]);
                float dadr2 = bq2 * (dk12dr2[iP] * c2 + k12 * dc2dr2[iP] + dk22dr2[iP] * d1 + k22 * dd1dr2[iP]);
                float dbdr1 = bq2 * bq2 * (dk12dr1[iP] * d1 + k12 * dd1dr1[iP]) - (dk22dr1[iP] * c2 + k22 * dc2dr1[iP]);
                float dbdr2 = bq2 * bq2 * (dk12dr2[iP] * d1 + k12 * dd1dr2[iP]) - (dk22dr2[iP] * c2 + k22 * dc2dr2[iP]);

                dS2dR1[0][iP] = 1 / bq2 * 1 / (b * b + a * a) * (dadr1 * b - dbdr1 * a);
                dS2dR2[0][iP] = 1 / bq2 * 1 / (b * b + a * a) * (dadr2 * b - dbdr2 * a);
            } else {
                dS2dR1[0][iP] = 0;
                dS2dR2[0][iP] = 0;
            }
        }

        a = bq2 * (k12 * c2 - k22 * d1);
        b = -bq2 * k12 * d1 * bq2 - k22 * c2;
        for (int iP = 0; iP < 6; ++iP) {
            if ((b * b + a * a) > 0) {
                float dadr1 = bq2 * (dk12dr1[iP] * c2 + k12 * dc2dr1[iP] - (dk22dr1[iP] * d1 + k22 * dd1dr1[iP]));
                float dadr2 = bq2 * (dk12dr2[iP] * c2 + k12 * dc2dr2[iP] - (dk22dr2[iP] * d1 + k22 * dd1dr2[iP]));
                float dbdr1 = -bq2 * bq2 * (dk12dr1[iP] * d1 + k12 * dd1dr1[iP]) - (dk22dr1[iP] * c2 + k22 * dc2dr1[iP]);
                float dbdr2 = -bq2 * bq2 * (dk12dr2[iP] * d1 + k12 * dd1dr2[iP]) - (dk22dr2[iP] * c2 + k22 * dc2dr2[iP]);

                dS2dR1[1][iP] = 1 / bq2 * 1 / (b * b + a * a) * (dadr1 * b - dbdr1 * a);
                dS2dR2[1][iP] = 1 / bq2 * 1 / (b * b + a * a) * (dadr2 * b - dbdr2 * a);
            } else {
                dS2dR1[1][iP] = 0;
                dS2dR2[1][iP] = 0;
            }
        }
    }
    if (isStraight1 && pt12 > 0.) {
        dS1[0] = (k11 * c1 + k21 * d1) / (-k21 * c1);
        dS1[1] = (k11 * c1 - k21 * d1) / (-k21 * c1);

        float a = k11 * c1 + k21 * d1;
        float b = -k21 * c1;

        for (int iP = 0; iP < 6; ++iP) {
            if (b * b > 0) {
                float dadr1 = (dk11dr1[iP] * c1 + k11 * dc1dr1[iP] + dk21dr1[iP] * d1 + k21 * dd1dr1[iP]);
                float dadr2 = (dk11dr2[iP] * c1 + k11 * dc1dr2[iP] + dk21dr2[iP] * d1 + k21 * dd1dr2[iP]);
                float dbdr1 = -(dk21dr1[iP] * c1 + k21 * dc1dr1[iP]);
                float dbdr2 = -(dk21dr2[iP] * c1 + k21 * dc1dr2[iP]);

                dS1dR1[0][iP] = dadr1 / b - dbdr1 * a / (b * b);
                dS1dR2[0][iP] = dadr2 / b - dbdr2 * a / (b * b);
            } else {
                dS1dR1[0][iP] = 0;
                dS1dR2[0][iP] = 0;
            }
        }

        a = k11 * c1 - k21 * d1;
        for (int iP = 0; iP < 6; ++iP) {
            if (b * b > 0) {
                float dadr1 = (dk11dr1[iP] * c1 + k11 * dc1dr1[iP] - dk21dr1[iP] * d1 - k21 * dd1dr1[iP]);
                float dadr2 = (dk11dr2[iP] * c1 + k11 * dc1dr2[iP] - dk21dr2[iP] * d1 - k21 * dd1dr2[iP]);
                float dbdr1 = -(dk21dr1[iP] * c1 + k21 * dc1dr1[iP]);
                float dbdr2 = -(dk21dr2[iP] * c1 + k21 * dc1dr2[iP]);

                dS1dR1[1][iP] = dadr1 / b - dbdr1 * a / (b * b);
                dS1dR2[1][iP] = dadr2 / b - dbdr2 * a / (b * b);
            } else {
                dS1dR1[1][iP] = 0;
                dS1dR2[1][iP] = 0;
            }
        }
    }
    if (isStraight2 && pt22 > 0.) {
        dS2[0] = (k12 * c2 + k22 * d1) / (-k22 * c2);
        dS2[1] = (k12 * c2 - k22 * d1) / (-k22 * c2);

        float a = k12 * c2 + k22 * d1;
        float b = -k22 * c2;

        for (int iP = 0; iP < 6; ++iP) {
            if (b * b > 0) {
                float dadr1 = (dk12dr1[iP] * c2 + k12 * dc2dr1[iP] + dk22dr1[iP] * d1 + k22 * dd1dr1[iP]);
                float dadr2 = (dk12dr2[iP] * c2 + k12 * dc2dr2[iP] + dk22dr2[iP] * d1 + k22 * dd1dr2[iP]);
                float dbdr1 = -(dk22dr1[iP] * c2 + k22 * dc2dr1[iP]);
                float dbdr2 = -(dk22dr2[iP] * c2 + k22 * dc2dr2[iP]);

                dS2dR1[0][iP] = dadr1 / b - dbdr1 * a / (b * b);
                dS2dR2[0][iP] = dadr2 / b - dbdr2 * a / (b * b);
            } else {
                dS2dR1[0][iP] = 0;
                dS2dR2[0][iP] = 0;
            }
        }

        a = k12 * c2 - k22 * d1;
        for (int iP = 0; iP < 6; ++iP) {
            if (b * b > 0) {
                float dadr1 = (dk12dr1[iP] * c2 + k12 * dc2dr1[iP] - dk22dr1[iP] * d1 - k22 * dd1dr1[iP]);
                float dadr2 = (dk12dr2[iP] * c2 + k12 * dc2dr2[iP] - dk22dr2[iP] * d1 - k22 * dd1dr2[iP]);
                float dbdr1 = -(dk22dr1[iP] * c2 + k22 * dc2dr1[iP]);
                float dbdr2 = -(dk22dr2[iP] * c2 + k22 * dc2dr2[iP]);

                dS2dR1[1][iP] = dadr1 / b - dbdr1 * a / (b * b);
                dS2dR2[1][iP] = dadr2 / b - dbdr2 * a / (b * b);
            } else {
                dS2dR1[1][iP] = 0;
                dS2dR2[1][iP] = 0;
            }
        }
    }

    // select a point which is close to the primary vertex (with the smallest r)

    float dr2[2];
    for (int iP = 0; iP < 2; ++iP) {
        float bs1 = bq1 * dS1[iP];
        float bs2 = bq2 * dS2[iP];
        float sss = std::sin(bs1);
        float ccc = std::cos(bs1);

        float sB = sss / bq1;
        float cB = (1. - ccc) / bq1;

        float x1 = fP[0] + sB * px1 + cB * py1;
        float y1 = fP[1] - cB * px1 + sB * py1;
        float z1 = fP[2] + dS1[iP] * fP[5];

        sss = std::sin(bs2);
        ccc = std::cos(bs2);

        sB = sss / bq2;
        cB = (1. - ccc) / bq2;

        float x2 = p.fP[0] + sB * px2 + cB * py2;
        float y2 = p.fP[1] - cB * px2 + sB * py2;
        float z2 = p.fP[2] + dS2[iP] * p.fP[5];

        float dx = x1 - x2;
        float dy = y1 - y2;
        float dz = z1 - z2;

        dr2[iP] = dx * dx + dy * dy + dz * dz;
    }

    bool isFirstRoot = dr2[0] < dr2[1];
    if (isFirstRoot) {
        dS[0] = dS1[0];
        dS[1] = dS2[0];

        for (int iP = 0; iP < 6; ++iP) {
            dsdr[0][iP] = dS1dR1[0][iP];
            dsdr[1][iP] = dS1dR2[0][iP];
            dsdr[2][iP] = dS2dR1[0][iP];
            dsdr[3][iP] = dS2dR2[0][iP];
        }
    } else {
        dS[0] = dS1[1];
        dS[1] = dS2[1];

        for (int iP = 0; iP < 6; ++iP) {
            dsdr[0][iP] = dS1dR1[1][iP];
            dsdr[1][iP] = dS1dR2[1][iP];
            dsdr[2][iP] = dS2dR1[1][iP];
            dsdr[3][iP] = dS2dR2[1][iP];
        }
    }

    // Line correction
    float bs1 = bq1 * dS[0];
    float bs2 = bq2 * dS[1];
    float sss = std::sin(bs1);
    float ccc = std::cos(bs1);

    float sB = sss / bq1;
    float cB = (1. - ccc) / bq1;

    float x1 = x01 + sB * px1 + cB * py1;
    float y1 = y01 - cB * px1 + sB * py1;
    float z1 = z01 + dS[0] * pz1;
    float ppx1 = ccc * px1 + sss * py1;
    float ppy1 = -sss * px1 + ccc * py1;
    float ppz1 = pz1;

    float sss1 = std::sin(bs2);
    float ccc1 = std::cos(bs2);

    float sB1 = sss1 / bq2;
    float cB1 = (1. - ccc1) / bq2;

    float x2 = x02 + sB1 * px2 + cB1 * py2;
    float y2 = y02 - cB1 * px2 + sB1 * py2;
    float z2 = z02 + dS[1] * pz2;
    float ppx2 = ccc1 * px2 + sss1 * py2;
    float ppy2 = -sss1 * px2 + ccc1 * py2;
    float ppz2 = pz2;

    float p12 = ppx1 * ppx1 + ppy1 * ppy1 + ppz1 * ppz1;
    float p22 = ppx2 * ppx2 + ppy2 * ppy2 + ppz2 * ppz2;
    float lp1p2 = ppx1 * ppx2 + ppy1 * ppy2 + ppz1 * ppz2;

    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;

    float ldrp1 = ppx1 * dx + ppy1 * dy + ppz1 * dz;
    float ldrp2 = ppx2 * dx + ppy2 * dy + ppz2 * dz;

    float detp = lp1p2 * lp1p2 - p12 * p22;
    if (std::abs(detp) < 1.e-4) detp = 1;  // TODO correct!!!

    // dsdr calculation
    float a1 = ldrp2 * lp1p2 - ldrp1 * p22;
    float a2 = ldrp2 * p12 - ldrp1 * lp1p2;
    float lp1p2_ds0 = bq1 * (ppx2 * ppy1 - ppy2 * ppx1);
    float lp1p2_ds1 = bq2 * (ppx1 * ppy2 - ppy1 * ppx2);
    float ldrp1_ds0 = -p12 + bq1 * (ppy1 * dx - ppx1 * dy);
    float ldrp1_ds1 = lp1p2;
    float ldrp2_ds0 = -lp1p2;
    float ldrp2_ds1 = p22 + bq2 * (ppy2 * dx - ppx2 * dy);
    float detp_ds0 = 2. * lp1p2 * lp1p2_ds0;
    float detp_ds1 = 2. * lp1p2 * lp1p2_ds1;
    float a1_ds0 = ldrp2_ds0 * lp1p2 + ldrp2 * lp1p2_ds0 - ldrp1_ds0 * p22;
    float a1_ds1 = ldrp2_ds1 * lp1p2 + ldrp2 * lp1p2_ds1 - ldrp1_ds1 * p22;
    float a2_ds0 = ldrp2_ds0 * p12 - ldrp1_ds0 * lp1p2 - ldrp1 * lp1p2_ds0;
    float a2_ds1 = ldrp2_ds1 * p12 - ldrp1_ds1 * lp1p2 - ldrp1 * lp1p2_ds1;

    float dsl1ds0 = a1_ds0 / detp - a1 * detp_ds0 / (detp * detp);
    float dsl1ds1 = a1_ds1 / detp - a1 * detp_ds1 / (detp * detp);
    float dsl2ds0 = a2_ds0 / detp - a2 * detp_ds0 / (detp * detp);
    float dsl2ds1 = a2_ds1 / detp - a2 * detp_ds1 / (detp * detp);

    float dsldr[4][6];
    for (int iP = 0; iP < 6; ++iP) {
        dsldr[0][iP] = dsl1ds0 * dsdr[0][iP] + dsl1ds1 * dsdr[2][iP];
        dsldr[1][iP] = dsl1ds0 * dsdr[1][iP] + dsl1ds1 * dsdr[3][iP];
        dsldr[2][iP] = dsl2ds0 * dsdr[0][iP] + dsl2ds1 * dsdr[2][iP];
        dsldr[3][iP] = dsl2ds0 * dsdr[1][iP] + dsl2ds1 * dsdr[3][iP];
    }

    for (int iDS = 0; iDS < 4; ++iDS) {
        for (int iP = 0; iP < 6; ++iP) dsdr[iDS][iP] += dsldr[iDS][iP];
    }

    float lp1p2_dr0[6]{0., 0., 0., ccc * ppx2 - ppy2 * sss, ccc * ppy2 + ppx2 * sss, pz2};
    float lp1p2_dr1[6]{0., 0., 0., ccc1 * ppx1 - ppy1 * sss1, ccc1 * ppy1 + ppx1 * sss1, pz1};
    float ldrp1_dr0[6]{
        -ppx1, -ppy1, -pz1, cB * ppy1 - ppx1 * sB + ccc * dx - sss * dy, -cB * ppx1 - ppy1 * sB + sss * dx + ccc * dy, -dS[0] * pz1 + dz};
    float ldrp1_dr1[6]{ppx1, ppy1, pz1, -cB1 * ppy1 + ppx1 * sB1, cB1 * ppx1 + ppy1 * sB1, dS[1] * pz1};
    float ldrp2_dr0[6]{-ppx2, -ppy2, -pz2, cB * ppy2 - ppx2 * sB, -cB * ppx2 - ppy2 * sB, -dS[0] * pz2};
    float ldrp2_dr1[6]{
        ppx2, ppy2, pz2, -cB1 * ppy2 + ppx2 * sB1 + ccc1 * dx - sss1 * dy, cB1 * ppx2 + ppy2 * sB1 + sss1 * dx + ccc1 * dy, dz + dS[1] * pz2};
    float p12_dr0[6]{0., 0., 0., 2 * px1, 2 * py1, 2 * pz1};
    float p22_dr1[6]{0., 0., 0., 2 * px2, 2 * py2, 2 * pz2};
    float a1_dr0[6], a1_dr1[6], a2_dr0[6], a2_dr1[6], detp_dr0[6], detp_dr1[6];
    for (int iP = 0; iP < 6; ++iP) {
        a1_dr0[iP] = ldrp2_dr0[iP] * lp1p2 + ldrp2 * lp1p2_dr0[iP] - ldrp1_dr0[iP] * p22;
        a1_dr1[iP] = ldrp2_dr1[iP] * lp1p2 + ldrp2 * lp1p2_dr1[iP] - ldrp1_dr1[iP] * p22 - ldrp1 * p22_dr1[iP];
        a2_dr0[iP] = ldrp2_dr0[iP] * p12 + ldrp2 * p12_dr0[iP] - ldrp1_dr0[iP] * lp1p2 - ldrp1 * lp1p2_dr0[iP];
        a2_dr1[iP] = ldrp2_dr1[iP] * p12 - ldrp1_dr1[iP] * lp1p2 - ldrp1 * lp1p2_dr1[iP];
        detp_dr0[iP] = 2 * lp1p2 * lp1p2_dr0[iP] - p12_dr0[iP] * p22;
        detp_dr1[iP] = 2 * lp1p2 * lp1p2_dr1[iP] - p12 * p22_dr1[iP];

        dsdr[0][iP] += a1_dr0[iP] / detp - a1 * detp_dr0[iP] / (detp * detp);
        dsdr[1][iP] += a1_dr1[iP] / detp - a1 * detp_dr1[iP] / (detp * detp);
        dsdr[2][iP] += a2_dr0[iP] / detp - a2 * detp_dr0[iP] / (detp * detp);
        dsdr[3][iP] += a2_dr1[iP] / detp - a2 * detp_dr1[iP] / (detp * detp);
    }

    dS[0] += a1 / detp;
    dS[1] += a2 / detp;
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
void KFParticle::GetDStoParticleLine(const KFParticle& p, float ds[2], float dsdr[4][6]) const {

    float p12 = fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5];
    float p22 = p.fP[3] * p.fP[3] + p.fP[4] * p.fP[4] + p.fP[5] * p.fP[5];
    float p1p2 = fP[3] * p.fP[3] + fP[4] * p.fP[4] + fP[5] * p.fP[5];

    float drp1 = fP[3] * (p.fP[0] - fP[0]) + fP[4] * (p.fP[1] - fP[1]) + fP[5] * (p.fP[2] - fP[2]);
    float drp2 = p.fP[3] * (p.fP[0] - fP[0]) + p.fP[4] * (p.fP[1] - fP[1]) + p.fP[5] * (p.fP[2] - fP[2]);

    float detp = p1p2 * p1p2 - p12 * p22;
    if (std::abs(detp) < 1.e-4) detp = 1;  // TODO correct!!!

    ds[0] = (drp2 * p1p2 - drp1 * p22) / detp;
    ds[1] = (drp2 * p12 - drp1 * p1p2) / detp;

    float x01 = fP[0];
    float y01 = fP[1];
    float z01 = fP[2];
    float px1 = fP[3];
    float py1 = fP[4];
    float pz1 = fP[5];

    float x02 = p.fP[0];
    float y02 = p.fP[1];
    float z02 = p.fP[2];
    float px2 = p.fP[3];
    float py2 = p.fP[4];
    float pz2 = p.fP[5];

    float drp1_dr1[6]{-px1, -py1, -pz1, -x01 + x02, -y01 + y02, -z01 + z02};
    float drp1_dr2[6]{px1, py1, pz1, 0., 0., 0.};
    float drp2_dr1[6]{-px2, -py2, -pz2, 0., 0., 0.};
    float drp2_dr2[6]{px2, py2, pz2, -x01 + x02, -y01 + y02, -z01 + z02};
    float dp1p2_dr1[6]{0., 0., 0., px2, py2, pz2};
    float dp1p2_dr2[6]{0., 0., 0., px1, py1, pz1};
    float dp12_dr1[6]{0., 0., 0., 2 * px1, 2 * py1, 2 * pz1};
    float dp12_dr2[6]{0., 0., 0., 0., 0., 0.};
    float dp22_dr1[6]{0., 0., 0., 0., 0., 0.};
    float dp22_dr2[6]{0., 0., 0., 2 * px2, 2 * py2, 2 * pz2};
    float ddetp_dr1[6]{0., 0., 0., -2 * p22 * px1 + 2 * p1p2 * px2, -2 * p22 * py1 + 2 * p1p2 * py2, -2 * p22 * pz1 + 2 * p1p2 * pz2};
    float ddetp_dr2[6]{0., 0., 0., 2 * p1p2 * px1 - 2 * p12 * px2, 2 * p1p2 * py1 - 2 * p12 * py2, 2 * p1p2 * pz1 - 2 * p12 * pz2};

    float da1_dr1[6];
    float da1_dr2[6];
    float da2_dr1[6];
    float da2_dr2[6];

    float a1 = drp2 * p1p2 - drp1 * p22;
    float a2 = drp2 * p12 - drp1 * p1p2;
    for (int i = 0; i < 6; ++i) {
        da1_dr1[i] = drp2_dr1[i] * p1p2 + drp2 * dp1p2_dr1[i] - drp1_dr1[i] * p22 - drp1 * dp22_dr1[i];
        da1_dr2[i] = drp2_dr2[i] * p1p2 + drp2 * dp1p2_dr2[i] - drp1_dr2[i] * p22 - drp1 * dp22_dr2[i];

        da2_dr1[i] = drp2_dr1[i] * p12 + drp2 * dp12_dr1[i] - drp1_dr1[i] * p1p2 - drp1 * dp1p2_dr1[i];
        da2_dr2[i] = drp2_dr2[i] * p12 + drp2 * dp12_dr2[i] - drp1_dr2[i] * p1p2 - drp1 * dp1p2_dr2[i];

        dsdr[0][i] = da1_dr1[i] / detp - a1 / (detp * detp) * ddetp_dr1[i];
        dsdr[1][i] = da1_dr2[i] / detp - a1 / (detp * detp) * ddetp_dr2[i];

        dsdr[2][i] = da2_dr1[i] / detp - a2 / (detp * detp) * ddetp_dr1[i];
        dsdr[3][i] = da2_dr2[i] / detp - a2 / (detp * detp) * ddetp_dr2[i];
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
// \param[out] F[36] - optional parameter, transport jacobian, 6x6 matrix F = d(fP new)/d(fP old)
// \param[out] F1[36] - optional parameter, corelation 6x6 matrix betweeen the current particle and particle or vertex
// with the state vector r1, to which the current particle is being transported, F1 = d(fP new)/d(r1)
void KFParticle::TransportBz(float bz, float dS, const float* dsdr, float P[], float C[], float* dsdr1, float* F, float* F1) const {

    float bq = bz * fQ * kCLight;
    float bs = bq * dS;
    float s = std::sin(bs);
    float c = std::cos(bs);
    float sB = s / bq;
    float cB = (1 - c) / bq;

    float px = fP[3];
    float py = fP[4];
    float pz = fP[5];

    P[0] = fP[0] + sB * px + cB * py;
    P[1] = fP[1] - cB * px + sB * py;
    P[2] = fP[2] + dS * pz;
    P[3] = c * px + s * py;
    P[4] = -s * px + c * py;
    P[5] = fP[5];
    P[6] = fP[6];
    P[7] = fP[7];

    float mJ[8][8];
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) mJ[i][j] = 0.;

    for (int i = 0; i < 8; ++i) mJ[i][i] = 1;
    mJ[0][3] = sB;
    mJ[0][4] = cB;
    mJ[1][3] = -cB;
    mJ[1][4] = sB;
    mJ[2][5] = dS;
    mJ[3][3] = c;
    mJ[3][4] = s;
    mJ[4][3] = -s;
    mJ[4][4] = c;

    float mJds[6][6];
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) mJds[i][j] = 0.;
    mJds[0][3] = c;
    mJds[0][4] = s;
    mJds[1][3] = -s;
    mJds[1][4] = c;
    mJds[2][5] = 1;
    mJds[3][3] = -bq * s;
    mJds[3][4] = bq * c;
    mJds[4][3] = -bq * c;
    mJds[4][4] = -bq * s;

    for (int i1 = 0; i1 < 6; ++i1) {
        for (int i2 = 0; i2 < 6; ++i2) mJ[i1][i2] += mJds[i1][3] * px * dsdr[i2] + mJds[i1][4] * py * dsdr[i2] + mJds[i1][5] * pz * dsdr[i2];
    }

    MultQSQt(mJ[0], fC, C, 8);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) F[i * 6 + j] = mJ[i][j];
    }

    for (int i1 = 0; i1 < 6; ++i1) {
        for (int i2 = 0; i2 < 6; ++i2) F1[i1 * 6 + i2] = mJds[i1][3] * px * dsdr1[i2] + mJds[i1][4] * py * dsdr1[i2] + mJds[i1][5] * pz * dsdr1[i2];
    }
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
void KFParticle::TransportLine(float ds, const float* dsdr, float P[], float C[], float* dsdr1, float F[], float F1[]) const {

    float mJ[8][8];
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) mJ[i][j] = 0.;
    }

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

    float px = fP[3], py = fP[4], pz = fP[5];

    P[0] = fP[0] + ds * fP[3];
    P[1] = fP[1] + ds * fP[4];
    P[2] = fP[2] + ds * fP[5];
    P[3] = fP[3];
    P[4] = fP[4];
    P[5] = fP[5];
    P[6] = fP[6];
    P[7] = fP[7];

    float mJds[6][6];
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) mJds[i][j] = 0.;
    }
    mJds[0][3] = 1.;
    mJds[1][4] = 1.;
    mJds[2][5] = 1.;

    for (int i1 = 0; i1 < 6; ++i1) {
        for (int i2 = 0; i2 < 6; ++i2) {
            mJ[i1][i2] += mJds[i1][3] * px * dsdr[i2] + mJds[i1][4] * py * dsdr[i2] + mJds[i1][5] * pz * dsdr[i2];
        }
    }
    MultQSQt(mJ[0], fC, C, 8);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) F[i * 6 + j] = mJ[i][j];
    }
    for (int i1 = 0; i1 < 6; ++i1) {
        for (int i2 = 0; i2 < 6; ++i2) {
            F1[i1 * 6 + i2] = mJds[i1][3] * px * dsdr1[i2] + mJds[i1][4] * py * dsdr1[i2] + mJds[i1][5] * pz * dsdr1[i2];
        }
    }
}

//  symmetric 3x3 matrix a using modified Choletsky decomposition. The result is stored to the same matrix a.
// \param[in,out] a - 3x3 symmetric matrix
void KFParticle::InvertCholetsky3(float a[6]) {

    float d[3];
    float uud;
    float u[3][3];
    for (int i = 0; i < 3; ++i) {
        d[i] = 0.;
        for (int j = 0; j < 3; ++j) u[i][j] = 0.;
    }

    for (int i = 0; i < 3; ++i) {
        uud = 0.;
        for (int j = 0; j < i; ++j) uud += u[j][i] * u[j][i] * d[j];
        uud = a[i * (i + 3) / 2] - uud;

        if (std::abs(uud) < 1.e-12f) uud = 1.e-12f;

        d[i] = uud / std::abs(uud);
        u[i][i] = std::sqrt(std::abs(uud));

        for (int j = i + 1; j < 3; ++j) {
            uud = 0;
            for (int k = 0; k < i; ++k) uud += u[k][i] * u[k][j] * d[k];
            uud = a[j * (j + 1) / 2 + i] - uud;
            u[i][j] = d[i] / u[i][i] * uud;
        }
    }

    float u1[3];

    for (int i = 0; i < 3; ++i) {
        u1[i] = u[i][i];
        u[i][i] = 1 / u[i][i];
    }
    for (int i = 0; i < 2; ++i) {
        u[i][i + 1] = -u[i][i + 1] * u[i][i] * u[i + 1][i + 1];
    }
    for (int i = 0; i < 1; ++i) {
        u[i][i + 2] = u[i][i + 1] * u1[i + 1] * u[i + 1][i + 2] - u[i][i + 2] * u[i][i] * u[i + 2][i + 2];
    }

    for (int i = 0; i < 3; ++i) a[i + 3] = u[i][2] * u[2][2] * d[2];
    for (int i = 0; i < 2; ++i) a[i + 1] = u[i][1] * u[1][1] * d[1] + u[i][2] * u[1][2] * d[2];
    a[0] = u[0][0] * u[0][0] * d[0] + u[0][1] * u[0][1] * d[1] + u[0][2] * u[0][2] * d[2];
}

// Matrix multiplication SOut = Q*S*Q^T, where Q - square matrix, S - symmetric matrix.
// \param[in] Q - square matrix
// \param[in] S - input symmetric matrix
// \param[out] SOut - output symmetric matrix
// \param[in] kN - dimensionality of the matrices
void KFParticle::MultQSQt(const float Q[], const float S[], float SOut[], const int kN) {

    float* mA = new float[kN * kN];

    for (int i{0}, ij{0}; i < kN; ++i) {
        for (int j{0}; j < kN; ++j, ++ij) {
            mA[ij] = 0;
            for (int k{0}; k < kN; ++k) mA[ij] += S[(k <= i) ? i * (i + 1) / 2 + k : k * (k + 1) / 2 + i] * Q[j * kN + k];
        }
    }

    for (int i{0}; i < kN; ++i) {
        for (int j{0}; j <= i; ++j) {
            int ij = (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
            SOut[ij] = 0;
            for (int k{0}; k < kN; ++k) SOut[ij] += Q[i * kN + k] * mA[k * kN + j];
        }
    }

    if (mA) delete[] mA;
}
