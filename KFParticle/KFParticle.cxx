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

#include <algorithm>
#include <limits>
#include <tuple>

#include "KFParticle.hxx"
#include "KFParticle_Const.hxx"
#include "KFParticle_Math.hxx"
#if KF_DEBUG
#include "KFParticle_Utils.hxx"
#endif

namespace KF {

// Obtain the measurements from the current particle and the daughter to be added for the Kalman filter
// mathematics.
// If these are two first daughters they are transported to the point of the closest approach,
// if the third or higher daughter is added it is transported to the DCA point of the already constructed
// vertex. The correlations are taken into account in the covariance matrices of both measurements,
// the correlation matrix of two measurements is also calculated.
// Input arguments:
// - `daughter` : the daughter `KF::Particle` to be added
// - `bz`       : z-component of homogeneous magnetic field
// Return: (packed in a single `Result::Measurement` struct)
// - `P` : the output parameters of the daughter particle at the DCA point
// - `V` : the output covariance matrix of the daughter parameters, takes into account the correlation
// - `D` : the correlation matrix between the current and daughter particles
Result::Measurement Particle::GetMeasurement(const Particle& daughter, double bz) const {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    if (fNDF == -1) {
        // case: second daughter to be added //
        auto [min1, min2] = Minimize(daughter, bz);
        auto tpr1 = Transport(min1, bz);
        auto tpr2 = daughter.Transport(min2, bz);

#if KF_DEBUG
        Utils::Print(__FUNCTION__, "F1", tpr1.jacob);
        Utils::Print(__FUNCTION__, "F2", tpr1.corr);
        Utils::Print(__FUNCTION__, "F3", tpr2.corr);
        Utils::Print(__FUNCTION__, "F4", tpr2.jacob);
#endif

        SymMatrix<6> V0Tmp{Math::MultiplyQSQT(tpr1.corr, Math::Slice<8, 6>(daughter.fC))};
        SymMatrix<6> V1Tmp{Math::MultiplyQSQT(tpr2.corr, Math::Slice<8, 6>(fC))};
#if KF_DEBUG
        Utils::Print(__FUNCTION__, "V0Tmp", V0Tmp);
        Utils::Print(__FUNCTION__, "V1Tmp", V1Tmp);
#endif

        Result::Measurement meas{tpr1.P, tpr2.P, tpr1.C, tpr2.C};

        for (size_t iC{0}; iC < 21; ++iC) {
            meas.C1[iC] += V0Tmp[iC];
            meas.C2[iC] += V1Tmp[iC];
        }

        Matrix<6, 6> C1F1T{Math::MultiplySymmWithNonSymm(Math::Slice<8, 6>(fC), Math::Transpose(tpr1.jacob))};  // = C1 x F1^T

        Matrix<6, 6> F3C1F1T{Math::MultiplyMatrices(tpr2.corr, C1F1T)};  // = F3 x C1 x F1^T

        Matrix<6, 6> C2F2T{Math::MultiplySymmWithNonSymm(Math::Slice<8, 6>(daughter.fC), Math::Transpose(tpr1.corr))};  // = C2 x F2^T

        Matrix<6, 6> F4C2F2T{Math::MultiplyMatrices(tpr2.jacob, C2F2T)};  // = F4 x C2 x F2^T

        meas.D = Math::Slice<6, 6, 3, 3>(Math::AddMatrices(F3C1F1T, F4C2F2T));

#if KF_DEBUG
        Utils::Print(__FUNCTION__, "meas.P1", meas.P1);
        Utils::Print(__FUNCTION__, "meas.P2", meas.P2);
        Utils::Print(__FUNCTION__, "meas.C1", meas.C1);
        Utils::Print(__FUNCTION__, "meas.C2", meas.C2);
        Utils::Print(__FUNCTION__, "C1F1T", C1F1T);
        Utils::Print(__FUNCTION__, "F3C1F1T", F3C1F1T);
        Utils::Print(__FUNCTION__, "C2F2T", C2F2T);
        Utils::Print(__FUNCTION__, "meas.D", meas.D);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

        return meas;
    }

    // case: third (or above) daughter to be added //
    // >> transport to vertex //
    auto min2 = daughter.Minimize({fP[0], fP[1], fP[2]}, bz);
    auto tpr2 = daughter.Transport(min2, bz);
    Result::Measurement meas{fP, tpr2.P, fC, tpr2.C};

    Matrix<3, 6> VFT{Math::MultiplySymmWithNonSymm(              //
        Math::Slice<8, 3>(fC),                                   //
        Math::Transpose(Math::Slice<6, 6, 6, 3>(tpr2.jacob)))};  // = V x F^T

    Matrix<6, 6> FVFT{Math::MultiplyMatrices(Math::Slice<6, 6, 6, 3>(tpr2.jacob), VFT)};  // = F x V x F^T

    meas.D = Math::Transpose(           //
        Math::MultiplySymmWithNonSymm(  //
            Math::Slice<8, 3>(fC),      //
            Math::Transpose(Math::Slice<6, 6, 3, 3>(tpr2.jacob))));

    meas.C2[0] += FVFT[0][0];
    meas.C2[1] += FVFT[1][0];
    meas.C2[2] += FVFT[1][1];
    meas.C2[3] += FVFT[2][0];
    meas.C2[4] += FVFT[2][1];
    meas.C2[5] += FVFT[2][2];
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "VFT", VFT);
    Utils::Print(__FUNCTION__, "FVFT", FVFT);
    Utils::Print(__FUNCTION__, "meas.P1", meas.P1);
    Utils::Print(__FUNCTION__, "meas.P2", meas.P2);
    Utils::Print(__FUNCTION__, "meas.C1", meas.C1);
    Utils::Print(__FUNCTION__, "meas.C2", meas.C2);
    Utils::Print(__FUNCTION__, "meas.D", meas.D);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

    return meas;
}

// Add daughter to the current particle. Uses simplified fast mathematics which consideres momentum
// and energy as independent variables and ignores constraint on the fixed mass.
// In this case the mass of the daughter particle can be corrupted when the constructed vertex
// is added as the measurement and the mass of the output short-lived particle can become
// unphysical - smaller then the threshold.
// Input arguments:
// - `daughter`       : the daughter `KF::Particle` to be added
// - `bz`             : z-component of magnetic field
// - `chi2_threshold` : do an early cut of chi2
// Note: it will modify the state of the current `KF::Particle`
void Particle::AddDaughterWithEnergyFit(const Particle& daughter, double bz, double chi2_threshold) {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    auto meas = GetMeasurement(daughter, bz);
    fPCAs.emplace_back(meas.P1[0], meas.P1[1], meas.P1[2], meas.P1[3], meas.P1[4], meas.P1[5]);
    fPCAs.emplace_back(meas.P2[0], meas.P2[1], meas.P2[2], meas.P2[3], meas.P2[4], meas.P2[5]);

    SymMatrix<3> mS_in{Math::AddMatrices(Math::Slice<8, 3>(meas.C1), Math::Slice<8, 3>(meas.C2))};
    SymMatrix<3> mS{Math::InvertCholesky3(mS_in)};

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "meas.P1", meas.P1);
    Utils::Print(__FUNCTION__, "meas.C1", meas.C1);
    Utils::Print(__FUNCTION__, "meas.P2", meas.P2);
    Utils::Print(__FUNCTION__, "meas.C2", meas.C2);
    Utils::Print(__FUNCTION__, "meas.D", meas.D);
    Utils::Print(__FUNCTION__, "mS_in", mS_in);
    Utils::Print(__FUNCTION__, "mS", mS);
#endif

    Vector<3> zeta{meas.P2[0] - meas.P1[0], meas.P2[1] - meas.P1[1], meas.P2[2] - meas.P1[2]};
    double dChi2{(mS[0] * zeta[0] + mS[1] * zeta[1] + mS[3] * zeta[2]) * zeta[0] + (mS[1] * zeta[0] + mS[2] * zeta[1] + mS[4] * zeta[2]) * zeta[1] +
                 (mS[3] * zeta[0] + mS[4] * zeta[1] + mS[5] * zeta[2]) * zeta[2]};
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "zeta", zeta);
    Utils::PrintDouble(__FUNCTION__, "dChi2", dChi2);
#endif
    if (dChi2 > chi2_threshold) return;

    // update current particle state //

    fP = meas.P1;
    for (size_t i{0}; i < 28; ++i) fC[i] = meas.C1[i];

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
    Utils::Print(__FUNCTION__, "fP (before Kalman gain)", fP);
    Utils::Print(__FUNCTION__, "fC (before Kalman gain)", fC);
#endif

    // CHt = CH' - D'
    // Kalman gain K = mCH'*S
    // New estimation of the vertex position r += K*zeta
    // New covariance matrix C -= K*(mCH')'

    // Build mCHt as a 7x3 matrix (each row i, column j)
    Matrix<7, 3> mCHt{};
    mCHt[0][0] = meas.C1[0];
    mCHt[0][1] = meas.C1[1];
    mCHt[0][2] = meas.C1[3];
    mCHt[1][0] = meas.C1[1];
    mCHt[1][1] = meas.C1[2];
    mCHt[1][2] = meas.C1[4];
    mCHt[2][0] = meas.C1[3];
    mCHt[2][1] = meas.C1[4];
    mCHt[2][2] = meas.C1[5];
    mCHt[3][0] = meas.C1[6] - meas.C2[6];
    mCHt[3][1] = meas.C1[7] - meas.C2[7];
    mCHt[3][2] = meas.C1[8] - meas.C2[8];
    mCHt[4][0] = meas.C1[10] - meas.C2[10];
    mCHt[4][1] = meas.C1[11] - meas.C2[11];
    mCHt[4][2] = meas.C1[12] - meas.C2[12];
    mCHt[5][0] = meas.C1[15] - meas.C2[15];
    mCHt[5][1] = meas.C1[16] - meas.C2[16];
    mCHt[5][2] = meas.C1[17] - meas.C2[17];
    mCHt[6][0] = meas.C1[21] - meas.C2[21];
    mCHt[6][1] = meas.C1[22] - meas.C2[22];
    mCHt[6][2] = meas.C1[23] - meas.C2[23];

    // K (7x3) = mCHt (7x3) × S (3x3 symmetric)
    Matrix<7, 3> mK{Math::MultiplyNonSymmWithSymm(mCHt, mS)};

    // fP += K × zeta
    for (size_t i{0}; i < 7; ++i) {
        for (size_t j{0}; j < 3; ++j) {
            fP[i] += mK[i][j] * zeta[j];
        }
    }

    // fC -= K × mCHt^T  (only lower triangle)
    for (size_t i{0}, idx{0}; i < 7; ++i) {
        for (size_t j{0}; j <= i; ++j, ++idx) {
            for (size_t k{0}; k < 3; ++k) {
                fC[idx] -= mK[i][k] * mCHt[j][k];
            }
        }
    }

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "mCHt", mCHt);
    Utils::Print(__FUNCTION__, "mK", mK);
    Utils::Print(__FUNCTION__, "fP (after Kalman gain)", fP);
    Utils::Print(__FUNCTION__, "fC (after Kalman gain)", fC);
#endif

    // do something else? //

    Matrix<3, 3> K{Math::MultiplySymmetricMatrices(Math::Slice<8, 3>(meas.C1), mS)};

    Matrix<3, 3> K2{Math::AddMatrices(Math::Identity<3>(), Math::Transpose(K), 1., -1.)};

    Matrix<3, 3> A{Math::MultiplyMatrices(meas.D, K2)};

    Matrix<3, 3> M{Math::MultiplyMatrices(K, A)};

    fC[0] += 2. * M[0][0];
    fC[1] += M[0][1] + M[1][0];
    fC[2] += 2. * M[1][1];
    fC[3] += M[0][2] + M[2][0];
    fC[4] += M[1][2] + M[2][1];
    fC[5] += 2. * M[2][2];
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "K", K);
    Utils::Print(__FUNCTION__, "K2", K2);
    Utils::Print(__FUNCTION__, "A", A);
    Utils::Print(__FUNCTION__, "M", M);
    Utils::Print(__FUNCTION__, "fC (the end)", fC);
#endif

    // update rest of properties //

    fNDF += 2;
    fQ += daughter.Charge();
    fChi2 += dChi2;
#if KF_DEBUG
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
}

// Set a topological constraint on the current particle.
// Input arguments:
// - `prod_vtx`       : assumed production vertex
// - `prod_cov`       : production vertex's covariance matrix
// - `bz`             : z-component of magnetic field
// - `chi2_threshold` : do an early cut in chi2
// Note: it will modify the state of the current `KF::Particle`
// Note: should be executed as final step, after the particle has been added all of its daughters!
void Particle::AddProductionVertex(const Vector<3>& prod_vtx, const SymMatrix<3>& prod_cov, double bz, double chi2_threshold) {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    // store current particle's vertex information //

    Vector<3> decay_vtx{fP[0], fP[1], fP[2]};
    SymMatrix<3> decay_cov{fC[0], fC[1], fC[2], fC[3], fC[4], fC[5]};

    // transport to production vertex //

    auto min = Minimize(prod_vtx, bz);
    fPCAs.emplace_back(min.pca);

    auto tpr = Transport(min, bz);

    SymMatrix<3> CTmp{Math::MultiplyQSQT(Math::Slice<6, 6, 3, 3>(tpr.corr), prod_cov)};

    SymMatrix<3> measC{Math::AddMatrices(Math::Slice<8, 3>(tpr.C), CTmp)};

    Matrix<3, 3> D{Math::Transpose(Math::MultiplySymmWithNonSymm(prod_cov, Math::Transpose(Math::Slice<6, 6, 3, 3>(tpr.corr))))};

    SymMatrix<3> mS_in{Math::AddMatrices(measC, prod_cov)};
    SymMatrix<3> mS{Math::InvertCholesky3(mS_in)};
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "CTmp", CTmp);
    Utils::Print(__FUNCTION__, "measC", measC);
    Utils::Print(__FUNCTION__, "D", D);
    Utils::Print(__FUNCTION__, "mS_in", mS_in);
    Utils::Print(__FUNCTION__, "mS", mS);
#endif

    Vector<3> res{prod_vtx[0] - tpr.P[0], prod_vtx[1] - tpr.P[1], prod_vtx[2] - tpr.P[2]};
    double dChi2{(mS[0] * res[0] + mS[1] * res[1] + mS[3] * res[2]) * res[0] + (mS[1] * res[0] + mS[2] * res[1] + mS[4] * res[2]) * res[1] +
                 (mS[3] * res[0] + mS[4] * res[1] + mS[5] * res[2]) * res[2]};
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "res", res);
    Utils::PrintDouble(__FUNCTION__, "dChi2", dChi2);
#endif
    if (dChi2 > chi2_threshold) return;

    // update current particle state //

    fP = tpr.P;
    for (size_t i{0}; i < 6; ++i) fC[i] = measC[i];
    for (size_t i{6}; i < 28; ++i) fC[i] = tpr.C[i];

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "fP (before Kalman gain)", fP);
    Utils::Print(__FUNCTION__, "fC (before Kalman gain)", fC);
#endif

    // Kalman gain calculation //

    // extract columns 0, 1, 2 from all 7 rows of the symm. matrix fC
    Matrix<7, 3> mCHt{};
    mCHt[0][0] = fC[0];
    mCHt[0][1] = fC[1];
    mCHt[0][2] = fC[3];
    mCHt[1][0] = fC[1];
    mCHt[1][1] = fC[2];
    mCHt[1][2] = fC[4];
    mCHt[2][0] = fC[3];
    mCHt[2][1] = fC[4];
    mCHt[2][2] = fC[5];
    mCHt[3][0] = fC[6];
    mCHt[3][1] = fC[7];
    mCHt[3][2] = fC[8];
    mCHt[4][0] = fC[10];
    mCHt[4][1] = fC[11];
    mCHt[4][2] = fC[12];
    mCHt[5][0] = fC[15];
    mCHt[5][1] = fC[16];
    mCHt[5][2] = fC[17];
    mCHt[6][0] = fC[21];
    mCHt[6][1] = fC[22];
    mCHt[6][2] = fC[23];

    // mK (7x3) = mCHt (7x3) × mS (3x3 symmetric)
    Matrix<7, 3> mK{Math::MultiplyNonSymmWithSymm(mCHt, mS)};

    // fP += K × res
    for (size_t i{0}; i < 7; ++i) {
        for (size_t j{0}; j < 3; ++j) {
            fP[i] += mK[i][j] * res[j];
        }
    }

    // fC -= K × mCHt^T (only lower triangle)
    for (size_t i{0}, idx{0}; i < 7; ++i) {
        for (size_t j{0}; j <= i; ++j, ++idx) {
            for (size_t k{0}; k < 3; ++k) {
                fC[idx] -= mK[i][k] * mCHt[j][k];
            }
        }
    }

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "mCHt", mCHt);
    Utils::Print(__FUNCTION__, "mK", mK);
    Utils::Print(__FUNCTION__, "fP (after Kalman gain)", fP);
    Utils::Print(__FUNCTION__, "fC (after Kalman gain)", fC);
#endif

    Matrix<3, 3> K{Math::MultiplySymmetricMatrices(measC, mS)};

    Matrix<3, 3> K2{Math::AddMatrices(Math::Identity<3>(), Math::Transpose(K), 1., -1.)};

    Matrix<3, 3> A{Math::MultiplyMatrices(Math::Transpose(D), K2)};

    Matrix<3, 3> M{Math::MultiplyMatrices(K, A)};

    fC[0] += 2. * M[0][0];
    fC[1] += M[0][1] + M[1][0];
    fC[2] += 2. * M[1][1];
    fC[3] += M[0][2] + M[2][0];
    fC[4] += M[1][2] + M[2][1];
    fC[5] += 2. * M[2][2];
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "K", K);
    Utils::Print(__FUNCTION__, "K2", K2);
    Utils::Print(__FUNCTION__, "A", A);
    Utils::Print(__FUNCTION__, "M", M);
    Utils::Print(__FUNCTION__, "fC", fC);
#endif

    // update chi2 and ndf //

    fChi2 += dChi2;
    fNDF += 2;

    // update ds to decay vertex and update Css //

    auto min2decay = Minimize(decay_vtx, bz);
    fP[7] = min2decay.ds;
    fC[35] = 0.;

    for (size_t iDsDr{0}; iDsDr < 6; ++iDsDr) {
        double dsdrC{0.};
        double dsdpV{0.};

        for (size_t k{0}; k < 6; ++k) dsdrC += min2decay.ds_dr[k] * fC[IJ(k, iDsDr)];

        fC[iDsDr + 28] = dsdrC;
        fC[35] += dsdrC * min2decay.ds_dr[iDsDr];
        if (iDsDr < 3) {
            for (size_t k{0}; k < 3; ++k) dsdpV -= min2decay.ds_dr[k] * decay_cov[IJ(k, iDsDr)];
            fC[35] -= dsdpV * min2decay.ds_dr[iDsDr];
        }
    }
#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "fChi2", fChi2);
    Utils::Print(__FUNCTION__, "fNDF", fNDF);
    Utils::Print(__FUNCTION__, "fP (the end)", fP);
    Utils::Print(__FUNCTION__, "fC (the end)", fC);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
}

// Set a mass constraint on the current particle.
// Constraint equation g(...) : E^2 - (Px^2 + Py^2 + Pz^2) - target_mass^2 = 0
// Input argument:
// - `mass` : the mass to be set on the state vector mP
// Note: it will modify the state of the current `KF::Particle`
void Particle::AddMassConstraint(double target_mass) {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    double m2{target_mass * target_mass};
    double p2{fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]};

    // jacobian = d(g)/dr //
    Vector<8> mH{0., 0., 0., -2 * fP[3], -2 * fP[4], -2 * fP[5], 2 * fP[6], 0.};

    // residual = target_mass^2 - current_mass^2 //
    double zeta{m2 - fP[6] * fP[6] + p2};

    double s2{0.};
    Vector<8> mCHt{};
    for (size_t i{0}; i < 8; ++i) {
        for (size_t j{0}; j < 8; ++j) mCHt[i] += fC[IJ(i, j)] * mH[j];
        s2 += mH[i] * mCHt[i];
    }

    if (std::abs(s2) < Const::AbsAlmostZero) return;  // protection

    // apply Kalman filter update //
    for (size_t i{0}, ii{0}; i < 8; ++i) {
        // i-th component of Kalman gain vector //
        double ki{mCHt[i] / s2};
        // update state //
        fP[i] += ki * zeta;
        // update cov matrix //
        for (size_t j{0}; j <= i; ++j, ++ii) fC[ii] -= ki * mCHt[j];
    }

    fChi2 += zeta * zeta / s2;
    fNDF += 1;  // one d.o.f. is added because a single independent constraint has been applied
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "fP", fP);
    Utils::Print(__FUNCTION__, "fC", fC);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
}

// Find point of closest approach (PCA) of this particle w.r.t. an arbitrary vertex.
// Input:
// - `v` : arbitrary vertex
// Return: (packed in a single `Result::Minimization` struct)
// - `ds_dr` : partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
Result::Minimization Particle::MinimizeLinePoint(const Vector<3>& v) const {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    double x0{fP[0]};
    double y0{fP[1]};
    double z0{fP[2]};
    double px0{fP[3]};
    double py0{fP[4]};
    double pz0{fP[5]};

    double dx{v[0] - x0};
    double dy{v[1] - y0};
    double dz{v[2] - z0};

    double p2{px0 * px0 + py0 * py0 + pz0 * pz0};
    double a{px0 * dx + py0 * dy + pz0 * dz};

    Result::Minimization min{};
    if (p2 < Const::AbsAlmostZero) return min;

    min.ds = a / p2;

    min.ds_dr[0] = -px0 / p2;
    min.ds_dr[1] = -py0 / p2;
    min.ds_dr[2] = -pz0 / p2;
    min.ds_dr[3] = (dx * p2 - 2. * px0 * a) / (p2 * p2);
    min.ds_dr[4] = (dy * p2 - 2. * py0 * a) / (p2 * p2);
    min.ds_dr[5] = (dz * p2 - 2. * pz0 * a) / (p2 * p2);

    min.ds_dr1[0] = -min.ds_dr[0];
    min.ds_dr1[1] = -min.ds_dr[1];
    min.ds_dr1[2] = -min.ds_dr[2];

    min.pca.xyz[0] = x0 + px0 * min.ds;
    min.pca.xyz[1] = y0 + py0 * min.ds;
    min.pca.xyz[2] = z0 + pz0 * min.ds;

    min.pca.dir[0] = px0;
    min.pca.dir[1] = py0;
    min.pca.dir[2] = pz0;

#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "min.ds", min.ds);
    Utils::Print(__FUNCTION__, "min.ds_dr", min.ds_dr);
    Utils::Print(__FUNCTION__, "min.ds_dr1", min.ds_dr1);
    Utils::Print(__FUNCTION__, "min.(x,y,z)", min.pca.xyz);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
    return min;
}

// Find point of closest approach (PCA) of this particle w.r.t. an arbitrary vertex.
// Input arguments:
// - `v`  : arbirtrary vertex
// - `bz` : z-component of homogeneouse magnetic field
// Return: (packed in a single `Result::Minimization` struct)
// - `ds`     : transport parameters
// - `ds_dr`  : partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
// - `ds_dr1` : partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds2)/dr1, d(ds1)/dr2
// - cache properties : dir, pca, theta, sin, cos, sB, cB
Result::Minimization Particle::MinimizeHelixPoint(const Vector<3>& v, double bz) const {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif
    Result::Minimization min{};

    // 1 -- find point of closest approach (PCA) in XY plane //

    double bq{bz * fQ * Const::Kappa};

    double px0{fP[3]};
    double py0{fP[4]};
    double pz0{fP[5]};
    double pt2{px0 * px0 + py0 * py0};
    double x0{fP[0]};
    double y0{fP[1]};
    double z0{fP[2]};

    double dx{v[0] - x0};
    double dy{v[1] - y0};
    double dz{v[2] - z0};
    double a{dx * px0 + dy * py0};

    double abq{bq * a};
    double bbq{bq * (dx * py0 - dy * px0) - pt2};

    // 1.a -- get solution and update cache properties //

    min.theta = std::atan2(abq, -bbq);
    std::tie(min.sin, min.cos) = Math::sincos(min.theta);
    min.sB = min.sin / bq;
    min.cB = (1. - min.cos) / bq;

    min.ds = min.theta / bq;

    min.pca.xyz[0] = x0 + min.sB * px0 + min.cB * py0;
    min.pca.xyz[1] = y0 - min.cB * px0 + min.sB * py0;
    min.pca.xyz[2] = z0 + min.ds * pz0;

    min.pca.dir[0] = min.cos * px0 + min.sin * py0;
    min.pca.dir[1] = -min.sin * px0 + min.cos * py0;
    min.pca.dir[2] = pz0;
#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "min.ds (no z-correction)", min.ds);
    Utils::Print(__FUNCTION__, "min.(x,y,z) (no z-correction)", min.pca.xyz);
    // Utils::Print(__FUNCTION__, "dx", dx); // PENDING
    // Utils::Print(__FUNCTION__, "dy", dy); // PENDING
    // Utils::Print(__FUNCTION__, "dz", dz); // PENDING
    // Utils::Print(__FUNCTION__, "dca", std::sqrt(dca_sq)); // PENDING
#endif

    // 1.b -- handle derivatives //

    double den{abq * abq + bbq * bbq};
    den = den < Const::AbsAlmostZero ? Const::AbsAlmostZero : den;

    min.ds_dr[0] = (px0 * bbq - py0 * abq) / den;
    min.ds_dr[1] = (px0 * abq + py0 * bbq) / den;
    min.ds_dr[2] = 0.;
    min.ds_dr[3] = -(dx * bbq + dy * abq + 2. * px0 * a) / den;
    min.ds_dr[4] = (dx * abq - dy * bbq - 2. * py0 * a) / den;
    min.ds_dr[5] = 0.;

    min.ds_dr1[0] = -min.ds_dr[0];
    min.ds_dr1[1] = -min.ds_dr[1];
    min.ds_dr1[2] = -min.ds_dr[2];
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "min.ds_dr (no z-correction)", min.ds_dr);
    Utils::Print(__FUNCTION__, "min.ds_dr1 (no z-correction)", min.ds_dr1);
#endif

    // 2 -- add z-component as small correction //

    double cbq{bbq * min.cos - abq * min.sin - pz0 * pz0};
    if (std::abs(cbq) < Const::AbsAlmostZero) return min;  // protection

    double sz{(min.ds * pz0 - dz) * pz0 / cbq};

    // 2.a -- update derivatives //

    Vector<6> dc_dr{-bq * py0 * min.cos - bbq * min.sin * bq * min.ds_dr[0] + px0 * bq * min.sin - abq * min.cos * bq * min.ds_dr[0],
                    bq * px0 * min.cos - bbq * min.sin * bq * min.ds_dr[1] + py0 * bq * min.sin - abq * min.cos * bq * min.ds_dr[1],
                    0.,
                    (-bq * dy - 2. * px0) * min.cos - bbq * min.sin * bq * min.ds_dr[3] - dx * bq * min.sin - abq * min.cos * bq * min.ds_dr[3],
                    (bq * dx - 2. * py0) * min.cos - bbq * min.sin * bq * min.ds_dr[4] - dy * bq * min.sin - abq * min.cos * bq * min.ds_dr[4],
                    -2. * pz0};

    for (size_t iP{0}; iP < 6; ++iP) min.ds_dr[iP] += pz0 * pz0 * min.ds_dr[iP] / cbq - sz / cbq * dc_dr[iP];
    min.ds_dr[2] += pz0 / cbq;
    min.ds_dr[5] += (2. * pz0 * min.ds - dz) / cbq;

    min.ds_dr1[0] = -min.ds_dr[0];
    min.ds_dr1[1] = -min.ds_dr[1];
    min.ds_dr1[2] = -min.ds_dr[2];
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "min.ds_dr (after z-correction)", min.ds_dr);
    Utils::Print(__FUNCTION__, "min.ds_dr1 (after z-correction)", min.ds_dr1);
#endif

    // 2.b -- update ds //

    min.ds += sz;
#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "min.ds (after z-correction)", min.ds);
#endif

    // 2.c -- update rest of cache properties //

    min.theta = bq * min.ds;
    std::tie(min.sin, min.cos) = Math::sincos(min.theta);
    min.sB = min.sin / bq;
    min.cB = (1. - min.cos) / bq;

    min.pca.xyz[0] = x0 + min.sB * px0 + min.cB * py0;
    min.pca.xyz[1] = y0 - min.cB * px0 + min.sB * py0;
    min.pca.xyz[2] = z0 + min.ds * pz0;

    min.pca.dir[0] = min.cos * px0 + min.sin * py0;
    min.pca.dir[1] = -min.sin * px0 + min.cos * py0;
    min.pca.dir[2] = pz0;

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
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "min.(x,y,z) (after z-correction)", min.pca.xyz);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

    return min;
}

// Find point of closest approach (PCA) of this particle w.r.t. another one, assuming both or one of them have a helical trajectory.
// Input arguments:
// - `p`  : second particle
// - `bz` : z-component of homogeneous magnetic field
// Return: (packed as a pair of `Result::Minimization` structs)
// - `ds`     : transport parameters
// - `ds_dr`  : partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1, d(ds2)/dr2
// - `ds_dr1` : partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds2)/dr1, d(ds1)/dr2
std::pair<Result::Minimization, Result::Minimization> Particle::MinimizeHelixHelix(const Particle& p, double bz) const {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif
    Result::Minimization min1{};
    Result::Minimization min2{};

    // 1 -- find points of closest approach (PCAs) in XY plane //

    double bq1{bz * fQ * Const::Kappa};
    double bq2{bz * p.fQ * Const::Kappa};

    bool isStraight1{std::abs(bq1) < Const::AbsAlmostZero};
    bool isStraight2{std::abs(bq2) < Const::AbsAlmostZero};

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
        Cache tmp1{};
        if (!isStraight1) {
            tmp1.theta = std::atan2(bq1 * (k11 * c1 + sign * k21 * d1), sign * bq1 * k11 * d1 * bq1 - k21 * c1);
            std::tie(tmp1.sin, tmp1.cos) = Math::sincos(tmp1.theta);
            tmp1.sB = tmp1.sin / bq1;
            tmp1.cB = (1. - tmp1.cos) / bq1;
            tmp1.ds = tmp1.theta / bq1;

            tmp1.pca.xyz[0] = x01 + tmp1.sB * px01 + tmp1.cB * py01;
            tmp1.pca.xyz[1] = y01 - tmp1.cB * px01 + tmp1.sB * py01;
            tmp1.pca.xyz[2] = z01 + tmp1.ds * pz01;
            tmp1.pca.dir[0] = tmp1.cos * px01 + tmp1.sin * py01;
            tmp1.pca.dir[1] = -tmp1.sin * px01 + tmp1.cos * py01;
            tmp1.pca.dir[2] = pz01;
        } else {
            tmp1.ds = (k11 * c1 + sign * k21 * d1) / (-k21 * c1);
            tmp1.pca.xyz[0] = x01 + px01 * tmp1.ds;
            tmp1.pca.xyz[1] = y01 + py01 * tmp1.ds;
            tmp1.pca.xyz[2] = z01 + pz01 * tmp1.ds;
            tmp1.pca.dir[0] = px01;
            tmp1.pca.dir[1] = py01;
            tmp1.pca.dir[2] = pz01;
        }

        // particle 2 //
        Cache tmp2{};
        if (!isStraight2) {
            tmp2.theta = std::atan2(bq2 * (k12 * c2 + sign * k22 * d1), sign * bq2 * k12 * d1 * bq2 - k22 * c2);
            std::tie(tmp2.sin, tmp2.cos) = Math::sincos(tmp2.theta);
            tmp2.sB = tmp2.sin / bq2;
            tmp2.cB = (1. - tmp2.cos) / bq2;
            tmp2.ds = tmp2.theta / bq2;

            tmp2.pca.xyz[0] = x02 + tmp2.sB * px02 + tmp2.cB * py02;
            tmp2.pca.xyz[1] = y02 - tmp2.cB * px02 + tmp2.sB * py02;
            tmp2.pca.xyz[2] = z02 + tmp2.ds * pz02;
            tmp2.pca.dir[0] = tmp2.cos * px02 + tmp2.sin * py02;
            tmp2.pca.dir[1] = -tmp2.sin * px02 + tmp2.cos * py02;
            tmp2.pca.dir[2] = pz02;
        } else {
            tmp2.ds = (k12 * c2 + sign * k22 * d1) / (-k22 * c2);
            tmp2.pca.xyz[0] = x02 + px02 * tmp2.ds;
            tmp2.pca.xyz[1] = y02 + py02 * tmp2.ds;
            tmp2.pca.xyz[2] = z02 + pz02 * tmp2.ds;
            tmp2.pca.dir[0] = px02;
            tmp2.pca.dir[1] = py02;
            tmp2.pca.dir[2] = pz02;
        }

        Vector<3> tmp_diff{tmp2.pca.xyz};
        for (size_t i{0}; i < 3; ++i) tmp_diff[i] -= tmp1.pca.xyz[i];
        double tmp_dca_sq{Math::SquaredNorm(tmp_diff)};

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
    Utils::PrintDouble(__FUNCTION__, "min1.ds (no z-correction)", min1.ds);
    Utils::Print(__FUNCTION__, "min1.(x,y,z)", min1.pca.xyz);
    Utils::PrintDouble(__FUNCTION__, "min2.ds (no z-correction)", min2.ds);
    Utils::Print(__FUNCTION__, "min2.(x,y,z)", min2.pca.xyz);
    Utils::PrintDouble(__FUNCTION__, "dx", dx);
    Utils::PrintDouble(__FUNCTION__, "dy", dy);
    Utils::PrintDouble(__FUNCTION__, "dz", dz);
    Utils::Print(__FUNCTION__, "w_sign", w_sign);
    Utils::PrintDouble(__FUNCTION__, "dca", std::sqrt(dca_sq));
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
        for (size_t i{0}; i < 6; ++i) {
            dd1dr1[i] = -kd * dkddr1[i] / d1;
            dd1dr2[i] = -kd * dkddr2[i] / d1;
        }
        dd1dr1[3] += px01 * pt22 / d1;
        dd1dr1[4] += py01 * pt22 / d1;
        dd1dr2[3] += px02 * pt12 / d1;
        dd1dr2[4] += py02 * pt12 / d1;
    }
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "dk11dr1", dk11dr1);
    Utils::Print(__FUNCTION__, "dk11dr2", dk11dr2);
    Utils::Print(__FUNCTION__, "dk12dr1", dk12dr1);
    Utils::Print(__FUNCTION__, "dk12dr2", dk12dr2);
    Utils::Print(__FUNCTION__, "dk21dr1", dk21dr1);
    Utils::Print(__FUNCTION__, "dk21dr2", dk21dr2);
    Utils::Print(__FUNCTION__, "dk22dr1", dk22dr1);
    Utils::Print(__FUNCTION__, "dk22dr2", dk22dr2);
    Utils::Print(__FUNCTION__, "dkddr1", dkddr1);
    Utils::Print(__FUNCTION__, "dkddr2", dkddr2);
    Utils::Print(__FUNCTION__, "dc1dr1", dc1dr1);
    Utils::Print(__FUNCTION__, "dc1dr2", dc1dr2);
    Utils::Print(__FUNCTION__, "dc2dr1", dc2dr1);
    Utils::Print(__FUNCTION__, "dc2dr2", dc2dr2);
    Utils::Print(__FUNCTION__, "dd1dr1", dd1dr1);
    Utils::Print(__FUNCTION__, "dd1dr2", dd1dr2);
#endif

    if (!isStraight1) {
        double a{bq1 * (k11 * c1 + w_sign * k21 * d1)};
        double b{w_sign * bq1 * k11 * d1 * bq1 - k21 * c1};
        double c{b * b + a * a};
        double d{c > 0. ? (1. / bq1 * 1. / c) : 0.};

        for (size_t iP{0}; iP < 6; ++iP) {
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

        for (size_t iP{0}; iP < 6; ++iP) {
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

        for (size_t iP{0}; iP < 6; ++iP) {
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

        for (size_t iP{0}; iP < 6; ++iP) {
            double dadr1{dk12dr1[iP] * c2 + k12 * dc2dr1[iP] + w_sign * dk22dr1[iP] * d1 + w_sign * k22 * dd1dr1[iP]};
            double dadr2{dk12dr2[iP] * c2 + k12 * dc2dr2[iP] + w_sign * dk22dr2[iP] * d1 + w_sign * k22 * dd1dr2[iP]};
            double dbdr1{-dk22dr1[iP] * c2 - k22 * dc2dr1[iP]};
            double dbdr2{-dk22dr2[iP] * c2 - k22 * dc2dr2[iP]};

            min2.ds_dr1[iP] = (std::abs(b) > Const::AbsAlmostZero ? dadr1 / b : 0.) - (b * b > Const::AbsAlmostZero ? dbdr1 * a / (b * b) : 0.);
            min2.ds_dr[iP] = (std::abs(b) > Const::AbsAlmostZero ? dadr2 / b : 0.) - (b * b > Const::AbsAlmostZero ? dbdr2 * a / (b * b) : 0.);
        }
    }

    double px1{min1.pca.dir[0]};
    double py1{min1.pca.dir[1]};
    double px2{min2.pca.dir[0]};
    double py2{min2.pca.dir[1]};
#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "px1", px1);
    Utils::PrintDouble(__FUNCTION__, "py1", py1);
    Utils::PrintDouble(__FUNCTION__, "px2", px2);
    Utils::PrintDouble(__FUNCTION__, "py2", py2);
#endif

    // 2 -- add z-component as small correction //

    double p12{px1 * px1 + py1 * py1 + pz01 * pz01};
    double p22{px2 * px2 + py2 * py2 + pz02 * pz02};
    double lp1p2{px1 * px2 + py1 * py2 + pz01 * pz02};

    double detp{lp1p2 * lp1p2 - p12 * p22};  // protection
    if (std::abs(detp) < Const::AbsAlmostZero || detp * detp < Const::AbsAlmostZero) {
#if KF_DEBUG
        std::println(stdout, "({}) early return has been called!", __FUNCTION__);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
        return {min1, min2};
    }

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
    Utils::PrintDouble(__FUNCTION__, "a1", a1);
    Utils::PrintDouble(__FUNCTION__, "a2", a2);
    Utils::PrintDouble(__FUNCTION__, "lp1p2_ds0", lp1p2_ds0);
    Utils::PrintDouble(__FUNCTION__, "lp1p2_ds1", lp1p2_ds1);
    Utils::PrintDouble(__FUNCTION__, "ldrp1_ds0", ldrp1_ds0);
    Utils::PrintDouble(__FUNCTION__, "ldrp1_ds1", ldrp1_ds1);
    Utils::PrintDouble(__FUNCTION__, "ldrp2_ds0", ldrp2_ds0);
    Utils::PrintDouble(__FUNCTION__, "ldrp2_ds1", ldrp2_ds1);
    Utils::PrintDouble(__FUNCTION__, "detp_ds0", detp_ds0);
    Utils::PrintDouble(__FUNCTION__, "detp_ds1", detp_ds1);
    Utils::PrintDouble(__FUNCTION__, "a1_ds0", a1_ds0);
    Utils::PrintDouble(__FUNCTION__, "a1_ds1", a1_ds1);
    Utils::PrintDouble(__FUNCTION__, "a2_ds0", a2_ds0);
    Utils::PrintDouble(__FUNCTION__, "a2_ds1", a2_ds1);
    Utils::PrintDouble(__FUNCTION__, "dsl1ds0", dsl1ds0);
    Utils::PrintDouble(__FUNCTION__, "dsl1ds1", dsl1ds1);
    Utils::PrintDouble(__FUNCTION__, "dsl2ds0", dsl2ds0);
    Utils::PrintDouble(__FUNCTION__, "dsl2ds1", dsl2ds1);
#endif

    Vector<6> dsldr0{};
    Vector<6> dsldr1{};
    Vector<6> dsldr2{};
    Vector<6> dsldr3{};
    for (size_t iP{0}; iP < 6; ++iP) {
        dsldr0[iP] = dsl1ds0 * min1.ds_dr[iP] + dsl1ds1 * min2.ds_dr1[iP];
        dsldr1[iP] = dsl1ds0 * min1.ds_dr1[iP] + dsl1ds1 * min2.ds_dr[iP];
        dsldr2[iP] = dsl2ds0 * min1.ds_dr[iP] + dsl2ds1 * min2.ds_dr1[iP];
        dsldr3[iP] = dsl2ds0 * min1.ds_dr1[iP] + dsl2ds1 * min2.ds_dr[iP];
    }

    for (size_t iP{0}; iP < 6; ++iP) {
        min1.ds_dr[iP] += dsldr0[iP];
        min1.ds_dr1[iP] += dsldr1[iP];
        min2.ds_dr1[iP] += dsldr2[iP];
        min2.ds_dr[iP] += dsldr3[iP];
    }
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "min1.ds_dr (after z-correction 1)", min1.ds_dr);
    Utils::Print(__FUNCTION__, "min1.ds_dr1 (after z-correction 1)", min1.ds_dr1);
    Utils::Print(__FUNCTION__, "min2.ds_dr1 (after z-correction 1)", min2.ds_dr1);
    Utils::Print(__FUNCTION__, "min2.ds_dr (after z-correction 1)", min2.ds_dr);
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
    Utils::Print(__FUNCTION__, "lp1p2_dr0", lp1p2_dr0);
    Utils::Print(__FUNCTION__, "lp1p2_dr1", lp1p2_dr1);
    Utils::Print(__FUNCTION__, "ldrp1_dr0", ldrp1_dr0);
    Utils::Print(__FUNCTION__, "ldrp1_dr1", ldrp1_dr1);
    Utils::Print(__FUNCTION__, "ldrp2_dr0", ldrp2_dr0);
    Utils::Print(__FUNCTION__, "ldrp2_dr1", ldrp2_dr1);
    Utils::Print(__FUNCTION__, "p12_dr0", p12_dr0);
    Utils::Print(__FUNCTION__, "p22_dr1", p22_dr1);
#endif

    for (size_t iP{0}; iP < 6; ++iP) {
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
    Utils::Print(__FUNCTION__, "min1.ds_dr (after z-correction 2)", min1.ds_dr);
    Utils::Print(__FUNCTION__, "min1.ds_dr1 (after z-correction 2)", min1.ds_dr1);
    Utils::Print(__FUNCTION__, "min2.ds_dr1 (after z-correction 2)", min2.ds_dr1);
    Utils::Print(__FUNCTION__, "min2.ds_dr (after z-correction 2)", min2.ds_dr);
#endif

    // 2.b -- update ds //

    min1.ds += a1 / detp;
    min2.ds += a2 / detp;
#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "min1.ds (after z-correction)", min1.ds);
    Utils::PrintDouble(__FUNCTION__, "min2.ds (after z-correction)", min2.ds);
#endif

    // 2.c -- update rest of cache properties //

    min1.theta = bq1 * min1.ds;
    std::tie(min1.sin, min1.cos) = Math::sincos(min1.theta);
    min1.sB = min1.sin / bq1;
    min1.cB = (1. - min1.cos) / bq1;

    min1.pca.xyz[0] = x01 + min1.sB * px01 + min1.cB * py01;
    min1.pca.xyz[1] = y01 - min1.cB * px01 + min1.sB * py01;
    min1.pca.xyz[2] = z01 + min1.ds * pz01;
    min1.pca.dir[0] = min1.cos * px01 + min1.sin * py01;
    min1.pca.dir[1] = -min1.sin * px01 + min1.cos * py01;
    min1.pca.dir[2] = pz01;

    min2.theta = bq2 * min2.ds;
    std::tie(min2.sin, min2.cos) = Math::sincos(min2.theta);
    min2.sB = min2.sin / bq2;
    min2.cB = (1. - min2.cos) / bq2;

    min2.pca.xyz[0] = x02 + min2.sB * px02 + min2.cB * py02;
    min2.pca.xyz[1] = y02 - min2.cB * px02 + min2.sB * py02;
    min2.pca.xyz[2] = z02 + min2.ds * pz02;
    min2.pca.dir[0] = min2.cos * px02 + min2.sin * py02;
    min2.pca.dir[1] = -min2.sin * px02 + min2.cos * py02;
    min2.pca.dir[2] = pz02;
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "min1.(x,y,z) (after z-correction)", min1.pca.xyz);
    Utils::Print(__FUNCTION__, "min2.(x,y,z) (after z-correction)", min2.pca.xyz);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

    return {min1, min2};
}

// Find point of closest approach (PCA) of this particle w.r.t. another one, assuming both of them have a straight line trajectory.
// Input arguments:
// - `p` : second particle
// Return: (packed as a pair of `Result::Minimization` structs)
// - `ds`     : transport parameters
// - `ds_dr`  : partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1, d(ds2)/dr2
// - `ds_dr1` : partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2, d(ds2)/dr1
std::pair<Result::Minimization, Result::Minimization> Particle::MinimizeLineLine(const Particle& p) const {
#if KF_DEBUG
    std::println(stdout, "-- started ({}) --", __FUNCTION__);
#endif
    Result::Minimization min1{};
    Result::Minimization min2{};

    // find points of closest approach (PCAs) in XY plane //

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
    if (std::abs(detp) < Const::AbsAlmostZero) return {MinimizeLinePoint({0., 0., 0.}), MinimizeLinePoint({0., 0., 0.})};  // protection

    min1.ds = (drp2 * p1p2 - drp1 * p22) / detp;
    min1.pca.xyz[0] = x01 + px01 * min1.ds;
    min1.pca.xyz[1] = y01 + py01 * min1.ds;
    min1.pca.xyz[2] = z01 + pz01 * min1.ds;
    min1.pca.dir[0] = px01;
    min1.pca.dir[1] = py01;
    min1.pca.dir[2] = pz01;
    min2.ds = (drp2 * p12 - drp1 * p1p2) / detp;
    min2.pca.xyz[0] = x02 + px02 * min2.ds;
    min2.pca.xyz[1] = y02 + py02 * min2.ds;
    min2.pca.xyz[2] = z02 + pz02 * min2.ds;
    min2.pca.dir[0] = px02;
    min2.pca.dir[1] = py02;
    min2.pca.dir[2] = pz02;
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "min1.ds", min1.ds);
    Utils::Print(__FUNCTION__, "min1.(x,y,z)", min1.pca.xyz);
    Utils::Print(__FUNCTION__, "min2.ds", min2.ds);
    Utils::Print(__FUNCTION__, "min2.(x,y,z)", min2.pca.xyz);
#endif

    // handle derivatives //

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

    for (size_t i{0}; i < 6; ++i) {
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
    Utils::Print(__FUNCTION__, "min1.ds_dr", min1.ds_dr);
    Utils::Print(__FUNCTION__, "min1.ds_dr1", min1.ds_dr1);
    Utils::Print(__FUNCTION__, "min2.ds_dr1", min2.ds_dr1);
    Utils::Print(__FUNCTION__, "min2.ds_dr", min2.ds_dr);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

    return {min1, min2};
}

// Transport current particle assuming a helical trajectory (charged and under homogeneous magnetic field `bz`)
// Input arguments: (packed in a single `Result::Minimization` struct)
// - `ds`     : transport parameter
// - `ds_dr`  : partial derivatives of ds w.r.t. state parameters = d(ds1)/dr1
// - `ds_dr1` : partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2
// - `bz`     : z-component of homogeneous magnetic field
// Return: (packed in a single `Result::Transport` struct)
// - `P`     : where transported parameters should be stored
// - `C`     : where transported covariance matrix (8x8) should be stored in the lower triangular form
// - `jacob` : transport jacobian = d(fP new)/d(fP old)
// - `corr`  : correlation matrix = d(fP new)/d(r1)
Result::Transport Particle::TransportBz(const Result::Minimization& min, double bz) const {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif
    Result::Transport tpr{};

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

    Matrix<8, 8> mJ{Math::Identity<8>()};
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

    for (size_t i1{0}; i1 < 6; ++i1) {
        for (size_t i2{0}; i2 < 6; ++i2) {
            mJ[i1][i2] += mJds[i1][3] * px0 * min.ds_dr[i2] + mJds[i1][4] * py0 * min.ds_dr[i2] + mJds[i1][5] * pz0 * min.ds_dr[i2];
        }
    }

    tpr.C = Math::MultiplyQSQT(mJ, fC);
    tpr.jacob = Math::Slice<8, 8, 6, 6>(mJ);

    for (size_t i1{0}; i1 < 6; ++i1) {
        for (size_t i2{0}; i2 < 6; ++i2) {
            tpr.corr[i1][i2] = mJds[i1][3] * px0 * min.ds_dr1[i2] + mJds[i1][4] * py0 * min.ds_dr1[i2] + mJds[i1][5] * pz0 * min.ds_dr1[i2];
        }
    }
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "State", tpr.P);
    Utils::Print(__FUNCTION__, "Cov", tpr.C);
    Utils::Print(__FUNCTION__, "Jacob", tpr.jacob);
    Utils::Print(__FUNCTION__, "Corr", tpr.corr);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

    return tpr;
}

// Transport current particle assuming an straight line trajectory (neutral or charged without magnetic field)
// Input arguments: (packed in a single `Result::Minimization` struct)
// - `ds`     : transport parameter
// - `ds_dr`  : partial derivatives of ds w.r.t. state parameters = d(ds1)/dr1
// - `ds_dr1` : partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2
// Return: (packed in a single `Result::Transport` struct)
// - `P`     : new transported state
// - `C`     : new transported covariance matrix
// - `jacob` : transport jacobian = d(fP new)/d(fP old)
// - `corr`  : correlation matrix = d(fP new)/d(r1)
Result::Transport Particle::TransportLine(const Result::Minimization& min) const {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif
    Result::Transport tpr{};

    Matrix<8, 8> mJ{Math::Identity<8>()};
    mJ[0][3] = min.ds;
    mJ[1][4] = min.ds;
    mJ[2][5] = min.ds;

    double px{fP[3]};
    double py{fP[4]};
    double pz{fP[5]};

    tpr.P[0] = fP[0] + min.ds * px;
    tpr.P[1] = fP[1] + min.ds * py;
    tpr.P[2] = fP[2] + min.ds * pz;
    tpr.P[3] = px;
    tpr.P[4] = py;
    tpr.P[5] = pz;
    tpr.P[6] = fP[6];
    tpr.P[7] = fP[7];

    Matrix<6, 6> mJds{};
    mJds[0][3] = 1.;
    mJds[1][4] = 1.;
    mJds[2][5] = 1.;

    for (size_t i1{0}; i1 < 6; ++i1) {
        for (size_t i2{0}; i2 < 6; ++i2) {
            mJ[i1][i2] += mJds[i1][3] * px * min.ds_dr[i2] + mJds[i1][4] * py * min.ds_dr[i2] + mJds[i1][5] * pz * min.ds_dr[i2];
        }
    }

    tpr.C = Math::MultiplyQSQT(mJ, fC);
    tpr.jacob = Math::Slice<8, 8, 6, 6>(mJ);

    for (size_t i1{0}; i1 < 6; ++i1) {
        for (size_t i2{0}; i2 < 6; ++i2) {
            tpr.corr[i1][i2] = mJds[i1][3] * px * min.ds_dr1[i2] + mJds[i1][4] * py * min.ds_dr1[i2] + mJds[i1][5] * pz * min.ds_dr1[i2];
        }
    }
#if KF_DEBUG
    Utils::Print(__FUNCTION__, "State", tpr.P);
    Utils::Print(__FUNCTION__, "Cov", tpr.C);
    Utils::Print(__FUNCTION__, "Jacob", tpr.jacob);
    Utils::Print(__FUNCTION__, "Corr", tpr.corr);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

    return tpr;
}

}  // namespace KF
