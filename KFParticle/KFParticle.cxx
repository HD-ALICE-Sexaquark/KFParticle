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

#include <armadillo>

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

        arma::mat::fixed<6, 6> V0Tmp = tpr1.corr * daughter.fC.submat(0, 0, arma::size(6, 6)) * tpr1.corr.t();
        arma::mat::fixed<6, 6> V1Tmp = tpr2.corr * fC.submat(0, 0, arma::size(6, 6)) * tpr2.corr.t();

#if KF_DEBUG
        Utils::Print(__FUNCTION__, "V0Tmp", V0Tmp);
        Utils::Print(__FUNCTION__, "V1Tmp", V1Tmp);
#endif

        Result::Measurement meas{tpr1.P, tpr2.P, tpr1.C, tpr2.C};

        meas.C1.submat(0, 0, arma::size(6, 6)) += V0Tmp;
        meas.C2.submat(0, 0, arma::size(6, 6)) += V1Tmp;

        arma::mat::fixed<6, 6> fC_sliced = fC.submat(0, 0, arma::size(6, 6));
        arma::mat::fixed<6, 6> d_fC_sliced = daughter.fC.submat(0, 0, arma::size(6, 6));

        arma::mat::fixed<6, 6> F3C1F1T = tpr2.corr * fC_sliced * tpr1.jacob.t();    // = F3 x C1 x F1^T
        arma::mat::fixed<6, 6> F4C2F2T = tpr2.jacob * d_fC_sliced * tpr1.corr.t();  // = F4 x C2 x F2^T

        arma::mat::fixed<6, 6> temp_addition = F3C1F1T + F4C2F2T;
        meas.D = temp_addition.submat(0, 0, arma::size(3, 3));

#if KF_DEBUG
        Utils::Print(__FUNCTION__, "meas.P1", meas.P1);
        Utils::Print(__FUNCTION__, "meas.P2", meas.P2);
        Utils::Print(__FUNCTION__, "meas.C1", meas.C1);
        Utils::Print(__FUNCTION__, "meas.C2", meas.C2);
        // Utils::Print(__FUNCTION__, "C1F1T", C1F1T);
        Utils::Print(__FUNCTION__, "F3C1F1T", F3C1F1T);
        // Utils::Print(__FUNCTION__, "C2F2T", C2F2T);
        Utils::Print(__FUNCTION__, "meas.D", meas.D);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

        return meas;
    }

    // case: third (or above) daughter to be added //
    // >> transport to vertex //
    auto min2 = daughter.Minimize(fP.head(3), bz);
    auto tpr2 = daughter.Transport(min2, bz);
    Result::Measurement meas{fP, tpr2.P, fC, tpr2.C};

    arma::mat::fixed<3, 3> sliced_fC = fC.submat(0, 0, arma::size(3, 3));
    arma::mat::fixed<6, 3> sliced_F_6x3 = tpr2.jacob.submat(0, 0, arma::size(6, 3));
    arma::mat::fixed<3, 3> sliced_F_3x3 = tpr2.jacob.submat(0, 0, arma::size(3, 3));

    arma::mat::fixed<6, 6> FVFT = sliced_F_6x3 * sliced_fC * sliced_F_6x3.t();  // = F x V x F^T

    arma::mat::fixed<3, 3> tmp = sliced_fC * sliced_F_3x3.t();
    meas.D = tmp.t();

    meas.C2.submat(0, 0, arma::size(3, 3)) += FVFT.submat(0, 0, arma::size(3, 3));

#if KF_DEBUG
    // Utils::Print(__FUNCTION__, "VFT", VFT);
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
// Return: true if successful, false if not.
// Note: it will modify the state of the current `KF::Particle`
bool Particle::AddDaughterWithEnergyFit(const Particle& daughter, double bz, double chi2_threshold) {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    auto meas = GetMeasurement(daughter, bz);
    fPCAs.emplace_back(meas.P1(0), meas.P1(1), meas.P1(2), meas.P1(3), meas.P1(4), meas.P1(5));
    fPCAs.emplace_back(meas.P2(0), meas.P2(1), meas.P2(2), meas.P2(3), meas.P2(4), meas.P2(5));

    arma::mat::fixed<3, 3> mS_in = meas.C1.submat(0, 0, arma::size(3, 3)) + meas.C2.submat(0, 0, arma::size(3, 3));

    arma::mat::fixed<3, 3> mS(arma::fill::zeros);
    bool inv_success = arma::inv_sympd(mS, mS_in);
    if (!inv_success) {
#if KF_DEBUG
        std::println(stdout, "({}) failed matrix inversion!", __FUNCTION__);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
        return false;
    }

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "meas.P1", meas.P1);
    Utils::Print(__FUNCTION__, "meas.C1", meas.C1);
    Utils::Print(__FUNCTION__, "meas.P2", meas.P2);
    Utils::Print(__FUNCTION__, "meas.C2", meas.C2);
    Utils::Print(__FUNCTION__, "meas.D", meas.D);
    Utils::Print(__FUNCTION__, "mS_in", mS_in);
    Utils::Print(__FUNCTION__, "mS", mS);
#endif

    arma::vec::fixed<3> zeta = meas.P2.head(3) - meas.P1.head(3);
    double dChi2{arma::as_scalar(zeta.t() * mS * zeta)};

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "zeta", zeta);
    Utils::PrintDouble(__FUNCTION__, "dChi2", dChi2);
#endif

    // optimization
    if (dChi2 > chi2_threshold) {
#if KF_DEBUG
        std::println(stdout, "({}) early return has been called!", __FUNCTION__);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
        return false;
    }

    // update current particle state //

    fP = meas.P1;
    fC.submat(0, 0, arma::size(7, 7)) = meas.C1.submat(0, 0, arma::size(7, 7));

    // add daughter's momentum to particle momentum //

    fP.subvec(3, 6) += meas.P2.subvec(3, 6);
    fC.submat(3, 3, arma::size(4, 4)) += meas.C2.submat(3, 3, arma::size(4, 4));
    fC = arma::symmatl(fC);  // force symmetry

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "fP (before Kalman gain)", fP);
    Utils::Print(__FUNCTION__, "fC (before Kalman gain)", fC);
#endif

    // CHt = CH' - D'
    arma::mat::fixed<7, 3> mCHt = meas.C1.submat(0, 0, arma::size(7, 3));
    mCHt.submat(3, 0, arma::size(4, 3)) -= meas.C2.submat(3, 0, arma::size(4, 3));

    // Kalman gain K = mCH'*S
    arma::mat::fixed<7, 3> mK = mCHt * mS;

    // New estimation of the vertex position r += K*zeta
    fP.head(7) += mK * zeta;

    // New covariance matrix C -= K*(mCH')'
    fC.submat(0, 0, arma::size(7, 7)) -= mK * mCHt.t();
    fC = arma::symmatl(fC);  // force symmetry

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "mCHt", mCHt);
    Utils::Print(__FUNCTION__, "mK", mK);
    Utils::Print(__FUNCTION__, "fP (after Kalman gain)", fP);
    Utils::Print(__FUNCTION__, "fC (after Kalman gain)", fC);
#endif

    // do something else? //

    arma::mat::fixed<3, 3> K = meas.C1.submat(0, 0, arma::size(3, 3)) * mS;
    arma::mat::fixed<3, 3> K2 = arma::eye(3, 3) - K.t();
    arma::mat::fixed<3, 3> A = meas.D * K2;
    arma::mat::fixed<3, 3> M = K * A;

    fC.submat(0, 0, arma::size(3, 3)) += M + M.t();

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

    return true;
}

// Set a topological constraint on the current particle.
// Input arguments:
// - `prod_vtx`       : assumed production vertex
// - `prod_cov`       : production vertex's covariance matrix -- 3x3 symm. matrix
// - `bz`             : z-component of magnetic field
// - `chi2_threshold` : do an early cut in chi2
// Return: true if successful, false if not.
// Note: it will modify the state of the current `KF::Particle`
// Note: should be executed as final step, after the particle has been added all of its daughters!
bool Particle::AddProductionVertex(const arma::vec::fixed<3>& prod_vtx, const arma::mat::fixed<3, 3>& prod_cov, double bz, double chi2_threshold) {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    // store current particle's vertex information //

    arma::vec::fixed<3> decay_vtx = fP.head(3);
    arma::mat::fixed<3, 3> decay_cov = fC.submat(0, 0, arma::size(3, 3));  // symm.

    // transport to production vertex //

    auto min = Minimize(prod_vtx, bz);
    fPCAs.emplace_back(min.pca);

    auto tpr = Transport(min, bz);

    arma::mat::fixed<3, 3> sliced_corr = tpr.corr.submat(0, 0, arma::size(3, 3));
    arma::mat::fixed<3, 3> CTmp = sliced_corr * prod_cov * sliced_corr.t();      // symm.
    arma::mat::fixed<3, 3> measC = tpr.C.submat(0, 0, arma::size(3, 3)) + CTmp;  // symm.
    arma::mat::fixed<3, 3> D = prod_cov * sliced_corr;

    arma::mat::fixed<3, 3> mS_in = measC + prod_cov;
    if (!mS_in.is_symmetric()) return false;

    arma::mat::fixed<3, 3> mS;
    bool inv_success = arma::inv_sympd(mS, mS_in);
    if (!inv_success) return false;

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "CTmp", CTmp);
    Utils::Print(__FUNCTION__, "measC", measC);
    Utils::Print(__FUNCTION__, "D", D);
    Utils::Print(__FUNCTION__, "mS_in", mS_in);
    Utils::Print(__FUNCTION__, "mS", mS);
#endif

    arma::vec::fixed<3> res = prod_vtx - tpr.P.head(3);
    double dChi2{arma::as_scalar(res.t() * mS * res)};

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "res", res);
    Utils::PrintDouble(__FUNCTION__, "dChi2", dChi2);
#endif

    // optimization
    if (dChi2 > chi2_threshold) {
#if KF_DEBUG
        std::println(stdout, "({}) early return has been called!", __FUNCTION__);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
        return false;
    }

    // update current particle state //

    fP = tpr.P;
    fC.submat(0, 0, arma::size(3, 3)) = measC;
    fC.submat(3, 3, arma::size(4, 4)) = tpr.C.submat(3, 3, arma::size(4, 4));

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "fP (before Kalman gain)", fP);
    Utils::Print(__FUNCTION__, "fC (before Kalman gain)", fC);
#endif

    // Kalman gain calculation //

    arma::mat::fixed<7, 3> mCHt = fC.submat(0, 0, arma::size(7, 3));
    arma::mat::fixed<7, 3> mK = mCHt * mS;

    // fP += K × res
    fP.head(7) += mK * res;

    // fC -= K × mCHt^T (only lower triangle)
    fC.submat(0, 0, arma::size(7, 7)) -= mK * mCHt.t();

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "mCHt", mCHt);
    Utils::Print(__FUNCTION__, "mK", mK);
    Utils::Print(__FUNCTION__, "fP (after Kalman gain)", fP);
    Utils::Print(__FUNCTION__, "fC (after Kalman gain)", fC);
#endif

    arma::mat::fixed<3, 3> K = measC * mS;
    arma::mat::fixed<3, 3> K2 = arma::eye(3, 3) - K.t();
    arma::mat::fixed<3, 3> A = D.t() * K2;
    arma::mat::fixed<3, 3> M = K * A;

    fC.submat(0, 0, arma::size(3, 3)) += M + M.t();

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

    // update ds to decay vertex //

    auto min2decay = Minimize(decay_vtx, bz);
    fP(7) = min2decay.ds;

    // update Css //

    fC(7, 7) = 0.;

    arma::vec::fixed<6> ds_dr_6 = min2decay.ds_dr;
    arma::rowvec::fixed<6> dsdrC = ds_dr_6.t() * fC.submat(0, 0, arma::size(6, 6));
    fC.row(7).head(6) = dsdrC;

    fC(7, 7) += arma::dot(dsdrC, ds_dr_6);

    arma::vec::fixed<3> ds_dr_3 = min2decay.ds_dr.head(3);
    arma::rowvec::fixed<3> dsdpV = ds_dr_3.t() * decay_cov;

    fC(7, 7) -= arma::dot(dsdpV, ds_dr_3);

#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "fChi2", fChi2);
    Utils::PrintDouble(__FUNCTION__, "fNDF", double(fNDF));
    Utils::Print(__FUNCTION__, "fP (the end)", fP);
    Utils::Print(__FUNCTION__, "fC (the end)", fC);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

    return true;
}

// Set a mass constraint on the current particle.
// Constraint equation g(...) : E^2 - (Px^2 + Py^2 + Pz^2) - target_mass^2 = 0
// Input argument:
// - `mass` : the mass to be set on the state vector mP
// Return: true if successful, false if not.
// Note: it will modify the state of the current `KF::Particle`
bool Particle::AddMassConstraint(double target_mass) {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    double m2{target_mass * target_mass};
    double p2{Px() * Px() + Py() * Py() + Pz() * Pz()};

    // jacobian = d(g)/dr //
    arma::vec::fixed<8> mH = {0., 0., 0., -2 * Px(), -2 * Py(), -2 * Pz(), 2 * E(), 0.};

    // residual = target_mass^2 - current_mass^2 //
    double zeta{m2 - E() * E() + p2};

    arma::vec::fixed<8> mCHt = fC * mH;
    double s2 = arma::dot(mH, mCHt);

    // protection
    if (s2 < Const::AbsAlmostZero) {
#if KF_DEBUG
        std::println(stdout, "({}) early return has been called!", __FUNCTION__);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
        return false;
    }

    arma::vec::fixed<8> K = mCHt / s2;  // Kalman gain
    fP += K * zeta;
    fC -= K * mCHt.t();
    fC = arma::symmatl(fC);  // force symmetry

    fChi2 += zeta * zeta / s2;
    fNDF += 1;  // one d.o.f. is added because a single independent constraint has been applied

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "fP", fP);
    Utils::Print(__FUNCTION__, "fC", fC);
    std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif

    return true;
}

// Find point of closest approach (PCA) of this particle w.r.t. an arbitrary vertex.
// Input:
// - `v` : arbitrary vertex
// Return: (packed in a single `Result::Minimization` struct)
// - `ds_dr` : partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
Result::Minimization Particle::MinimizeLinePoint(const arma::vec::fixed<3>& v) const {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    double x0{X()};
    double y0{Y()};
    double z0{Z()};
    double px0{Px()};
    double py0{Py()};
    double pz0{Pz()};

    double dx{v(0) - x0};
    double dy{v(1) - y0};
    double dz{v(2) - z0};

    double p2{px0 * px0 + py0 * py0 + pz0 * pz0};
    double a{px0 * dx + py0 * dy + pz0 * dz};

    Result::Minimization min{};

    // protection
    if (p2 < Const::AbsAlmostZero || p2 * p2 < Const::AbsAlmostZero) {
#if KF_DEBUG
        std::println(stdout, "({}) early return has been called!", __FUNCTION__);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
        return min;
    }

    min.ds = a / p2;

    min.ds_dr(0) = -px0 / p2;
    min.ds_dr(1) = -py0 / p2;
    min.ds_dr(2) = -pz0 / p2;
    min.ds_dr(3) = (dx * p2 - 2. * px0 * a) / (p2 * p2);
    min.ds_dr(4) = (dy * p2 - 2. * py0 * a) / (p2 * p2);
    min.ds_dr(5) = (dz * p2 - 2. * pz0 * a) / (p2 * p2);

    min.ds_dr1(0) = -min.ds_dr(0);
    min.ds_dr1(1) = -min.ds_dr(1);
    min.ds_dr1(2) = -min.ds_dr(2);

    min.pca.xyz = {x0 + px0 * min.ds, y0 + py0 * min.ds, z0 + pz0 * min.ds};
    min.pca.dir = {px0, py0, pz0};

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
Result::Minimization Particle::MinimizeHelixPoint(const arma::vec::fixed<3>& v, double bz) const {
#if KF_DEBUG
    std::println(stdout, "-- starting ({}) --", __FUNCTION__);
#endif

    Result::Minimization min{};

    // 1 -- find point of closest approach (PCA) in XY plane //

    double bq{bz * fQ * Const::Kappa};

    double px0{fP(3)};
    double py0{fP(4)};
    double pz0{fP(5)};
    double pt2{px0 * px0 + py0 * py0};
    double x0{fP(0)};
    double y0{fP(1)};
    double z0{fP(2)};

    double dx{v(0) - x0};
    double dy{v(1) - y0};
    double dz{v(2) - z0};
    double a{dx * px0 + dy * py0};

    double abq{bq * a};
    double bbq{bq * (dx * py0 - dy * px0) - pt2};

    // 1.a -- get solution and update cache properties //

    min.theta = std::atan2(abq, -bbq);
    std::tie(min.sin, min.cos) = Math::sincos(min.theta);
    min.sB = min.sin / bq;
    min.cB = (1. - min.cos) / bq;

    min.ds = min.theta / bq;
    min.pca.xyz = {x0 + min.sB * px0 + min.cB * py0, y0 - min.cB * px0 + min.sB * py0, z0 + min.ds * pz0};
    min.pca.dir = {min.cos * px0 + min.sin * py0, -min.sin * px0 + min.cos * py0, pz0};

#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "min.ds (no z-correction)", min.ds);
    Utils::Print(__FUNCTION__, "min.(x,y,z) (no z-correction)", min.pca.xyz);
#endif

    // 1.b -- handle derivatives //

    double den{abq * abq + bbq * bbq};
    den = den < Const::AbsAlmostZero ? Const::AbsAlmostZero : den;

    min.ds_dr(0) = (px0 * bbq - py0 * abq) / den;
    min.ds_dr(1) = (px0 * abq + py0 * bbq) / den;
    min.ds_dr(2) = 0.;
    min.ds_dr(3) = -(dx * bbq + dy * abq + 2. * px0 * a) / den;
    min.ds_dr(4) = (dx * abq - dy * bbq - 2. * py0 * a) / den;
    min.ds_dr(5) = 0.;

    min.ds_dr1(0) = -min.ds_dr(0);
    min.ds_dr1(1) = -min.ds_dr(1);
    min.ds_dr1(2) = -min.ds_dr(2);

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "min.ds_dr (no z-correction)", min.ds_dr);
    Utils::Print(__FUNCTION__, "min.ds_dr1 (no z-correction)", min.ds_dr1);
#endif

    // 2 -- add z-component as small correction //

    double cbq{bbq * min.cos - abq * min.sin - pz0 * pz0};

    // protection
    if (std::abs(cbq) < Const::AbsAlmostZero) {
#if KF_DEBUG
        std::println(stdout, "({}) early return has been called!", __FUNCTION__);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
        return min;
    }

    double sz{(min.ds * pz0 - dz) * pz0 / cbq};

    // 2.a -- update derivatives //

    arma::vec::fixed<6> dc_dr = {
        -bq * py0 * min.cos - bbq * min.sin * bq * min.ds_dr(0) + px0 * bq * min.sin - abq * min.cos * bq * min.ds_dr(0),
        bq * px0 * min.cos - bbq * min.sin * bq * min.ds_dr(1) + py0 * bq * min.sin - abq * min.cos * bq * min.ds_dr(1),
        0.,
        (-bq * dy - 2. * px0) * min.cos - bbq * min.sin * bq * min.ds_dr(3) - dx * bq * min.sin - abq * min.cos * bq * min.ds_dr(3),
        (bq * dx - 2. * py0) * min.cos - bbq * min.sin * bq * min.ds_dr(4) - dy * bq * min.sin - abq * min.cos * bq * min.ds_dr(4),
        -2. * pz0};

    min.ds_dr += pz0 * pz0 * min.ds_dr / cbq - sz / cbq * dc_dr;
    min.ds_dr(2) += pz0 / cbq;
    min.ds_dr(5) += (2. * pz0 * min.ds - dz) / cbq;

    min.ds_dr1(0) = -min.ds_dr(0);
    min.ds_dr1(1) = -min.ds_dr(1);
    min.ds_dr1(2) = -min.ds_dr(2);

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

    min.pca.xyz = {x0 + min.sB * px0 + min.cB * py0, y0 - min.cB * px0 + min.sB * py0, z0 + min.ds * pz0};
    min.pca.dir = {min.cos * px0 + min.sin * py0, -min.sin * px0 + min.cos * py0, pz0};

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
            tmp1.pca.xyz = {x01 + tmp1.sB * px01 + tmp1.cB * py01, y01 - tmp1.cB * px01 + tmp1.sB * py01, z01 + tmp1.ds * pz01};
            tmp1.pca.dir = {tmp1.cos * px01 + tmp1.sin * py01, -tmp1.sin * px01 + tmp1.cos * py01, pz01};
        } else {
            tmp1.ds = (k11 * c1 + sign * k21 * d1) / (-k21 * c1);
            tmp1.pca.xyz = {x01 + px01 * tmp1.ds, y01 + py01 * tmp1.ds, z01 + pz01 * tmp1.ds};
            tmp1.pca.dir = {px01, py01, pz01};
        }

        // particle 2 //
        Cache tmp2{};
        if (!isStraight2) {
            tmp2.theta = std::atan2(bq2 * (k12 * c2 + sign * k22 * d1), sign * bq2 * k12 * d1 * bq2 - k22 * c2);
            std::tie(tmp2.sin, tmp2.cos) = Math::sincos(tmp2.theta);
            tmp2.sB = tmp2.sin / bq2;
            tmp2.cB = (1. - tmp2.cos) / bq2;
            tmp2.ds = tmp2.theta / bq2;
            tmp2.pca.xyz = {x02 + tmp2.sB * px02 + tmp2.cB * py02, y02 - tmp2.cB * px02 + tmp2.sB * py02, z02 + tmp2.ds * pz02};
            tmp2.pca.dir = {tmp2.cos * px02 + tmp2.sin * py02, -tmp2.sin * px02 + tmp2.cos * py02, pz02};
        } else {
            tmp2.ds = (k12 * c2 + sign * k22 * d1) / (-k22 * c2);
            tmp2.pca.xyz = {x02 + px02 * tmp2.ds, y02 + py02 * tmp2.ds, z02 + pz02 * tmp2.ds};
            tmp2.pca.dir = {px02, py02, pz02};
        }

        arma::vec::fixed<3> tmp_diff = tmp2.pca.xyz - tmp1.pca.xyz;
        double tmp_dca_sq{arma::dot(tmp_diff, tmp_diff)};

        // store //
        if (tmp_dca_sq < dca_sq) {
            static_cast<Cache&>(min1) = tmp1;
            static_cast<Cache&>(min2) = tmp2;

            dx = tmp_diff(0);
            dy = tmp_diff(1);
            dz = tmp_diff(2);

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
    Utils::PrintDouble(__FUNCTION__, "w_sign", double(w_sign));
    Utils::PrintDouble(__FUNCTION__, "dca", std::sqrt(dca_sq));
#endif

    // 1.b -- handle derivatives //

    arma::vec::fixed<6> dk11dr1 = {bq2 * px01, bq2 * py01, 0., bq2 * dx0 - py02, bq2 * dy0 + px02, 0.};
    arma::vec::fixed<6> dk11dr2 = {-bq2 * px01, -bq2 * py01, 0., py01, -px01, 0.};
    arma::vec::fixed<6> dk12dr1 = {bq1 * px02, bq1 * py02, 0., -py02, px02, 0.};
    arma::vec::fixed<6> dk12dr2 = {-bq1 * px02, -bq1 * py02, 0., bq1 * dx0 + py01, bq1 * dy0 - px01, 0.};
    arma::vec::fixed<6> dk21dr1 = {
        bq1 * bq2 * py01, -bq1 * bq2 * px01, 0., 2. * bq2 * px01 + bq1 * (-(bq2 * dy0) - px02), 2. * bq2 * py01 + bq1 * (bq2 * dx0 - py02), 0.};
    arma::vec::fixed<6> dk21dr2 = {-(bq1 * bq2 * py01), bq1 * bq2 * px01, 0., -(bq1 * px01), -(bq1 * py01), 0.};
    arma::vec::fixed<6> dk22dr1 = {bq1 * bq2 * py02, -(bq1 * bq2 * px02), 0., bq2 * px02, bq2 * py02, 0.};
    arma::vec::fixed<6> dk22dr2 = {
        -(bq1 * bq2 * py02), bq1 * bq2 * px02, 0., bq2 * (-(bq1 * dy0) + px01) - 2. * bq1 * px02, bq2 * (bq1 * dx0 + py01) - 2. * bq1 * py02, 0.};

    arma::vec::fixed<6> dkddr1 = {
        bq1 * bq2 * dx0 + bq2 * py01 - bq1 * py02, bq1 * bq2 * dy0 - bq2 * px01 + bq1 * px02, 0., -bq2 * dy0 - px02, bq2 * dx0 - py02, 0.};
    arma::vec::fixed<6> dkddr2 = {
        -bq1 * bq2 * dx0 - bq2 * py01 + bq1 * py02, -bq1 * bq2 * dy0 + bq2 * px01 - bq1 * px02, 0., bq1 * dy0 - px01, -bq1 * dx0 - py01, 0.};

    arma::vec::fixed<6> dc1dr1 = {-(bq1 * (bq1 * bq2 * dx0 + bq2 * py01 - bq1 * py02)), -(bq1 * (bq1 * bq2 * dy0 - bq2 * px01 + bq1 * px02)), 0.,
                                  -2. * bq2 * px01 - bq1 * (-bq2 * dy0 - px02),         -2. * bq2 * py01 - bq1 * (bq2 * dx0 - py02),          0.};
    arma::vec::fixed<6> dc1dr2 = {-bq1 * (-bq1 * bq2 * dx0 - bq2 * py01 + bq1 * py02),
                                  -bq1 * (-bq1 * bq2 * dy0 + bq2 * px01 - bq1 * px02),
                                  0.,
                                  -bq1 * (bq1 * dy0 - px01),
                                  -bq1 * (-bq1 * dx0 - py01),
                                  0.};
    arma::vec::fixed<6> dc2dr1 = {bq2 * (bq1 * bq2 * dx0 + bq2 * py01 - bq1 * py02),
                                  bq2 * (bq1 * bq2 * dy0 - bq2 * px01 + bq1 * px02),
                                  0.,
                                  bq2 * (-bq2 * dy0 - px02),
                                  bq2 * (bq2 * dx0 - py02),
                                  0.};
    arma::vec::fixed<6> dc2dr2 = {bq2 * (-bq1 * bq2 * dx0 - bq2 * py01 + bq1 * py02), bq2 * (-bq1 * bq2 * dy0 + bq2 * px01 - bq1 * px02), 0.,
                                  bq2 * (bq1 * dy0 - px01) + 2 * bq1 * px02,          bq2 * (-(bq1 * dx0) - py01) + 2 * bq1 * py02,       0.};
    arma::vec::fixed<6> dd1dr1(arma::fill::zeros);
    arma::vec::fixed<6> dd1dr2(arma::fill::zeros);

    if (d1 > 0.) {
        dd1dr1 = -kd * dkddr1 / d1;
        dd1dr2 = -kd * dkddr2 / d1;
        dd1dr1(3) += px01 * pt22 / d1;
        dd1dr1(4) += py01 * pt22 / d1;
        dd1dr2(3) += px02 * pt12 / d1;
        dd1dr2(4) += py02 * pt12 / d1;
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
            double dadr1{bq1 * (dk11dr1(iP) * c1 + k11 * dc1dr1(iP) + w_sign * dk21dr1(iP) * d1 + w_sign * k21 * dd1dr1(iP))};
            double dadr2{bq1 * (dk11dr2(iP) * c1 + k11 * dc1dr2(iP) + w_sign * dk21dr2(iP) * d1 + w_sign * k21 * dd1dr2(iP))};
            double dbdr1{w_sign * bq1 * bq1 * (dk11dr1(iP) * d1 + k11 * dd1dr1(iP)) - (dk21dr1(iP) * c1 + k21 * dc1dr1(iP))};
            double dbdr2{w_sign * bq1 * bq1 * (dk11dr2(iP) * d1 + k11 * dd1dr2(iP)) - (dk21dr2(iP) * c1 + k21 * dc1dr2(iP))};

            min1.ds_dr(iP) = d * (dadr1 * b - dbdr1 * a);
            min1.ds_dr1(iP) = d * (dadr2 * b - dbdr2 * a);
        }
    } else {
        double a{k11 * c1 + w_sign * k21 * d1};
        double b{-k21 * c1};

        for (size_t iP{0}; iP < 6; ++iP) {
            double dadr1{dk11dr1(iP) * c1 + k11 * dc1dr1(iP) + w_sign * dk21dr1(iP) * d1 + w_sign * k21 * dd1dr1(iP)};
            double dadr2{dk11dr2(iP) * c1 + k11 * dc1dr2(iP) + w_sign * dk21dr2(iP) * d1 + w_sign * k21 * dd1dr2(iP)};
            double dbdr1{-dk21dr1(iP) * c1 - k21 * dc1dr1(iP)};
            double dbdr2{-dk21dr2(iP) * c1 - k21 * dc1dr2(iP)};

            min1.ds_dr(iP) = dadr1 / b - dbdr1 * a / (b * b);
            min1.ds_dr1(iP) = dadr2 / b - dbdr2 * a / (b * b);
        }
    }

    if (!isStraight2) {
        double a{bq2 * (k12 * c2 + w_sign * k22 * d1)};
        double b{w_sign * bq2 * k12 * d1 * bq2 - k22 * c2};
        double c{b * b + a * a};
        double d{c > 0. ? (1. / bq2 * 1. / c) : 0.};

        for (size_t iP{0}; iP < 6; ++iP) {
            double dadr1{bq2 * (dk12dr1(iP) * c2 + k12 * dc2dr1(iP) + w_sign * dk22dr1(iP) * d1 + w_sign * k22 * dd1dr1(iP))};
            double dadr2{bq2 * (dk12dr2(iP) * c2 + k12 * dc2dr2(iP) + w_sign * dk22dr2(iP) * d1 + w_sign * k22 * dd1dr2(iP))};
            double dbdr1{w_sign * bq2 * bq2 * (dk12dr1(iP) * d1 + k12 * dd1dr1(iP)) - (dk22dr1(iP) * c2 + k22 * dc2dr1(iP))};
            double dbdr2{w_sign * bq2 * bq2 * (dk12dr2(iP) * d1 + k12 * dd1dr2(iP)) - (dk22dr2(iP) * c2 + k22 * dc2dr2(iP))};

            min2.ds_dr1(iP) = d * (dadr1 * b - dbdr1 * a);
            min2.ds_dr(iP) = d * (dadr2 * b - dbdr2 * a);
        }
    } else {
        double a{k12 * c2 + w_sign * k22 * d1};
        double b{-k22 * c2};

        for (size_t iP{0}; iP < 6; ++iP) {
            double dadr1{dk12dr1(iP) * c2 + k12 * dc2dr1(iP) + w_sign * dk22dr1(iP) * d1 + w_sign * k22 * dd1dr1(iP)};
            double dadr2{dk12dr2(iP) * c2 + k12 * dc2dr2(iP) + w_sign * dk22dr2(iP) * d1 + w_sign * k22 * dd1dr2(iP)};
            double dbdr1{-dk22dr1(iP) * c2 - k22 * dc2dr1(iP)};
            double dbdr2{-dk22dr2(iP) * c2 - k22 * dc2dr2(iP)};

            min2.ds_dr1(iP) = (std::abs(b) > Const::AbsAlmostZero ? dadr1 / b : 0.) - (b * b > Const::AbsAlmostZero ? dbdr1 * a / (b * b) : 0.);
            min2.ds_dr(iP) = (std::abs(b) > Const::AbsAlmostZero ? dadr2 / b : 0.) - (b * b > Const::AbsAlmostZero ? dbdr2 * a / (b * b) : 0.);
        }
    }

    double px1{min1.pca.dir(0)};
    double py1{min1.pca.dir(1)};
    double px2{min2.pca.dir(0)};
    double py2{min2.pca.dir(1)};

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

    double detp{lp1p2 * lp1p2 - p12 * p22};

    // protection
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

    arma::vec::fixed<6> dsldr0 = dsl1ds0 * min1.ds_dr + dsl1ds1 * min2.ds_dr1;
    arma::vec::fixed<6> dsldr1 = dsl1ds0 * min1.ds_dr1 + dsl1ds1 * min2.ds_dr;
    arma::vec::fixed<6> dsldr2 = dsl2ds0 * min1.ds_dr + dsl2ds1 * min2.ds_dr1;
    arma::vec::fixed<6> dsldr3 = dsl2ds0 * min1.ds_dr1 + dsl2ds1 * min2.ds_dr;

    min1.ds_dr += dsldr0;
    min1.ds_dr1 += dsldr1;
    min2.ds_dr1 += dsldr2;
    min2.ds_dr += dsldr3;

#if KF_DEBUG
    Utils::Print(__FUNCTION__, "min1.ds_dr (after z-correction 1)", min1.ds_dr);
    Utils::Print(__FUNCTION__, "min1.ds_dr1 (after z-correction 1)", min1.ds_dr1);
    Utils::Print(__FUNCTION__, "min2.ds_dr1 (after z-correction 1)", min2.ds_dr1);
    Utils::Print(__FUNCTION__, "min2.ds_dr (after z-correction 1)", min2.ds_dr);
#endif

    arma::vec::fixed<6> lp1p2_dr0 = {0., 0., 0., min1.cos * px2 - py2 * min1.sin, min1.cos * py2 + px2 * min1.sin, pz02};
    arma::vec::fixed<6> lp1p2_dr1 = {0., 0., 0., min2.cos * px1 - py1 * min2.sin, min2.cos * py1 + px1 * min2.sin, pz01};
    arma::vec::fixed<6> ldrp1_dr0 = {-px1,
                                     -py1,
                                     -pz01,
                                     min1.cB * py1 - px1 * min1.sB + min1.cos * dx - min1.sin * dy,
                                     -min1.cB * px1 - py1 * min1.sB + min1.sin * dx + min1.cos * dy,
                                     -min1.ds * pz01 + dz};
    arma::vec::fixed<6> ldrp1_dr1 = {px1, py1, pz01, -min2.cB * py1 + px1 * min2.sB, min2.cB * px1 + py1 * min2.sB, min2.ds * pz01};
    arma::vec::fixed<6> ldrp2_dr0 = {-px2, -py2, -pz02, min1.cB * py2 - px2 * min1.sB, -min1.cB * px2 - py2 * min1.sB, -min1.ds * pz02};
    arma::vec::fixed<6> ldrp2_dr1 = {px2,
                                     py2,
                                     pz02,
                                     -min2.cB * py2 + px2 * min2.sB + min2.cos * dx - min2.sin * dy,
                                     min2.cB * px2 + py2 * min2.sB + min2.sin * dx + min2.cos * dy,
                                     dz + min2.ds * pz02};
    arma::vec::fixed<6> p12_dr0 = {0., 0., 0., 2. * px01, 2. * py01, 2. * pz01};
    arma::vec::fixed<6> p22_dr1 = {0., 0., 0., 2. * px02, 2. * py02, 2. * pz02};

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

    arma::vec::fixed<6> a1_dr0 = ldrp2_dr0 * lp1p2 + ldrp2 * lp1p2_dr0 - ldrp1_dr0 * p22;
    arma::vec::fixed<6> a1_dr1 = ldrp2_dr1 * lp1p2 + ldrp2 * lp1p2_dr1 - ldrp1_dr1 * p22 - ldrp1 * p22_dr1;
    arma::vec::fixed<6> a2_dr0 = ldrp2_dr0 * p12 + ldrp2 * p12_dr0 - ldrp1_dr0 * lp1p2 - ldrp1 * lp1p2_dr0;
    arma::vec::fixed<6> a2_dr1 = ldrp2_dr1 * p12 - ldrp1_dr1 * lp1p2 - ldrp1 * lp1p2_dr1;

    arma::vec::fixed<6> detp_dr0 = 2. * lp1p2 * lp1p2_dr0 - p12_dr0 * p22;
    arma::vec::fixed<6> detp_dr1 = 2. * lp1p2 * lp1p2_dr1 - p12 * p22_dr1;

    min1.ds_dr += a1_dr0 / detp - a1 * detp_dr0 / (detp * detp);
    min1.ds_dr1 += a1_dr1 / detp - a1 * detp_dr1 / (detp * detp);
    min2.ds_dr1 += a2_dr0 / detp - a2 * detp_dr0 / (detp * detp);
    min2.ds_dr += a2_dr1 / detp - a2 * detp_dr1 / (detp * detp);

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

    min1.pca.xyz = {x01 + min1.sB * px01 + min1.cB * py01, y01 - min1.cB * px01 + min1.sB * py01, z01 + min1.ds * pz01};
    min1.pca.dir = {min1.cos * px01 + min1.sin * py01, -min1.sin * px01 + min1.cos * py01, pz01};

    min2.theta = bq2 * min2.ds;
    std::tie(min2.sin, min2.cos) = Math::sincos(min2.theta);
    min2.sB = min2.sin / bq2;
    min2.cB = (1. - min2.cos) / bq2;

    min2.pca.xyz = {x02 + min2.sB * px02 + min2.cB * py02, y02 - min2.cB * px02 + min2.sB * py02, z02 + min2.ds * pz02};
    min2.pca.dir = {min2.cos * px02 + min2.sin * py02, -min2.sin * px02 + min2.cos * py02, pz02};

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

    double x01{X()};
    double y01{Y()};
    double z01{Z()};
    double px01{Px()};
    double py01{Py()};
    double pz01{Pz()};

    double x02{p.X()};
    double y02{p.Y()};
    double z02{p.Z()};
    double px02{p.Px()};
    double py02{p.Py()};
    double pz02{p.Pz()};

    double p12{px01 * px01 + py01 * py01 + pz01 * pz01};
    double p22{px02 * px02 + py02 * py02 + pz02 * pz02};
    double p1p2{px01 * px02 + py01 * py02 + pz01 * pz02};

    double drp1{px01 * (x02 - x01) + py01 * (y02 - y01) + pz01 * (z02 - z01)};
    double drp2{px02 * (x02 - x01) + py02 * (y02 - y01) + pz02 * (z02 - z01)};

    double detp{p1p2 * p1p2 - p12 * p22};

    // protection
    if (std::abs(detp) < Const::AbsAlmostZero || std::abs(detp * detp) < Const::AbsAlmostZero) {
#if KF_DEBUG
        std::println(stdout, "({}) early return has been called!", __FUNCTION__);
        std::println(stdout, "-- finished ({}) --", __FUNCTION__);
#endif
        return {MinimizeLinePoint({0., 0., 0.}), MinimizeLinePoint({0., 0., 0.})};
    }

    min1.ds = (drp2 * p1p2 - drp1 * p22) / detp;
    min1.pca.xyz = {x01 + px01 * min1.ds, y01 + py01 * min1.ds, z01 + pz01 * min1.ds};
    min1.pca.dir = {px01, py01, pz01};

    min2.ds = (drp2 * p12 - drp1 * p1p2) / detp;
    min2.pca.xyz = {x02 + px02 * min2.ds, y02 + py02 * min2.ds, z02 + pz02 * min2.ds};
    min2.pca.dir = {px02, py02, pz02};

#if KF_DEBUG
    Utils::PrintDouble(__FUNCTION__, "min1.ds", min1.ds);
    Utils::Print(__FUNCTION__, "min1.(x,y,z)", min1.pca.xyz);
    Utils::PrintDouble(__FUNCTION__, "min2.ds", min2.ds);
    Utils::Print(__FUNCTION__, "min2.(x,y,z)", min2.pca.xyz);
#endif

    // handle derivatives //

    arma::vec::fixed<6> drp1_dr1 = {-px01, -py01, -pz01, -x01 + x02, -y01 + y02, -z01 + z02};
    arma::vec::fixed<6> drp1_dr2 = {px01, py01, pz01, 0., 0., 0.};
    arma::vec::fixed<6> drp2_dr1 = {-px02, -py02, -pz02, 0., 0., 0.};
    arma::vec::fixed<6> drp2_dr2 = {px02, py02, pz02, -x01 + x02, -y01 + y02, -z01 + z02};
    arma::vec::fixed<6> dp1p2_dr1 = {0., 0., 0., px02, py02, pz02};
    arma::vec::fixed<6> dp1p2_dr2 = {0., 0., 0., px01, py01, pz01};
    arma::vec::fixed<6> dp12_dr1 = {0., 0., 0., 2. * px01, 2. * py01, 2. * pz01};
    arma::vec::fixed<6> dp12_dr2 = {0., 0., 0., 0., 0., 0.};
    arma::vec::fixed<6> dp22_dr1 = {0., 0., 0., 0., 0., 0.};
    arma::vec::fixed<6> dp22_dr2 = {0., 0., 0., 2. * px02, 2. * py02, 2. * pz02};
    arma::vec::fixed<6> ddetp_dr1 = {
        0., 0., 0., -2 * p22 * px01 + 2. * p1p2 * px02, -2 * p22 * py01 + 2. * p1p2 * py02, -2 * p22 * pz01 + 2. * p1p2 * pz02};
    arma::vec::fixed<6> ddetp_dr2 = {
        0., 0., 0., 2. * p1p2 * px01 - 2. * p12 * px02, 2. * p1p2 * py01 - 2. * p12 * py02, 2. * p1p2 * pz01 - 2. * p12 * pz02};

    double a1{drp2 * p1p2 - drp1 * p22};
    double a2{drp2 * p12 - drp1 * p1p2};

    arma::vec::fixed<6> da1_dr1 = drp2_dr1 * p1p2 + drp2 * dp1p2_dr1 - drp1_dr1 * p22 - drp1 * dp22_dr1;
    arma::vec::fixed<6> da1_dr2 = drp2_dr2 * p1p2 + drp2 * dp1p2_dr2 - drp1_dr2 * p22 - drp1 * dp22_dr2;
    arma::vec::fixed<6> da2_dr1 = drp2_dr1 * p12 + drp2 * dp12_dr1 - drp1_dr1 * p1p2 - drp1 * dp1p2_dr1;
    arma::vec::fixed<6> da2_dr2 = drp2_dr2 * p12 + drp2 * dp12_dr2 - drp1_dr2 * p1p2 - drp1 * dp1p2_dr2;

    min1.ds_dr = da1_dr1 / detp - a1 * ddetp_dr1 / (detp * detp);
    min1.ds_dr1 = da1_dr2 / detp - a1 * ddetp_dr2 / (detp * detp);
    min2.ds_dr1 = da2_dr1 / detp - a2 * ddetp_dr1 / (detp * detp);
    min2.ds_dr = da2_dr2 / detp - a2 * ddetp_dr2 / (detp * detp);

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

    double px0{fP(3)};
    double py0{fP(4)};
    double pz0{fP(5)};

    tpr.P(0) = fP(0) + min.sB * px0 + min.cB * py0;
    tpr.P(1) = fP(1) - min.cB * px0 + min.sB * py0;
    tpr.P(2) = fP(2) + min.ds * pz0;
    tpr.P(3) = min.cos * px0 + min.sin * py0;
    tpr.P(4) = -min.sin * px0 + min.cos * py0;
    tpr.P(5) = fP(5);
    tpr.P(6) = fP(6);
    tpr.P(7) = fP(7);

    arma::mat::fixed<8, 8> mJ(arma::fill::eye);
    mJ(0, 3) = min.sB;
    mJ(0, 4) = min.cB;
    mJ(1, 3) = -min.cB;
    mJ(1, 4) = min.sB;
    mJ(2, 5) = min.ds;
    mJ(3, 3) = min.cos;
    mJ(3, 4) = min.sin;
    mJ(4, 3) = -min.sin;
    mJ(4, 4) = min.cos;

    arma::mat::fixed<6, 6> mJds(arma::fill::zeros);
    mJds(0, 3) = min.cos;
    mJds(0, 4) = min.sin;
    mJds(1, 3) = -min.sin;
    mJds(1, 4) = min.cos;
    mJds(2, 5) = 1.;
    mJds(3, 3) = -bq * min.sin;
    mJds(3, 4) = bq * min.cos;
    mJds(4, 3) = -bq * min.cos;
    mJds(4, 4) = -bq * min.sin;

    arma::vec::fixed<6> mJds_p = mJds.col(3) * px0 + mJds.col(4) * py0 + mJds.col(5) * pz0;
    mJ.submat(0, 0, arma::size(6, 6)) += mJds_p * min.ds_dr.t();

    tpr.C = mJ * fC * mJ.t();
    tpr.jacob = mJ.submat(0, 0, arma::size(6, 6));

    tpr.corr = mJds_p * min.ds_dr1.t();

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

    arma::mat::fixed<8, 8> mJ(arma::fill::eye);
    mJ(0, 3) = min.ds;
    mJ(1, 4) = min.ds;
    mJ(2, 5) = min.ds;

    double px{Px()};
    double py{Py()};
    double pz{Pz()};

    tpr.P(0) = fP(0) + min.ds * px;
    tpr.P(1) = fP(1) + min.ds * py;
    tpr.P(2) = fP(2) + min.ds * pz;
    tpr.P(3) = px;
    tpr.P(4) = py;
    tpr.P(5) = pz;
    tpr.P(6) = fP(6);
    tpr.P(7) = fP(7);

    arma::mat::fixed<6, 6> mJds(arma::fill::zeros);
    mJds(0, 3) = 1.;
    mJds(1, 4) = 1.;
    mJds(2, 5) = 1.;

    arma::vec::fixed<6> mJds_p = mJds.col(3) * px + mJds.col(4) * py + mJds.col(5) * pz;

    mJ.submat(0, 0, arma::size(6, 6)) += mJds_p * min.ds_dr.t();

    tpr.C = mJ * fC * mJ.t();
    tpr.jacob = mJ.submat(0, 0, arma::size(6, 6));
    tpr.corr = mJds_p * min.ds_dr1.t();

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
