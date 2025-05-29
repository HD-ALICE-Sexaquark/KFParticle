/*
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
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

#ifndef KFPTrack_H
#define KFPTrack_H

#include <cmath>

// @class KFPTrack
// @brief A scalar class for storage of the track in the cartesian parametrisation.
// @author  M.Zyzak, I.Kisel
// @date 05.02.2019
// @version 1.0
//
// A track is described with the state vector { X, Y, Z, Px, Py, Pz }
// and the corresponding covariance matrix. Also contains charge of the
// track, chi2 of the track fit, the corresponding number of degrees of freedom,
// the unique Id of the track and the field approximation along the track trajectory.
class KFPTrack {
   public:
    KFPTrack() : fChi2{-1.}, fQ{0}, fNDF{-1} {}
    ~KFPTrack() = default;

    // Fill an array p with the parameters of the track.
    // \param[out] p - array where { X, Y, Z, Px, Py, Pz } are copied
    bool GetXYZPxPyPz(float *p) const {
        for (int i = 0; i < 6; ++i) p[i] = fP[i];
        return 1;
    }
    // Copy the covariance matrix of the track to the array of floats.
    // \param[out] cv[21] - the output array, where the covariance matrix is copied
    bool GetCovarianceXYZPxPyPz(float cv[21]) const {
        for (int i = 0; i < 21; ++i) cv[i] = fC[i];
        return 1;
    }
    // Copy the covariance matrix of the track to the array of doubles.
    // \param[out] cv[21] - the output array, where the covariance matrix is copied
    bool GetCovarianceXYZPxPyPz(double cv[21]) const {
        for (int i = 0; i < 21; ++i) cv[i] = fC[i];
        return 1;
    }

    // Copy position of the track to the output array of floats.
    // \param[out] position - the output array with the position of the track
    void GetXYZ(float *position) const {
        position[0] = fP[0];
        position[1] = fP[1];
        position[2] = fP[2];
    }
    // Copy 3 momentum components of the track to the output array of floats.
    // \param[out] position - the output array with the momentum of the track
    void GetPxPyPz(float *position) const {
        position[0] = fP[3];
        position[1] = fP[4];
        position[2] = fP[5];
    }
    // Copy position of the track to the output array of floats.
    // \param[out] position - the output array with the position of the track
    void XvYvZv(float *position) const {
        position[0] = fP[0];
        position[1] = fP[1];
        position[2] = fP[2];
    }
    // Copy 3 momentum components of the track to the output array of floats.
    // \param[out] position - the output array with the momentum of the track
    void PxPyPz(float *position) const {
        position[0] = fP[3];
        position[1] = fP[4];
        position[2] = fP[5];
    }
    // Copy position of the track to the output array of doubles.
    // \param[out] position - the output array with the position of the track
    void XvYvZv(double *position) const {
        position[0] = fP[0];
        position[1] = fP[1];
        position[2] = fP[2];
    }
    // Copy 3 momentum components of the track to the output array of doubles.
    // \param[out] position - the output array with the momentum of the track
    void PxPyPz(double *position) const {
        position[0] = fP[3];
        position[1] = fP[4];
        position[2] = fP[5];
    }

    float GetX() const { return fP[0]; }   // return X coordinate of the track
    float GetY() const { return fP[1]; }   // return Y coordinate of the track
    float GetZ() const { return fP[2]; }   // return Z coordinate of the track
    float GetPx() const { return fP[3]; }  // return Px component of the momentum of the track
    float GetPy() const { return fP[4]; }  // return Py component of the momentum of the track
    float GetPz() const { return fP[5]; }  // return Pz component of the momentum of the track

    float GetPt() const { return std::sqrt(fP[3] * fP[3] + fP[4] * fP[4]); }                 // return Pt - transverse momentum of the track
    float GetP() const { return std::sqrt(fP[3] * fP[3] + fP[4] * fP[4] + fP[5] * fP[5]); }  // return P - momentum of the track

    // Copy the covariance matrix of the track to the array of floats
    // \param[out] covmatrix[21] - the output array, where the covariance matrix is copied
    void GetCovarianceMatrix(float *covmatrix) {
        for (int i = 0; i < 21; ++i) covmatrix[i] = fC[i];
    }
    float GetParameter(int i) const { return fP[i]; }  // return parameter "i" of the track. \param[in] i - index of the parameter to be returned

    // return element of the covariance matrix "i" of the track. \param[in] i - index of the element to be returned
    float GetCovariance(int i) const { return fC[i]; }

    int Charge() const { return fQ; }                     // return charge of the track
    float GetChi2perNDF() const { return fChi2 / fNDF; }  // return Chi2/NDF of the track, NDF is a number of degrees of freedom
    float GetChi2() const { return fChi2; }               // return Chi2 of the track
    int GetNDF() const { return fNDF; }                   // return number of degrees of freedom of the track

    // return a pointer to the array of track parameters
    const float *GetTrack() const { return fP; }
    // return a pointer to the array of the covariance matrix elements stored in a lower triangular form
    const float *GetCovMatrix() const { return fC; }

    // set parameters { X, Y, Z, Px, Py, Pz } of the track from the input array of floats
    // \param[in] position - input array with the track parameters
    void SetParameters(const float *position) {
        for (int i = 0; i < 6; ++i) fP[i] = position[i];
    }
    // set parameters { X, Y, Z, Px, Py, Pz } of the track from the input array of doubles
    // \param[in] position - input array with the track parameters
    void SetParameters(double *position) {
        for (int i = 0; i < 6; ++i) fP[i] = position[i];
    }

    // set parameters { X, Y, Z, Px, Py, Pz } of the track.
    // \param[in] x - X coordinate to be set
    // \param[in] y - Y coordinate to be set
    // \param[in] z - Z coordinate to be set
    // \param[in] Px - Px momentum component to be set
    // \param[in] Py - Py momentum component to be set
    // \param[in] Pz - Pz momentum component to be set
    void SetParameters(float x, float y, float z, float px, float py, float pz) {
        fP[0] = x;
        fP[1] = y;
        fP[2] = z;
        fP[3] = px;
        fP[4] = py;
        fP[5] = pz;
    }
    // set position { X, Y, Z } of the track.
    // \param[in] x - X coordinate to be set
    // \param[in] y - Y coordinate to be set
    // \param[in] z - Z coordinate to be set
    void SetXYZ(float x, float y, float z) {
        fP[0] = x;
        fP[1] = y;
        fP[2] = z;
    }
    // set momentum { Px, Py, Pz } of the track.
    // \param[in] Px - Px momentum component to be set
    // \param[in] Py - Py momentum component to be set
    // \param[in] Pz - Pz momentum component to be set
    void SetPxPyPz(float px, float py, float pz) {
        fP[3] = px;
        fP[4] = py;
        fP[5] = pz;
    }

    void SetX(float x) { fP[0] = x; }         // set X coordinate of the track
    void SetY(float y) { fP[1] = y; }         // set Y coordinate of the track
    void SetZ(float z) { fP[2] = z; }         // set Z coordinate of the track
    void SetPx(float px) { fP[3] = px; }      // set Px component of the track momentum
    void SetPy(float py) { fP[4] = py; }      // set Py component of the track momentum
    void SetPz(float pz) { fP[5] = pz; }      // set Pz component of the track momentum
    void SetCharge(int q) { fQ = q; }         // set charge of the track
    void SetChi2(float chi) { fChi2 = chi; }  // set a value of the track Chi2
    void SetNDF(int ndf) { fNDF = ndf; }      // set a value of the number of degrees of freedom

    // set the covariance matrix from the input array of floats.
    // \param[in] C[21] - array with the input elements of the covariance matrix stored in the lower triangular form
    void SetCovarianceMatrix(const float *C) {
        for (int i = 0; i < 21; ++i) fC[i] = C[i];
    }
    // set the covariance matrix from the input array of doubles.
    // \param[in] C[21] - array with the input elements of the covariance matrix stored in the lower triangular form
    void SetCovarianceMatrix(const double *C) {
        for (int i = 0; i < 21; ++i) fC[i] = C[i];
    }

    // set an element of the covariance matrix with index "i".
    // \param[in] c - value to be set
    // \param[in] i - index of the element
    void SetCovariance(const int i, const float c) { fC[i] = c; }

   private:
    float fP[6];   // Parameters of the track: { X, Y, Z, Px, Py, Pz }
    float fC[21];  // Covariance matrix of the track parameters Stored in the lower triangular form
    float fChi2;   // Chi-square of the track fit
    char fQ;       // Charge of the track
    short fNDF;    // Number of degree of freedom of the fit
};

#endif
