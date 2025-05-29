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

#ifndef KFPVertex_H
#define KFPVertex_H

// @class KFPVertex
// @brief A scalar class for storage of the vertex in the cartesian parametrisation.
// @author M.Zyzak, I.Kisel
// @date 05.02.2019
// @version 1.0
//
// A vertex is described with the state vector { X, Y, Z }
// and the corresponding covariance matrix. Also contains chi2 of the fit,
// corresponding number of degrees of freedom,
// and number of tracks which were used to construct current vertex.
// The class is used to provide an external vertex through the interfaces
// to the KF Particle package.
class KFPVertex {
   public:
    KFPVertex() : fChi2{-1.}, fNContributors{0}, fNDF{-1} {
        for (int iP = 0; iP < 3; ++iP) fP[iP] = 0;
        for (int iC = 0; iC < 6; ++iC) fC[iC] = 0;
    }
    ~KFPVertex() = default;

    float GetX() const { return fP[0]; }  // return X coordinate of the vertex
    float GetY() const { return fP[1]; }  // return Y coordinate of the vertex
    float GetZ() const { return fP[2]; }  // return Z coordinate of the vertex

    // Copy position of the vertex to the output array of floats.
    // \param[out] position - the output array with the position of the vertex
    void GetXYZ(float *position) const {
        position[0] = fP[0];
        position[1] = fP[1];
        position[2] = fP[2];
    }
    // Copy position of the vertex to the output array of doubles.
    // \param[out] position - the output array with the position of the vertex
    void GetXYZ(double *position) const {
        position[0] = fP[0];
        position[1] = fP[1];
        position[2] = fP[2];
    }
    // Copy the covariance matrix of the vertex to the array of floats.
    // \param[out] covmatrix[6] - the output array, where the covariance matrix is copied
    void GetCovarianceMatrix(float *covmatrix) const {
        for (int i = 0; i < 6; ++i) covmatrix[i] = fC[i];
    }
    // Copy the covariance matrix of the vertex to the array of doubles.
    // \param[out] covmatrix[6] - the output array, where the covariance matrix is copied
    void GetCovarianceMatrix(double *covmatrix) const {
        for (int i = 0; i < 6; ++i) covmatrix[i] = fC[i];
    }

    float GetChi2perNDF() const { return fChi2 / fNDF; }     // return Chi2/NDF of the vertex, NDF is a number of degrees of freedom
    float GetChi2() const { return fChi2; }                  // return Chi2 of the vertex fit
    int GetNDF() const { return fNDF; }                      // return number of degrees of freedom of the vertex
    int GetNContributors() const { return fNContributors; }  // return number of tracks which were used for construction of the vertex

    // return parameter "i" of the vertex.
    // \param[in] i - index of the parameter to be returned
    float GetParameter(int i) const { return fP[i]; }
    // return element of the covariance matrix "i" of the vertex.
    // \param[in] i - index of the element to be returned
    float GetCovariance(int i) const { return fC[i]; }

    // set position { X, Y, Z } of the vertex from the input array of doubles.
    // \param[in] position - input array with the vertex parameters
    void SetXYZ(float *position) {
        fP[0] = position[0];
        fP[1] = position[1];
        fP[2] = position[2];
    }
    // set position { X, Y, Z } of the vertex.
    // \param[in] x - X coordinate to be set
    // \param[in] y - Y coordinate to be set
    // \param[in] z - Z coordinate to be set
    void SetXYZ(float x, float y, float z) {
        fP[0] = x;
        fP[1] = y;
        fP[2] = z;
    }
    void SetX(float x) { fP[0] = x; }                       // set X coordinate of the vertex
    void SetY(float y) { fP[1] = y; }                       // set Y coordinate of the vertex
    void SetZ(float z) { fP[2] = z; }                       // set Z coordinate of the vertex
    void SetChi2(float chi) { fChi2 = chi; }                // set Chi2 of the vertex
    void SetNDF(int ndf) { fNDF = ndf; }                    // set number of degrees of freedom of the vertex
    void SetNContributors(int nc) { fNContributors = nc; }  // set number of tracks which were used for construction of the vertex

    // set the covariance matrix from the input array of floats.
    // \param[in] C[6] - array with the input elements of the covariance matrix stored in the lower triangular form
    void SetCovarianceMatrix(float *C) {
        for (int i = 0; i < 6; ++i) fC[i] = C[i];
    }

    // set the covariance matrix from the input array of floats.
    // \param[in] C00 - Cxx
    // \param[in] C10 - Cxy = Cyx
    // \param[in] C11 - Cyy
    // \param[in] C20 - Cxz = Czx
    // \param[in] C21 - Cyz = Czy
    // \param[in] C22 - Czz
    void SetCovarianceMatrix(float C00, float C10, float C11, float C20, float C21, float C22) {
        fC[0] = C00;
        fC[1] = C10;
        fC[2] = C11;
        fC[3] = C20;
        fC[4] = C21;
        fC[5] = C22;
    }

   private:
    float fP[3];         // coordinates of the vertex
    float fC[6];         // covariance matrix of the vertex parameters
    float fChi2;         // chi-square of the vertex fit
    int fNContributors;  // number of tracks, from which the vertex was built
    int fNDF;            // number of degrees of freedom of the vertex fit
};

#endif
