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

#include "KFParticleBase.h"
#include <cmath>

class KFPTrack;
class KFPVertex;

// @class KFParticle
// @brief The main scalar class of KF Particle package, describes particle objects.
// @author  S.Gorbunov, I.Kisel, M.Zyzak
// @date 05.02.2019
// @version 1.0
//
// The main scalar class of KF Particle pacakge, describes particle objects.
// The particle is described with the state vector { X, Y, Z, Px, Py, Pz, E }
// and the corresponding covariance matrix.
// It contains functionality to create particle-object from track, to construct
// short-lived particles from other tracks or particles. The mathematics is
// based on the Kalman filter method. It also allows to subtract particles from
// the already constructed object,
// to transport particles, get parameters together with their errors, get distance
// to other particles and vertices, get deviations from them in terms of errors, etc.
class KFParticle : public KFParticleBase {

   public:
    //  INITIALIZATION

    // Set magnetic field for all particles
#ifdef HomogeneousField
    static void SetField(float Bz);
#endif
    // Constructor (empty)

    KFParticle() : KFParticleBase() { ; }

    // Destructor (empty)

    virtual ~KFParticle() { ; }

    // Construction of mother particle by its 2-3-4 daughters

    KFParticle(const KFParticle &d1, const KFParticle &d2);

    KFParticle(const KFParticle &d1, const KFParticle &d2, const KFParticle &d3);

    KFParticle(const KFParticle &d1, const KFParticle &d2, const KFParticle &d3, const KFParticle &d4);

    // Initialization from "cartesian" coordinates ( X Y Z Px Py Pz )
    // Parameters, covariance matrix, charge and PID hypothesis should be provided

    void Create(const float Param[], const float Cov[], int Charge, float mass /*int PID*/);
    void Create(const double Param[], const double Cov[], int Charge, float mass /*int PID*/);

    // Initialization from ALICE track, PID hypothesis shoould be provided

    KFParticle(const KFPTrack &track, const float mass);

    // Initialization from VVertex

    KFParticle(const KFPVertex &vertex);

    // Initialize covariance matrix and set current parameters to 0.0

    void Initialize();

    //  ACCESSORS

    // Simple accessors

    float GetX() const;     // Return X coordinate of the particle, fP[0]
    float GetY() const;     // Return Y coordinate of the particle, fP[1]
    float GetZ() const;     // Return Z coordinate of the particle, fP[2]
    float GetPx() const;    // Return X component of the momentum, fP[3]
    float GetPy() const;    // Return Y component of the momentum, fP[4]
    float GetPz() const;    // Return Z component of the momentum, fP[5]
    float GetE() const;     // return energy of the particle, fP[6]
    float GetS() const;     // return dS=l/p, l - decay length, fP[7], defined if production vertex is set
    char GetQ() const;      // return charge of the particle
    float GetChi2() const;  // return Chi2 of the fit
    int GetNDF() const;     // return number of decrease of freedom

    bool GetAtProductionVertex() const {
        return fAtProductionVertex;
    }  // return a flag which shows if the particle is located at the production point
    void SetAtProductionVertex(bool b) { fAtProductionVertex = b; }  // Set a flag that particle is at the production point

#ifdef NonhomogeneousField
    const float *GetFieldCoeff() const { return fieldRegion; }  // return the field approximation for the current particle
    void SetFieldCoeff(float c, int i) { fieldRegion[i] = c; }  // set the field coefficient with index i
#endif

    const float &X() const { return fP[0]; }     // Return X coordinate of the particle, fP[0]
    const float &Y() const { return fP[1]; }     // Return Y coordinate of the particle, fP[1]
    const float &Z() const { return fP[2]; }     // Return Z coordinate of the particle, fP[2]
    const float &Px() const { return fP[3]; }    // Return X component of the momentum, fP[3]
    const float &Py() const { return fP[4]; }    // Return Y component of the momentum, fP[4]
    const float &Pz() const { return fP[5]; }    // Return Z component of the momentum, fP[5]
    const float &E() const { return fP[6]; }     // return energy of the particle, fP[6]
    const float &S() const { return fP[7]; }     // return dS=l/p, l - decay length, fP[7], defined if production vertex is set
    const char &Q() const { return fQ; }         // return charge of the particle
    const float &Chi2() const { return fChi2; }  // return Chi2 of the fit
    const int &NDF() const { return fNDF; }      // return number of decrease of freedom

    float GetParameter(int i) const;          // return P[i] parameter
    float GetCovariance(int i) const;         // return C[i] element of the covariance matrix in the lower triangular form
    float GetCovariance(int i, int j) const;  // return C[i,j] element of the covariance matrix

    // Accessors with calculations, value returned w/o error flag

    float GetP() const;              // return momentum
    float GetPt() const;             // return transverse momentum
    float GetEta() const;            // return pseudorapidity
    float GetPhi() const;            // return the azimuthal angle phi
    float GetMomentum() const;       // return momentum
    float GetMass() const;           // return mass
    float GetDecayLength() const;    // return decay length
    float GetDecayLengthXY() const;  // return decay length in XY
    float GetLifeTime() const;       // return life time ctau [cm]
    float GetR() const;              // return distance to the origin of the coordinate system {0,0,0}

    // Accessors to estimated errors

    float GetErrX() const;              // return the error of X of current position
    float GetErrY() const;              // return the error of Y of current position
    float GetErrZ() const;              // return the error of Z of current position
    float GetErrPx() const;             // return the error of X-compoment of the particle momentum
    float GetErrPy() const;             // return the error of Y-compoment of the particle momentum
    float GetErrPz() const;             // return the error of Z-compoment of the particle momentum
    float GetErrE() const;              // return the error of energy
    float GetErrS() const;              // return the error of decay length / momentum
    float GetErrP() const;              // return the error of momentum
    float GetErrPt() const;             // return the error of transverse momentum
    float GetErrEta() const;            // return the error of pseudorapidity
    float GetErrPhi() const;            // return the error of the azimuthal angle phi
    float GetErrMomentum() const;       // return the error of momentum
    float GetErrMass() const;           // return the error of mass
    float GetErrDecayLength() const;    // return the error of decay length
    float GetErrDecayLengthXY() const;  // return the error of decay length in XY
    float GetErrLifeTime() const;       // return the error of life time
    float GetErrR() const;              // return the error of distance to the origin of the coordinate system {0,0,0}

    // Accessors with calculations( &value, &estimated sigma )
    // error flag returned (0 means no error during calculations)

    int GetP(float &P, float &SigmaP) const;                                            // momentum
    int GetPt(float &Pt, float &SigmaPt) const;                                         // transverse momentum
    int GetEta(float &Eta, float &SigmaEta) const;                                      // pseudorapidity
    int GetPhi(float &Phi, float &SigmaPhi) const;                                      // phi
    int GetMomentum(float &P, float &SigmaP) const;                                     // momentum
    int GetMass(float &M, float &SigmaM) const;                                         // mass
    int GetDecayLength(float &L, float &SigmaL) const;                                  // decay length
    int GetDecayLengthXY(float &L, float &SigmaL) const;                                // decay length in XY
    int GetLifeTime(float &T, float &SigmaT) const;                                     // life time
    int GetR(float &R, float &SigmaR) const;                                            // R
    float GetRapidity() const { return 0.5 * log((fP[6] + fP[5]) / (fP[6] - fP[5])); }  // return rapidity of the particle
    float GetTheta() const { return atan2(GetPt(), fP[5]); }                            // return the polar angle in RZ

    //  MODIFIERS

    float &X();     // Modifier of X coordinate of the particle, fP[0].
    float &Y();     // Modifier of Y coordinate of the particle, fP[1].
    float &Z();     // Modifier of Z coordinate of the particle, fP[2].
    float &Px();    // Modifier of X component of the momentum, fP[3].
    float &Py();    // Modifier of Y component of the momentum, fP[4].
    float &Pz();    // Modifier of Z component of the momentum, fP[5].
    float &E();     // Modifier of energy of the particle, fP[6].
    float &S();     // Modifier of dS=l/p, l - decay length, fP[7], defined if production vertex is set.
    char &Q();      // Modifier of charge of the particle.
    float &Chi2();  // Modifier of Chi2 of the fit.
    int &NDF();     // Modifier of number of decrease of freedom.

    float &Parameter(int i);          // Modifier of P[i] parameter.
    float &Covariance(int i);         // Modifier of C[i] element of the covariance matrix in the lower triangular form.
    float &Covariance(int i, int j);  // Modifier of C[i,j] element of the covariance matrix.
    float *Parameters();              // return pointer to the parameters fP
    float *CovarianceMatrix();        // return pointer to the covariance matrix fC

    // CONSTRUCTION OF THE PARTICLE BY ITS DAUGHTERS AND MOTHER USING THE KALMAN FILTER METHOD

    // Add daughter to the particle
    void AddDaughter(const KFParticle &Daughter);

    // Add daughter via += operator: ex.{ D0; D0+=Pion; D0+= Kaon; }
    void operator+=(const KFParticle &Daughter);

    // Everything in one go
    void Construct(const KFParticle *vDaughters[], int nDaughters, const KFParticle *ProdVtx = 0, float Mass = -1);

    //                   TRANSPORT
    //  ( main transportation parameter is S = SignedPath/Momentum )
    //  ( parameters of decay & production vertices are stored locally )

    // Transport the particle close to xyz[] point

    void TransportToPoint(const float xyz[]);

    // Transport the particle close to VVertex
#ifdef HomogeneousField
    void TransportToVertex(const KFPVertex &v);
#endif
    // Transport the particle close to another particle p
    void TransportToParticle(const KFParticle &p);

    // Get dS to a certain space point
    float GetDStoPoint(const float xyz[3], float dsdr[6]) const;

    // Get dS to other particle p (dSp for particle p also returned)
    void GetDStoParticle(const KFParticleBase &p, float dS[2], float dsdr[4][6]) const;

    // OTHER UTILITIES

    // Calculate distance from another object [cm] in XY-plane

    bool GetDistanceFromVertexXY(const float vtx[], float &val, float &err) const;
    bool GetDistanceFromVertexXY(const float vtx[], const float Cv[], float &val, float &err) const;
    bool GetDistanceFromVertexXY(const KFParticle &Vtx, float &val, float &err) const;
#ifdef HomogeneousField
    bool GetDistanceFromVertexXY(const KFPVertex &Vtx, float &val, float &err) const;
#endif

    float GetDistanceFromVertexXY(const float vtx[]) const;
    float GetDistanceFromVertexXY(const KFParticle &Vtx) const;
#ifdef HomogeneousField
    float GetDistanceFromVertexXY(const KFPVertex &Vtx) const;
#endif
    float GetDistanceFromParticleXY(const KFParticle &p) const;

    // Calculate sqrt(Chi2/ndf) deviation from another object in XY plane
    // ( v = [xyz]-vertex, Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix )

    float GetDeviationFromVertexXY(const float v[], const float Cv[] = 0) const;
    float GetDeviationFromVertexXY(const KFParticle &Vtx) const;
#ifdef HomogeneousField
    float GetDeviationFromVertexXY(const KFPVertex &Vtx) const;
#endif
    float GetDeviationFromParticleXY(const KFParticle &p) const;

    // Get parameters at an arbitrary reconstructed point taking into account its errors
    void GetParametersAtPoint(const float *point, const float *pointCov, float *m, float *mV);

    // Calculate opennig angle between two particles

    float GetAngle(const KFParticle &p) const;
    float GetAngleXY(const KFParticle &p) const;
    float GetAngleRZ(const KFParticle &p) const;

    float GetPseudoProperDecayTime(const KFParticle &primVertex, const float &mass, float *timeErr2 = 0) const;

    void GetFieldValue(const float xyz[], float B[]) const;

    void Transport(float dS, const float *dsdr, float P[], float C[], float *dsdr1 = 0, float *F = 0, float *F1 = 0) const;

   protected:
    //  INTERNAL STUFF

    // Method to access ALICE field
#ifdef HomogeneousField
    static float GetFieldAlice();
#endif

   private:
#ifdef HomogeneousField
    static float fgBz;  // Bz compoment of the magnetic field (is defined in case of #ifdef HomogeneousField)
#endif
#ifdef NonhomogeneousField
    // \brief Approximation of the magnetic field along the track trajectory.
    // Each component (Bx, By, Bz) is approximated with the parabola depending on Z coordinate. Is defined in case of #ifdef NonhomogeneousField.
    //
    float fieldRegion[10];
#endif
};

// Inline implementation of the KFParticle methods

#ifdef HomogeneousField
// set the constant homogemeous one-component magnetic field Bz (is defined in case of #ifdef HomogeneousField).
// \param[in] Bz - Z-component of the magnetic field
inline void KFParticle::SetField(float Bz) { fgBz = Bz; }
#endif

// Constructs a particle from three input daughter particles
// \param[in] d1 - the first daughter particle
// \param[in] d2 - the second daughter particle
// \param[in] d3 - the third daughter particle
inline KFParticle::KFParticle(const KFParticle &d1, const KFParticle &d2, const KFParticle &d3) {
    KFParticle mother;
    mother += d1;
    mother += d2;
    mother += d3;
    *this = mother;
}

// Constructs a particle from four input daughter particles
// \param[in] d1 - the first daughter particle
// \param[in] d2 - the second daughter particle
// \param[in] d3 - the third daughter particle
// \param[in] d4 - the fourth daughter particle
inline KFParticle::KFParticle(const KFParticle &d1, const KFParticle &d2, const KFParticle &d3, const KFParticle &d4) {
    KFParticle mother;
    mother += d1;
    mother += d2;
    mother += d3;
    mother += d4;
    *this = mother;
}

inline void KFParticle::Initialize() { KFParticleBase::Initialize(); }

inline float KFParticle::GetX() const { return KFParticleBase::GetX(); }
inline float KFParticle::GetY() const { return KFParticleBase::GetY(); }
inline float KFParticle::GetZ() const { return KFParticleBase::GetZ(); }
inline float KFParticle::GetPx() const { return KFParticleBase::GetPx(); }
inline float KFParticle::GetPy() const { return KFParticleBase::GetPy(); }
inline float KFParticle::GetPz() const { return KFParticleBase::GetPz(); }
inline float KFParticle::GetE() const { return KFParticleBase::GetE(); }
inline float KFParticle::GetS() const { return KFParticleBase::GetS(); }
inline char KFParticle::GetQ() const { return KFParticleBase::GetQ(); }
inline float KFParticle::GetChi2() const { return KFParticleBase::GetChi2(); }
inline int KFParticle::GetNDF() const { return KFParticleBase::GetNDF(); }
inline float KFParticle::GetParameter(int i) const { return KFParticleBase::GetParameter(i); }
inline float KFParticle::GetCovariance(int i) const { return KFParticleBase::GetCovariance(i); }
inline float KFParticle::GetCovariance(int i, int j) const { return KFParticleBase::GetCovariance(i, j); }
inline float KFParticle::GetP() const {
    float par, err;
    if (KFParticleBase::GetMomentum(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetPt() const {
    float par, err;
    if (KFParticleBase::GetPt(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetEta() const {
    float par, err;
    if (KFParticleBase::GetEta(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetPhi() const {
    float par, err;
    if (KFParticleBase::GetPhi(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetMomentum() const {
    float par, err;
    if (KFParticleBase::GetMomentum(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetMass() const {
    float par, err;
    if (KFParticleBase::GetMass(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetDecayLength() const {
    float par, err;
    if (KFParticleBase::GetDecayLength(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetDecayLengthXY() const {
    float par, err;
    if (KFParticleBase::GetDecayLengthXY(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetLifeTime() const {
    float par, err;
    if (KFParticleBase::GetLifeTime(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetR() const {
    float par, err;
    if (KFParticleBase::GetR(par, err))
        return 0;
    else
        return par;
}

inline float KFParticle::GetErrX() const { return std::sqrt(fabs(GetCovariance(0, 0))); }
inline float KFParticle::GetErrY() const { return std::sqrt(fabs(GetCovariance(1, 1))); }
inline float KFParticle::GetErrZ() const { return std::sqrt(fabs(GetCovariance(2, 2))); }
inline float KFParticle::GetErrPx() const { return std::sqrt(fabs(GetCovariance(3, 3))); }
inline float KFParticle::GetErrPy() const { return std::sqrt(fabs(GetCovariance(4, 4))); }
inline float KFParticle::GetErrPz() const { return std::sqrt(fabs(GetCovariance(5, 5))); }
inline float KFParticle::GetErrE() const { return std::sqrt(fabs(GetCovariance(6, 6))); }
inline float KFParticle::GetErrS() const { return std::sqrt(fabs(GetCovariance(7, 7))); }
inline float KFParticle::GetErrP() const {
    float par, err;
    if (KFParticleBase::GetMomentum(par, err))
        return 1.e10;
    else
        return err;
}
inline float KFParticle::GetErrPt() const {
    float par, err;
    if (KFParticleBase::GetPt(par, err))
        return 1.e10;
    else
        return err;
}
inline float KFParticle::GetErrEta() const {
    float par, err;
    if (KFParticleBase::GetEta(par, err))
        return 1.e10;
    else
        return err;
}
inline float KFParticle::GetErrPhi() const {
    float par, err;
    if (KFParticleBase::GetPhi(par, err))
        return 1.e10;
    else
        return err;
}
inline float KFParticle::GetErrMomentum() const {
    float par, err;
    if (KFParticleBase::GetMomentum(par, err))
        return 1.e10;
    else
        return err;
}
inline float KFParticle::GetErrMass() const {
    float par, err;
    if (KFParticleBase::GetMass(par, err))
        return 1.e10;
    else
        return err;
}

inline float KFParticle::GetErrDecayLength() const {
    float par, err;
    if (KFParticleBase::GetDecayLength(par, err))
        return 1.e10;
    else
        return err;
}

inline float KFParticle::GetErrDecayLengthXY() const {
    float par, err;
    if (KFParticleBase::GetDecayLengthXY(par, err))
        return 1.e10;
    else
        return err;
}

inline float KFParticle::GetErrLifeTime() const {
    float par, err;
    if (KFParticleBase::GetLifeTime(par, err))
        return 1.e10;
    else
        return err;
}

inline float KFParticle::GetErrR() const {
    float par, err;
    if (KFParticleBase::GetR(par, err))
        return 1.e10;
    else
        return err;
}

// Calculate particle momentum and its error. If they are well defined return 0, otherwise 1.
// \param[out] P - momentum of the particle
// \param[out] SigmaP - its error
inline int KFParticle::GetP(float &P, float &SigmaP) const { return KFParticleBase::GetMomentum(P, SigmaP); }

// Calculate particle transverse  momentum and its error. If they are well defined return 0, otherwise 1.
// \param[out] Pt - transverse momentum of the particle
// \param[out] SigmaPt - its error
inline int KFParticle::GetPt(float &Pt, float &SigmaPt) const { return KFParticleBase::GetPt(Pt, SigmaPt); }

// Calculate particle pseudorapidity and its error. If they are well defined return 0, otherwise 1.
// \param[out] Eta - pseudorapidity of the particle
// \param[out] SigmaEta - its error
inline int KFParticle::GetEta(float &Eta, float &SigmaEta) const { return KFParticleBase::GetEta(Eta, SigmaEta); }

// Calculate particle polar angle at the current point and its error. If they are well defined return 0, otherwise 1.
// \param[out] Phi - polar angle of the particle
// \param[out] SigmaPhi - its error
inline int KFParticle::GetPhi(float &Phi, float &SigmaPhi) const { return KFParticleBase::GetPhi(Phi, SigmaPhi); }

// Calculate particle momentum and its error. If they are well defined return 0, otherwise 1.
// \param[out] P - momentum of the particle
// \param[out] SigmaP - its error
inline int KFParticle::GetMomentum(float &P, float &SigmaP) const { return KFParticleBase::GetMomentum(P, SigmaP); }

// Calculate the mass of the particle and its error. If they are well defined return 0, otherwise 1.
// \param[out] M - mass of the particle
// \param[out] SigmaM - its error

inline int KFParticle::GetMass(float &M, float &SigmaM) const { return KFParticleBase::GetMass(M, SigmaM); }

// Calculate the decay length of the particle in the laboratory system and its error. If they are well defined return 0, otherwise 1.
// The production point should be set before calling this function.
// \param[out] L - the decay length
// \param[out] SigmaL - its error
inline int KFParticle::GetDecayLength(float &L, float &SigmaL) const { return KFParticleBase::GetDecayLength(L, SigmaL); }

// Calculate the projection in the XY plane of the decay length of the particle in the laboratory
//  system and its error. If they are well defined return 0, otherwise 1.
//  The production point should be set before calling this function.
//  \param[out] L - the decay length
//  \param[out] SigmaL - its error
inline int KFParticle::GetDecayLengthXY(float &L, float &SigmaL) const { return KFParticleBase::GetDecayLengthXY(L, SigmaL); }

// Calculate the lifetime times speed of life (ctau) [cm] of the particle in the
// center of mass frame and its error. If they are well defined return 0, otherwise 1.
// The production point should be set before calling this function.
// \param[out] T - lifetime of the particle [cm]
// \param[out] SigmaT - its error
inline int KFParticle::GetLifeTime(float &T, float &SigmaT) const { return KFParticleBase::GetLifeTime(T, SigmaT); }

// Calculate the distance to the point {0,0,0} and its error. If they are well defined return 0, otherwise 1.
// \param[out] R - polar angle of the particle
// \param[out] SigmaR - its error
inline int KFParticle::GetR(float &R, float &SigmaR) const { return KFParticleBase::GetR(R, SigmaR); }

inline float &KFParticle::X() { return KFParticleBase::X(); }
inline float &KFParticle::Y() { return KFParticleBase::Y(); }
inline float &KFParticle::Z() { return KFParticleBase::Z(); }
inline float &KFParticle::Px() { return KFParticleBase::Px(); }
inline float &KFParticle::Py() { return KFParticleBase::Py(); }
inline float &KFParticle::Pz() { return KFParticleBase::Pz(); }
inline float &KFParticle::E() { return KFParticleBase::E(); }
inline float &KFParticle::S() { return KFParticleBase::S(); }
inline char &KFParticle::Q() { return KFParticleBase::Q(); }
inline float &KFParticle::Chi2() { return KFParticleBase::Chi2(); }
inline int &KFParticle::NDF() { return KFParticleBase::NDF(); }
inline float &KFParticle::Parameter(int i) { return KFParticleBase::Parameter(i); }
inline float &KFParticle::Covariance(int i) { return KFParticleBase::Covariance(i); }
inline float &KFParticle::Covariance(int i, int j) { return KFParticleBase::Covariance(i, j); }
inline float *KFParticle::Parameters() { return fP; }
inline float *KFParticle::CovarianceMatrix() { return fC; }

// Operator to add daughter to the current particle. Calls AddDaughter() function.
// \param[in] Daughter - the daughter particle
inline void KFParticle::operator+=(const KFParticle &Daughter) {
#ifdef NonhomogeneousField
    for (int i = 0; i < 10; i++) SetFieldCoeff(Daughter.GetFieldCoeff()[i], i);
#endif
    KFParticleBase::operator+=(Daughter);
}

// Adds daughter to the current particle. Depending on the selected construction method uses: \n
// 1) Either simplifyed fast mathematics which consideres momentum and energy as
// independent variables and thus ignores constraint on the fixed mass (fConstructMethod = 0).
// In this case the mass of the daughter particle can be corrupted when the constructed vertex
// is added as the measurement and the mass of the output short-lived particle can become
// unphysical - smaller then the threshold. Implemented in the
// AddDaughterWithEnergyFit() function \n
// 2) Or slower but correct mathematics which requires that the masses of daughter particles
// stays fixed in the construction process (fConstructMethod = 2). Implemented in the
// AddDaughterWithEnergyFitMC() function.
// \param[in] Daughter - the daughter particle
inline void KFParticle::AddDaughter(const KFParticle &Daughter) {
#ifdef NonhomogeneousField
    for (int i = 0; i < 10; i++) SetFieldCoeff(Daughter.GetFieldCoeff()[i], i);
#endif
    KFParticleBase::AddDaughter(Daughter);
}

// Constructs a short-lived particle from a set of daughter particles:
// 1) all parameters of the "this" objects are initialised;
// 2) daughters are added one after another;
// 3) if Parent pointer is not null, the production vertex is set to it;
// 4) if Mass hypothesis >=0 the mass constraint is set.
// \param[in] vDaughters - array of daughter particles
// \param[in] nDaughters - number of daughter particles in the input array
// \param[in] Parent - optional parrent particle
// \param[in] Mass - optional mass hypothesis
inline void KFParticle::Construct(const KFParticle *vDaughters[], int nDaughters, const KFParticle *ProdVtx, float Mass) {
#ifdef NonhomogeneousField
    for (int i = 0; i < 10; i++) SetFieldCoeff(vDaughters[0]->GetFieldCoeff()[i], i);
#endif
    KFParticleBase::Construct((const KFParticleBase **)vDaughters, nDaughters, (const KFParticleBase *)ProdVtx, Mass);
}

// Transports particle to the distance of closest approach to the point xyz.
// \param[in] xyz[3] - point, where particle should be transported
inline void KFParticle::TransportToPoint(const float xyz[]) {
    float dsdr[6] = {0.f};
    float dS = GetDStoPoint(xyz, dsdr);
    TransportToDS(dS, dsdr);
}
#ifdef HomogeneousField
// Transports particle to the distance of closest approach to the vertex v.
// \param[in] v - vertex, where particle should be transported
inline void KFParticle::TransportToVertex(const KFPVertex &v) { TransportToPoint(KFParticle(v).fP); }
#endif
// Transports particle to the distance of closest approach to the particle p.
// \param[in] p - particle, to which the current particle should be transported.
inline void KFParticle::TransportToParticle(const KFParticle &p) {
    float dsdr[4][6];
    float dS[2];
    GetDStoParticle(p, dS, dsdr);
    TransportToDS(dS[0], dsdr[0]);
}

// return dS = l/p parameter, where
// 1) l - signed distance to the DCA point with the input xyz point;
// 2) p - momentum of the particle;
// Also calculate partial derivatives dsdr of the parameter dS over the state vector of the current particle.
// If "HomogeneousField" is defined KFParticleBase::GetDStoPointBz() is called,
// if "NonhomogeneousField" is defined - KFParticleBase::GetDStoPointCBM()
// \param[in] xyz[3] - point, to which particle should be transported
// \param[out] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
// \param[in] param - optional parameter, is used in case if the parameters of the particle are rotated
// to other coordinate system (see GetDStoPointBy() function), otherwise fP are used
inline float KFParticle::GetDStoPoint(const float xyz[], float *dsdr) const {
#ifdef HomogeneousField
    return KFParticleBase::GetDStoPointBz(GetFieldAlice(), xyz, dsdr);
#endif
#ifdef NonhomogeneousField
    return KFParticleBase::GetDStoPointCBM(xyz, dsdr);
#endif
}

#ifdef HomogeneousField
// return value of the constant homogemeous one-component magnetic field Bz, (is defined in case of #ifdef HomogeneousField).
inline float KFParticle::GetFieldAlice() { return fgBz; }
#endif

#ifdef HomogeneousField
// Calculate the Bx, By, Bz components at the point xyz using approximation of the
// magnetic field along the particle trajectory.
// \param[in] xyz[3] - X, Y, Z coordiantes of the point where the magnetic field should be calculated
// \param[out] B[3] - value of X, Y, Z components of the calculated magnetic field at the given point
inline void KFParticle::GetFieldValue(const float * /*xyz*/, float B[]) const {

    B[0] = B[1] = 0;
    B[2] = GetFieldAlice();
}
#endif

#ifdef NonhomogeneousField
// Calculate the Bx, By, Bz components at the point xyz using approximation of the
// magnetic field along the particle trajectory.
// \param[in] xyz[3] - X, Y, Z coordiantes of the point where the magnetic field should be calculated
// \param[out] B[3] - value of X, Y, Z components of the calculated magnetic field at the given point
inline void KFParticle::GetFieldValue(const float xyz[], float B[]) const {

    const float dz = (xyz[2] - fieldRegion[9]);
    const float dz2 = dz * dz;

    B[0] = fieldRegion[0] + fieldRegion[1] * dz + fieldRegion[2] * dz2;
    B[1] = fieldRegion[3] + fieldRegion[4] * dz + fieldRegion[5] * dz2;
    B[2] = fieldRegion[6] + fieldRegion[7] * dz + fieldRegion[8] * dz2;
}
#endif

// Calculate dS = l/p parameters for two particles, where
// 1) l - signed distance to the DCA point with the other particle;
// 2) p - momentum of the particle
// dS[0] is the transport parameter for the current particle, dS[1] - for the particle "p".
// Also calculate partial derivatives dsdr of the parameters dS[0] and dS[1] over the state vectors of the particles:
// 1) dsdr[0][6] = d(dS[0])/d(param1);
// 2) dsdr[1][6] = d(dS[0])/d(param2);
// 3) dsdr[2][6] = d(dS[1])/d(param1);
// 4) dsdr[3][6] = d(dS[1])/d(param2);
// where param1 are parameters of the current particle fP and
// param2 are parameters of the second particle p.fP. If "HomogeneousField" is defined KFParticleBase::GetDStoParticleBz() is called,
// if "NonhomogeneousField" is defined - KFParticleBase::GetDStoParticleCBM()
// \param[in] p - second particle
// \param[out] dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
// \param[out] dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
inline void KFParticle::GetDStoParticle(const KFParticleBase &p, float dS[2], float dsdr[4][6]) const {
#ifdef HomogeneousField
    KFParticleBase::GetDStoParticleBz(GetFieldAlice(), p, dS, dsdr);
#endif
#ifdef NonhomogeneousField
    KFParticleBase::GetDStoParticleCBM(p, dS, dsdr);
#endif
}

// Transports the parameters and their covariance matrix of the current particle
// on a length defined by the transport parameter dS = l/p, where l is the signed distance and p is
// the momentum of the current particle. If "HomogeneousField" is defined KFParticleBase::TransportBz()
// is called, if "NonhomogeneousField" - KFParticleBase::TransportCBM().
// The obtained parameters and covariance matrix are stored to the arrays P and
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
inline void KFParticle::Transport(float dS, const float *dsdr, float P[], float C[], float *dsdr1, float *F, float *F1) const {
#ifdef HomogeneousField
    KFParticleBase::TransportBz(GetFieldAlice(), dS, dsdr, P, C, dsdr1, F, F1);
#endif
#ifdef NonhomogeneousField
    KFParticleBase::TransportCBM(dS, dsdr, P, C, dsdr1, F, F1);
#endif
}

#endif
