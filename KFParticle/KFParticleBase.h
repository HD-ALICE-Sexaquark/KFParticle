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

#ifdef __ROOT__  // for the STAR experiment
#define HomogeneousField
#endif

#include <cmath>
#include <iostream>
#include <vector>

// @class KFParticleBase
// @brief The base of KFParticle class, describes particle objects.
// @author  S.Gorbunov, I.Kisel, M.Zyzak
// @date 05.02.2019
// @version 1.0
//
// Contains the main mathematics of the KFParticle.
class KFParticleBase {
   public:
    // Abstract methods are defined in the KFParticle class
    // Virtual method to access the magnetic field
    virtual void GetFieldValue(const float xyz[], float B[]) const = 0;

    // Virtual method to get extrapolation parameter dS=l/p to. Is defined in KFParticle.
    virtual float GetDStoPoint(const float xyz[3], float dsdr[6]) const = 0;

    float GetDStoPointLine(const float xyz[3], float dsdr[6]) const;
    float GetDStoPointBz(float B, const float xyz[3], float dsdr[6], const float *param = 0) const;
    float GetDStoPointBy(float By, const float xyz[3], float dsdr[6]) const;
    float GetDStoPointB(const float *B, const float xyz[3], float dsdr[6]) const;
    float GetDStoPointCBM(const float xyz[3], float dsdr[6]) const;

    // Virtual method to get extrapolation parameter dS=l/p to another particle. Is defined in KFParticle.
    virtual void GetDStoParticle(const KFParticleBase &p, float dS[2], float dsdr[4][6]) const = 0;

    void GetDStoParticleLine(const KFParticleBase &p, float dS[2], float dsdr[4][6]) const;
    void GetDStoParticleBz(float Bz, const KFParticleBase &p, float dS[2], float dsdr[4][6], const float *param1 = 0, const float *param2 = 0) const;
    void GetDStoParticleBy(float B, const KFParticleBase &p, float dS[2], float dsdr[4][6]) const;
    void GetDStoParticleCBM(const KFParticleBase &p, float dS[2], float dsdr[4][6]) const;

    // Virtual method to transport a particle on a certain distance along the trajectory. Is defined in KFParticle.
    virtual void Transport(float dS, const float dsdr[6], float P[], float C[], float *dsdr1 = 0, float *F = 0, float *F1 = 0) const = 0;

    KFParticleBase();
    virtual ~KFParticleBase() { ; }  // The default destructor

    void Initialize(const float Param[], const float Cov[], int Charge, float Mass);
    void Initialize();

    // define the construction method for the current particle (see description of fConstructMethod)
    void SetConstructMethod(int m) { fConstructMethod = m; }
    void SetMassHypo(float m) { fMassHypo = m; }                         // set the mass hypothesis to the particle, is used when fConstructMethod = 2
    const float &GetMassHypo() const { return fMassHypo; }               // return the mass hypothesis
    const float &GetSumDaughterMass() const { return SumDaughterMass; }  // return the sum of masses of the daughters

    //  ACCESSORS

    // Simple accessors

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

    // Accessors with calculations( &value, &estimated sigma )
    // error flag returned (0 means no error during calculations)

    int GetMomentum(float &p, float &error) const;
    int GetPt(float &pt, float &error) const;
    int GetEta(float &eta, float &error) const;
    int GetPhi(float &phi, float &error) const;
    int GetMass(float &m, float &error) const;
    int GetDecayLength(float &l, float &error) const;
    int GetDecayLengthXY(float &l, float &error) const;
    int GetLifeTime(float &ctau, float &error) const;
    int GetR(float &r, float &error) const;

    //  MODIFIERS
    float &X() { return fP[0]; }     // Modifier of X coordinate of the particle, fP[0]
    float &Y() { return fP[1]; }     // Modifier of Y coordinate of the particle, fP[1]
    float &Z() { return fP[2]; }     // Modifier of Z coordinate of the particle, fP[2]
    float &Px() { return fP[3]; }    // Modifier of X component of the momentum, fP[3]
    float &Py() { return fP[4]; }    // Modifier of Y component of the momentum, fP[4]
    float &Pz() { return fP[5]; }    // Modifier of Z component of the momentum, fP[5]
    float &E() { return fP[6]; }     // Modifier of energy of the particle, fP[6]
    float &S() { return fP[7]; }     // Modifier of dS=l/p, l - decay length, fP[7], defined if production vertex is set
    char &Q() { return fQ; }         // Modifier of charge of the particle
    float &Chi2() { return fChi2; }  // Modifier of Chi2 of the fit
    int &NDF() { return fNDF; }      // Modifier of number of decrease of freedom

    float &Parameter(int i) { return fP[i]; }                 // Modifier of P[i] parameter
    float &Covariance(int i) { return fC[i]; }                // Modifier of C[i] element of the covariance matrix in the lower triangular form
    float &Covariance(int i, int j) { return fC[IJ(i, j)]; }  // Modifier of C[i,j] element of the covariance matrix

    // CONSTRUCTION OF THE PARTICLE BY ITS DAUGHTERS AND MOTHER USING THE KALMAN FILTER METHOD

    // Simple way to add daughter ex. D0+= Pion;
    void operator+=(const KFParticleBase &Daughter);

    // Add daughter track to the particle
    void AddDaughter(const KFParticleBase &Daughter);
    void SubtractDaughter(const KFParticleBase &Daughter);

    void AddDaughterWithEnergyFit(const KFParticleBase &Daughter);
    void AddDaughterWithEnergyFitMC(const KFParticleBase &Daughter);

    // Set production vertex
    void SetProductionVertex(const KFParticleBase &Vtx);

    // Set mass constraint
    void SetNonlinearMassConstraint(float Mass);
    void SetMassConstraint(float Mass, float SigmaMass = 0);

    // Set no decay length for resonances
    void SetNoDecayLength();

    // Everything in one go
    void Construct(const KFParticleBase *vDaughters[], int nDaughters, const KFParticleBase *ProdVtx = 0, float Mass = -1);

    // Transport functions
    void TransportToDecayVertex();
    void TransportToProductionVertex();
    void TransportToDS(float dS, const float *dsdr);
    void TransportBz(float Bz, float dS, const float *dsdr, float P[], float C[], float *dsdr1 = 0, float *F = 0, float *F1 = 0) const;
    void TransportCBM(float dS, const float *dsdr, float P[], float C[], float *dsdr1 = 0, float *F = 0, float *F1 = 0) const;

    // OTHER UTILITIES

    // Calculate distance from another object [cm]
    float GetDistanceFromVertex(const float vtx[]) const;
    float GetDistanceFromVertex(const KFParticleBase &Vtx) const;
    float GetDistanceFromParticle(const KFParticleBase &p) const;

    // Calculate sqrt(Chi2/ndf) deviation from vertex
    // v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix
    float GetDeviationFromVertex(const float v[], const float Cv[] = 0) const;
    float GetDeviationFromVertex(const KFParticleBase &Vtx) const;
    float GetDeviationFromParticle(const KFParticleBase &p) const;

    void SubtractFromVertex(KFParticleBase &Vtx) const;
    void SubtractFromParticle(KFParticleBase &Vtx) const;

    static void GetArmenterosPodolanski(KFParticleBase &positive, KFParticleBase &negative, float QtAlfa[2]);
    void RotateXY(float angle, float Vtx[3]);

    int Id() const { return fId; }                                         // return Id of the particle
    int NDaughters() const { return fDaughtersIds.size(); }                // return number of daughter particles
    const std::vector<int> &DaughterIds() const { return fDaughtersIds; }  // return the vector with the indices of daughter particles
    void CleanDaughtersId() { fDaughtersIds.clear(); }                     // clean the vector with the indices of daughter particles

    void SetId(int id) { fId = id; }  // set the Id of the particle. After the construction of a particle should be set by user.
    void AddDaughterId(int id) { fDaughtersIds.push_back(id); }  // add index of the daughter particle

    void SetPDG(int pdg) { fPDG = pdg; }  // set the PDG hypothesis
    int GetPDG() const { return fPDG; }   // return the PDG hypothesis

#ifdef __ROOT__  // for the STAR experiment
    virtual void Print(Option_t *opt = "") const;
    int IdTruth() const { return fIdTruth; }
    int QaTruth() const { return fQuality; }
    int IdParentMcVx() const { return fIdParentMcVx; }
    int IdParentVx() const { return IdParentMcVx(); }
    void SetParentID(int id = 0) { fParentID = id; }
    int GetParentID() const { return fParentID; }
    void SetIdParentMcVx(int id) { fIdParentMcVx = id; }
    void SetIdTruth(int idtru, int qatru = 0) {
        fIdTruth = (unsigned short)idtru;
        fQuality = (unsigned short)qatru;
    }
    virtual void Clear(Option_t * /*option*/ = "");
#endif

    static void InvertCholetsky3(float a[6]);
    static void MultQSQt(const float Q[], const float S[], float SOut[], const int kN);

    template <typename D, size_t N>
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

    template <typename D, size_t N>
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

    template <typename D, size_t N>
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
    template <typename D, size_t N, size_t M>
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

    template <typename D, size_t N>
    void PrintJoinedMatrix(std::string_view name, D arr1[N], D arr2[N], D arr3[N]) const {
        std::cout << "  " << name << " =\n";
        for (size_t i = 0; i < N; ++i) {
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
    bool GetMeasurement(const KFParticleBase &daughter, float m[], float V[], float D[3][3]);
    void SetMassConstraint(float *mP, float *mC, float mJ[7][7], float mass);

    float fP[8];            // Particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}
    float fC[36];           // Low-triangle covariance matrix of fP
    float fChi2;            // Chi^2
    float fSFromDecay;      // Distance from the decay vertex to the current position
    float SumDaughterMass;  // Sum of the daughter particles masses Needed to set the constraint on the minimum mass during particle construction
    float fMassHypo;        // The mass hypothesis, used for the constraints during particle construction
    int fNDF;               // Number of degrees of freedom
    int fId;                // Id of the particle
#ifdef __ROOT__             // for the STAR experiment
    short fParentID;        // Id of the parent particle
    short fIdTruth;         // MC track id
    short fQuality;         // quality of this information (percentage of hits coming from the above MC track)
    short fIdParentMcVx;    // for track and McTrack for vertex
#endif
    bool fAtProductionVertex;  // Flag shows if particle is at the production point
    char fQ;                   // The charge of the particle in the units of the elementary charge

    // \brief Determines the method for the particle construction.
    // 0 - Energy considered as an independent veriable, fitted independently from momentum, without any constraints on mass
    // 2 - Energy considered as an independent variable, fitted independently from momentum, with constraints on mass of daughter particle
    char fConstructMethod;
    int fPDG;  // The PDG hypothesis assigned to the particle.

    // \brief A vector with ids of the daughter particles:
    // 1) if particle is created from a track - the index of the track, in this case the size of the vector is always equal to one;
    // 2) if particle is constructed from other particles - indices of these particles in the same array.
    std::vector<int> fDaughtersIds;
};

#ifdef __ROOT__  // for the STAR experiment
std::ostream &operator<<(std::ostream &os, KFParticleBase const &particle);
#endif

#endif
