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

#ifndef KFVERTEX_H
#define KFVERTEX_H

#include "KFPVertex.h"
#include "KFParticleBase.h"

// @class KFVertex
// @brief Mathematics for reconstruction of primary vertices based on KFParticle.
// @author S.Gorbunov, I.Kisel, M.Zyzak
// @date 05.02.2019
// @version 1.0
//
// The class is inherited from KFParticle, adds functionality for reconstruction of primary vertices.
class KFVertex : public KFParticleBase {
   public:
    KFVertex() : KFParticleBase(), fIsConstrained(0) {}
    KFVertex(const KFParticleBase &particle)
        : KFParticleBase(particle), fIsConstrained(0) {}  // Vertex is constructed from the current position of a given particle
    KFVertex(const KFPVertex &vertex);
    ~KFVertex() = default;

    // Return number of particles used for construction of the vertex
    int GetNContributors() const { return fIsConstrained ? fNDF / 2 : (fNDF + 3) / 2; }

    void SetBeamConstraint(float X, float Y, float Z, float ErrX, float ErrY, float ErrZ);
    void SetBeamConstraintOff();

    void ConstructPrimaryVertex(float bz, const KFParticleBase *vDaughters[], int nDaughters, bool vtxFlag[], float ChiCut = 3.5);

   protected:
    bool fIsConstrained;  // Flag showing if the the beam constraint is set
};

#endif
