/*

Copyright (c) 2005-2022, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "FixedDurationCellCycleModel.hpp"

#include <iostream>

#include "RandomNumberGenerator.hpp"

template<class Archive>
void FixedDurationCellCycleModel::serialize(Archive & archive, const unsigned int version)
{
    // Archive cell-cycle model using serialization code from AbstractSimplePhaseBasedCellCycleModel
    archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);
}

FixedDurationCellCycleModel::FixedDurationCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModel()
{
  SetPhaseDurations();
  mPhaseTimer = 0.0;
  mLastCellAge = 0.0;
}

void FixedDurationCellCycleModel::SetPhaseDurations()
{
    SetStemCellG1Duration(7.0);
    SetTransitCellG1Duration(7.0);
    SetSDuration(6.0);
    SetG2Duration(3.0);
    SetMDuration(2.0);
    SetMinimumGapDuration(3.0);
}

void FixedDurationCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);  // Make sure cell exists

    mG1Duration = 7.0;
}

AbstractCellCycleModel* FixedDurationCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    FixedDurationCellCycleModel* pCellCycleModel = new FixedDurationCellCycleModel();
    return pCellCycleModel;
}


void FixedDurationCellCycleModel::SetPhaseTimer(const double value) {
    mPhaseTimer = value;
}

bool FixedDurationCellCycleModel::ReadyToDivide() {
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();
        if ((mCurrentCellCyclePhase != G_ZERO_PHASE) &&
            (mPhaseTimer >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()) &&
            !mpCell->GetCellData()->GetItem("growth inhibited")) 
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}

void FixedDurationCellCycleModel::UpdateCellCyclePhase()
{
    double timeSinceBirth = GetAge();
    assert(timeSinceBirth >= 0);
    
    // Update cell growth phase timer
    double change_in_cell_age = timeSinceBirth - mLastCellAge;
    mPhaseTimer += change_in_cell_age;
    mLastCellAge = timeSinceBirth;

    // Select the correct phase
    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (mPhaseTimer < GetG1Duration())
    {
        if (!mpCell->GetCellData()->GetItem("growth inhibited")) {
            mCurrentCellCyclePhase = G_ONE_PHASE;
        } else {
            std::cout << "Cell inhibited\n";
            if (mCurrentCellCyclePhase != G_ONE_PHASE) {
                mPhaseTimer -= change_in_cell_age;
            }
        }
    }
    else if (mPhaseTimer <  GetG1Duration() + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if (mPhaseTimer < GetG1Duration() + GetSDuration() + GetG2Duration())
    {
        if (!mpCell->GetCellData()->GetItem("growth inhibited")) {
            mCurrentCellCyclePhase = G_TWO_PHASE;
        } else {
            std::cout << "Cell inhibited\n";
            if (mCurrentCellCyclePhase != G_TWO_PHASE) {
                mPhaseTimer -= change_in_cell_age;
            }
        }
    }
    else if (mPhaseTimer < GetG1Duration() + GetSDuration() + GetG2Duration() + GetMDuration())
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedDurationCellCycleModel)
