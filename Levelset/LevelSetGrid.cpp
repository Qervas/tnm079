#include <Levelset/LevelSetGrid.h>

LevelSetGrid::Iterator::Iterator(const BitMask3D* mask, size_t i, size_t j, size_t k) {
    // Initialize counters
    this->mask = mask;
    this->iMax = mask->GetDimX();
    this->jMax = mask->GetDimY();
    this->kMax = mask->GetDimZ();
    this->i = glm::clamp(i, size_t{0}, iMax);
    this->j = glm::clamp(j, size_t{0}, jMax);
    this->k = glm::clamp(k, size_t{0}, kMax);

    // Check if we're within range
    if (i == iMax && j == jMax && k == kMax) {
        endState = true;
    } else {
        endState = false;
    }

    // Go to the first valid entry
    if (!endState && !mask->GetValue(i, j, k)) (*this)++;
}

void LevelSetGrid::Dilate() {
    BitMask3D newMask(mMask);
    Iterator it = BeginNarrowBand();
    Iterator iend = EndNarrowBand();
    while (it != iend) {
        size_t i = it.GetI();
        size_t j = it.GetJ();
        size_t k = it.GetK();
        newMask.SetValue(i, j, k, true);
        if (k < GetDimZ() - 1) {
            newMask.SetValue(i, j, k + 1, true);
        }
        if (k > 0) {
            newMask.SetValue(i, j, k - 1, true);
        }
        if (j < GetDimY() - 1) {
            newMask.SetValue(i, j + 1, k, true);
        }
        if (j > 0) {
            newMask.SetValue(i, j - 1, k, true);
        }
        if (i < GetDimX() - 1) {
            newMask.SetValue(i + 1, j, k, true);
        }
        if (i > 0) {
            newMask.SetValue(i - 1, j, k, true);
        }
        it++;
    }
    mMask = newMask;
}

void LevelSetGrid::Rebuild() {
    // Initialize a new mask to track cells in the narrow band
    BitMask3D newMask(mMask.GetDimX(), mMask.GetDimY(), mMask.GetDimZ());
    
    // Check if narrow band is disabled (extreme constant values)
    bool narrowBandDisabled = (mInsideConstant <= -std::numeric_limits<float>::max() / 2.0f) && 
                             (mOutsideConstant >= std::numeric_limits<float>::max() / 2.0f);
    
    if (narrowBandDisabled) {
        // If narrow band is disabled, include all cells
        for (size_t i = 0; i < mPhi.GetDimX(); i++) {
            for (size_t j = 0; j < mPhi.GetDimY(); j++) {
                for (size_t k = 0; k < mPhi.GetDimZ(); k++) {
                    newMask.SetValue(i, j, k, true);
                }
            }
        }
    } else {
        // First pass: identify cells that should be in the narrow band
        for (size_t i = 0; i < mPhi.GetDimX(); i++) {
            for (size_t j = 0; j < mPhi.GetDimY(); j++) {
                for (size_t k = 0; k < mPhi.GetDimZ(); k++) {
                    float value = mPhi.GetValue(i, j, k);
                    
                    // Check if the value is within the narrow band range
                    if (value >= mInsideConstant && value <= mOutsideConstant) {
                        newMask.SetValue(i, j, k, true);
                    }
                }
            }
        }
        
        // Second pass: clamp values outside the narrow band
        for (size_t i = 0; i < mPhi.GetDimX(); i++) {
            for (size_t j = 0; j < mPhi.GetDimY(); j++) {
                for (size_t k = 0; k < mPhi.GetDimZ(); k++) {
                    if (!newMask.GetValue(i, j, k)) {
                        // If outside the narrow band, set to appropriate constant
                        if (mPhi.GetValue(i, j, k) > 0) {
                            mPhi.SetValue(i, j, k, mOutsideConstant);
                        } else {
                            mPhi.SetValue(i, j, k, mInsideConstant);
                        }
                    }
                }
            }
        }
        
        // Dilate the narrow band to ensure we have enough cells for gradient calculations
        mMask = newMask;
        Dilate();
    }
    
    // Update the mask
    if (narrowBandDisabled) {
        mMask = newMask;
    }
}

glm::ivec3 LevelSetGrid::GetDimensions() { return glm::ivec3(GetDimX(), GetDimY(), GetDimZ()); }
