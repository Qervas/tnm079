/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#pragma once

#include "Levelset/LevelSetOperator.h"
#include "Util/Stopwatch.h"

/*! \brief A level set operator that does morphing
 */
//! \lab5 Implement morphing
class OperatorMorph : public LevelSetOperator {
protected:
    const Implicit* mTarget;
    //! Stopwatch for performance measurements
    Stopwatch mStopwatch;

public:
    OperatorMorph(LevelSet* LS, const Implicit* target) : LevelSetOperator(LS), mTarget(target) {}

    virtual float ComputeTimestep() {
        float dx = mLS->GetDx();
        float maxSpeed = 0.0f;
        
        // Find maximum speed (difference between target and current) in narrow band
        LevelSetGrid::Iterator iter = mLS->BeginNarrowBand();
        LevelSetGrid::Iterator iend = mLS->EndNarrowBand();
        while (iter != iend) {
            size_t i = iter.GetI();
            size_t j = iter.GetJ();
            size_t k = iter.GetK();
            float x = i, y = j, z = k;
            mLS->TransformGridToWorld(x, y, z);
            float currentPhi = mLS->GetValue(x, y, z);
            float targetPhi = mTarget->GetValue(x, y, z);
            float speed = std::abs(currentPhi - targetPhi);
            if (speed > maxSpeed) maxSpeed = speed;
            iter++;
        }
        
        // Ensure stable timestep
        if (maxSpeed < 1e-6f) maxSpeed = 1.0f;
        return 0.5f * dx / maxSpeed;  // Conservative factor
    }

    virtual void Propagate(float time) {
        // Start timing
        mStopwatch.start();

        // Propagate level set with stable timestep dt
        // until requested time is reached
        for (float elapsed = 0.f; elapsed < time;) {
            // Determine timestep for stability
            float dt = ComputeTimestep();

            if (dt > time - elapsed) {
                dt = time - elapsed;
            }
            elapsed += dt;

            // Integrate level set function in time using Euler integration
            IntegrateEuler(dt);
        }
        
        // Stop timing and report
        double elapsedTime = mStopwatch.stop();
        
        // Get narrow band information
        int width = mLS->GetNarrowBandWidth();
        std::string bandStatus = (width > 0) ? "enabled" : "disabled";
        
        // Report performance
        std::cout << "Morphing Operation Performance:" << std::endl;
        std::cout << "  Narrow band: " << bandStatus << std::endl;
        if (width > 0) {
            std::cout << "  Band width: " << width << std::endl;
        }
        std::cout << "  Elapsed time: " << elapsedTime << " seconds" << std::endl;
    }

    virtual float Evaluate(size_t i, size_t j, size_t k) {
        // Convert grid coordinates to world coordinates
        float x = i, y = j, z = k;
        mLS->TransformGridToWorld(x, y, z);
        
        // Get current and target level set values
        float currentPhi = mLS->GetValue(x, y, z);  // Use world coordinates
        float targetPhi = mTarget->GetValue(x, y, z);
        
        // Speed function: difference between current and target (flipped)
        // This makes the surface move away from target toward current shape
        float speed = currentPhi - targetPhi;
        
        // Use existing Godunov method for upwind differencing
        float ddx2, ddy2, ddz2;
        Godunov(i, j, k, speed, ddx2, ddy2, ddz2);
        
        // Compute gradient magnitude
        float gradientMagnitude = std::sqrt(ddx2 + ddy2 + ddz2);
        
        // Return rate of change: ∂ϕ/∂t = -F|∇ϕ|
        // where F is the target distance (speed function)
        return -speed * gradientMagnitude;
    }
};
