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
        // Compute and return a stable timestep
        return 1.f;
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
        // Compute the rate of change (dphi/dt)
        
        // 1. Get the current position in world coordinates
        float x = i;
        float y = j;
        float z = k;
        mLS->TransformGridToWorld(x, y, z);
        
        // 2. Get the current level set value (phi) at this position
        float currentPhi = mLS->GetValue(i, j, k);
        
        // 3. Get the target implicit value at this position
        float targetPhi = mTarget->GetValue(x, y, z);
        
        // 4. Calculate the gradient magnitude of the level set
        float dx = mLS->DiffXpm(i, j, k);
        float dy = mLS->DiffYpm(i, j, k);
        float dz = mLS->DiffZpm(i, j, k);
        float gradientMagnitude = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        // 5. Implement morphing as a speed function based on the difference
        // between the current level set and the target implicit
        // The speed is proportional to the difference between the two surfaces
        float speed = targetPhi - currentPhi;
        
        // 6. Return the rate of change: speed * gradient magnitude
        // This follows the level set equation: dϕ/dt = speed * |∇ϕ|
        return speed * gradientMagnitude;
    }
};
