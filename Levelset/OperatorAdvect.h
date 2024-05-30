#pragma once

#include "Levelset/LevelSetOperator.h"
#include "Math/Function3D.h"

/*! \brief A level set operator that does external advection
 *
 * This class implements level set advectionr in an external vector field by the
 * PDE
 *
 *  \f$
 *  \dfrac{\partial \phi}{\partial t} + \mathbf{V}(\mathbf{x})\cdot \nabla \phi
 * = 0 \f$
 */
//! \lab4 Implement advection in external vector field
class OperatorAdvect : public LevelSetOperator {
protected:
    Function3D<glm::vec3>* mVectorField;

public:
    OperatorAdvect(LevelSet* LS, Function3D<glm::vec3>* vf)
        : LevelSetOperator(LS), mVectorField(vf) {}

    virtual float ComputeTimestep() {
        // Compute and return a stable timestep
        // (Hint: Function3D::GetMaxValue())
        float maxVelocity = glm::length(mVectorField->GetMaxValue());
        float dx = mLS->GetDx();
        return dx / maxVelocity;
    }

    virtual void Propagate(float time) {
        // Determine timestep for stability
        float dt = ComputeTimestep();

        // Propagate level set with stable timestep dt
        // until requested time is reached
        for (float elapsed = 0.f; elapsed < time;) {
            if (dt > time - elapsed) {
                dt = time - elapsed;
            }
            elapsed += dt;

            // IntegrateEuler(dt);
            IntegrateRungeKutta(dt);
        }
    }

    virtual float Evaluate(size_t i, size_t j, size_t k) {
        // Compute the rate of change (dphi/dt)

        // Remember that the point (i,j,k) is given in grid coordinates, while
        // the velocity field used for advection needs to be sampled in
        // world coordinates (x,y,z). You can use LevelSet::TransformGridToWorld()
        // for this task.
            // Transform grid coordinates to world coordinates
        float x = static_cast<float>(i) * mLS->GetDx();
        float y = static_cast<float>(j) * mLS->GetDx();
        float z = static_cast<float>(k) * mLS->GetDx();
        mLS->TransformGridToWorld(x, y, z);

        // Sample the external velocity field
        glm::vec3 velocity = mVectorField->GetValue(x, y, z);

        // Compute the rate of change using the dot product of the velocity field and the gradient
        float phi_x = mLS->DiffXp(i, j, k) - mLS->DiffXm(i, j, k);
        float phi_y = mLS->DiffYp(i, j, k) - mLS->DiffYm(i, j, k);
        float phi_z = mLS->DiffZp(i, j, k) - mLS->DiffZm(i, j, k);

        float rateOfChange = velocity.x * phi_x + velocity.y * phi_y + velocity.z * phi_z;

        return -rateOfChange;
    }
};