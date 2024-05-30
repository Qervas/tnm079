/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 * Acknowledgements for original code base:
 * - Gunnar Johansson
 * - Ken Museth
 * - Michael Bang Nielsen
 * - Ola Nilsson
 * - Andreas Soderstrom
 *
 * Code updated in the period 2017-2018 by Jochen Jankowai
 *
 *************************************************************************************************/
#include <Geometry/Quadric.h>

Quadric::Quadric(const glm::mat4& q) : mQuadric(q) {}

Quadric::~Quadric() {}

/*!
 * Evaluation of world coordinates are done through either transformation
 * of the world-coordinates by mWorld2Obj, or transformation of the quadric
 * coefficient matrix by GetTransform() ONCE (see Section 2.2 in lab text).
 */
float Quadric::GetValue(float x, float y, float z) const {
    TransformW2O(x, y, z);
    glm::vec4 p(x, y, z, 1.0f);
    float value = glm::dot(p, mQuadric * p);
    return value;
}

/*!
 * Use the quadric matrix to evaluate the gradient.
 */
glm::vec3 Quadric::GetGradient(float x, float y, float z) const {
    TransformW2O(x, y, z);
    glm::mat3 Qsub(mQuadric);
    glm::vec3 translation(mQuadric[3]);
    glm::vec3 p(x, y, z);
    glm::vec3 grad = 2.0f *( Qsub * p + translation);
    return grad;
}
