#include "QuadricDecimationMesh.h"
#include <VC++/glm/gtx/euler_angles.hpp>

const QuadricDecimationMesh::VisualizationMode QuadricDecimationMesh::QuadricIsoSurfaces =
    NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize() {
    // Allocate memory for the quadric array
    size_t numVerts = mVerts.size();
    mQuadrics.reserve(numVerts);
    std::streamsize width = std::cerr.precision();  // store stream precision
    for (size_t i = 0; i < numVerts; i++) {

        // Compute quadric for vertex i here
        mQuadrics.push_back(createQuadricForVert(i));

        // Calculate initial error, should be numerically close to 0

        glm::vec3 v0 = mVerts[i].pos;
        glm::vec4 v(v0[0], v0[1], v0[2], 1);
        glm::mat4 m = mQuadrics.back();

        // TODO CHECK
        float error = glm::dot(v, (m * v));
        // std::cerr << std::scientific << std::setprecision(2) << error << " ";
    }
    std::cerr << std::setprecision(width) << std::fixed;  // reset stream precision

    // Run the initialize for the parent class to initialize the edge collapses
    DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute,
 * DecimationMesh::EdgeCollapse
 */
// Existing includes and class definition...

#include "QuadricDecimationMesh.h"
#include <VC++/glm/gtc/matrix_inverse.hpp>
#include <VC++/glm/glm.hpp>
#include <cmath>

void QuadricDecimationMesh::computeCollapse(EdgeCollapse* collapse) {
    // Retrieve the vertices at the endpoints of the edge to be collapsed
    size_t v1 = e(collapse->halfEdge).vert;
    size_t v2 = e(e(collapse->halfEdge).next).vert;

    glm::vec4 posV1(v(v1).pos, 1.0f);
    glm::vec4 posV2(v(v2).pos, 1.0f);
    glm::vec4 between = 0.5f * (posV1 + posV2);

    // Create and combine quadrics for the vertices
    auto Q1 = createQuadricForVert(v1);
    auto Q2 = createQuadricForVert(v2);
    auto Q = Q1 + Q2;

    // Modify Q for solving the position
    auto Qhat = Q;
    Qhat[0][3] = Qhat[1][3] = Qhat[2][3] = 0.0f;
    Qhat[3][3] = 1.0f;

    glm::vec4 zero(0.0f, 0.0f, 0.0f, 1.0f);

    float EPSILON = 1e-9f;
    bool notInvertible = abs(glm::determinant(Qhat)) < EPSILON;

    float costV1 = glm::dot(posV1, Q * posV1);
    float costV2 = glm::dot(posV2, Q * posV2);
    float costBetween = glm::dot(between, Q * between);

    if (!notInvertible) {
        glm::vec4 v = glm::inverse(Qhat) * zero;
        collapse->position = glm::vec3(v);
    } else {
        if (costV1 < costV2 && costV1 < costBetween) {
            collapse->position = glm::vec3(posV1);
        } else if (costV2 < costV1 && costV2 < costBetween) {
            collapse->position = glm::vec3(posV2);
        } else {
            collapse->position = glm::vec3(between);
        }
    }

    glm::vec4 position(collapse->position, 1.0f);
    float deltaV = glm::dot(position, Q * position);
    collapse->cost = deltaV;
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(size_t ind) {
    DecimationMesh::updateVertexProperties(ind);
    mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
glm::mat4 QuadricDecimationMesh::createQuadricForVert(size_t indx) const {
    glm::mat4 Q(0.0f);
    auto vertex = v(indx);
    auto neighboring_faces = FindNeighborFaces(indx);
    for(auto & F_indx : neighboring_faces){
        Q += createQuadricForFace(F_indx);
    } 
    // The quadric for a vertex is the sum of all the quadrics for the adjacent
    // faces Tip: Matrix4x4 has an operator +=
    return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
glm::mat4 QuadricDecimationMesh::createQuadricForFace(size_t indx) const {
    const Face& F = f(indx);
    const glm::vec3& normal = F.normal;
    const Vertex& v0 = v(e(F.edge).vert); // any vertex in the triangle

    float a = normal.x;
    float b = normal.y;
    float c = normal.z;
    float d = -glm::dot(normal, v0.pos);

    glm::mat4 K(
        {a * a, a * b, a * c, a * d},
        {a * b, b * b, b * c, b * d},
        {a * c, b * c, c * c, c * d},
        {a * d, b * d, c * d, d * d}
    );
    
    return K;
}

void QuadricDecimationMesh::Render() {
    DecimationMesh::Render();

    glEnable(GL_LIGHTING);
    glMatrixMode(GL_MODELVIEW);

    if (mVisualizationMode == QuadricIsoSurfaces) {
        // Apply transform
        glPushMatrix();  // Push modelview matrix onto stack

        // Implement the quadric visualization here
        std::cout << "Quadric visualization not implemented" << std::endl;

        // Restore modelview matrix
        glPopMatrix();
    }
}
