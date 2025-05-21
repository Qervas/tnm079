#include "QuadricDecimationMesh.h"
#include <VC++/glm/gtx/euler_angles.hpp>
#include <GUI/GLViewer.h>
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

    // Retrieve the camera position
    glm::vec3 cameraPos = GLViewer::GetCamera().GetPosition();

    // Calculate distance from camera to vertices
    float distV1 = glm::distance(glm::vec3(posV1), cameraPos);
    float distV2 = glm::distance(glm::vec3(posV2), cameraPos);
    float distBetween = glm::distance(glm::vec3(between), cameraPos);

    // Introduce a weighting factor for distance
    float distWeight = 0.1f;  
    costV1 += distWeight * distV1;
    costV2 += distWeight * distV2;
    costBetween += distWeight * distBetween;

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
// void QuadricDecimationMesh::computeCollapse(EdgeCollapse* collapse) {
//     // Retrieve the vertices at the endpoints of the edge to be collapsed
//     size_t v1 = e(collapse->halfEdge).vert;
//     size_t v2 = e(e(collapse->halfEdge).next).vert;

//     glm::vec4 posV1(v(v1).pos, 1.0f);
//     glm::vec4 posV2(v(v2).pos, 1.0f);
//     glm::vec4 between = 0.5f * (posV1 + posV2);

//     // Create and combine quadrics for the vertices
//     auto Q1 = createQuadricForVert(v1);
//     auto Q2 = createQuadricForVert(v2);
//     auto Q = Q1 + Q2;

//     // Modify Q for solving the position
//     auto Qhat = Q;
//     Qhat[0][3] = Qhat[1][3] = Qhat[2][3] = 0.0f;
//     Qhat[3][3] = 1.0f;

//     glm::vec4 zero(0.0f, 0.0f, 0.0f, 1.0f);

//     float EPSILON = 1e-9f;
//     bool notInvertible = abs(glm::determinant(Qhat)) < EPSILON;

//     float costV1 = glm::dot(posV1, Q * posV1);
//     float costV2 = glm::dot(posV2, Q * posV2);
//     float costBetween = glm::dot(between, Q * between);

//     if (!notInvertible) {
//         glm::vec4 v = glm::inverse(Qhat) * zero;
//         collapse->position = glm::vec3(v);
//     } else {
//         if (costV1 < costV2 && costV1 < costBetween) {
//             collapse->position = glm::vec3(posV1);
//         } else if (costV2 < costV1 && costV2 < costBetween) {
//             collapse->position = glm::vec3(posV2);
//         } else {
//             collapse->position = glm::vec3(between);
//         }
//     }

//     glm::vec4 position(collapse->position, 1.0f);
//     float deltaV = glm::dot(position, Q * position);
//     collapse->cost = deltaV;
// }

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
    
    // If QuadricIsoSurfaces mode is active, render the iso-surfaces
    if (mVisualizationMode == QuadricIsoSurfaces) {
        // Set up OpenGL state for transparent ellipsoids
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(0.0f, 0.8f, 0.8f, 0.4f); // Cyan, semi-transparent
        glDisable(GL_LIGHTING);
        
        // Save current OpenGL state
        glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
        
        // Set polygon mode to wireframe for better visualization
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        RenderQuadricIsoSurfaces();
        
        // Restore OpenGL state
        glPopAttrib();
        glEnable(GL_LIGHTING);
        glDisable(GL_BLEND);
    }
}

// Render quadric error ellipsoids for all vertices
void QuadricDecimationMesh::RenderQuadricIsoSurfaces() {
    // Set up OpenGL state for rendering ellipsoids
    glEnable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    
    // Use wireframe mode for better visualization
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
    // Set line width for better visibility
    glLineWidth(1.0f);
    
    // Enable blending for transparency
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Set material properties
    float diffuse[] = {0.0f, 0.8f, 0.0f, 0.7f}; // Green color like in Garland's paper
    float ambient[] = {0.0f, 0.2f, 0.0f, 0.7f};
    float specular[] = {0.5f, 1.0f, 0.5f, 0.7f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0f);
    
    // Use a smaller stride to show more ellipsoids on important features
    // but not too many to avoid clutter
    size_t stride = std::max(size_t(1), mVerts.size() / 100);
    float visualEpsilon = mQuadricEpsilon * 0.1f; // Smaller epsilon for more precise ellipsoids
    int displayedCount = 0;
    
    // Display ellipsoids for selected vertices
    for (size_t i = 0; i < mVerts.size(); i += stride) {
        // Get the quadric for this vertex
        glm::mat4 Q = mQuadrics[i];
        
        // Get the vertex position
        glm::vec3 vertexPos = mVerts[i].pos;
        
        // Visualize the quadric error ellipsoid
        if (VisualizeQuadricEllipsoid(Q, vertexPos, visualEpsilon)) {
            displayedCount++;
        }
    }
    
    // Reset OpenGL state
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable(GL_BLEND);
    
    // Debug output
    std::cout << "Displayed " << displayedCount << " quadric ellipsoids" << std::endl;
}

// Visualize a single quadric error ellipsoid
// Following Garland's thesis section 4.1.2
bool QuadricDecimationMesh::VisualizeQuadricEllipsoid(const glm::mat4& Q, 
                                                     const glm::vec3& center, 
                                                     float epsilon) {
    // 1. Extract the upper 3x3 submatrix of Q (A) and the vector part (b)
    glm::mat3 A;
    glm::vec3 b;
    float c = Q[3][3];
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A[i][j] = Q[i][j];
        }
        b[i] = Q[i][3];
    }
    
    // 2. Check if A is invertible
    float detA = glm::determinant(A);
    if (std::abs(detA) < 1e-10) {
        return false; // Skip if not invertible
    }
    
    // 3. Compute the center of the ellipsoid: v0 = -A^(-1)b
    glm::vec3 v0 = -glm::inverse(A) * b;
    
    // 4. Compute the offset quadric: K = c - b^T * A^(-1) * b
    float K = c - glm::dot(b, v0);
    
    // 5. If K is zero or negative, the ellipsoid doesn't exist
    if (K <= 1e-10) { // Use a small threshold instead of exactly 0
        return false;
    }
    
    // 6. Perform eigendecomposition of A to get principal directions
    // For simplicity, we'll use Cholesky factorization as an approximation
    glm::mat4 Atemp(0.0f);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Atemp[i][j] = A[i][j];
        }
        Atemp[i][3] = 0.0f;
        Atemp[3][i] = 0.0f;
    }
    Atemp[3][3] = 1.0f;
    
    glm::mat4 Rtemp(0.0f);
    if (!CholeskyFactorization(Atemp, Rtemp)) {
        return false; // Skip if factorization fails
    }
    
    // Extract the 3x3 part
    glm::mat3 R3;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R3[i][j] = Rtemp[i][j];
        }
    }
    
    // 7. Compute the transformation matrix
    // Use K to scale properly according to Garland's method
    float scaleFactor = std::sqrt(epsilon / K);
    glm::mat3 T = glm::inverse(R3) * scaleFactor;
    
    // 8. Render the ellipsoid
    glPushMatrix();
    
    // Place ellipsoid on the surface
    // Use a small fraction of v0 to keep ellipsoids close to the surface
    // but still showing the error direction
    glm::vec3 ellipsoidCenter = center + v0 * 0.1f;
    glTranslatef(ellipsoidCenter.x, ellipsoidCenter.y, ellipsoidCenter.z);
    
    // Apply the transformation matrix
    float glMat[16] = {0};
    // Convert glm::mat3 to OpenGL matrix (column-major)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            glMat[i*4+j] = T[j][i]; // Note the transpose here
        }
    }
    glMat[15] = 1.0f; // Set the homogeneous coordinate
    
    glMultMatrixf(glMat);
    
    // Draw a unit sphere that will be transformed into an ellipsoid
    GLUquadric* quad = gluNewQuadric();
    gluQuadricDrawStyle(quad, GLU_LINE); // Use wireframe for better visualization
    gluQuadricNormals(quad, GLU_SMOOTH);
    gluSphere(quad, 1.0, 12, 12); // Fewer slices for cleaner look
    gluDeleteQuadric(quad);
    
    glPopMatrix();
    return true;
}
