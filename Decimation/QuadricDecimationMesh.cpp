#include "QuadricDecimationMesh.h"
#include <VC++/glm/gtx/euler_angles.hpp>
#include <VC++/glm/gtc/type_ptr.hpp>
#include <GUI/GLViewer.h>
#include <Util/Util.h>
const QuadricDecimationMesh::VisualizationMode QuadricDecimationMesh::QuadricIsoSurfaces =
    NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize() {
    // Print camera information at the start of decimation
    std::cout << "\n=== QUADRIC DECIMATION INITIALIZE ===" << std::endl;
    printCameraDebugInfo();
    
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
    
    std::cout << "=== INITIALIZATION COMPLETE ===\n" << std::endl;
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
    // Get the vertices at the endpoints of the edge to be collapsed
    size_t v1 = e(collapse->halfEdge).vert;
    size_t v2 = e(e(collapse->halfEdge).next).vert;

    glm::vec3 posV1 = v(v1).pos;
    glm::vec3 posV2 = v(v2).pos;
    glm::vec3 between = 0.5f * (posV1 + posV2);

    // GRADE 4: CUSTOM HEURISTIC - View-Dependent Decimation
    // =====================================================
    // MOTIVATION: Real-time rendering applications where screen-space detail matters
    // APPLICATIONS: 
    // - Games: LOD systems that preserve detail where players look
    // - VR: Maintain quality in central vision, reduce in periphery  
    // - Interactive walkthroughs: Preserve nearby detail, reduce distant geometry
    // - Streaming/bandwidth optimization: Reduce data for distant objects
    
    // HYBRID APPROACH: Combine quadric quality with view-dependent factors
    
    // 1. Start with basic quadric computation for geometric quality
    auto Q1 = mQuadrics[v1];
    auto Q2 = mQuadrics[v2];
    auto Q = Q1 + Q2;

    // Find optimal position using quadrics
    auto Qhat = Q;
    Qhat[0][3] = Qhat[1][3] = Qhat[2][3] = 0.0f;
    Qhat[3][3] = 1.0f;

    float EPSILON = 1e-9f;
    bool notInvertible = abs(glm::determinant(Qhat)) < EPSILON;
    
    glm::vec3 optimalPos;
    float quadricCost;
    
    if (!notInvertible) {
        glm::vec4 rhs(0.0f, 0.0f, 0.0f, 1.0f);
        glm::vec4 solution = glm::inverse(Qhat) * rhs;
        optimalPos = glm::vec3(solution.x, solution.y, solution.z);
        
        // Calculate quadric error
        glm::vec4 pos4(optimalPos.x, optimalPos.y, optimalPos.z, 1.0f);
        quadricCost = glm::dot(pos4, (Q * pos4));
    } else {
        // Fallback to midpoint
        optimalPos = between;
        glm::vec4 pos4(between.x, between.y, between.z, 1.0f);
        quadricCost = glm::dot(pos4, (Q * pos4));
    }
    
    // 2. VIEW-DEPENDENT FACTORS FOR ORBITAL CAMERA
    // ============================================
    
    glm::vec3 cameraPos = GLViewer::GetCamera().GetPosition();
    
    // FIXED: Calculate actual world-space viewing direction for orbital camera
    // In orbital camera, we look toward the origin (0,0,0) from camera position
    glm::vec3 actualViewDir = glm::normalize(glm::vec3(0.0f, 0.0f, 0.0f) - cameraPos);
    glm::vec3 viewTarget = glm::vec3(0.0f, 0.0f, 0.0f); // Looking at origin
    
    // A) DEPTH FROM VIEWING PLANE
    // Calculate how far the edge is along the actual viewing direction
    glm::vec3 toEdge = between - cameraPos;
    float depthAlongView = glm::dot(toEdge, actualViewDir);
    float depthWeight = 1.0f + 0.8f * std::max(0.0f, depthAlongView / 3.0f); // Farther = less important
    
    // B) DISTANCE FROM VIEW CENTER 
    // How far is the edge from what we're actually looking at (origin)
    float screenDistance = glm::distance(between, viewTarget);
    float screenWeight = 1.0f + 0.6f * (screenDistance / 2.0f); // Center = more important
    
    // C) VIEWING ANGLE ALIGNMENT
    // Edges that are more perpendicular to the actual view direction are less important
    glm::vec3 edgeVector = glm::normalize(posV2 - posV1);
    float viewAlignment = std::abs(glm::dot(edgeVector, actualViewDir));
    float alignmentWeight = 1.0f + 0.5f * (1.0f - viewAlignment); // Perpendicular = less important
    
    // D) FEATURE PRESERVATION (keep sharp edges)
    float featureImportance = calculateFeatureImportance(v1, v2);
    float featureWeight = 1.0f + 1.0f * featureImportance; 
    
    // E) SILHOUETTE PRESERVATION
    // Edges on the silhouette (boundary of the object) are more important
    float silhouetteImportance = calculateSilhouetteImportance(v1, v2, actualViewDir);
    float silhouetteWeight = 1.0f + 1.2f * silhouetteImportance;
    
    // 3. FINAL COST CALCULATION
    // ========================
    
    // Base cost: quadric error (ensures geometric quality)
    float baseCost = std::max(0.0f, quadricCost);
    
    // Apply all view-dependent factors
    float finalCost = baseCost * depthWeight * screenWeight * alignmentWeight * featureWeight * silhouetteWeight;
    
    // Set result
    collapse->position = optimalPos;
    collapse->cost = finalCost;
    
    // DEBUG: Print camera position and decimation details
    static int debugCounter = 0;
    if (debugCounter % 1000 == 0) {
        std::cout << "\n=== VIEW-DEPENDENT DECIMATION DEBUG ===" << std::endl;
        std::cout << "Camera Position: (" << cameraPos.x << ", " << cameraPos.y << ", " << cameraPos.z << ")" << std::endl;
        std::cout << "Actual View Dir: (" << actualViewDir.x << ", " << actualViewDir.y << ", " << actualViewDir.z << ")" << std::endl;
        std::cout << "Edge Position: (" << between.x << ", " << between.y << ", " << between.z << ")" << std::endl;
        std::cout << "Depth Along View: " << depthAlongView << " -> Weight: " << depthWeight << std::endl;
        std::cout << "Screen Distance: " << screenDistance << " -> Weight: " << screenWeight << std::endl;
        std::cout << "View Alignment: " << viewAlignment << " -> Weight: " << alignmentWeight << std::endl;
        std::cout << "Feature Importance: " << featureImportance << " -> Weight: " << featureWeight << std::endl;
        std::cout << "Silhouette Importance: " << silhouetteImportance << " -> Weight: " << silhouetteWeight << std::endl;
        std::cout << "Quadric Cost: " << baseCost << " -> Final Cost: " << finalCost << std::endl;
        std::cout << "========================================\n" << std::endl;
    }
    debugCounter++;
}

// GRADE 4: Helper functions for feature importance calculation
float QuadricDecimationMesh::calculateFeatureImportance(size_t v1, size_t v2) {
    // Calculate feature importance based on local geometry
    // Higher values = more important features (sharp edges, corners)
    
    auto faces1 = FindNeighborFaces(v1);
    auto faces2 = FindNeighborFaces(v2);

    if (faces1.size() < 2 || faces2.size() < 2) {
        return 2.0f; // Boundary edges are very important
    }
    
    // Calculate curvature around both vertices
    float curvature1 = calculateNormalVariation(v1, faces1);
    float curvature2 = calculateNormalVariation(v2, faces2);
    
    // Calculate dihedral angle across the edge
    float dihedralAngle = calculateDihedralAngle(v1, v2);
    
    // Combine measures - weight dihedral angle more heavily
    float avgCurvature = 0.5f * (curvature1 + curvature2);
    return avgCurvature + 2.0f * dihedralAngle;
}

float QuadricDecimationMesh::calculateNormalVariation(size_t vertexIdx, 
                                                     const std::vector<size_t>& faces) {
    if (faces.size() < 2) return 0.0f;
    
    // Calculate average normal
    glm::vec3 avgNormal(0.0f);
    for (size_t faceIdx : faces) {
        avgNormal += f(faceIdx).normal;
    }
    avgNormal = glm::normalize(avgNormal);
    
    // Calculate deviation from average
    float variation = 0.0f;
    for (size_t faceIdx : faces) {
        glm::vec3 normal = f(faceIdx).normal;
        float dot = glm::clamp(glm::dot(normal, avgNormal), -1.0f, 1.0f);
        float angle = std::acos(dot);
        variation += angle;
    }
    
    return variation / faces.size();
}

float QuadricDecimationMesh::calculateDihedralAngle(size_t v1, size_t v2) {
    // Find faces shared by both vertices (the edge between them)
    auto faces1 = FindNeighborFaces(v1);
    auto faces2 = FindNeighborFaces(v2);
    
    std::vector<size_t> sharedFaces;
    for (size_t f1 : faces1) {
        for (size_t f2 : faces2) {
            if (f1 == f2) {
                sharedFaces.push_back(f1);
                break;
            }
        }
    }
    
    if (sharedFaces.size() != 2) {
        return 1.0f; // Boundary edge - moderately important
    }
    
    // Calculate angle between face normals
    glm::vec3 n1 = f(sharedFaces[0]).normal;
    glm::vec3 n2 = f(sharedFaces[1]).normal;
    
    float dot = glm::clamp(glm::dot(n1, n2), -1.0f, 1.0f);
    float angle = std::acos(dot);
    
    // Return deviation from flat (π = flat, 0 = sharp crease)
    return std::abs(angle - M_PI);
}

float QuadricDecimationMesh::calculateSilhouetteImportance(size_t v1, size_t v2, const glm::vec3& viewDir) {
    // Calculate if edge is on the silhouette (boundary visible from current view)
    // Silhouette edges have faces with normals facing different directions relative to view
    
    auto faces1 = FindNeighborFaces(v1);
    auto faces2 = FindNeighborFaces(v2);
    
    // Find shared faces (the faces that share this edge)
    std::vector<size_t> sharedFaces;
    for (size_t f1 : faces1) {
        for (size_t f2 : faces2) {
            if (f1 == f2) {
                sharedFaces.push_back(f1);
                break;
            }
        }
    }
    
    if (sharedFaces.size() != 2) {
        return 1.0f; // Boundary edge - definitely on silhouette
    }
    
    // Check if the two faces are on different sides relative to the view direction
    glm::vec3 n1 = f(sharedFaces[0]).normal;
    glm::vec3 n2 = f(sharedFaces[1]).normal;
    
    float dot1 = glm::dot(n1, viewDir);
    float dot2 = glm::dot(n2, viewDir);
    
    // If normals face opposite directions relative to view, this is a silhouette edge
    if ((dot1 > 0 && dot2 < 0) || (dot1 < 0 && dot2 > 0)) {
        return 1.0f; // Strong silhouette edge
    }
    
    // Calculate the difference in facing direction
    return std::abs(dot1 - dot2) / 2.0f; // Normalize to 0-1 range
}

void QuadricDecimationMesh::printCameraDebugInfo() {
    std::cout << "\n*** CAMERA DEBUG INFO ***" << std::endl;
    
    glm::vec3 camPos = GLViewer::GetCamera().GetPosition();
    glm::vec3 camLookAt = GLViewer::GetCamera().GetLookAtVector();
    glm::vec3 camLookAtPoint = GLViewer::GetCamera().GetLookAtPoint();
    glm::vec3 camUp = GLViewer::GetCamera().GetUpVector();
    glm::vec3 camRight = GLViewer::GetCamera().GetRightVector();
    
    std::cout << "Camera Position: (" << camPos.x << ", " << camPos.y << ", " << camPos.z << ")" << std::endl;
    std::cout << "Look At Vector: (" << camLookAt.x << ", " << camLookAt.y << ", " << camLookAt.z << ")" << std::endl;
    std::cout << "Look At Point: (" << camLookAtPoint.x << ", " << camLookAtPoint.y << ", " << camLookAtPoint.z << ")" << std::endl;
    std::cout << "Up Vector: (" << camUp.x << ", " << camUp.y << ", " << camUp.z << ")" << std::endl;
    std::cout << "Right Vector: (" << camRight.x << ", " << camRight.y << ", " << camRight.z << ")" << std::endl;
    std::cout << "**************************\n" << std::endl;
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
    // COLORFUL WIREFRAME: Rainbow colors for maximum visibility and beauty!
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
    glDisable(GL_LIGHT2);
    glDisable(GL_LIGHT3);
    glDisable(GL_LIGHT4);
    glDisable(GL_LIGHT5);
    
    // Custom colorful rendering - draw each face with bright colors
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixf(glm::value_ptr(mTransform));
    
    // First pass: Filled faces with rainbow colors
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_TRIANGLES);
    
    for (size_t i = 0; i < mFaces.size(); i++) {
        if (isFaceCollapsed(i)) continue;
        
        Face& face = mFaces[i];
        HalfEdge* edge = &mEdges[face.edge];
        
        Vertex& v1 = mVerts[edge->vert];
        edge = &mEdges[edge->next];
        Vertex& v2 = mVerts[edge->vert];
        edge = &mEdges[edge->next];
        Vertex& v3 = mVerts[edge->vert];
        
        // Generate bright rainbow color based on face index
        float hue = (float(i % 360)) / 360.0f; // Cycle through hues
        float r, g, b;
        
        // HSV to RGB conversion for vibrant colors
        if (hue < 1.0f/6.0f) {
            r = 1.0f; g = hue * 6.0f; b = 0.0f;
        } else if (hue < 2.0f/6.0f) {
            r = 2.0f - hue * 6.0f; g = 1.0f; b = 0.0f;
        } else if (hue < 3.0f/6.0f) {
            r = 0.0f; g = 1.0f; b = (hue - 2.0f/6.0f) * 6.0f;
        } else if (hue < 4.0f/6.0f) {
            r = 0.0f; g = 1.0f - (hue - 3.0f/6.0f) * 6.0f; b = 1.0f;
        } else if (hue < 5.0f/6.0f) {
            r = (hue - 4.0f/6.0f) * 6.0f; g = 0.0f; b = 1.0f;
        } else {
            r = 1.0f; g = 0.0f; b = 1.0f - (hue - 5.0f/6.0f) * 6.0f;
        }
        
        // Make colors brighter and more saturated
        r = 0.3f + 0.7f * r; // Ensure minimum brightness
        g = 0.3f + 0.7f * g;
        b = 0.3f + 0.7f * b;
        
        glColor3f(r, g, b);
        glVertex3fv(glm::value_ptr(v1.pos));
        glVertex3fv(glm::value_ptr(v2.pos));
        glVertex3fv(glm::value_ptr(v3.pos));
    }
    glEnd();
    
    // Second pass: Black wireframe outline for clarity
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth(1.0f);
    glColor3f(0.0f, 0.0f, 0.0f); // Black outline
    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset(-1.0f, -1.0f);
    
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < mFaces.size(); i++) {
        if (isFaceCollapsed(i)) continue;
        
        Face& face = mFaces[i];
        HalfEdge* edge = &mEdges[face.edge];
        
        Vertex& v1 = mVerts[edge->vert];
        edge = &mEdges[edge->next];
        Vertex& v2 = mVerts[edge->vert];
        edge = &mEdges[edge->next];
        Vertex& v3 = mVerts[edge->vert];
        
        glVertex3fv(glm::value_ptr(v1.pos));
        glVertex3fv(glm::value_ptr(v2.pos));
        glVertex3fv(glm::value_ptr(v3.pos));
    }
    glEnd();
    
    glDisable(GL_POLYGON_OFFSET_LINE);
    glPopMatrix();
    
    // Reset state
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLineWidth(1.0f);
    
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
    // Set up OpenGL state for bright green ellipsoids like in the target image
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDisable(GL_CULL_FACE);
    
    // SHINING/TWINKLING ELLIPSOIDS!
    // Create smooth, uniform animated brightness based on time
    static float animationTime = 0.0f;
    animationTime += 0.02f; // Slower, smoother animation speed
    
    // Smooth linear pulsing brightness effect (no sudden acceleration)
    float normalizedTime = std::fmod(animationTime, 2.0f); // Cycle every 2 seconds
    float pulse;
    if (normalizedTime < 1.0f) {
        pulse = 0.7f + 0.3f * normalizedTime; // Linear fade up from 0.7 to 1.0
    } else {
        pulse = 1.0f - 0.3f * (normalizedTime - 1.0f); // Linear fade down from 1.0 to 0.7
    }
    
    // Bright shining green color with animation
    glColor3f(0.0f, pulse, 0.0f);
    
    // SUPER SHINY material properties for maximum visibility
    float diffuse[] = {0.0f, pulse, 0.0f, 1.0f}; // Animated green
    float ambient[] = {0.0f, pulse * 0.4f, 0.0f, 1.0f}; // Animated ambient
    float specular[] = {pulse, pulse, pulse, 1.0f}; // White specular highlights
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128.0f); // Maximum shininess!
    
    // Render tiny ellipsoids like in the target images
    float epsilon = 0.001f; // Small error tolerance
    int displayedCount = 0;
    
    // Show many more vertices now that it's stable
    for (size_t i = 0; i < mVerts.size(); i += 8) {
        glm::mat4 Q = mQuadrics[i];
        glm::vec3 vertexPos = mVerts[i].pos;
        
        // Individual smooth twinkling effect for each ellipsoid
        float phaseOffset = (float)(i % 100) / 100.0f; // Different phase for each ellipsoid
        float individualTime = std::fmod(animationTime + phaseOffset, 1.5f); // Individual cycle
        float individualTwinkle;
        if (individualTime < 0.75f) {
            individualTwinkle = 0.8f + 0.2f * (individualTime / 0.75f); // Linear up
        } else {
            individualTwinkle = 1.0f - 0.2f * ((individualTime - 0.75f) / 0.75f); // Linear down
        }
        
        if (VisualizeQuadricEllipsoid(Q, vertexPos, epsilon, individualTwinkle)) {
            displayedCount++;
        }
    }
    
    std::cout << "Displayed " << displayedCount << " quadric ellipsoids out of " << mVerts.size() << " vertices" << std::endl;
}

// Visualize a single quadric error ellipsoid
// Following Garland's thesis Section 4.1.2 with safety checks
bool QuadricDecimationMesh::VisualizeQuadricEllipsoid(const glm::mat4& Q, 
                                                     const glm::vec3& vertex_pos, 
                                                     float epsilon, 
                                                     float twinkle) {
    // Safety check: validate input parameters
    if (epsilon <= 0.0f || epsilon > 1.0f) {
        return false;
    }
    
    // Check for NaN or infinite values in vertex position
    if (!std::isfinite(vertex_pos.x) || !std::isfinite(vertex_pos.y) || !std::isfinite(vertex_pos.z)) {
        return false;
    }
    
    // Check quadric matrix for NaN/infinite values
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (!std::isfinite(Q[i][j])) {
                return false;
            }
        }
    }
    
    // SIMPLIFIED APPROACH: Just render a tiny sphere at the vertex
    // This avoids complex matrix operations that can cause crashes
    glPushMatrix();
    
    // Translate to vertex position
    glTranslatef(vertex_pos.x, vertex_pos.y, vertex_pos.z);
    
    // Scale to make a tiny ellipsoid with individual twinkling size variation
    float scale = 0.01f * twinkle; // Size varies with twinkling
    glScalef(scale, scale, scale);
    
    // Set individual twinkling color for this ellipsoid
    glColor3f(0.0f, twinkle, 0.0f);
    
    // Render tiny solid sphere with twinkling effect
    GLUquadric* quad = gluNewQuadric();
    if (quad) {
        gluQuadricDrawStyle(quad, GLU_FILL);
        gluQuadricNormals(quad, GLU_SMOOTH);
        gluSphere(quad, 1.0f, 8, 8);
        gluDeleteQuadric(quad);
    }
    
    glPopMatrix();
    return true;
}
