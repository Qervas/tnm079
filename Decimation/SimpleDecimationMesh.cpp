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
#include "SimpleDecimationMesh.h"

void SimpleDecimationMesh::computeCollapse(EdgeCollapse* collapse) {
    // CUSTOM HEURISTIC: Curvature-based decimation metric
    // Motivation: Preserve high-curvature features (edges, corners, details)
    // while aggressively decimating flat regions
    // Applications: Architectural models, character models, CAD objects
    
    size_t v0_idx = mEdges[collapse->halfEdge].vert;
    size_t v1_idx = mEdges[mEdges[collapse->halfEdge].pair].vert;
    
    const glm::vec3& v0 = mVerts[v0_idx].pos;
    const glm::vec3& v1 = mVerts[v1_idx].pos;

    // Position new vertex at midpoint
    collapse->position = 0.5f * (v0 + v1);
    
    // Calculate curvature-based cost
    float curvature0 = calculateVertexCurvature(v0_idx);
    float curvature1 = calculateVertexCurvature(v1_idx);
    
    // Average curvature of the two vertices
    float avgCurvature = 0.5f * (curvature0 + curvature1);
    
    // Edge length component (basic geometric cost)
    float edgeLength = glm::distance(v0, v1);
    
    // Curvature weight: high curvature = high cost (preserve detail)
    // Low curvature = low cost (decimate flat areas first)
    float curvatureWeight = 2.0f;
    
    // Final cost: combine edge length with curvature penalty
    // High curvature areas get much higher cost (preserved longer)
    // Flat areas get lower cost (decimated first)
    collapse->cost = edgeLength * (1.0f + curvatureWeight * avgCurvature);
}

float SimpleDecimationMesh::calculateVertexCurvature(size_t vertexIndex) {
    // Calculate discrete mean curvature at vertex using normal variation
    // Method: Measure how much normals change around the vertex
    
    const Vertex& vertex = mVerts[vertexIndex];
    std::vector<size_t> neighborFaces = FindNeighborFaces(vertexIndex);
    
    if (neighborFaces.size() < 2) {
        return 0.0f; // Boundary vertex or isolated
    }
    
    // Calculate average normal of neighboring faces
    glm::vec3 avgNormal(0.0f);
    for (size_t faceIdx : neighborFaces) {
        avgNormal += mFaces[faceIdx].normal;
    }
    avgNormal = glm::normalize(avgNormal);
    
    // Measure curvature as variance in normal directions
    float curvature = 0.0f;
    for (size_t faceIdx : neighborFaces) {
        glm::vec3 faceNormal = mFaces[faceIdx].normal;
        // Angle between face normal and average normal
        float dot = glm::clamp(glm::dot(faceNormal, avgNormal), -1.0f, 1.0f);
        float angle = std::acos(dot);
        curvature += angle * angle; // Square for more emphasis
    }
    
    // Normalize by number of faces
    curvature /= neighborFaces.size();
    
    // Add edge-based curvature component
    std::vector<size_t> neighborVerts = FindNeighborVertices(vertexIndex);
    if (neighborVerts.size() > 2) {
        float edgeCurvature = calculateEdgeCurvature(vertexIndex, neighborVerts);
        curvature += 0.5f * edgeCurvature;
    }
    
    return curvature;
}

float SimpleDecimationMesh::calculateEdgeCurvature(size_t vertexIndex, 
                                                  const std::vector<size_t>& neighbors) {
    // Calculate curvature based on how much the vertex deviates 
    // from the plane of its neighbors
    
    if (neighbors.size() < 3) return 0.0f;
    
    const glm::vec3& centerPos = mVerts[vertexIndex].pos;
    
    // Calculate centroid of neighbor positions
    glm::vec3 centroid(0.0f);
    for (size_t neighborIdx : neighbors) {
        centroid += mVerts[neighborIdx].pos;
    }
    centroid /= float(neighbors.size());
    
    // Calculate average distance from center to neighbors
    float avgDistance = 0.0f;
    for (size_t neighborIdx : neighbors) {
        avgDistance += glm::distance(centerPos, mVerts[neighborIdx].pos);
    }
    avgDistance /= float(neighbors.size());
    
    // Curvature is how much the center deviates from the neighbor centroid
    // relative to the average edge length
    float deviation = glm::distance(centerPos, centroid);
    float relativeCurvature = avgDistance > 0.0f ? deviation / avgDistance : 0.0f;
    
    return relativeCurvature;
}
