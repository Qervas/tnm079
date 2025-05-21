#ifndef _strange_dubdivmesh_
#define _strange_dubdivmesh_

#include "AdaptiveLoopSubdivisionMesh.h"
#include <algorithm>
#include <cmath>

class StrangeSubdivisionMesh : public AdaptiveLoopSubdivisionMesh {
public:
    StrangeSubdivisionMesh() : AdaptiveLoopSubdivisionMesh() {
        // Initialize with a moderate threshold that will adapt over time
        mCurvatureThreshold = 0.5f;
        // Start with a moderate reduction factor
        mThresholdReductionFactor = 0.75f;
    }

    StrangeSubdivisionMesh(const HalfEdgeMesh& m, size_t s = 0) 
        : AdaptiveLoopSubdivisionMesh(m, s) {
        // Initialize with a moderate threshold that will adapt over time
        mCurvatureThreshold = 0.5f;
        // Start with a moderate reduction factor
        mThresholdReductionFactor = 0.75f;
    }

    virtual void Subdivide() {
        // Reduce the threshold with each subdivision step to focus on finer details
        mCurvatureThreshold *= mThresholdReductionFactor;
        
        // Compute face curvatures before subdivision
        ComputeFaceCurvatures();
        
        // Call the parent class subdivision method
        AdaptiveLoopSubdivisionMesh::Subdivide();
    }

protected:
    // Curvature-based adaptive subdivision
    bool Subdividable(size_t fi) {
        // If we haven't computed curvatures yet or the face doesn't exist, use default behavior
        if (mFaceCurvatures.empty() || fi >= mFaceCurvatures.size()) {
            return true;
        }
        
        // Subdivide faces with curvature above the current threshold
        return mFaceCurvatures[fi] > mCurvatureThreshold;
    }
    
    // Compute curvature for each face in the mesh
    void ComputeFaceCurvatures() {
        size_t numFaces = GetNumFaces();
        mFaceCurvatures.resize(numFaces);
        
        // Find the maximum curvature to normalize values
        float maxCurvature = 0.0f;
        
        // Compute the curvature for each face
        for (size_t i = 0; i < numFaces; i++) {
            mFaceCurvatures[i] = ComputeFaceCurvature(i);
            maxCurvature = std::max(maxCurvature, mFaceCurvatures[i]);
        }
        
        // Normalize curvature values if we have a non-zero maximum
        if (maxCurvature > 0.0f) {
            for (size_t i = 0; i < numFaces; i++) {
                mFaceCurvatures[i] /= maxCurvature;
            }
        }
    }
    
    // Compute curvature for a single face based on normal variation
    float ComputeFaceCurvature(size_t faceIndex) {
        // Get the face normal
        glm::vec3 faceNormal = f(faceIndex).normal;
        
        // Get the three vertices of the face
        EdgeIterator eit = GetEdgeIterator(f(faceIndex).edge);
        size_t v0 = eit.GetEdgeVertexIndex();
        size_t v1 = eit.Next().GetEdgeVertexIndex();
        size_t v2 = eit.Next().GetEdgeVertexIndex();
        
        // Find neighboring faces
        std::vector<size_t> neighborFaces;
        
        // Add neighbors from each edge
        eit = GetEdgeIterator(f(faceIndex).edge);
        for (int i = 0; i < 3; i++) {
            size_t pairFace = eit.Pair().GetEdgeFaceIndex();
            if (pairFace != faceIndex) {
                neighborFaces.push_back(pairFace);
            }
            eit.Pair().Next();
        }
        
        // Compute curvature as average deviation of normals
        float curvature = 0.0f;
        for (size_t neighborFace : neighborFaces) {
            glm::vec3 neighborNormal = f(neighborFace).normal;
            // Use 1 - dot product as a measure of normal deviation (0 for parallel, 2 for opposite)
            curvature += 1.0f - glm::dot(faceNormal, neighborNormal);
        }
        
        // Average the curvature if we have neighbors
        if (!neighborFaces.empty()) {
            curvature /= static_cast<float>(neighborFaces.size());
        }
        
        return curvature;
    }

private:
    // Threshold for determining which faces to subdivide based on curvature
    float mCurvatureThreshold;
    
    // Factor to reduce the threshold with each subdivision step
    float mThresholdReductionFactor;
    
    // Store curvature values for each face
    std::vector<float> mFaceCurvatures;
};

#endif
