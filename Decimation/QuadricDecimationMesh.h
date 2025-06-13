#pragma once

#include "Decimation/DecimationMesh.h"
#include <iomanip>

#ifdef __APPLE__
#include "GLUT/glut.h"
#else
#include "GL/glut.h"
#endif

class QuadricDecimationMesh : public virtual DecimationMesh {
public:
    static const VisualizationMode QuadricIsoSurfaces;

    virtual std::list<VisualizationMode> GetVisualizationModes() {
        std::list<VisualizationMode> L = DecimationMesh::GetVisualizationModes();
        L.push_back(QuadricIsoSurfaces);
        return L;
    }

    QuadricDecimationMesh() {}
    virtual ~QuadricDecimationMesh() {}

    //! Initialize member data (error quadrics)
    virtual void Initialize();

protected:
    //! Compute the cost and new position for an edge collapse
    virtual void computeCollapse(EdgeCollapse* collapse);
    //! Update vertex properties. Used after an edge collapse
    virtual void updateVertexProperties(size_t ind);
    //! Compute the quadric for a vertex
    glm::mat4 createQuadricForVert(size_t indx) const;
    //! Copmute the quadric for a face
    glm::mat4 createQuadricForFace(size_t indx) const;
    
    //! Grade 4: Custom heuristic helper functions for view-dependent decimation
    float calculateFeatureImportance(size_t v1, size_t v2);
    float calculateNormalVariation(size_t vertexIdx, const std::vector<size_t>& faces);
    float calculateDihedralAngle(size_t v1, size_t v2);
    float calculateSilhouetteImportance(size_t v1, size_t v2, const glm::vec3& viewDir);
    
    //! Debug function to print camera information
    void printCameraDebugInfo();
    
    //! Render (redefined)
    virtual void Render();
    
    //! Render quadric error ellipsoids
    void RenderQuadricIsoSurfaces();
    
    //! Visualize a single quadric error ellipsoid
    //! Returns true if the ellipsoid was successfully visualized
    bool VisualizeQuadricEllipsoid(const glm::mat4& Q, const glm::vec3& center, float epsilon, float twinkle = 1.0f);

    //! The quadrics used in the decimation
    std::vector<glm::mat4> mQuadrics;
    
    //! Error tolerance for quadric visualization
    float mQuadricEpsilon = 0.01f;
};
