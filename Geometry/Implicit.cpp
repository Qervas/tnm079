#include <Geometry/Implicit.h>
#include <gtc/type_ptr.hpp>
#include <cmath>

#ifdef __APPLE__
#include "GLUT/glut.h"
#else
#include "GL/glut.h"
#endif

const Implicit::VisualizationMode Implicit::Gradients = NewVisualizationMode("Gradients");
const Implicit::VisualizationMode Implicit::Curvature = NewVisualizationMode("Curvature");

Implicit::Implicit() : mMesh(NULL), mMeshSampling(0.1f), mDelta(0.1f) {}

Implicit::~Implicit() {
    if (mMesh != NULL) {
        delete mMesh;
        mMesh = NULL;
    }
}

void Implicit::Update() {
    if (mVisualizationMode == Curvature) {
        if (typeid(*mMesh) == typeid(SimpleMesh)) {
            SimpleMesh* ptr = static_cast<SimpleMesh*>(mMesh);
            std::vector<SimpleMesh::Vertex>& verts = ptr->GetVerts();

            glm::mat4 M = glm::transpose(GetTransform());

            // Compute curvature of implicit geometry and assign to the vertex
            // property
            for (size_t i = 0; i < verts.size(); i++) {
                const glm::vec3 vObject = verts.at(i).pos;

                // Transform vertex position to world space
                glm::vec4 vWorld =
                    GetTransform() * glm::vec4(vObject[0], vObject[1], vObject[2], 1);

                // Get curvature in world space
                verts.at(i).curvature = GetCurvature(vWorld[0], vWorld[1], vWorld[2]);

                // Get gradient in world space (used for lighting)
                glm::vec3 nWorld = GetGradient(vWorld[0], vWorld[1], vWorld[2]);

                // Transform gradient to object space
                glm::vec4 nObject = M * glm::vec4(nWorld[0], nWorld[1], nWorld[2], 0);
                verts.at(i).normal = glm::normalize(glm::vec3(nObject[0], nObject[1], nObject[2]));
            }

            ptr->mAutoMinMax = mAutoMinMax;
            ptr->mMinCMap = mMinCMap;
            ptr->mMaxCMap = mMaxCMap;
            ptr->SetVisualizationMode(Mesh::CurvatureVertex);
            ptr->Update();
        } else {
            std::cerr << "No Curvature visualization mode implemented for mesh type: "
                      << typeid(*mMesh).name() << std::endl;
        }
    }
}

void Implicit::Initialize() {
    Geometry* mesh = dynamic_cast<Geometry*>(mMesh);
    if (mesh == NULL)
        std::cerr << "Error: implicit geometry not triangulated, add call to "
                     "triangulate()"
                  << std::endl;
    else {
        std::cerr << "Computing normals etc... ";
        mesh->Initialize();
        std::cerr << " done" << std::endl;
    }
}

/*!
 * Evaluates gradient at (x,y,z) through discrete finite difference scheme.
 */
glm::vec3 Implicit::GetGradient(float x, float y, float z) const {
    // Implement finite difference evaluation of gradient at world coordinates
    // (x,y,z) Use mDelta variable as epsilon in eqn. 16 in lab text
    // Central difference for gradient calculation
    float epsilon = mDelta;

    float phi_x = (GetValue(x + epsilon, y, z) - GetValue(x - epsilon, y, z)) / (2 * epsilon);
    float phi_y = (GetValue(x, y + epsilon, z) - GetValue(x, y - epsilon, z)) / (2 * epsilon);
    float phi_z = (GetValue(x, y, z + epsilon) - GetValue(x, y, z - epsilon)) / (2 * epsilon);

    return glm::vec3(phi_x, phi_y, phi_z);
}

/*!
 * Evaluates curvature at (x,y,z) through discrete finite difference scheme.
 */
float Implicit::GetCurvature(float x, float y, float z) const {
    // Central difference for second-order partial derivatives
    float epsilon = mDelta;
    
    // First-order partial derivatives
    float phi_x = (GetValue(x + epsilon, y, z) - GetValue(x - epsilon, y, z)) / (2 * epsilon);
    float phi_y = (GetValue(x, y + epsilon, z) - GetValue(x, y - epsilon, z)) / (2 * epsilon);
    float phi_z = (GetValue(x, y, z + epsilon) - GetValue(x, y, z - epsilon)) / (2 * epsilon);

    // Second-order partial derivatives
    float phi_xx = (GetValue(x + epsilon, y, z) - 2 * GetValue(x, y, z) + GetValue(x - epsilon, y, z)) / (epsilon * epsilon);
    float phi_yy = (GetValue(x, y + epsilon, z) - 2 * GetValue(x, y, z) + GetValue(x, y - epsilon, z)) / (epsilon * epsilon);
    float phi_zz = (GetValue(x, y, z + epsilon) - 2 * GetValue(x, y, z) + GetValue(x, y, z - epsilon)) / (epsilon * epsilon);

    float phi_xy = (GetValue(x + epsilon, y + epsilon, z) - GetValue(x + epsilon, y - epsilon, z) -
                    GetValue(x - epsilon, y + epsilon, z) + GetValue(x - epsilon, y - epsilon, z)) / (4 * epsilon * epsilon);

    float phi_yz = (GetValue(x, y + epsilon, z + epsilon) - GetValue(x, y + epsilon, z - epsilon) -
                    GetValue(x, y - epsilon, z + epsilon) + GetValue(x, y - epsilon, z - epsilon)) / (4 * epsilon * epsilon);

    float phi_zx = (GetValue(x + epsilon, y, z + epsilon) - GetValue(x + epsilon, y, z - epsilon) -
                    GetValue(x - epsilon, y, z + epsilon) + GetValue(x - epsilon, y, z - epsilon)) / (4 * epsilon * epsilon);

    // Compute mean curvature using the derived formula
    float gradient_magnitude = glm::sqrt(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z);
    float curvature = ((phi_xx * (phi_y * phi_y + phi_z * phi_z) - 2 * phi_x * phi_y * phi_xy - 2 * phi_x * phi_z * phi_zx) +
                       (phi_yy * (phi_x * phi_x + phi_z * phi_z) - 2 * phi_y * phi_z * phi_yz) +
                       (phi_zz * (phi_x * phi_x + phi_y * phi_y))) / glm::pow(gradient_magnitude, 3);

    return curvature;
}

float Implicit::ComputeArea(float dx) const { return 0; }

float Implicit::ComputeVolume(float dx) const {
    Bbox box = GetBoundingBox();
    float volume = 0;

    float H;
    for (float x = box.pMin[0]; x <= box.pMax[0] + 0.5f * dx; x += dx) {
        for (float y = box.pMin[1]; y <= box.pMax[1] + 0.5f * dx; y += dx) {
            for (float z = box.pMin[2]; z <= box.pMax[2] + 0.5f * dx; z += dx) {
                float val = GetValue(x, y, z);
                if (val < -dx) {
                    H = 1;
                } else if (val > dx) {
                    H = 0;
                } else {
                    H = 0.5f *
                        (1.f + val / dx +
                         sin(-static_cast<float>(M_PI) * val / dx) / static_cast<float>(M_PI));
                }

                volume += H;
            }
        }
    }

    return volume * static_cast<float>(std::pow(dx, 3.0));
}

Bbox Implicit::GetBoundingBox() const {
    // transform returns a copy
    return mBox.Transform(GetTransform());
}

void Implicit::SetBoundingBox(const Bbox& b) { mBox = b.Transform(mWorld2Obj); }

void Implicit::SetTransform(const glm::mat4& transform) {
    Geometry::SetTransform(transform);
    mWorld2Obj = glm::inverse(GetTransform());
}

void Implicit::TransformW2O(float& x, float& y, float& z) const {
    glm::vec4 vprim, v = glm::vec4(x, y, z, 1.f);
    vprim = mWorld2Obj * v;
    x = vprim[0];
    y = vprim[1];
    z = vprim[2];
}

void Implicit::Render() {
    // Draw bounding box for debugging
    Bbox b = GetBoundingBox();

    glm::vec3& v0 = b.pMin;
    glm::vec3& v1 = b.pMax;

    if (mSelected) {
        glLineWidth(2.f);
        glColor3f(0.8f, 0.8f, 0.8f);
    } else {
        glColor3f(0.1f, 0.1f, 0.1f);
    }

    glBegin(GL_LINE_STRIP);
    glVertex3f(v0[0], v0[1], v0[2]);
    glVertex3f(v1[0], v0[1], v0[2]);
    glVertex3f(v1[0], v1[1], v0[2]);
    glVertex3f(v0[0], v1[1], v0[2]);
    glVertex3f(v0[0], v0[1], v0[2]);
    glEnd();

    glBegin(GL_LINE_STRIP);
    glVertex3f(v0[0], v0[1], v1[2]);
    glVertex3f(v1[0], v0[1], v1[2]);
    glVertex3f(v1[0], v1[1], v1[2]);
    glVertex3f(v0[0], v1[1], v1[2]);
    glVertex3f(v0[0], v0[1], v1[2]);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(v0[0], v0[1], v0[2]);
    glVertex3f(v0[0], v0[1], v1[2]);

    glVertex3f(v1[0], v0[1], v0[2]);
    glVertex3f(v1[0], v0[1], v1[2]);

    glVertex3f(v0[0], v1[1], v0[2]);
    glVertex3f(v0[0], v1[1], v1[2]);

    glVertex3f(v1[0], v1[1], v0[2]);
    glVertex3f(v1[0], v1[1], v1[2]);
    glEnd();
    glLineWidth(1.f);
    glPushMatrix();
    glMultMatrixf(glm::value_ptr(mTransform));

    Geometry* mesh = dynamic_cast<Geometry*>(mMesh);
    if (mesh == NULL) {
        std::cerr << "Error: implicit geometry not triangulated, add call to "
                     "triangulate()"
                  << std::endl;
    } else {
        mesh->SetShowNormals(mShowNormals);
        mesh->SetWireframe(mWireframe);
        mesh->SetOpacity(mOpacity);

        mesh->Render();
    }

    if (mVisualizationMode == Gradients) {
        if (typeid(*mMesh) == typeid(SimpleMesh)) {
            SimpleMesh* ptr = static_cast<SimpleMesh*>(mMesh);
            const std::vector<SimpleMesh::Vertex>& verts = ptr->GetVerts();

            glDisable(GL_LIGHTING);

            glm::mat4 M = glm::transpose(GetTransform());

            glColor3f(0, 0, 1);
            glBegin(GL_LINES);
            for (size_t i = 0; i < verts.size(); i++) {
                const glm::vec3 vObject = verts.at(i).pos;

                // Transform vertex position to world space
                glm::vec4 vWorld =
                    GetTransform() * glm::vec4(vObject[0], vObject[1], vObject[2], 1);

                // Get gradient in world space
                glm::vec3 nWorld = GetGradient(vWorld[0], vWorld[1], vWorld[2]);

                // Transform gradient to object space
                glm::vec4 nObject = M * glm::vec4(nWorld[0], nWorld[1], nWorld[2], 0);
                glm::vec3 n = glm::vec3(nObject[0], nObject[1], nObject[2]);

                glVertex3fv(glm::value_ptr(vObject));
                glVertex3fv(glm::value_ptr(vObject + n * 0.1f));
            }
            glEnd();
        } else {
            std::cerr << "No Gradient visualization mode implemented for mesh type: "
                      << typeid(*mMesh).name() << std::endl;
        }
    }

    glPopMatrix();

    GLObject::Render();
}

void Implicit::SetVisualizationMode(const VisualizationMode& source) {
    Geometry::SetVisualizationMode(source);
    Update();
}

void Implicit::SetColorMap(ColorMap* colormap) {
    Geometry* mesh = dynamic_cast<Geometry*>(mMesh);
    if (mesh != NULL) mesh->SetColorMap(colormap);

    Geometry::SetColorMap(colormap);
    Update();
}
