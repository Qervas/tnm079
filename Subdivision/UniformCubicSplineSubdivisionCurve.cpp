#include <Subdivision/UniformCubicSplineSubdivisionCurve.h>
#include <glm.hpp>
#include <gtc/type_ptr.hpp>

UniformCubicSplineSubdivisionCurve::UniformCubicSplineSubdivisionCurve(
    const std::vector<glm::vec3>& joints, glm::vec3 lineColor, float lineWidth)
    : mCoefficients(joints), mControlPolygon(joints) {
    this->mLineColor = lineColor;
    this->mLineWidth = lineWidth;
}

void UniformCubicSplineSubdivisionCurve::Subdivide() {
    std::vector<glm::vec3> newc;

    assert(mCoefficients.size() > 4 && "Need at least 5 points to subdivide");

    // Handle the first boundary point
    newc.push_back(mCoefficients.front());

    // Apply the subdivision rules to the internal coefficients
    for (int i = 1; i < mCoefficients.size() - 1; ++i) {
        glm::vec3 newCoefficient1 = 0.125f * (mCoefficients[i - 1] + 6.0f * mCoefficients[i] + mCoefficients[i + 1]);
        glm::vec3 newCoefficient2 = 0.125f * (4.0f * mCoefficients[i] + 4.0f * mCoefficients[i + 1]);

        newc.push_back(newCoefficient1);
        newc.push_back(newCoefficient2);
    }

    // Handle the last boundary point
    newc.push_back(0.125f * (mCoefficients[mCoefficients.size() - 2] + 6.0f * mCoefficients.back() + mCoefficients[mCoefficients.size() - 1]));
    newc.push_back(mCoefficients.back());

    // Verify that the new number of coefficients is correct
    assert(mCoefficients.size() * 2 == newc.size() + 1 && "Incorrect number of new coefficients!");

    mCoefficients = newc;
}


void UniformCubicSplineSubdivisionCurve::Render() {
    // Apply transform
    glPushMatrix();  // Push modelview matrix onto stack

    // Convert transform-matrix to format matching GL matrix format
    // Load transform into modelview matrix
    glMultMatrixf(glm::value_ptr(mTransform));

    mControlPolygon.Render();

    // save line point and color states
    glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_CURRENT_BIT);

    // draw segments
    glLineWidth(mLineWidth);
    glColor3fv(glm::value_ptr(mLineColor));
    glBegin(GL_LINE_STRIP);
    // just draw the spline as a series of connected linear segments
    for (size_t i = 0; i < mCoefficients.size(); i++) {
        glVertex3fv(glm::value_ptr(mCoefficients.at(i)));
    }
    glEnd();

    // restore attribs
    glPopAttrib();

    glPopMatrix();

    GLObject::Render();
}
