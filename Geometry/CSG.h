#pragma once

#include <Geometry/Implicit.h>
#include <cmath>

/*! \brief CSG Operator base class */
class CSG_Operator : public Implicit {
protected:
    //! Constructor
    CSG_Operator(Implicit* l, Implicit* r) : left(l), right(r) {}

    //! Pointers to left and right child nodes
    Implicit *left, *right;
};

/*! \brief Union boolean operation */
class Union : public CSG_Operator {
public:
    Union(Implicit* l, Implicit* r) : CSG_Operator(l, r) {
        // Compute the resulting (axis aligned) bounding box from
        // the left and right children
        mBox = Bbox::BoxUnion(l->GetBoundingBox(), r->GetBoundingBox());
    }

    virtual float GetValue(float x, float y, float z) const {
        // The coordinates (x,y,z) are passed in from world space,
        // remember to transform them into object space
        // (Hint: Implicit::TransformW2O()). This
        // is needed because the CSG operators are also implicit geometry
        // and can be transformed like all implicit surfaces.
        // Then, get values from left and right children and perform the
        // boolean operation.
        TransformW2O(x, y, z);
        float leftValue = left->GetValue(x, y, z);
        float rightValue = right->GetValue(x, y, z);
        return std::min(leftValue, rightValue);
    }
};

/*! \brief Intersection boolean operation */
class Intersection : public CSG_Operator {
public:
    Intersection(Implicit* l, Implicit* r) : CSG_Operator(l, r) {
        mBox = Bbox::BoxIntersection(l->GetBoundingBox(), r->GetBoundingBox());
    }

    virtual float GetValue(float x, float y, float z) const override{
        TransformW2O(x, y, z);
        float leftValue = left->GetValue(x, y, z);
        float rightValue = right->GetValue(x, y, z);
        return std::max(leftValue, rightValue);
    }
};

/*! \brief Difference boolean operation */
class Difference : public CSG_Operator {
public:
    Difference(Implicit* l, Implicit* r) : CSG_Operator(l, r) { mBox = l->GetBoundingBox(); }

    virtual float GetValue(float x, float y, float z) const override{
        TransformW2O(x, y, z);
        float leftValue = left->GetValue(x, y, z);
        float rightValue = right->GetValue(x, y, z);
        return std::max(leftValue, -rightValue);
    }
};

/*! \brief BlendedUnion boolean operation */
class BlendedUnion : public CSG_Operator {
public:
    BlendedUnion(Implicit* l, Implicit* r, int blend) : CSG_Operator(l, r), mBlend(blend) {
        mBox = Bbox::BoxUnion(l->GetBoundingBox(), r->GetBoundingBox());
    }

    virtual float GetValue(float x, float y, float z) const override {
        TransformW2O(x, y, z);
        float leftValue = left->GetValue(x, y, z);
        float rightValue = right->GetValue(x, y, z);
        float k = static_cast<float>(mBlend);
        if (k <= 0.0f) k = 0.0001f;
        return -std::log(std::exp(-k * leftValue) + std::exp(-k * rightValue)) / k;
    }

protected:
    int mBlend;
};

/*! \brief BlendedIntersection boolean operation */
class BlendedIntersection : public CSG_Operator {
public:
    BlendedIntersection(Implicit* l, Implicit* r, int blend) : CSG_Operator(l, r), mBlend(blend) {
        mBox = Bbox::BoxIntersection(l->GetBoundingBox(), r->GetBoundingBox());
    }

    virtual float GetValue(float x, float y, float z) const override {
        TransformW2O(x, y, z);
        float leftValue = left->GetValue(x, y, z);
        float rightValue = right->GetValue(x, y, z);
        float k = static_cast<float>(mBlend);
        if (k <= 0.0f) k = 0.0001f;
        
        // Super-elliptic blending for intersection
        // Using the formula: log(exp(k*a) + exp(k*b))/k
        return std::log(std::exp(k * leftValue) + std::exp(k * rightValue)) / k;
    }

protected:
    int mBlend;
};

/*! \brief BlendedDifference boolean operation */
class BlendedDifference : public CSG_Operator {
public:
    BlendedDifference(Implicit* l, Implicit* r, int blend) : CSG_Operator(l, r), mBlend(blend) {
        mBox = l->GetBoundingBox();
    }

    virtual float GetValue(float x, float y, float z) const override {
        TransformW2O(x, y, z);
        float leftValue = left->GetValue(x, y, z);
        float rightValue = right->GetValue(x, y, z);
        float k = static_cast<float>(mBlend);
        if (k <= 0.0f) k = 0.0001f;
        
        // Super-elliptic blending for difference
        // Using the formula: log(exp(k*a) + exp(-k*b))/k
        return std::log(std::exp(k * leftValue) + std::exp(-k * rightValue)) / k;
    }

protected:
    int mBlend;
};
