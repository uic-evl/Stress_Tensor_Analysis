// Vector.cpp -- implementation of Vector class

#include "math.h"
#include "vector.hh"

// Note that much of the implementation of this class is in-line in Vector.h


// Convert this Vector from degrees to radians
Vector& Vector::Degrees_to_Radians(void)
{
    this->ScalarMult((float)(M_PI / 180.0));
    return *this;
}

// Convert this Vector from radians to degrees
Vector& Vector::Radians_to_Degrees(void)
{
    this->ScalarMult((float)(180.0 / M_PI));
    return *this;
}
