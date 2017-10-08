#include <array>
#include <initializer_list>
#include <algorithm>
#include <cmath>

#include <TVector3.h>

// Temporary
#include <typeinfo>
#include <iostream>

#ifndef THREEVECTOR_H
#define THREEVECTOR_H

template<typename ValueType>
class ThreeVector {
private:
    // using std::array<float,3>::array to store data in vec
    std::array<ValueType, 3> vec;

public:
    ThreeVector(); // Default constructor
    ThreeVector(std::array<ValueType, 3>); // Construction by std::array
    ThreeVector(ValueType[3]);

    ThreeVector(ValueType, ValueType, ValueType);

    ThreeVector(std::initializer_list<ValueType>); // Construction by {x,y,z}
    ThreeVector(const TVector3 &InputVector);

    // += and -= operators add entries to each other
    ThreeVector<ValueType> &operator+=(const ThreeVector<ValueType> &);

    ThreeVector<ValueType> &operator-=(const ThreeVector<ValueType> &);

    ThreeVector<ValueType> &operator*=(const ValueType &);

    ThreeVector<ValueType> &operator/=(const ValueType &);

    const bool operator==(const ThreeVector<ValueType> &);

    const bool operator!=(const ThreeVector<ValueType> &);

    // [] operator accesses the stored data in vec
    ValueType &at(const unsigned long);

    ValueType &operator[](const unsigned long);

    unsigned long size();

    void print();

    // Friend functions for + and - operators derived from += and -=
    template<typename ValueType_Op>
    friend ThreeVector<ValueType_Op> operator+(ThreeVector<ValueType_Op>, const ThreeVector<ValueType_Op> &);

    template<typename ValueType_Op>
    friend ThreeVector<ValueType_Op> operator-(ThreeVector<ValueType_Op>, const ThreeVector<ValueType_Op> &);

    template<typename ValueType_Op>
    friend ThreeVector<ValueType_Op> operator*(const ValueType_Op &, ThreeVector<ValueType_Op>);

    template<typename ValueType_Op>
    friend ThreeVector<ValueType_Op> operator*(ThreeVector<ValueType_Op>, const ValueType_Op &);

    template<typename ValueType_Op>
    friend ThreeVector<ValueType_Op> operator/(ThreeVector<ValueType_Op>, const ValueType_Op &);

    // Vector operations non-static
    ValueType GetNorm() const;

    ThreeVector<ValueType> GetUnitVector() const;

    // Vector operations as static functions
    static ValueType VectorNorm(const ThreeVector<ValueType> &InputVector);

    static ThreeVector<ValueType> UnitVector(const ThreeVector<ValueType> &InputVector);

    static ThreeVector<ValueType> VectorProduct(const ThreeVector<ValueType> &, const ThreeVector<ValueType> &);

    static ValueType DotProduct(const ThreeVector<ValueType> &, const ThreeVector<ValueType> &);

    std::array<ValueType, 3> GetArray();

    std::vector<ValueType> GetStdVector();
};

// Inline implimentation of all members (due to template!)
#include "../source/ThreeVector.cpp"

#endif