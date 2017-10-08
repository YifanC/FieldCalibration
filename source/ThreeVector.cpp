#include "../include/ThreeVector.hpp"

template<typename ValueType>
ThreeVector<ValueType>::ThreeVector() {}

template<typename ValueType>
ThreeVector<ValueType>::ThreeVector(std::array<ValueType, 3> input_array) {
    this->vec = input_array;
}

template<typename ValueType>
ThreeVector<ValueType>::ThreeVector(ValueType input_array[3]) {
    for (unsigned coord = 0; coord < vec.size(); coord++)
        this->vec[coord] = input_array[coord];
}

template<typename ValueType>
ThreeVector<ValueType>::ThreeVector(ValueType x_input, ValueType y_input, ValueType z_input) {
    this->vec = {x_input, y_input, z_input};
}

template<typename ValueType>
ThreeVector<ValueType>::ThreeVector(std::initializer_list<ValueType> input_list) {
    std::copy(input_list.begin(), input_list.end(), this->vec.begin());
}

// From three variable input constructor inherited
template<typename ValueType>
ThreeVector<ValueType>::ThreeVector(const TVector3 &InputVector) : ThreeVector(InputVector.X(), InputVector.Y(),
                                                                               InputVector.Z()) {

}

template<typename ValueType>
ThreeVector<ValueType> &ThreeVector<ValueType>::operator+=(const ThreeVector<ValueType> &vec_right) {
    for (unsigned coord = 0; coord < this->vec.size(); coord++)
        this->vec[coord] += vec_right.vec[coord];
    return *this;
}

template<typename ValueType>
ThreeVector<ValueType> &ThreeVector<ValueType>::operator-=(const ThreeVector<ValueType> &vec_right) {
    for (unsigned coord = 0; coord < this->vec.size(); coord++)
        this->vec[coord] -= vec_right.vec[coord];
    return *this;
}

template<typename ValueType>
ThreeVector<ValueType> &ThreeVector<ValueType>::operator*=(const ValueType &scal_right) {
    for (unsigned coord = 0; coord < this->vec.size(); coord++)
        this->vec[coord] *= scal_right;
    return *this;
}

template<typename ValueType>
ThreeVector<ValueType> &ThreeVector<ValueType>::operator/=(const ValueType &scal_right) {
    for (unsigned coord = 0; coord < this->vec.size(); coord++)
        this->vec[coord] /= scal_right;
    return *this;
}

template<typename ValueType>
ThreeVector<ValueType> operator+(ThreeVector<ValueType> vec_left, const ThreeVector<ValueType> &vec_right) {
    return vec_left += vec_right;
}

template<typename ValueType>
ThreeVector<ValueType> operator-(ThreeVector<ValueType> vec_left, const ThreeVector<ValueType> &vec_right) {
    return vec_left -= vec_right;
}

template<typename ValueType>
ThreeVector<ValueType> operator*(ThreeVector<ValueType> vec_left, const ValueType &scal_right) {
    return vec_left *= scal_right;
}

template<typename ValueType>
ThreeVector<ValueType> operator*(const ValueType &scal_left, ThreeVector<ValueType> vec_right) {
    return vec_right *= scal_left;
}

template<typename ValueType>
ThreeVector<ValueType> operator/(ThreeVector<ValueType> vec_left, const ValueType &scal_right) {
    return vec_left /= scal_right;
}

template<typename ValueType>
ValueType &ThreeVector<ValueType>::at(const unsigned long index) {
    return vec.at(index);
}

template<typename ValueType>
ValueType &ThreeVector<ValueType>::operator[](const unsigned long index) {
    return this->vec[index];
}

template<typename ValueType>
const bool ThreeVector<ValueType>::operator==(const ThreeVector<ValueType> &vec_right) {
    return ((this->vec[0] == vec_right.vec[0]) && (this->vec[1] == vec_right.vec[1]) &&
            (this->vec[2] == vec_right.vec[2]));
}

template<typename ValueType>
const bool ThreeVector<ValueType>::operator!=(const ThreeVector<ValueType> &vec_right) {
    return !(*this == vec_right);
}

template<typename ValueType>
unsigned long ThreeVector<ValueType>::size() {
    return vec.size();
}

template<typename ValueType>
void ThreeVector<ValueType>::print() {
    std::cout << "(x,y,z) = (" << vec[0] << "," << vec[1] << "," << vec[2] << ")" << std::endl;
}

template<typename ValueType>
ValueType ThreeVector<ValueType>::GetNorm() const {
    ValueType norm = 0;
    for (const auto &vec_iterator : vec)
        norm += std::pow(vec_iterator, 2);
    return std::sqrt(norm);
}

template<typename ValueType>
ThreeVector<ValueType> ThreeVector<ValueType>::GetUnitVector() const {
    return *this / GetNorm();
}

template<typename ValueType>
ValueType ThreeVector<ValueType>::VectorNorm(const ThreeVector<ValueType> &vec_input) {
    return vec_input.GetNorm();
}

template<typename ValueType>
ThreeVector<ValueType> ThreeVector<ValueType>::UnitVector(const ThreeVector<ValueType> &vec_input) {
    return vec_input.GetUnitVector();
}

template<typename ValueType>
ValueType
ThreeVector<ValueType>::DotProduct(const ThreeVector<ValueType> &vec_left, const ThreeVector<ValueType> &vec_right) {
    ValueType scal_res = 0.0;
    for (unsigned coord = 0; coord < vec_left.vec.size(); coord++)
        scal_res += vec_left.vec[coord] * vec_right.vec[coord];

    return scal_res;
}

template<typename ValueType>
ThreeVector<ValueType>
ThreeVector<ValueType>::VectorProduct(const ThreeVector<ValueType> &vec_left, const ThreeVector<ValueType> &vec_right) {
    ThreeVector<ValueType> vec_res;
    for (unsigned coord = 0; coord < vec_res.vec.size(); coord++)
        vec_res.vec[coord] =
                vec_left.vec[(1 + coord) % vec_res.vec.size()] * vec_right.vec[(2 + coord) % vec_res.vec.size()] -
                vec_right.vec[(1 + coord) % vec_res.vec.size()] * vec_left.vec[(2 + coord) % vec_res.vec.size()];
    return vec_res;
}

template<typename ValueType>
std::array<ValueType, 3> ThreeVector<ValueType>::GetArray() {
    return this->vec;
}

template<typename ValueType>
std::vector<ValueType> ThreeVector<ValueType>::GetStdVector() {
    std::vector<ValueType> vec_res;
    for (unsigned coord = 0; coord < this->vec.size(); coord++) {
        vec_res.push_back(this->vec[coord]);
    }
    return vec_res;
}