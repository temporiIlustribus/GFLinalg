#ifndef GFLINALG_GFSTORAGE_H
#define GFLINALG_GFSTORAGE_H

#include "GFSPlinalg.hpp"

/**
 * @defgroup Storage
 */

namespace GFlinalg {
template <typename T>
struct GFElemRef<BasicGFElem<T>> {
    using TBase = BasicGFElem<T>;

    constexpr GFElemRef() = default;
    constexpr GFElemRef(T& ref, const GFElemState<T> st) : mState(st), mValue(ref) {}
    explicit constexpr GFElemRef(TBase& base) : mState(base.getState()), mValue(base.val()) {}
    constexpr GFElemRef(GFElemRef& other): mState(other.mState), mValue(other.mValue) {}

    operator TBase() {
        return TBase(mValue, mState);
    }

    inline T& val() {
        return mValue;
    }

    [[nodiscard]] inline size_t gfDegree() const {
        return BasicGFElem<T>(mValue, mState).gfDegree();
    }

    [[nodiscard]] inline size_t gfOrder() const {
        return BasicGFElem<T>(mValue, mState).gfOrder();
    }

    inline T getMod() const {
        return BasicGFElem<T>(mValue, mState).getMod();
    }

    inline T reduce() {
        return BasicGFElem<T>(mValue, mState).reduce();
    }

    inline TBase getInverse() const {
        return BasicGFElem<T>(mValue, mState).getInverse();
    }

    inline GFElemRef& invert() {
        auto res = BasicGFElem<T>(mValue, mState).invert();
        mValue = res.val();
        mState = res.getState();
        return *this;
    }

    [[nodiscard]] inline size_t degree(size_t startPos = 1) const {
        return BasicGFElem<T>(mValue, mState).degree(startPos);
    }

    inline GFElemRef& operator=(const GFElemRef& other) {
        mState = other.mState;
        mValue = other.mValue;
        return *this;
    }

    inline GFElemRef& operator=(const TBase& other) {
        mState = other.getState();
        mValue = other.val();
        return *this;
    }

    inline friend TBase operator+(const GFElemRef& a, const GFElemRef& b) {
        return BasicGFElem<T>(a.mValue, a.mState) + BasicGFElem<T>(b.mValue, b.mState);
    }

    inline friend TBase operator*(const GFElemRef& a, const GFElemRef& b) {
        return BasicGFElem<T>(a.mValue, a.mState) * BasicGFElem<T>(b.mValue, b.mState);
    }

    inline friend TBase operator/(const GFElemRef& a, const GFElemRef& b) {
        return BasicGFElem<T>(a.mValue, a.mState) / BasicGFElem<T>(b.mValue, b.mState);
    }

    inline GFElemRef& operator+=(const GFElemRef& a) {
        *this = *this + a;
        return *this;
    }

    inline GFElemRef& operator*=(const GFElemRef& a) {
        *this = *this * a;
        return *this;
    }

    inline GFElemRef& operator/=(const GFElemRef& a) {
        *this = *this / a;
        return *this;
    }

    inline GFElemRef& operator+=(const TBase& a) {
        *this = *this + a;
        return *this;
    }

    inline GFElemRef& operator*=(const TBase& a) {
        *this = *this * a;
        return *this;
    }

    inline GFElemRef& operator/=(const TBase& a) {
        *this = *this / a;
        return *this;
    }

    template <class T1>
    inline friend bool operator==(const GFElemRef<T1>& lhs, const GFElemRef<T1>& rhs) {
        return lhs.mValue == rhs.mValue && lhs.mState == rhs.mState;
    }

    template <class T1>
    inline friend bool operator==(const GFElemRef<T1>& lhs, const TBase & rhs) {
        return lhs.mValue == rhs.val() && lhs.mState == rhs.getState();
    }

    template <class T1>
    inline friend bool operator!=(const GFElemRef<T1>& lhs, const TBase & rhs) {
        return !(lhs == rhs);
    }

    template <class T1>
    inline friend bool operator!=(const GFElemRef<T1>& lhs, const GFElemRef<T1>& rhs) {
        return !(lhs == rhs);
    }

    template <class T1>
    inline friend bool operator<(const GFElemRef<T1>& lhs, const GFElemRef<T1>& rhs) {
        return lhs.mValue < rhs.mValue;
    }

    template <class T1>
    inline friend bool operator>(const GFElemRef<T1>& lhs, const GFElemRef<T1>& rhs) {
        return rhs.mValue < lhs.mValue;
    }

    template <class T1>
    inline friend bool operator<=(const GFElemRef<T1>& lhs, const GFElemRef<T1>& rhs) {
        return lhs.mValue <= rhs.mValue;
    }

    template <class T1>
    inline friend bool operator>=(const GFElemRef<T1>& lhs, const GFElemRef<T1>& rhs) {
        return rhs <= lhs;
    }

protected:
    T& mValue;
    const GFElemState<T>& mState;
};

template <typename T>
class GFElemPtr<BasicGFElem<T>> {
public:
    using Elem = BasicGFElem<T>;

    constexpr GFElemPtr() = default;
    constexpr GFElemPtr(const GFElemPtr& other): mRef(other.mRef) {}
    constexpr GFElemPtr(GFElemRef<Elem>& other): mRef(other) {}

    inline GFElemPtr& operator=(const GFElemPtr& other) {
        mRef = other.mRef;
        return *this;
    }

    GFElemRef<Elem> operator*() {
        return mRef;
    }

    GFElemRef<Elem>* operator->() {
        return &mRef;
    }

protected:
    GFElemRef<Elem> mRef;
};

template<typename T, size_t R, size_t C>
class MatrixEngine<BasicGFElem<T>, R, C> {
public:
    constexpr MatrixEngine() = default;

    constexpr explicit MatrixEngine(GFElemState<T>& state): mState(state) {}

    constexpr MatrixEngine& operator =(MatrixEngine&&) noexcept = default;
    constexpr MatrixEngine& operator =(MatrixEngine const&) = default;

    constexpr GFElemRef<BasicGFElem<T>> operator()(size_t i, size_t j) {
        T& ref = mData.at(C * i + j);
        return GFElemRef<BasicGFElem<T>>(ref, mState);
    }

    [[nodiscard]] constexpr size_t columns() const noexcept {
        return C;
    }

    [[nodiscard]] constexpr size_t rows() const noexcept {
        return R;
    }

    [[nodiscard]] constexpr size_t size() const noexcept {
        return R * C;
    }

    constexpr void swap(MatrixEngine& rhs) noexcept {
        std::swap(mData, rhs.mData);
        std::swap(mState, rhs.mState);
    }

private:
    GFElemState<T> mState;
    std::array<T, R * C> mData;
};

struct AccessorBasic {
    template <typename T, size_t R>
    struct Accessor {
        using pointer      = MatrixEngine<T, R, 1>*;
        using reference    = GFElemRef<T>;

        constexpr reference operator()(pointer p, ptrdiff_t i) const noexcept {
            return (*p)(i, 0);
        }
    };
};
}

#endif // GFLINALG_GFSTORAGE_H
