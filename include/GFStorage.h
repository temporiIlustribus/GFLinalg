#ifndef GFLINALG_GFSTORAGE_H
#define GFLINALG_GFSTORAGE_H

#include "GFSPlinalg.hpp"

/**
 * @defgroup Storage
 */

namespace GFlinalg {

template <typename T>
class GFElemRef;

template <typename T>
class GFElemPtr;

template <typename T>
struct GFElemRef<BasicGFElem<T>> {
    using TBase = BasicGFElem<T>;

    constexpr GFElemRef() = default;
    constexpr GFElemRef(T& ref, const GFElemState<T> st) : mState(st), mValue(ref) {}
    explicit constexpr GFElemRef(const TBase& base) : mState(base.getState()), mValue(base.val()) {}
    constexpr GFElemRef(const GFElemRef& other): mState(other.mState), mValue(other.mValue) {}

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

    inline TBase getInverse() const;

    inline GFElemRef& invert();

    [[nodiscard]] inline size_t degree(size_t startPos = 1) const {
        return BasicGFElem<T>(mValue, mState).degree(startPos);
    }

    inline GFElemRef& operator=(const GFElemRef& other) {
        mState = other.mState;
        mValue = other.mValue;
        return *this;
    }

    inline friend GFElemRef operator+(const GFElemRef& a, const GFElemRef& b) {
        return BasicGFElem<T>(a.mValue, a.mState) + BasicGFElem<T>(b.mValue, b.mState);
    }

    inline friend GFElemRef operator*(const GFElemRef& a, const GFElemRef& b) {
        return BasicGFElem<T>(a.mValue, a.mState) * BasicGFElem<T>(b.mValue, b.mState);
    }

    inline friend GFElemRef operator/(const GFElemRef& a, const GFElemRef& b) {
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

    inline  GFElemRef& operator/=(const GFElemRef& a) {
        *this = *this / a;
        return *this;
    }

    template <class T1>
    inline friend bool operator==(const GFElemRef<T1>& lhs, const GFElemRef<T1>& rhs) {
        return lhs.mValue == rhs.mValue && lhs.mState == rhs.mState;
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
        return (lhs < rhs) || (lhs == rhs);
    }

    template <class T1>
    inline friend bool operator>=(const GFElemRef<T1>& lhs, const GFElemRef<T1>& rhs) {
        return (rhs < lhs) || (lhs == rhs);
    }

protected:
    T& mValue;
    const GFElemState<T>& mState;
};

template <typename T>
class GFElemPtr<BasicGFElem<T>> {
public:
    constexpr GFElemPtr() = default;
    constexpr GFElemPtr(const GFElemPtr& other): mRef(other.mRef) {}

    inline GFElemPtr& operator=(const GFElemPtr& other) {
        mRef = other.mRef;
        return *this;
    }

    GFElemRef<T> operator*() {
        return mRef;
    }

    GFElemRef<T>* operator->() {
        return &mRef;
    }

protected:
    GFElemRef<T>& mRef;
};

/*template <class T>
class MDSpanIter {};

template <typename Accessor>
class MDSpan {
    using TVal  = typename Accessor::TVal;
    using TRef  = GFElemRef<TVal>;
    using TIter = MDSpanIter<TVal>;

    template <class... IndexType>
    constexpr TRef operator()(IndexType... indices) const noexcept;

    template <class IndexType, size_t N>
    constexpr TRef operator()(const std::array<IndexType, N>& indices) const noexcept;

    constexpr TIter begin() noexcept;
    constexpr TIter end() noexcept;
};*/
}

#endif // GFLINALG_GFSTORAGE_H
