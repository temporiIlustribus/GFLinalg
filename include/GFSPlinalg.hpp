#pragma once

#include "GFbase.hpp"

/**
 * @ingroup linalg
 *
 * Below are listed single parameter version.
 *
 * @note Comments to these classes are identical to their analogues with
 * two template parameters.
 *
 */
namespace GFlinalg {

template <typename MPT>
struct GFElemState {
    size_t SZ, order;
    MPT modPol;

    GFElemState() : SZ(), order(), modPol() {};
    GFElemState(size_t size, size_t order, const MPT& modPol) : SZ(size), order(order), modPol(modPol) {}
    GFElemState(size_t size, size_t order) : SZ(size), order(order) {}

    GFElemState(const GFElemState& other):
        SZ(other.SZ), order(other.order), modPol(other.modPol) {}

    bool operator != (const GFElemState<MPT>& other) {
        return !(*this == other);
    }

    GFElemState& operator=(const GFElemState& other) {
        SZ = other.SZ;
        order = other.order;
        modPol = other.modPol;
        return *this;
    }
};

template <typename T>
bool operator == (const GFElemState<T>& one, const GFElemState<T>& other) {
    return one.SZ == other.SZ && one.order == other.order && one.modPol == other.modPol;
}

template <typename T>
class IGFElem {
public:
    using State = GFElemState<T>;

    State& getState() { return mState; }

    const State& getState() const { return mState; }

    virtual ~IGFElem() = default;

protected:
    State mState;
};

/**
 * Polynomial based single template parameter GF element class
 *
 * All math operations are performed directly in polynomial form.
 *
 * Time complexity:
 * <ul>
 *   <li>"+" - O(1)</li>
 *   <li>"*" - O(n^2)</li>
 *   <li>"/" - O(log^2(order) + n^2)</li>
 * </ul>
 *
 *  @note Primitive modulus polynomial is stored in class instances.
 */
template <class T>
class BasicGFElem : public IGFElem<T> {
protected:
    T value;

    using State = typename IGFElem<T>::State;
    using IGFElem<T>::mState;

public:
    explicit constexpr BasicGFElem() : value() {
        this->mState(0, 0);
    }

    explicit constexpr BasicGFElem(const T& value, const T& modulus, bool doReduce = true):
        value(value) {

        mState = State();
        mState.modPol = modulus;
        mState.SZ = op::modPolDegree<T>(modulus);
        mState.order = 1 << mState.SZ;

        if (doReduce)
            this->reduce();
    }

    explicit constexpr BasicGFElem(const T& value, const GFElemState<T>& state):
        value(value) {
        this->mState = state;
    }

    /**
     * Construct %BasicGFElem using a pair of \c Iter s that provide target polynomial's coefficients.
     * Coefficients are passed left to right.
     *
     * \c Iter s must meet the requirements of \c LegacyRandomAccessIterator<T>.
     *
     * @example <tt>{1,0,1,0,0} -> x^4 + x^2 (10100)</tt>
     * @example <tt>{1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)</tt>
     *
     * @note Polynomial is automatically reduced after initialization.
     */
    template <typename Iter>
    explicit constexpr BasicGFElem(Iter begin, Iter end, const T& modulus): value(0) {
        static_assert(std::is_convertible_v<decltype(*begin), T>);

        mState = State(0, 0, modulus);

        while (begin++ != end) {
            value |= (static_cast<T>(*begin) & 1);
            value <<= 1;
        }

        mState.SZ = op::modPolDegree<T>();
        mState.order = 1 << mState.SZ;

        this->reduce();
    }

    /**
     * Construct %BasicGFElem using a pair of \c Iter s that provide target polynomial's coefficients and a pair
     * of \c Iter s that provide target modulus polynomial's coefficients.
     *
     * @see BasicGFElem::BasicGFElem(Iter, Iter , const T&)
     */
    template <typename Iter>
    explicit constexpr BasicGFElem(Iter begin, Iter end, Iter beginMod, Iter endMod): value(0) {
        static_assert(std::is_convertible_v<decltype(*begin), T>);

        mState = State(0, 0, 0);

        while (begin++ != end) {
            value |= (static_cast<T>(*begin) & 1);
            value <<= 1;
        }

        while (beginMod++ != endMod) {
            mState.modPol |= (static_cast<T>(*beginMod) & 1);
            mState.modPol <<= 1;
        }

        mState.SZ = op::modPolDegree<T>();
        mState.order = 1 << mState.SZ;

        reduce();
    }

    /**
     * @return Copy of the polynomial in the contained form.
     */
    T val() const { return value; }

    /**
     * @warning Do not use this function to alter the contained value. Use the \c reduce() method instead.
     *
     * @see BasicGFElem::reduce()
     *
     * @return Polynomial in the contained form.
     */
    T& val() { return value; }

    /**
     * @return \c n For \c GF(2^n).
     */
    [[nodiscard]] size_t gfDegree() const { return mState.SZ; }

    /**
     * @return \c 2^n For \c GF(2^n).
     */
    [[nodiscard]] size_t gfOrder() const { return mState.order; }

    T getMod() const { return mState.modPol; }

    /**
     * Reduce the element by modulus polynomial.
     * @return \c value reduced by \c mState.modPol.
     */
    T reduce() {
        auto pos  = mState.order - mState.SZ;
        uint8_t i = 1;

        while (value >= 1U << mState.SZ) {
            if ((value >> (mState.order - i)) & 1)
                value ^= mState.modPol << (pos - i);
            ++i;
        }

        return value;
    }

    /**
     * @return For polynomial \c pol return \c pol^(-1), where <tt>pol * pol^(-1) = pol^(-1) * pol = 1</tt>
     */
    BasicGFElem getInverse() const { return pow<BasicGFElem>(*this, mState.order - 2); }

    /**
     * Set \c *this = \c getInverse();
     */
    BasicGFElem& invert() {
        *this = pow<BasicGFElem>(*this, mState.order - 2);
        return *this;
    }

    [[nodiscard]] size_t degree(size_t startPos = 1) const { return mState.order - leadElemPos<T>(value, startPos); }

    friend BasicGFElem operator+(const BasicGFElem& a, const BasicGFElem& b) {
        if (a.mState.modPol != b.mState.modPol)
            throw std::runtime_error("Cannot perform addition for elements of different fields");

        return BasicGFElem(a.val() ^ b.val(), a.mState.modPol);
    }

    BasicGFElem& operator+=(const BasicGFElem& other) {
        *this = BasicGFElem(val() ^ other.val(), mState.modPol);
        return *this;
    }

    friend BasicGFElem operator*(const BasicGFElem& a, const BasicGFElem& b) {
        if (a.mState.modPol != b.mState.modPol)
            throw std::runtime_error("Cannot perform multiplication for elements of different fields");

        return {op::polMul<BasicGFElem>(a, b)};
    }

    BasicGFElem& operator*=(const BasicGFElem& other) {
        if (mState.modPol != other.mState.modPol)
            throw std::runtime_error("Cannot perform multiplication for elements of different fields");

        *this = *this * other;

        return *this;
    }

    friend BasicGFElem operator/ (const BasicGFElem& a, const BasicGFElem& b) {
        if (a.mState.modPol != b.mState.modPol)
            throw std::runtime_error("Cannot perform division for elements of different fields");

        return {op::polDiv<BasicGFElem>(a, b)};
    }

    BasicGFElem& operator/=(const BasicGFElem& other) {
        *this = *this / other;
        return *this;
    }

    template <class T1>
    friend bool operator==(const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
    template <class T1>
    friend bool operator!=(const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);

    template <class T1>
    friend bool operator<(const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
    template <class T1>
    friend bool operator>(const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
    template <class T1>
    friend bool operator<=(const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
    template <class T1>
    friend bool operator>=(const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
};

template <class T>
bool operator<(const BasicGFElem<T>& lhs, const BasicGFElem<T>& rhs) {
    if (lhs.mState.modPol != rhs.mState.modPol)
        throw std::runtime_error("Cannot compare elements of different fields");

    return lhs.value < rhs.value;
}

template <class T>
bool operator>(const BasicGFElem<T>& lhs, const BasicGFElem<T>& rhs) {
    if (lhs.mState.modPol != rhs.mState.modPol)
        throw std::runtime_error("Cannot compare elements of different fields");

    return lhs.value > rhs.value;
}
template <class T>
bool operator<=(const BasicGFElem<T>& lhs, const BasicGFElem<T>& rhs) {
    if (lhs.mState.modPol != rhs.mState.modPol)
        throw std::runtime_error("Cannot compare elements of different fields");

    return lhs.value <= rhs.value;
}

template <class T>
bool operator>=(const BasicGFElem<T>& lhs, const BasicGFElem<T>& rhs) {
    if (lhs.mState.modPol != rhs.mState.modPol)
        throw std::runtime_error("Cannot compare elements of different fields");

    return lhs.value >= rhs.value;
}

template <class T>
bool operator==(const BasicGFElem<T>& lhs, const BasicGFElem<T>& rhs) {
    return (lhs.value == rhs.value) && (lhs.mState.modPol == rhs.mState.modPol);
}

template <class T>
bool operator!=(const BasicGFElem<T>& lhs, const BasicGFElem<T>& rhs) {
    return !(lhs == rhs);
}

/**
 * LUT based single template parameter GF element class. Uses conversion to power of the
 * primitive element through a LUT internally.
 *
 * LUT stands for "look-up table".
 *
 * Time complexity:
 * <ul>
 *  <li>"+" - O(1)</li>
 *  <li>"*" - O(1)</li>
 *  <li>"/" - O(1)</li>
 * </ul>
 *
 * Memory complexity: \c O(2^n)
 *
 * Primitive modulus polynomial is stored in individual instances.
 *
 * Uses a LUTPair pointer internally. Use \c LUTVectPair<T>(params...) to make the LUTPair.
 *
 * The LUT is generated for the specific field. Use LUT-s with proper parameters.
 *
 */
template <class T>
class PowGFElem : public BasicGFElem<T> {
    using BasicGFElem<T>::BasicGFElem;

protected:
    LUTVectPair<T> const* alphaToIndex = nullptr;

public:
    explicit PowGFElem(const BasicGFElem<T>& pol) : BasicGFElem<T>(pol) {}

    explicit PowGFElem(const T& val, const T& modPol, LUTVectPair<T> const* lut)
        : BasicGFElem<T>(val, modPol), alphaToIndex(lut) {}

    explicit PowGFElem(const BasicGFElem<T>& pol, LUTVectPair<T> const* lut)
        : BasicGFElem<T>(pol), alphaToIndex(lut) { }

    LUTVectPair<T> const* lut() const { return alphaToIndex; }

    /**
     * Multiply elements using LUTs.
     *
     * Algorithm:
     * <ol>
     *  <li>\c a is the primitive element</li>
     *  <li><tt>x * y = z</tt></li>
     *  <li><tt>x-LUT-> a^u</tt></li>
     *  <li><tt>y-LUT-> a^v</tt></li>
     *  <li><tt>x * y = a^u * a^v = a^(u+v) -LUT-> z</tt></li>
     * </ol>
     */
    PowGFElem operator*(const PowGFElem& other) const {
        if (other.mState.modPol != this->mState.modPol)
            throw std::runtime_error("Cannot perform multiplication for elements of different fields");

        if (this->value == 0 || other.value == 0)
            return PowGFElem(0, this->mState.modPol);

        return PowGFElem(alphaToIndex->indToPol[
            alphaToIndex->polToInd[this->value] +
            alphaToIndex->polToInd[other.value]],
                         this->mState.modPol, alphaToIndex);
    }

    PowGFElem& operator*=(const PowGFElem& other) {
        this->val() = alphaToIndex->indToPol[
            alphaToIndex->polToInd[this->value] +
            alphaToIndex->polToInd[other.value]];

        return *this;
    }

    /**
     * Divide elements using LUTs.
     *
     * Algorithm:
     * <ol>
     *  <li>\c a is the primitive element</li>
     *  <li><tt>x / y = z</tt></li>
     *  <li><tt>x-LUT-> a^u</tt></li>
     *  <li><tt>y-LUT-> a^v</tt></li>
     *  <li><tt> x / y = a^u * a^v = a^(u-v) -LUT-> z</tt></li>
     * </ol>
     */
    PowGFElem operator/(const PowGFElem& other) const {
        if (other.mState.modPol != this->mState.modPol)
            throw std::runtime_error("Cannot perform division for elements of different fields");

        if (this->value == 0)
            return PowGFElem(0, this->mState.modPol);

        if (other.value == 0)
            throw std::out_of_range("Division by zero");

        auto temp = alphaToIndex->polToInd[this->value];

        if (temp < alphaToIndex->polToInd[other.value])
            temp += this->mState.order - 1;

        return PowGFElem(alphaToIndex->indToPol[temp - alphaToIndex->polToInd[other.value]],
            this->mState.modPol, alphaToIndex);
    }

    PowGFElem& operator/=(const PowGFElem& other) {
        if (this->value == 0)
            return PowGFElem(0);

        auto temp = alphaToIndex->polToInd[this->value];

        if (temp < alphaToIndex->polToInd[other.value])
            temp += this->order - 1;

        *this->val = alphaToIndex->indToPol[temp - alphaToIndex->polToInd[other.value]];

        return *this;
    }

    PowGFElem operator+(const PowGFElem& other) const {
        if (other.mState.modPol != this->mState.modPol)
            throw std::runtime_error("Cannot perform addition for elements of different fields");

        return PowGFElem(this->val() ^ other.val(), this->mState.modPol, alphaToIndex);
    }

    PowGFElem operator+=(const PowGFElem& other) {
        *this = *this + other;
        return *this;
    }

    template <class T1>
    friend PowGFElem<T1> pow(const PowGFElem<T1>& val, size_t power);
};

template <class T>
PowGFElem<T> pow(const PowGFElem<T>& val, size_t power) {
    size_t index = (val.alphaToIndex->polToInd[val.value] * power) % (val.gfOrder() - 1);
    return PowGFElem<T>(val.alphaToIndex->indToPol[index], val.mState.modPol, val.alphaToIndex);
}

/**
 * Table based GF element class. All math operations are performed via the multiplication table
 * and the division table.
 *
 * Time complexity:
 * <ul>
 *  <li>"+" - O(1)</li>
 *  <li>"*" - O(1)</li>
 *  <li>"/" - O(1)</li>
 * </ul>
 *
 * Memory complexity: \c O(2^n)
 *
 * Primitive modulus polynomial is stored in individual instances.
 *
 * %TableGFElem uses a \c std::vector<T> pointer internally to access the tables.
 * Use \c makeMulTable(params...) to make the tables.
 *
 * The LUT is generated for the specific field.
 * Use LUT-s with proper parameters.
 *
 * @note This documentation will imply <tt>Table[a][b] = Table[a * order +  b]</tt>.
 */
template <class T>
class TableGFElem : public BasicGFElem<T> {
    using BasicGFElem<T>::BasicGFElem;
private:
    std::vector<T> const* mulTable = nullptr;
    std::vector<T> const* divTable = nullptr;

public:
    explicit TableGFElem(const BasicGFElem<T>& pol) : BasicGFElem<T>(pol) {}

    explicit TableGFElem(const T& val, const T& modPol, std::vector<T> const* mulTable_,
        std::vector<T> const* divTable_) : BasicGFElem<T>(val, modPol) {
        mulTable = mulTable_;
        divTable = divTable_;
    }

    explicit TableGFElem(const BasicGFElem<T>& pol, std::vector<T> const* mulTable_,
        std::vector<T> const* divTable_): BasicGFElem<T>(pol) {
        mulTable = mulTable_;
        divTable = divTable_;
    }

    /**
     * Create the multiplication table (<tt>table[p1][p2] = p1 * p2</tt>) for given \c modPol:
     */
    static std::vector<T> makeMulTable(const T& modPol) {
        size_t order = 1u << op::modPolDegree<T>(modPol);

        std::vector<T> out(order * order);
        size_t res;

        for (size_t i = 0; i < order; ++i) {
            for (size_t j = i; j < order; ++j) {
                BasicGFElem<T> a(i, modPol);
                BasicGFElem<T> b(j, modPol);

                res = (a * b).val();

                out[a.val() * order + b.val()] = res;
                out[b.val() * order + a.val()] = res;
            }
        }

        return out;
    }

    /**
     * Create the inverted multiplication table using the multiplication table.
     *
     * Inverted multiplication table is a table such as:
     * <ol>
     *  <li><tt>table[p1 * p2][p2] = p1</tt></li>
     *  <li><tt>table[p1 * p2][p1] = p2</tt></li>
     * </ol>
     */
    static std::vector<T> makeInvMulTable(std::vector<T> const* mulTable, const T& modPol) {
        size_t order = 1U << op::modPolDegree<T>(modPol);

        std::vector<T> out(order * order);

        size_t res = 0;

        for (size_t i = 0; i < order; ++i) {
            for (size_t j = i; j < order; ++j) {
                BasicGFElem<T> a(i, modPol);
                BasicGFElem<T> b(j, modPol);

                res = mulTable->at(a.val() * order + b.val()) * order;

                out[res + a.val()] = b.val();
                out[res + b.val()] = a.val();
            }
        }

        return out;
    }

    /**
     * Create the division table (<tt>table[po1][p2] = p1 / p2</tt>) for given \c modPol.
     */
    static std::vector<T> makeDivTable(const T& modPol) {
        size_t order = 1U << op::modPolDegree<T>(modPol);

        std::vector<T> out(order * order);

        for (size_t i = 0; i < order; ++i) {
            for (size_t j = 0; j < order; ++j) {
                BasicGFElem<T> a(i, modPol);
                BasicGFElem<T> b(j, modPol);

                out[a.value * order + b.value] = (a / b).val();
            }
        }

        return out;
    }

    std::vector<T>* multable() { return mulTable; }
    std::vector<T>* multable() const { return mulTable; }
    std::vector<T>* divtable() { return divTable; }
    std::vector<T>* divtable() const { return divTable; }

    TableGFElem operator+(const TableGFElem& other) const {
        if (other.mState.modPol != this->mState.modPol)
            throw std::runtime_error("Cannot perform addition for elements of different fields");

        return TableGFElem(this->val() ^ other.val(), this->mState.modPol, mulTable, divTable);
    }

    TableGFElem& operator+=(const TableGFElem& other) {
        *this = *this + other;
        return *this;
    }

    TableGFElem operator*(const TableGFElem& other) const {
        if (other.mState.modPol != this->mState.modPol)
            throw std::runtime_error("Cannot perform multiplication for elements of different fields");

        return TableGFElem(mulTable->at(this->val() * this->mState.order + other.val()),
            this->mState.modPol, mulTable, divTable);
    }

    TableGFElem& operator*=(const TableGFElem& other) {
        *this = *this * other;
        return *this;
    }

    TableGFElem operator/(const TableGFElem& other) const {
        if (other.mState.modPol != this->mState.modPol)
            throw std::runtime_error("Cannot perform division for elements of different fields");

        if (other.value == 0)
            throw std::out_of_range("Division by zero");

        return TableGFElem(divTable->at(this->val() * this->mState.order + other.val()),
            this->mState.modPol, mulTable, divTable);
    }

    TableGFElem& operator/=(const TableGFElem& other) {
        *this = *this / other;
        return *this;
    }
};
}
