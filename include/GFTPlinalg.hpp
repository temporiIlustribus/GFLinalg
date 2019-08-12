#pragma once

#include "GFbase.hpp"

namespace GFlinalg {

    //
    //
    // Two parameter templated versions
    //
    //

    //! Polynomial based GF element class
    /**
     *
     * Basic GF element class. All math operations are done directly in polynomial form.
     *
     *  Operation complexity:
     *  * " + " - O(1)
     *  * " * " - O(n^2)
     *  * " / " - O(log^2(order) + n^2)
     *
     *  Memory complexity: O(1)
     *  
     */
    template <class T, T modPol>
    class  BasicBinPolynomial {
    protected:
        T value;
        constexpr static size_t SZ = op::modPolDegree<T>(modPol);
        constexpr static size_t order = 1 << SZ;

    public:
        //! Default constructor
        explicit BasicBinPolynomial() : value(0) {}
        //! Simple constructor.
        /*!
         * \param val polynomial written in binary form e.g. x^4 + x^2 -> 10100 (20)
         * \param doReduce by default is true. Determines if the stored value should be reduced.
         */
        explicit BasicBinPolynomial(const T& val, bool doReduce = true) : value(val) {
            if (doReduce) reduce();
        }
        //! Iterator based constructor
        /*!
         * \param first iterator to the start of the container
         * \param second iterator to the end of the container
         *
         * Construct polynomial from a container (coefficients are passed in left to right)
         *
         *    Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
         *
         *    Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
         *
         */
        template<typename Iter>
        explicit constexpr BasicBinPolynomial(Iter first, Iter last) : value(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
                ++first;
            }
            reduce();
        }

        //! Returns copy of the polynomial in the contained form
        T val() const noexcept { return value; }
        /*!
         * Returns polynomial in the contained form
         *
         * Note: using this operator to alter contained value violates internal
         *      invariant. To resolve, use "reduce()" method.
         */
        T& val() noexcept { return value; }

        //! Primitive modulus polynomial (modPol)
        static constexpr T modpol = modPol;

        static constexpr T getMod() { return modPol; }

        //! For GF(2^n) returns n
        static size_t gfDegree() { return SZ; }
        //! For GF(2^n) returns 2^n
        static size_t  gfOrder() { return order; }
        //! Calculates the degree of the polynomial
        size_t degree(size_t startPos = 1) {
            return order - op::leadElemPos<T>(value, startPos);
        }
        //! Reduces polynomial by modulus polynomial (modPol)
        T reduce() {
            auto pos = order - SZ;
            // Reduce by modulus
            uint8_t i = 1;
            while (value >= 1U << SZ) {
                if ((value >> (order - i)) & 1)
                    value ^= modPol << (pos - i);
                ++i;
            }
            return value;
        }
        /*!
         * Returns inverse of polynomial
         *
         *    e.g. For polynomial "a" returns "a^(-1)" such that:
         *
         *    a * a^(-1) == a^(-1) * a == 1
         *
         */
        BasicBinPolynomial getInverse() {
            return op::pow<BasicBinPolynomial>(*this, order - 2);
        }
        /*!
         *Inverts polynomial
         *
         *    e.g. For polynomial "a" calculates a^(-1) such that:
         *
         *    a * a^(-1) == a^(-1) * a == 1
         *
         *    and sets a = a^(-1)
         *
         */
        BasicBinPolynomial& invert() {
            *this = op::pow<BasicBinPolynomial>(*this, order - 2);
            return *this;
        }
        //! Cast stored value to T1
        template <class T1>
        explicit operator T1() { return static_cast<T1>(value); }
        //! Returns the stored value
        friend const T& operator + (const BasicBinPolynomial& a) {
            return a.value;
        }

        friend BasicBinPolynomial operator + (const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            return op::polSum<BasicBinPolynomial>(a, b);
        }

        BasicBinPolynomial& operator += (const BasicBinPolynomial& other) {
            *this = op::polSum<BasicBinPolynomial>(*this, other);
            return *this;
        }
        //! Multiplies elements in Galois field as polynomials
        friend BasicBinPolynomial operator * (const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            return op::polMul<BasicBinPolynomial>(a, b);
        }
        //! Multiplies elements in Galois Field as polynomials
        BasicBinPolynomial& operator *= (const BasicBinPolynomial& other) {
            *this = op::polMul<BasicBinPolynomial>(*this, other);
            return *this;
        }
        //! Divides elements in Galois field as polynomials
        friend BasicBinPolynomial operator / (const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            return op::polDiv<BasicBinPolynomial>(a, b);
        }
        //! Divides elements in Galois field as polynomials
        BasicBinPolynomial& operator /= (const BasicBinPolynomial& other) {
            *this = op::polDiv<BasicBinPolynomial>(*this, other);
            return *this;
        }

        template <class T1, T1 modPol1>
        friend bool operator == (const BasicBinPolynomial<T1, modPol1>& lhs, const BasicBinPolynomial<T1, modPol1>& rhs);
        template <class T1, T1 modPol1>
        friend bool operator != (const BasicBinPolynomial<T1, modPol1>& lhs, const BasicBinPolynomial<T1, modPol1>& rhs);

        template <class T1, T1 modPol1>
        friend bool operator < (const BasicBinPolynomial<T1, modPol1>& lhs, const BasicBinPolynomial<T1, modPol1>& rhs);
        template <class T1, T1 modPol1>
        friend bool operator > (const BasicBinPolynomial<T1, modPol1>& lhs, const BasicBinPolynomial<T1, modPol1>& rhs);
        template <class T1, T1 modPol1>
        friend bool operator <= (const BasicBinPolynomial<T1, modPol1>& lhs, const BasicBinPolynomial<T1, modPol1>& rhs);
        template <class T1, T1 modPol1>
        friend bool operator >= (const BasicBinPolynomial<T1, modPol1>& lhs, const BasicBinPolynomial<T1, modPol1>& rhs);
    };
    template <class T, T modPol>
    bool operator == (const BasicBinPolynomial<T, modPol>&  lhs, const BasicBinPolynomial<T, modPol>& rhs) {
        return lhs.value == rhs.value;
    }

    template <class T, T modPol>
    bool operator != (const BasicBinPolynomial<T, modPol>&  lhs, const BasicBinPolynomial<T, modPol>& rhs) {
        return lhs.value != rhs.value;
    }

    template <class T, T modPol>
    bool operator < (const BasicBinPolynomial<T, modPol>&  lhs, const BasicBinPolynomial<T, modPol>& rhs) {
        return lhs.value < rhs.value;
    }

    template <class T, T modPol>
    bool operator > (const BasicBinPolynomial<T, modPol>&  lhs, const BasicBinPolynomial<T, modPol>& rhs) {
        return lhs.value > rhs.value;
    }

    template <class T, T modPol>
    bool operator <= (const BasicBinPolynomial<T, modPol>&  lhs, const BasicBinPolynomial<T, modPol>& rhs) {
        return lhs.value <= rhs.value;
    }

    template <class T, T modPol>
    bool operator >= (const BasicBinPolynomial<T, modPol>&  lhs, const BasicBinPolynomial<T, modPol>& rhs) {
        return lhs.value >= rhs.value;
    }



    //! Look-up table based GF element class
    /**
     * GF element class; uses convertion to power of primitive element Through a LUT internaly
     *
     *  Operation complexity:
     *  *  " + " - O(1)
     *  *  " * " - O(1)
     *  *  " / " - O(1)
     *
     *  Memory complexity: O(2^n)
     *
     */
    template <class T, T modPol>
    class PowBinPolynomial : public BasicBinPolynomial<T, modPol> {
    private:
        using BasicBinPolynomial<T, modPol>::BasicBinPolynomial;
        using BasicBinPolynomial<T, modPol>::order;
        using BasicBinPolynomial<T, modPol>::value;
        

        const static LUTArrPair<T, modPol> alphaToIndex;  /*!<Internal pair of look-up arrays (LUTArrPair)*/

    public:
        //! Default constructor
        explicit PowBinPolynomial(const BasicBinPolynomial<T, modPol>& pol) : BasicBinPolynomial<T, modPol>(pol) {}
        /*!
         *Returns a pair of look-up vectors (LUTArrPair):
         * * indToPol: power of primitive element -> polynomial
         * * polToInd: polynomial -> power of primitive element
         */
        static LUTArrPair<T, modPol> getAlphaToIndex() {
            return alphaToIndex;
        }
        //! Multiplies elements in Galois field using LUTs
        /*!
         *Theory:
         *
         *    a is the primitive element
         *
         *    x * y = z;
         *
         *    x -LUT-> a^u; y-LUT-> a^v;
         *
         *    x * y = a^u * a^v = a^(u+v) -LUT-> z
         *
         */
        PowBinPolynomial operator * (const PowBinPolynomial& other) const {
            if (this->value == 0 || other.value == 0)
                return PowBinPolynomial(0, false);
            return PowBinPolynomial(alphaToIndex.indToPol[alphaToIndex.polToInd[this->value] + alphaToIndex.polToInd[other.value]]);
        }

        //! Multiplies elements in Galois field using LUTs
        /*!
         *Theory:
         *
         *    a is the primitive element
         *
         *    x * y = z;
         *
         *    x -LUT-> a^u; y-LUT-> a^v;
         *
         *    x * y = a^u * a^v = a^(u+v) -LUT-> z
         *
         */
        PowBinPolynomial& operator *= (const PowBinPolynomial& other) {
            this->val() = alphaToIndex.indToPol[alphaToIndex.polToInd[this->val()] +
                alphaToIndex.polToInd[other.value]];
            return (*this);
        }
        //! Divides elements in Galois field using LUTs
        /*!
         *Theory:
         *
         *    a is the primitive element
         *
         *    x / y = z;
         *
         *    x -LUT-> a^u; y-LUT-> a^v;
         *
         *    x / y = a^u * a^v = a^(u-v) -LUT-> z
         *
         */
        PowBinPolynomial operator / (const PowBinPolynomial& other) const {
            if (value == 0)
                return PowBinPolynomial(0);
            if (other.value == 0)
                throw std::out_of_range("Division by zero");
            auto temp(alphaToIndex.polToInd[this->value]);
            if (temp < alphaToIndex.polToInd[other.value])
                temp += order - 1;
            return PowBinPolynomial(alphaToIndex.indToPol[temp - alphaToIndex.polToInd[other.value]]);
        }

        //! Divides elements in Galois field using LUTs
        /*!
         *Theory:
         *
         *    a is the primitive element
         *
         *    x / y = z;
         *
         *    x -LUT-> a^u; y-LUT-> a^v;
         *
         *    x / y = a^u * a^v = a^(u-v) -LUT-> z
         *
         */
        PowBinPolynomial& operator /= (const PowBinPolynomial& other) {
            if (value == 0)
                return PowBinPolynomial(0);
            auto temp(alphaToIndex.polToInd[this->value]);
            if (temp < alphaToIndex.polToInd[other.value])
                temp += order - 1;
            *this->val = alphaToIndex.indToPol[temp - alphaToIndex.polToInd[other.value]];
            return *this;
        }

        PowBinPolynomial operator + (const PowBinPolynomial& other) const {
            return PowBinPolynomial(op::polSum<PowBinPolynomial>(*this, other));
        }

        PowBinPolynomial operator += (const PowBinPolynomial& other) {
            *this = *this + other;
            return *this;
        }
        template <class T1, T1 modPol1>
        friend PowBinPolynomial<T1, modPol1> pow(const PowBinPolynomial<T1, modPol1>& val, size_t power);
        friend bool operator == (const PowBinPolynomial<T, modPol>& a, const PowBinPolynomial<T, modPol>& b) {

            return PowBinPolynomial<T, modPol>(a).value == PowBinPolynomial<T, modPol>(b).value;
        }
        friend bool operator != (const PowBinPolynomial<T, modPol>& a, const PowBinPolynomial<T, modPol>& b) {
            return PowBinPolynomial<T, modPol>(a).value != PowBinPolynomial<T, modPol>(b).value;
        }
    };
    template <class T, T modPol>
    PowBinPolynomial<T, modPol> pow(const PowBinPolynomial<T, modPol>& val, size_t power) {

        return PowBinPolynomial<T, modPol>(
            PowBinPolynomial<T, modPol>::alphaToIndex.indToPol[
                (PowBinPolynomial<T, modPol>::alphaToIndex.polToInd[val.value] * power) % (val.gfOrder() - 1)]);
    }


    //! Table based GF element class
    /**
     * Table based GF element class. All math operations are done using multiplication table and division table
     *
     *   Operation complexity:
     *   * " + " - O(1)
     *   * " * " - O(1)
     *   * " / " - O(1)
     *
     *   Memory complexity: O(4^n)
     *
     */
    template <class T, T modPol>
    class TableBinPolynomial : public BasicBinPolynomial<T, modPol> {
    public:
        using BasicBinPolynomial<T, modPol>::order;
        using BasicBinPolynomial<T, modPol>::value;
        using BasicBinPolynomial<T, modPol>::BasicBinPolynomial;
        using GFtable = std::array<T, order * order>;
    private:
        static const GFtable mulTable; /*<Internal multiplication table */
        static const GFtable divTable; /*<Internal Division table */
    public:
        //! Simple constructor
        explicit TableBinPolynomial(const BasicBinPolynomial<T, modPol>& pol) : BasicBinPolynomial<T, modPol>(pol) {}


        /*!
         *Creates multiplication table (GFtable):
         *
         *    Table[pol1][pol2] = pol1 * pol2;
         *
        */
        static constexpr GFtable makeMulTable() {
            GFtable temp{0};
            // Note: requires checking whether we have traversed through all elements
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicBinPolynomial<T, modPol> a(i);
                    BasicBinPolynomial<T, modPol> b(j);
                    temp[a.val() * order + b.val()] = (a * b).val();
                    temp[b.val() * order + a.val()] = temp[a.val() * order + b.val()];
                }
            }
            return temp;
        }

        /*!
         *Creates division table (GFtable) using multiplication table:
         *
         *    Table[pol1 * pol2][pol2] = pol1;
         *    Table[pol1 * pol2][pol1] = pol2;
         *
         *Note: multiplication table is required to create this table
        */
        static constexpr GFtable makeInvMulTable() {
            GFtable temp{0};
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicBinPolynomial<T, modPol> a(i);
                    BasicBinPolynomial<T, modPol> b(j);
                    temp[mulTable[a.val() * order + b.val()] * order + a.val()] = b.val();
                    temp[mulTable[a.val() * order + b.val()] * order + b.val()] = a.val();
                }
            }
            return temp;
        }
        /*!
         *Creates division table (GFtable) :
         *
         *    Table[pol1][pol2] = pol1 / pol2;
         *
         */
        static constexpr GFtable makeDivTable() {
            GFtable temp;
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = 0; j < order; ++j) {
                    BasicBinPolynomial<T, modPol> a(i);
                    BasicBinPolynomial<T, modPol> b(j);
                    temp[a.value * order + b.value] = (a / b).val();
                }
            }
            divTable = temp;
            return temp;
        }


        TableBinPolynomial operator + (const TableBinPolynomial& other) const {
            return (op::polSum<TableBinPolynomial>(*this, other));
        }

        TableBinPolynomial& operator += (const TableBinPolynomial& other) {
            *this = *this + other;
            return *this;
        }
        //! Multiplies elements in Galois field using multiplication table
        TableBinPolynomial operator * (const TableBinPolynomial& other) const {
            return TableBinPolynomial(mulTable[this->val() * order + other.val()], false);
        }

        //! Multiplies elements in Galois field using multiplication table
        TableBinPolynomial& operator *= (const TableBinPolynomial& other) {
            *this = *this * other;
            return *this;
        }
        //! Divides elements in Galois field using division table
        TableBinPolynomial operator / (const TableBinPolynomial& other) const {
            if (other.value == 0)
                throw std::out_of_range("Division by zero");
            return TableBinPolynomial(divTable[this->val() * order + other.val()], false);
        }
        //! Divides elements in Galois field using division table
        TableBinPolynomial& operator /= (const TableBinPolynomial& other) {
            *this = *this / other;
            return *this;
        }
    };
}