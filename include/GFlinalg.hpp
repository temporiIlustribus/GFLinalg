#pragma once
#include <any>
#include <array>
#include <stdexcept>
#include <iostream>
#include <string>
#include <memory>


namespace GFlinalg {
    template <class T, T modPol>
    class  BasicBinPolynomial;
    template <class T, T modPol>
    class  PowBinPolynomial;
    template <class T, T modPol>
    class  TableBinPolynomial;
    template <class T>
    class BasicGFElem;
    template <class T>
    class PowGFElem;
    template <class T>
    class TableGFElem;
    namespace op {
        //! Get the position of the leading 1 in the polynomial (from the left)
        template<class T>
        uint8_t leadElemPos(const T& pol, uint8_t startPos = 1, size_t order = sizeof(T) << 3) {
            uint8_t pos = 0;
            for (uint8_t i = startPos; i < order + 1; ++i) {
                pos = order - i;
                if ((pol >> pos) & 1) break;
            }
            return order - pos;
        }
        //! Get the degree of the modulus polynomial
        template<class T>
        constexpr uint8_t modPolDegree(const T& modPol) {
            uint8_t pos = 0;
            uint8_t containerSize = sizeof(T) << 3;
            for (uint8_t i = 1; i < containerSize + 1; ++i) {
                pos = containerSize - i;
                if ((modPol >> pos) & 1) break;
            }
            return pos;
        }
        //! Internal polynomial addition (is equivalent to XOR)
        template<class Polynomial>
        Polynomial polSum(const Polynomial& a, const Polynomial& b) {
            Polynomial res(a.val() ^ b.val());
            return res;
        }
        //! Internal multiplication version 1
        template<class Polynomial>
        Polynomial polMulOld(const Polynomial& a, const Polynomial& b) {
            Polynomial res(0);
            for (size_t i = 0; i < a.gfOrder(); ++i) {
                if ((b.value >> i) & 1) {
                    res.val() ^= a.val() << i;
                }
            }
            res.reduce();
            return res;
        }
        //! Internal multiplication version 2
        template<class Polynomial>
        Polynomial polMul(const Polynomial &a, const Polynomial &b) {
            Polynomial res{a};
            res.val() = 0;
            auto av = a.val();
            auto bv = b.val();
            while (bv > 0) {
                if (bv & 1)
                    res.val() ^= av;
                bv >>= 1;
                av <<= 1;
                if (av & a.gfOrder())
                    av ^= a.modpol;
            }
            return res;
        }
        //! Binary power function in Galois field, performing operations in polynomial form
        template <class Polynomial>
        Polynomial pow(Polynomial val, size_t power) {
            Polynomial res(val);
            res.val() = 1;
            while (power) {
                if (power & 1) {
                    res = res * val;
                    --power;
                } else {
                    val *= val;
                    power >>= 1;
                }
            }
            return res;
        }
        //! Internal division
        template<class Polynomial>
        Polynomial polDiv(const Polynomial& a, const Polynomial& b) {
            // Get b^-1: a / b = a * b^-1
            if (b.val() == 0)
                throw std::out_of_range("Division by zero");
            Polynomial invB = op::pow<Polynomial>(b, b.gfOrder() - 2);
            return polMul<Polynomial>(a, invB);
        }
        //! Writes Galois Field element to std::ostream in polynomial form
        /*!
         * * \param out ostream that the elemnt will be written to
         * * \param pol GF element to write to the ostream
         *
         *  GF element will be written to the ostream in polynomial form with higher powers of x
         *  written first.
         *
         *  Example:
         *
         *  Let a be a GF element with coefficients {1,0,1,1};
         *
         *  out << a; - Results in x^3 + x + 1 being written to "out";
         */
        template <class Polynomial>
        std::ostream& operator << (std::ostream& out, Polynomial pol) {
            size_t deg = pol.degree();
            bool flag = false;
            while (pol.val()) {
                if (deg > 0) {
                    out.operator<<('x');
                    if (deg > 1) {
                        out.operator<<('^');
                        out.operator<<(deg);
                    }
                } else {
                    out.operator<<('1');
                }
                flag = true;
                pol.val() ^= (1 << pol.degree());
                deg = pol.degree(deg - 1);
                if (pol.val())
                    out.operator<<('+');
            }
            if (!flag)
                out.operator<<('0');
            return out;
        }

        //! Pair of arrays, used as LUTs
        template <class T, T modPol>
        struct  LUTArrPair {
            static constexpr uint64_t order = (1U << modPolDegree<T>(modPol));
            std::array<T, (order - 1) << 1> indToPol; /*!<Look-up table, converts powers of primitive element to polynomials */
            std::array<size_t, order> polToInd; /*!<Look-up table, converts polynomials to powers of primitive element*/
            /*!
             * Creates a pair of look-up arrays (LUTPair):
             * * indToPol: power of primitive element -> polynomial
             * * polToInd: polynomial -> power of primitive element
             */
            LUTArrPair() : polToInd(), indToPol() {
                BasicBinPolynomial<T, modPol> counter{1};
                BasicBinPolynomial<T, modPol> modifier{2};
                for (size_t i = 0; i < order - 1; ++i) {
                    indToPol[i] = counter.val();
                    polToInd[indToPol[i]] = i;
                    counter *= modifier;
                }
                // This is to avoid % operations in math operators
                for (size_t i = order - 1; i < indToPol.size(); ++i)
                    indToPol[i] = indToPol[i - order + 1];
            }
            LUTArrPair(const std::array<T, (order - 1) * 2>& alph, const std::array<size_t, order>& ind) :
                indToPol(alph), polToInd(ind) {
            }
            LUTArrPair(const std::pair<std::array<T, (order - 1) * 2>, std::array<size_t, order>>& val) :
                indToPol(val.first), polToInd(val.second) {
            }
        };
        template <class T>
        struct LUTVectPair {
            std::vector<T> indToPol;
            std::vector<T> polToInd;
            size_t order;
            LUTVectPair(const T& modPol) : order(1U << modPolDegree<T>(modPol)) {
                polToInd.resize(order); 
                indToPol.resize((order - 1) << 1);
                BasicGFElem<T> counter{1, modPol};
                BasicGFElem<T> modifier{2, modPol};
                for (size_t i = 0; i < order - 1; ++i) {
                    indToPol[i] = counter.val();
                    polToInd[indToPol[i]] = i;
                    counter *= modifier;
                }
                // This is to avoid % operations in math operators
                for (size_t i = order - 1; i < indToPol.size(); ++i) {
                    indToPol[i] = indToPol[i - order + 1];
                }
            }
            LUTVectPair(const LUTVectPair&& lut) : polToInd(std::move(lut.polToInd)),
                indToPol(std::move(lut.indToPol)) {
            }
        };
    }
    using op::pow;
    using op::leadElemPos;
    using op::operator<<;
    using op::LUTArrPair;
    using op::LUTVectPair;
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
         *   Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
         *
         *   Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
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
        //! For GF(2^n) returns n
        static size_t gfDegree() { return SZ; }
        //! For GF(2^n) returns 2^n
        static size_t  gfOrder() { return order; }
        //! Ñalculates the degree of the polynomial
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
         */
        BasicBinPolynomial getInverse() {
            return op::pow<BasicBinPolynomial>(*this, order - 2);
        }
        /*!
         *Inverts polynomial
         *
         *   e.g. For polynomial "a" calculates a^(-1) such that:
         *
         *   a * a^(-1) == a^(-1) * a == 1
         *
         *   and sets a = a^(-1)
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
         *  a is the primitive element
         *
         *  x * y = z;
         *
         *  x -LUT-> a^u; y-LUT-> a^v;
         *
         *  x * y = a^u * a^v = a^(u+v) -LUT-> z
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
         *  a is the primitive element
         *
         *  x * y = z;
         *
         *  x -LUT-> a^u; y-LUT-> a^v;
         *
         *  x * y = a^u * a^v = a^(u+v) -LUT-> z
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
         *  a is the primitive element
         *
         *  x / y = z;
         *
         *  x -LUT-> a^u; y-LUT-> a^v;
         *
         *  x / y = a^u * a^v = a^(u-v) -LUT-> z
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
         *  a is the primitive element
         *
         *  x / y = z;
         *
         *  x -LUT-> a^u; y-LUT-> a^v;
         *
         *  x / y = a^u * a^v = a^(u-v) -LUT-> z
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
         *   Table[pol1][pol2] = pol1 * pol2;
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
         *
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


    //
    //
    // Single parameter templated versions
    //
    //


    //
    // Comments for Single parameter classes are indentical 
    // to their two template parameter analogs
    //

    //! Polynomial based GF element class
    /**
     *
     *Basic GF element class. All math operations are done directly in polynomial form.
     *  Operation complexity:
     *  " + " - O(1)
     *  " * " - O(n^2)
     *  " / " - O(log^2(order) + n^2)
     *
     *  Note: Primitive modulus polynomial is stored in class instances.
     */
    template <class T>
    class BasicGFElem {
    protected:
        T value;
        T modPol;
        size_t SZ;
        size_t order;
    public:
        //! Default constructor
        explicit constexpr BasicGFElem() : value(), modPol(), SZ(0), order(0) {}
        //! Simple constructor.
        /*!
          \param value unsigned integer type template class. Used to store the value internaly
          \param modulus unsigned integer type template class. Modulus polynomial.
          \param doReduce bool parameter, by default is true. Determines reduction of the stored value.
         */
        explicit constexpr BasicGFElem(const T& value, const T& modulus, bool doReduce = true) : value(value), 
                                                    modPol(modulus), SZ(op::modPolDegree<T>(modPol)), order(1 << SZ) {
            if (doReduce) reduce();
        }
        //! Iterator based constructor
        /*!
         *\param first iterator to the start of the container
         *\param second iterator to the end of the container
         *\param modulus modulus polynomial
         *
         *Construct polynomial from a container (coefficients are passed in left to right)
         *
         *  Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
         *
         *  Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
         *
         */
        template<typename Iter>
        explicit constexpr BasicGFElem(Iter first, Iter last, const T& modulus) : value(0), modPol(modulus), SZ(0), order(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
                ++first;
            }
            SZ = op::modPolDegree<T>();
            order = 1 << SZ;
            reduce();
        }
        //! Iterator based constructor
        /*!
         *\param first iterator to the start of the container with GF element
         *\param second iterator to the end of the container with GF element
         *\param firstMod iterator to the start of the container with modulus polynomial
         *\param secondMod iterator to the end of the container with modulus polynomial
         *
         *Construct polynomial from a container (coefficients are passed in left to right)
         *
         *  Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
         *
         *  Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
         *
         */
        template<typename Iter>
        explicit constexpr BasicGFElem(Iter first, Iter last, Iter firstMod, Iter lastMod) : value(0), modPol(0), SZ(0), order(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
                ++first;
            }
            while (firstMod != lastMod) {
                modPol |= (static_cast<T>(*firstMod) & 1);
                modPol <<= 1;
                ++firstMod;
            }
            SZ = op::modPolDegree<T>();
            order = 1 << SZ;
            reduce();
        }
        //! Returns copy of the polynomial in the contained form
        T val() const { return value; }
        /*!
         * Returns polynomial in the contained form
         *
         * Note: using this operator to alter contained value violates internal
         *       invariant. To resolve, use "reduce()" method.
         *
         */
        T& val() { return value; }
        T modpol = modPol;
        //! For GF(2^n) returns n
        size_t gfDegree() const { return SZ; }
        //! For GF(2^n) returns 2^n
        size_t  gfOrder() const { return order; }
        //! Returns the primitive modulus polynomial (modPol) 
        T getMod() const { return modPol; }
        //! Reduces the element by modulus polynomial (modPol)
        T reduce() {
            auto pos = order - SZ;
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
         *  e.g. For polynomial "a" returns "a^(-1)" such that:
         *
         *  a * a^(-1) == a^(-1) * a == 1
         *
         */
        BasicGFElem getInverse() {
            return pow<BasicGFElem>(*this, order - 2);
        }
        /*!
         * Inverts polynomial
         *   e.g. For polynomial "a" calculates a^(-1) such that:
         *
         *   a * a^(-1) == a^(-1) * a == 1
         *
         *   and sets a = a^(-1)
         *
         */
        BasicGFElem& invert() {
            *this = pow<BasicGFElem>(*this, order - 2);
            return *this;
        }
        //! Cast stored value to T1
        template <class T1>
        explicit operator T1() { return static_cast<T1>(value); }

        //! Calculates the degree of the polynomial
        size_t degree(size_t startPos = 1) const {
            return order - leadElemPos<T>(value, startPos);
        }

        friend BasicGFElem operator + (const BasicGFElem& a, const BasicGFElem& b) {
            if (a.modPol != b.modPol)
                throw std::runtime_error("Cannot perform addition for elements of different fields");
            return BasicGFElem(a.val() ^ b.val(), a.modPol);
        }

        BasicGFElem& operator += (const BasicGFElem& other) {
            *this = BasicGFElem(val() ^ other.val(), modPol);
            return *this;
        }
        //! Multiplies elements in Galois field as polynomials
        friend BasicGFElem operator * (const BasicGFElem& a, const BasicGFElem& b) {
            if (a.modPol != b.modPol)
                throw std::runtime_error("Cannot perform multiplication for elements of different fields");
            BasicGFElem res(op::polMul<BasicGFElem>(a, b));
            return res;
        }
        //! Multiplies elements in Galois field as polynomials
        BasicGFElem& operator *= (const BasicGFElem& other) {
            if (modPol != other.modPol)
                throw std::runtime_error("Cannot perform multiplication for elements of different fields");
            *this = *this * other;
            return *this;
        }
        //! Divides elements in Galois field as polynomials
        friend BasicGFElem operator / (const BasicGFElem& a, const BasicGFElem& b) {
            if (a.modPol != b.modPol)
                throw std::runtime_error("Cannot perform division for elements of different fields");
            BasicGFElem res(op::polDiv<BasicGFElem>(a, b));
            return res;
        }
        //! Divides elements in Galois field as polynomials
        BasicGFElem& operator /= (const BasicGFElem& other) {
            *this = *this / other;
            return *this;
        }

        template <class T1>
        friend bool operator == (const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
        template <class T1>
        friend bool operator != (const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);

        template <class T1>
        friend bool operator < (const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
        template <class T1>
        friend bool operator > (const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
        template <class T1>
        friend bool operator <= (const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);
        template <class T1>
        friend bool operator >= (const BasicGFElem<T1>& lhs, const BasicGFElem<T1>& rhs);

    };

    template <class T>
    bool operator < (const BasicGFElem<T>&  lhs, const BasicGFElem<T>& rhs) {
        if (lhs.modPol != rhs.modPol)
            throw std::runtime_error("Cannot compare elements of different fields");
        return (lhs.value < rhs.value);
    }

    template <class T>
    bool operator > (const BasicGFElem<T>& lhs, const BasicGFElem<T>& rhs) {
        if (lhs.modPol != rhs.modPol)
            throw std::runtime_error("Cannot compare elements of different fields");
        return (lhs.value > rhs.value);
    }
    template <class T>
    bool operator <= (const BasicGFElem<T>&  lhs, const BasicGFElem<T>& rhs) {
        if (lhs.modPol != rhs.modPol)
            throw std::runtime_error("Cannot compare elements of different fields");
        return lhs.value <= rhs.value;
    }

    template <class T>
    bool operator >= (const BasicGFElem<T>&  lhs, const BasicGFElem<T>& rhs) {
        if (lhs.modPol != rhs.modPol)
            throw std::runtime_error("Cannot compare elements of different fields");
        return lhs.value >= rhs.value;
    }
    template <class T>
    bool operator == (const BasicGFElem<T>&  lhs, const BasicGFElem<T>& rhs) {
        return (lhs.value == rhs.value) && (lhs.modPol == rhs.modPol);
    }

    template <class T>
    bool operator != (const BasicGFElem<T>&  lhs, const BasicGFElem<T>& rhs) {
        return (lhs.value != rhs.value) || (lhs.modPol != rhs.modPol);
    }


    /*!
     * GF element class; uses convertion to power of primitive element through a LUT internaly
     *
     *   Operation complexity:
     *   " + " - O(1)
     *   " * " - O(1)
     *   " / " - O(1)
     *
     *   Memory complexity: O(2^n)
     *
     *   Note 1: Primitive modulus polynomial is stored in individual instances
     *
     *   Note 2: Uses a LUTPair pointer internally; To make the LUTPair use LUTVectPair<T>(params...)
     *           the LUT is generated for the specific field; use LUT-s with propper parameters.
     *
     */
    template<class T>
    class PowGFElem : public BasicGFElem<T> {
    protected:
        using BasicGFElem<T>::value;
        using BasicGFElem<T>::order;
        using BasicGFElem<T>::modPol;
        LUTVectPair<T> const* alphaToIndex = nullptr;

    public:
        using BasicGFElem<T>::BasicGFElem;
        explicit PowGFElem(const BasicGFElem<T>& pol) : BasicGFElem<T>(pol) {}
        explicit PowGFElem(const T& val, const T& modPol, LUTVectPair<T> const* lut) : BasicGFElem<T>(val, modPol) {
            alphaToIndex = lut;
        }
        explicit PowGFElem(const BasicGFElem<T>& pol, LUTVectPair<T> const* lut) : BasicGFElem<T>(pol) {
            alphaToIndex = lut;
        }

        //! Access the internal look-up table
        LUTVectPair<T> const* lut() const {
            return alphaToIndex;
        }


        PowGFElem operator * (const PowGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform multiplication for elements of different fields");
            if (this->value == 0 || other.value == 0)
                return PowGFElem(0, modPol);
            // x * y = z
            // Let a be the primitive element: x * y -> a^u * a^v = a^(u+v) -> z
            return PowGFElem(alphaToIndex->indToPol[alphaToIndex->polToInd[this->value] +
                             alphaToIndex->polToInd[other.value]], modPol, alphaToIndex);
        }

        PowGFElem& operator *= (const PowGFElem& other) {
            this->val() = alphaToIndex->indToPol[alphaToIndex->polToInd[this->value] +
                alphaToIndex->polToInd[other.value]];
            return *this;
        }

        PowGFElem operator / (const PowGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform division for elements of different fields");
            if (value == 0)
                return PowGFElem(0, modPol);
            if (other.value == 0)
                throw std::out_of_range("Division by zero");
            auto temp = alphaToIndex->polToInd[this->value];
            if (temp < alphaToIndex->polToInd[other.value])
                temp += order - 1;
            // x / y = z
            // Let a be the primitive element: x / y -> a^u / a^v = a^(u-v) -> z
            return PowGFElem(alphaToIndex->indToPol[temp - alphaToIndex->polToInd[other.value]], modPol, alphaToIndex);
        }
        
        PowGFElem& operator /= (const PowGFElem& other) {
            if (value == 0)
                return PowGFElem(0);
            auto temp(alphaToIndex->polToInd[this->value]);
            if (temp < alphaToIndex->polToInd[other.value])
                temp += order - 1;
            *this->val = alphaToIndex->indToPol[temp - alphaToIndex->polToInd[other.value]];
            return *this;
        }

        PowGFElem operator + (const PowGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform addition for elements of different fields");
            return  PowGFElem(this->val() ^ other.val(),modPol, alphaToIndex);
        }

        PowGFElem operator += (const PowGFElem& other) {
            *this = *this + other;
            return *this;
        }
        template<class T1>
        friend PowGFElem<T1> pow(const PowGFElem<T1>& val, size_t power);
    };
    template<class T>
    PowGFElem<T> pow(const PowGFElem<T>& val, size_t power) {
        return PowGFElem<T>(
            val.alphaToIndex->indToPol[
                (val.alphaToIndex->polToInd[val.value] * power) % (val.gfOrder() - 1)], val.modPol, val.alphaToIndex);
    }


    /*!
     * Table based GF element class. All math operations are done using multiplication table and division table
     *
     *   Operation complexity:
     *   * " + " - O(1)
     *   * " * " - O(1)
     *   * " / " - O(1)
     *
     *   Memory complexity: O(4^n)
     *
     *   Note 1: Primitive modulus polynomial is stored in individual instances
     *
     *   Note 2: Uses a std::vector<T> pointer internally to access the tables; To make the tables use makeMulTable(params...)
     *           the LUT is generated for the specific field; use LUT-s with propper parameters.
     *
     *   Note 3: Data access whithin the tables is done using the following method:
     *
     *           Table[a][b] -> Table[a * order +  b];
     *
     *           This documentation will imply Table[a * order +  b] when using Table[a][b] notation.
     *
     */
    template <class T>
    class TableGFElem : public BasicGFElem<T> {
    public:
    private:
        std::vector<T> const* mulTable = nullptr;
        std::vector<T> const* divTable = nullptr;
        using BasicGFElem<T>::order;
        using BasicGFElem<T>::value;
        using BasicGFElem<T>::modPol;
        using BasicGFElem<T>::BasicGFElem;
    public:

        explicit TableGFElem(const BasicGFElem<T>& pol) : BasicGFElem<T>(pol) {}
        explicit TableGFElem(const T& val, const T& modPol,  std::vector<T> const* mulTable_,  std::vector<T> const* divTable_) : BasicGFElem<T>(val, modPol) {
            mulTable = mulTable_;
            divTable = divTable_;
        }
        explicit TableGFElem(const BasicGFElem<T>& pol,  std::vector<T> const* mulTable_,  std::vector<T> const* divTable_) : BasicGFElem<T>(pol) {
            mulTable = mulTable_;
            divTable = divTable_;
        }

        /*!
         * Creates multiplication table (std::vector<T>):
         *
         *    Table[pol1][pol2] = pol1 * pol2;
         *
         */
        static std::vector<T> makeMulTable(const T& modPol) {
            size_t order = 1U << op::modPolDegree<T>(modPol);
            std::vector<T> temp(order * order);
            size_t tempVal = 0;
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicGFElem<T> a(i, modPol);
                    BasicGFElem<T> b(j, modPol);
                    tempVal = a.val() * order + b.val();
                    temp[tempVal] = (a * b).val();
                    temp[b.val() * order + a.val()] = temp[tempVal];
                }
            }
            return temp;
        }

        /*
         * Creates division table ( std::vector<T>) using multiplication table:
         *
         *    Table[pol1 * pol2][pol2] = pol1;
         *
         *    Table[pol1 * pol2][pol1] = pol2;
         *
         *
         * Note: multiplication table is required to create this table
         *
         */
        static std::vector<T> makeInvMulTable(std::vector<T> const* mulTable, const T& modPol) {
            size_t order = 1U << op::modPolDegree<T>(modPol);
            std::vector<T> temp(order * order);
            size_t tempVal = 0;
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicGFElem<T> a(i, modPol);
                    BasicGFElem<T> b(j, modPol);
                    tempVal = a.val() * order + b.val();
                    temp[mulTable->at(tempVal) * order + a.val()] = b.val();
                    temp[mulTable->at(tempVal) * order + b.val()] = a.val();
                }
            }
            return temp;
        }
        /*!
         * Creates division table (std::vector<T>) :
         *
         *    Table[pol1][pol2] = pol1 / pol2;
         *
         */
        static std::vector<T> makeDivTable(const T& modPol) {
            size_t order = 1U << op::modPolDegree<T>(modPol);
            std::vector<T> temp(order * order);
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = 0; j < order; ++j) {
                    BasicGFElem<T> a(i, modPol);
                    BasicGFElem<T> b(j, modPol);
                    temp[a.value * order + b.value] = (a / b).val();
                }
            }
            return temp;
        }
         std::vector<T>* multable() {
            return mulTable;
        }
         std::vector<T>* multable() const {
            return mulTable;
        }
         std::vector<T>* divtable() {
            return divTable;
        }
         std::vector<T>* divtable() const {
            return divTable;
        }

        TableGFElem operator + (const TableGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform addition for elements of different fields");
            return TableGFElem(this->val() ^ other.val(), modPol, mulTable, divTable);
        }

        TableGFElem& operator += (const TableGFElem& other) {
            *this = *this + other;
            return *this;
        }

        TableGFElem operator * (const TableGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform multiplication for elements of different fields");
            return TableGFElem(mulTable->at(this->val() * order + other.val()), modPol, mulTable, divTable);
        }


        TableGFElem& operator *= (const TableGFElem& other) {
            *this = *this * other;
            return *this;
        }

        TableGFElem operator / (const TableGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform division for elements of different fields");
            if (other.value == 0)
                throw std::out_of_range("Division by zero");
            return TableGFElem(divTable->at(this->val() * order + other.val()), modPol, mulTable, divTable);
        }

        TableGFElem& operator /= (const TableGFElem& other) {
            *this = *this / other;
            return *this;
        }
    };

}