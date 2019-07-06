#pragma once
#include <iostream>
#include <memory>
#include <array>

namespace GFlinalg {

    /*
   Basic GF element class. All math operations are done directly in polynomial form
       Operation complexity:
       " + " - O(1)
       " * " - O(n^2)
       " / " - O(log^2(order) + n^2)
   */
    template <class T, T modPol>
    class  BasicBinPolynomial {
    protected:
        constexpr static uint8_t modPolDegree() {
            uint8_t pos = 0;
            uint8_t containerSize = sizeof(T) << 3;
            for (uint8_t i = 1; i < containerSize + 1; ++i) {
                pos = containerSize - i;
                if ((modPol >> pos) & 1) break;
            }
            return pos;
        }

        T value;
        constexpr static size_t SZ = modPolDegree();
        uint8_t deg;
        constexpr static size_t order = 1 << SZ;
    protected:
        // Get the position of the leading 1 in the polynomial
        static uint8_t leadElemPos(const T& pol, uint8_t startPos = 1) {
            uint8_t pos = 0;
            for (uint8_t i = startPos; i < order + 1; ++i) {
                pos = order - i;
                if ((pol >> pos) & 1) break;
            }
            return order - pos;
        }
        // Internal multiplication version 1
        static BasicBinPolynomial polMul(const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            BasicBinPolynomial res(0);
            for (size_t i = 0; i < order; ++i) {
                if ((b.value >> i) & 1) {
                    res.value ^= a.getVal() << i;
                }
            }
            res.reduce();
            return res;
        }

        // Internal addition
        static BasicBinPolynomial polSum(const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            BasicBinPolynomial res(a.value ^ b.value);
            return res;
        }

        // Internal division
        static BasicBinPolynomial polDiv(const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            // Get b^-1: a / b = a * b^-1
            auto invB = pow(b, order - 2);
            return polMul(a, invB);
        }

        void reduce(const T& modulus) {
            auto pos = leadElemPos(modulus);
            // Reduce by modulus
            uint8_t i = 1;
            while (value >= 1U << (order - pos)) {
                if ((value >> (order - i)) & 1)
                    value ^= modulus << (pos - i);
                ++i;
            }
            updateDegree();
        }

    public:
        explicit BasicBinPolynomial() : value(0), deg(0) {}
        explicit BasicBinPolynomial(const T& val) : value(val), deg(0) {
            reduce();
        }
        /*
        Construct polynomial from a container (coefficients are passed in left to right)
            Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
            Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
        */
        template<typename Iter>
        explicit BasicBinPolynomial(Iter first, Iter last) : value(0), deg(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
            }
            reduce();
        }

        // Returns copy of the polynomial in the contained form
        T getVal() const { return value; }
        /*
        Returns polynomial in the contained form
            Note: using this operator to alter contained value violates internal
                  invariant. To resolve, use "reduce()" method.

        */
        T& val() { return value; }
        // For GF(2^n) returns n
        static size_t gfSize() { return SZ; }
        // For GF(2^n) returns 2^n
        static size_t  gfOrder() { return order; }
        // Returns the degree of the polynomial
        size_t degree() const { return deg; }
        // Recalculates the degree of the polynomial
        size_t updateDegree(size_t startPos = 1) {
            deg = order - leadElemPos(value, startPos);
            return  deg;
        }
        // Reduces polynomial by modulus polynomial (modPol)
        void reduce() { reduce(modPol); }
        /*
        Returns inverse of polynomial
            e.g. For polynomial "a" returns "a^(-1)" such that:
            a * a^(-1) == a^(-1) * a == 1
        */
        BasicBinPolynomial getInverse() {
            return pow(*this, order - 2);
        }
        /*
        Inverts polynomial
            e.g. For polynomial "a" calculates a^(-1) such that:
            a * a^(-1) == a^(-1) * a == 1
            and sets a = a^(-1)
        */
        BasicBinPolynomial& invert() {
            *this = pow(*this, order - 2);
            return *this;
        }
        template <class T1>
        operator T1() { return static_cast<T1>(value); }

        BasicBinPolynomial operator + (const BasicBinPolynomial& other) const {
            return this->polSum(*this, other);
        }

        BasicBinPolynomial& operator += (const BasicBinPolynomial& other) {
            *this = this->polSum(*this, other);
            return *this;
        }

        BasicBinPolynomial operator * (const BasicBinPolynomial& other) const {
            return this->polMul(*this, other);
        }

        BasicBinPolynomial& operator *= (const BasicBinPolynomial& other) {
            *this = this->polMul(*this, other);
            return *this;
        }

        BasicBinPolynomial operator / (const BasicBinPolynomial& other) const {
            return this->polDiv(*this, other);
        }

        BasicBinPolynomial& operator /= (const BasicBinPolynomial& other) {
            *this = this->polDiv(*this, other);
            return *this;
        }

        template <class T1, T1 modPol1>
        friend BasicBinPolynomial<T1, modPol1> pow(BasicBinPolynomial<T1, modPol1> val, size_t power);

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

    template <class T, T modPol>
    std::ostream& operator << (std::ostream& out, BasicBinPolynomial<T, modPol> pol) {
        size_t deg = pol.degree();
        bool flag = false;
        while (pol.getVal()) {
            if (deg > 0) {
                out << 'x';
                if (deg > 1)
                    out << '^' << deg;
            } else {
                out << '1';
            }
            flag = true;
            pol.val() ^= (1 << pol.degree());
            deg = pol.updateDegree(deg - 1);
            if (pol.val())
                out << '+';
        }
        if (!flag)
            out << '0';
        return out;
    }

    template <class T, T modPol>
    BasicBinPolynomial<T, modPol> pow(BasicBinPolynomial<T, modPol> val, size_t power) {
        BasicBinPolynomial<T, modPol> res{1};
        while (power) {
            if (power & 1) {
                res = BasicBinPolynomial<T, modPol>::polMul(res, val);
                --power;
            } else {
                val = BasicBinPolynomial<T, modPol>::polMul(val, val);
                power >>= 1;
            }
        }
        return res;
    }


    /*
    GF element class; uses convertion to power of primitive element internaly
        Operation complexity:
        " + " - O(1)
        " * " - O(1)
        " / " - O(1)
    */
    template <class T, T modPol>
    class PowBinPolynomial : public BasicBinPolynomial<T, modPol> {
    private:
        using BasicBinPolynomial<T, modPol>::order;
        using BasicBinPolynomial<T, modPol>::value;
        using BasicBinPolynomial<T, modPol>::polSum;
        struct ArrayPair {
            std::array<T, (order - 1) << 1> indToPol;
            std::array<size_t, order> polToInd;
            ArrayPair() : polToInd(), indToPol() {}
            ArrayPair(const std::array<T, (order - 1) * 2>& alph, const std::array<size_t, order>& ind) : indToPol(alph), polToInd(ind) {}
            ArrayPair(const std::pair<std::array<T, (order - 1) * 2>, std::array<size_t, order>>& val) : indToPol(val.first), polToInd(val.second) {}
        };
        static ArrayPair alphaToIndex;
    public:
        using BasicBinPolynomial<T, modPol>::BasicBinPolynomial;
        explicit PowBinPolynomial(const BasicBinPolynomial<T, modPol>& pol) : BasicBinPolynomial<T, modPol>(pol) {}
        /*
        Creates a pair of vectors (ArrayPair):
            indToPol: power of primitive element -> polynomial
            polToInd: polynomial -> power of primitive element
        */
        static constexpr ArrayPair makeAlphaToIndex() {
            ArrayPair temp;
            T counter = 1;
            for (size_t i = 0; i < order - 1; ++i) {
                temp.indToPol[i] = BasicBinPolynomial<T, modPol>(counter).getVal();
                temp.polToInd[temp.indToPol[i]] = i;
                counter <<= 1;
            }
            // This is to avoid % operations in math operators
            for (size_t i = order - 1; i < temp.indToPol.size(); ++i)
                temp.indToPol[i] = temp.indToPol[i - order + 1];

            return temp;
        }

        /*
        Returns a pair of vectors (ArrayPair):
            indToPol: power of primitive element -> polynomial
            polToInd: polynomial -> power of primitive element
        */
        static ArrayPair getAlphaToIndex() {
            return alphaToIndex;
        }

        PowBinPolynomial operator * (const PowBinPolynomial& other) const {
            if (this->value == 0 || other.value == 0)
                return PowBinPolynomial(0);
            return PowBinPolynomial(alphaToIndex.indToPol[alphaToIndex.polToInd[this->value] +
                                    alphaToIndex.polToInd[other.value]]);
        }

        PowBinPolynomial& operator *= (const PowBinPolynomial& other) {
            this->val() = alphaToIndex.indToPol[alphaToIndex.polToInd[this->value] +
                alphaToIndex.polToInd[other.value]];
            return *this;
        }

        PowBinPolynomial operator / (const PowBinPolynomial& other) const {
            if (value == 0)
                return PowBinPolynomial(0);
            auto temp(alphaToIndex.polToInd[this->value]);
            if (temp < alphaToIndex.polToInd[other.value])
                temp += order - 1;
            return PowBinPolynomial(alphaToIndex.indToPol[temp - alphaToIndex.polToInd[other.value]]);
        }

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
            return PowBinPolynomial(polSum(*this, other));
        }

        PowBinPolynomial operator += (const PowBinPolynomial& other) {
            *this = *this + other;
            return *this;
        }
    };

    /*
    Table based GF element class. All math operations are done using multiplication table and division table
        Operation complexity:
        " + " - O(1)
        " * " - O(1)
        " / " - O(1)
    */
    template <class T, T modPol>
    class TableBinPolynomial : public BasicBinPolynomial<T, modPol> {
    public:
        using BasicBinPolynomial<T, modPol>::order;
        using BasicBinPolynomial<T, modPol>::value;
        using BasicBinPolynomial<T, modPol>::polSum;
        using BasicBinPolynomial<T, modPol>::BasicBinPolynomial;
        using GFtable = std::array<std::array<T, order>, order>;
    private:
        static GFtable mulTable;
        static GFtable divTable;
    public:

        explicit TableBinPolynomial(const BasicBinPolynomial<T, modPol>& pol) : BasicBinPolynomial<T, modPol>(pol) {}


        /*
        Creates multiplication table (GFtable):
            Table[pol1][pol2] = pol1 * pol2;
        */
        static constexpr GFtable makeMulTable() {
            GFtable temp;
            // Note: requires checking whether we have traversed through all elements
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicBinPolynomial<T, modPol> a(i);
                    BasicBinPolynomial<T, modPol> b(j);
                    temp[a.val()][b.val()] = (a * b).val();
                    temp[b.val()][a.val()] = temp[a.val()][b.val()];
                }
            }
            //mulTable = temp;
            return temp;
        }

        /*
        Creates division table (GFtable) using multiplication table:
            Table[pol1 * pol2][pol2] = pol1;
            Table[pol1 * pol2][pol1] = pol2;

        Note: multiplication table is required to create this table
        */
        static constexpr GFtable makeInvMulTable() {
            GFtable temp;
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicBinPolynomial<T, modPol> a(i);
                    BasicBinPolynomial<T, modPol> b(j);
                    temp[mulTable[a.val()][b.val()]][a.val()] = b.val();
                    temp[mulTable[a.val()][b.val()]][b.val()] = a.val();
                }
            }
            //divTable = temp;
            return temp;
        }
        /*
        Creates division table (GFtable) :
            Table[pol1][pol2] = pol1 / pol2;
        */
        static constexpr GFtable makeDivTable() {
            GFtable temp;
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = 0; j < order; ++j) {
                    BasicBinPolynomial<T, modPol> a(i);
                    BasicBinPolynomial<T, modPol> b(j);
                    temp[a.value][b.value] = (a / b).getVal();
                }
            }
            divTable = temp;
            return temp;
        }


        TableBinPolynomial operator + (const TableBinPolynomial& other) const {
            return TableBinPolynomial(polSum(*this, other));
        }

        TableBinPolynomial& operator += (const TableBinPolynomial& other) {
            *this = *this + other;
            return *this;
        }

        TableBinPolynomial operator * (const TableBinPolynomial& other) const {
            return TableBinPolynomial(mulTable[this->getVal()][other.getVal()]);
        }

        TableBinPolynomial& operator *= (const TableBinPolynomial& other) {
            *this = *this * other;
            return *this;
        }

        TableBinPolynomial operator / (const TableBinPolynomial& other) const {
            return TableBinPolynomial(divTable[this->getVal()][other.getVal()]);
        }

        TableBinPolynomial& operator /= (const TableBinPolynomial& other) {
            *this = *this / other;
            return *this;
        }
    };
}
