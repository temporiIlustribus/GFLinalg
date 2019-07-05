#pragma once
#include <iostream>
#include <memory>
#include <array>

namespace GFlinalg {

    template <class T, size_t SZ, T modPol>
    class  BaseBinPolynomial {
    public:
        constexpr static size_t order = 1 << SZ;
    protected:
        T value;
        uint8_t deg;

        // Internal multiplication version 1
        static BaseBinPolynomial polMul(const BaseBinPolynomial& a, const BaseBinPolynomial& b) {
            BaseBinPolynomial res(0);
            uint8_t leadPos = a.degree();
            uint8_t modLeadPos = b.degree();
            for (size_t i = 0; i < order; ++i) {
                if ((b.value >> i) & 1) {
                    res.value ^= a.getVal() << i;
                }
            }
            res.reduce();
            return res;
        }

        // Internal addition
        static BaseBinPolynomial polSum(const BaseBinPolynomial& a, const BaseBinPolynomial& b) {
            BaseBinPolynomial res(a.value ^ b.value);
            return res;
        }

        // Internal division
        static BaseBinPolynomial polDiv(const BaseBinPolynomial& a, const BaseBinPolynomial& b) {
            // Get b^-1: a / b = a * b^-1
            auto invB = galoisPow(b, order - 2);
            return polMul(a, invB);
        }

        // Get the position of the leading 1 in the polynomial
        static uint8_t leadElemPos(const T& pol, size_t startPos=1) {
            uint8_t pos = 0;
            for (uint8_t i = startPos; i < order+1; ++i) {
                if (pol >> (order - i) & 1) {
                    pos = i;
                    break;
                }
            }
            return pos;
        }

        void reduce(const T& modulus) {
            auto pos = leadElemPos(modulus);
            // Reduce by modulus
            uint8_t i = 1;
            while (value >= 1 << (order - pos)) {
                if ((value >> (order - i)) & 1)
                    value ^= modulus << (pos - i);
                ++i;
            }
            deg = order - leadElemPos(value, pos);
        }

    public:

        BaseBinPolynomial() : value(0), deg(0) {}
        BaseBinPolynomial(const T& val) : value(val), deg(0) {
            reduce();
        }
        /*
        Construct polynomial from a container (coefficients are passed in left to right)
            Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
            Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
        */
        template<typename Iter>
        BaseBinPolynomial(Iter first, Iter last): value(0), deg(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
            }
            reduce();
        }

        T getVal() const { return value; }
        T& val() { return value; }
        size_t gfSize() const { return SZ; }
        size_t  gfOrder() const { return order; }
        size_t degree() const { return deg; }
        size_t updateDegree(size_t startPos = 1) { 
            deg = order - leadElemPos(value, startPos);
            return  deg;
        }
        void reduce() { reduce(modPol); }
        template <class T1>
        operator T1() { return static_cast<T1>(value); }

        BaseBinPolynomial& operator ++() {
            ++value;
            reduce();
            return *this;
        }
        BaseBinPolynomial& operator --() {
            --value;
            if (value == T(-1)) {
                reduce();
            }
            return *this;
        }

        template <class T1, size_t SZ1, T1 modPol1>
        friend BaseBinPolynomial<T1, SZ1, modPol1> galoisPow(BaseBinPolynomial<T1, SZ1, modPol1> val, size_t power);

        template <class T1, size_t SZ1, T1 modPol1>
        friend bool operator == (const BaseBinPolynomial<T1, SZ1, modPol1>&  lhs, const BaseBinPolynomial<T1, SZ1, modPol1>& rhs);
        template <class T1, size_t SZ1, T1 modPol1>
        friend bool operator != (const BaseBinPolynomial<T1, SZ1, modPol1>&  lhs, const BaseBinPolynomial<T1, SZ1, modPol1>& rhs);

        template <class T1, size_t SZ1, T1 modPol1>
        friend bool operator < (const BaseBinPolynomial<T1, SZ1, modPol1>&  lhs, const BaseBinPolynomial<T1, SZ1, modPol1>& rhs);
        template <class T1, size_t SZ1, T1 modPol1>
        friend bool operator > (const BaseBinPolynomial<T1, SZ1, modPol1>&  lhs, const BaseBinPolynomial<T1, SZ1, modPol1>& rhs);
        template <class T1, size_t SZ1, T1 modPol1>
        friend bool operator <= (const BaseBinPolynomial<T1, SZ1, modPol1>&  lhs, const BaseBinPolynomial<T1, SZ1, modPol1>& rhs);
        template <class T1, size_t SZ1, T1 modPol1>
        friend bool operator >= (const BaseBinPolynomial<T1, SZ1, modPol1>&  lhs, const BaseBinPolynomial<T1, SZ1, modPol1>& rhs);
    };
    template <class T, size_t SZ, T modPol>
    bool operator == (const BaseBinPolynomial<T, SZ, modPol>&  lhs, const BaseBinPolynomial<T, SZ, modPol>& rhs) {
        return lhs.value == rhs.value;
    }

    template <class T, size_t SZ, T modPol>
    bool operator != (const BaseBinPolynomial<T, SZ, modPol>&  lhs, const BaseBinPolynomial<T, SZ, modPol>& rhs) {
        return lhs.value != rhs.value;
    }

    template <class T, size_t SZ, T modPol>
    bool operator < (const BaseBinPolynomial<T, SZ, modPol>&  lhs, const BaseBinPolynomial<T, SZ, modPol>& rhs) {
        return lhs.value < rhs.value;
    }

    template <class T, size_t SZ, T modPol>
    bool operator > (const BaseBinPolynomial<T, SZ, modPol>&  lhs, const BaseBinPolynomial<T, SZ, modPol>& rhs) {
        return lhs.value > rhs.value;
    }

    template <class T, size_t SZ, T modPol>
    bool operator <= (const BaseBinPolynomial<T, SZ, modPol>&  lhs, const BaseBinPolynomial<T, SZ, modPol>& rhs) {
        return lhs.value <= rhs.value;
    }

    template <class T, size_t SZ, T modPol>
    bool operator >= (const BaseBinPolynomial<T, SZ, modPol>&  lhs, const BaseBinPolynomial<T, SZ, modPol>& rhs) {
        return lhs.value >= rhs.value;
    }

    template <class T, size_t SZ, T modPol>
    std::ostream& operator << (std::ostream& out, BaseBinPolynomial<T, SZ, modPol> pol) {
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

    template <class T, size_t SZ, T modPol>
    BaseBinPolynomial<T, SZ, modPol> galoisPow(BaseBinPolynomial<T, SZ, modPol> val, size_t power) {
        BaseBinPolynomial<T, SZ, modPol> res{1};
        while (power) {
            if (power & 1) {
                res = BaseBinPolynomial<T, SZ, modPol>::polMul(res, val);
                --power;
            } else {
                val = BaseBinPolynomial<T, SZ, modPol>::polMul(val, val);
                power >>= 1;
            }
        }
        return res;
    }

    /*
    Basic GF element class. All math operations are done directly in polynomial form
        Operation complexity:
        " + " - O(1)
        " * " - O(n^2)
        " / " - O(log^2(order) + n^2)
    */
    template <class T, size_t SZ, T modPol>
    class BasicBinPolynomial : public BaseBinPolynomial<T, SZ, modPol> {
    public:
        using BaseBinPolynomial<T, SZ, modPol>::BaseBinPolynomial;
        BasicBinPolynomial(const BaseBinPolynomial<T, SZ, modPol>& pol) : BaseBinPolynomial<T, SZ, modPol>(pol) {}
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

        BasicBinPolynomial operator + (const BasicBinPolynomial& other) const {
            return this->polSum(*this, other);
        }

        BasicBinPolynomial& operator += (const BasicBinPolynomial& other) {
            *this = this->polSum(*this, other);
            return *this;
        }
    };



    /*
    GF element class; uses convertion to power of primitive element internaly
        Operation complexity:
        " + " - O(1)
        " * " - O(1)
        " / " - O(1)
    */
    template <class T, size_t SZ, T modPol>
    class PowBinPolynomial : public BaseBinPolynomial<T, SZ, modPol> {
    public:
        using BaseBinPolynomial<T, SZ, modPol>::value;
        using BaseBinPolynomial<T, SZ, modPol>::order;
        using BaseBinPolynomial<T, SZ, modPol>::polSum;
    private:
        using arrayPair = std::pair<std::array<T, (order - 1) << 1>, std::array<size_t, order>>;
        static arrayPair alphaToIndex;
    public:
        PowBinPolynomial(const BaseBinPolynomial<T, SZ, modPol>& pol) : BaseBinPolynomial<T, SZ, modPol>(pol) {}
        /*
        Creates a pair of vectors (arrayPair):
            first: power of primitive element -> polynomial
            second: polynomial -> power of primitive element
        */
        static constexpr auto makeAlphaToIndex() {
            std::array<T, (order - 1) * 2> alpha;
            std::array<size_t, order> index;
            T counter = 1;
            for (size_t i = 0; i < order - 1; ++i) {
                alpha[i] = BasicBinPolynomial<T, SZ, modPol>(counter).getVal();
                index[alpha[i]] = i;
                counter <<= 1;
            }
            // This is to avoid % operations in math operators
            for (size_t i = order - 1; i < alpha.size(); ++i)
                alpha[i] = alpha[i - order + 1];

            return std::make_pair(alpha, index);
        }

        // Debug function
        static void printAlpha() {
            for (size_t i = 0; i < order - 1; ++i) {
                std::cout << static_cast<size_t>(alphaToIndex.first[i]);
                std::cout << " : " << static_cast<size_t>(alphaToIndex.second[alphaToIndex.first[i]]) << "\n";
            }
        }

        /*
        Returns a pair of vectors (arrayPair):
            first: power of primitive element -> polynomial
            second: polynomial -> power of primitive element
        */
        static arrayPair getAlphaToIndex() {
            return alphaToIndex;
        }

        PowBinPolynomial operator * (const PowBinPolynomial& other) const {
            if (this->value == 0 || other.value == 0)
                return PowBinPolynomial(0);
            return PowBinPolynomial(alphaToIndex.first[alphaToIndex.second[this->value] +
                                    alphaToIndex.second[other.value]]);
        }

        PowBinPolynomial& operator *= (const PowBinPolynomial& other) {
            this->val() = alphaToIndex.first[alphaToIndex.second[this->value] +
                alphaToIndex.second[other.value]];
            return *this;
        }

        PowBinPolynomial operator / (const PowBinPolynomial& other) const {
            if (value == 0)
                return PowBinPolynomial(0);
            auto temp(alphaToIndex.second[this->value]);
            if (temp < alphaToIndex.second[other.value])
                temp += order - 1;
            return PowBinPolynomial(alphaToIndex.first[temp - alphaToIndex.second[other.value]]);
        }

        PowBinPolynomial& operator /= (const PowBinPolynomial& other) {
            if (value == 0)
                return PowBinPolynomial(0);
            auto temp(alphaToIndex.second[this->value]);
            if (temp < alphaToIndex.second[other.value])
                temp += order - 1;
            *this->val = alphaToIndex.first[temp - alphaToIndex.second[other.value]];
            return *this;
        }

        PowBinPolynomial operator + (const PowBinPolynomial& other) const {
            return PowBinPolynomial(polSum(this->value, other.value));
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
    template <class T, size_t SZ, T modPol>
    class TableBinPolynomial : public BasicBinPolynomial<T, SZ, modPol> {
    public:
        using BaseBinPolynomial<T, SZ, modPol>::order;
        using BaseBinPolynomial<T, SZ, modPol>::getVal;
        using BaseBinPolynomial<T, SZ, modPol>::val;
        using BasicBinPolynomial<T, SZ, modPol>::operator+;
        using GFtable = std::array<std::array<T, order>, order>;
    private:
        static GFtable mulTable;
        static GFtable divTable;
    public:

        TableBinPolynomial(const BasicBinPolynomial<T, SZ, modPol>& pol) : BasicBinPolynomial<T, SZ, modPol>(pol) {}
        TableBinPolynomial(const PowBinPolynomial<T, SZ, modPol>& pol) : BasicBinPolynomial<T, SZ, modPol>(pol.getVal()) {}


        /*
        Creates multiplication table (GFtable):
            Table[pol1][pol2] = pol1 * pol2;
        */
        static constexpr GFtable makeMulTable() {
            GFtable temp;
            // Note: requires checking whether we have traversed through all elements
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicBinPolynomial<T, SZ, modPol> a(i);
                    BasicBinPolynomial<T, SZ, modPol> b(j);
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
                    BasicBinPolynomial<T, SZ, modPol> a(i);
                    BasicBinPolynomial<T, SZ, modPol> b(j);
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
                    BasicBinPolynomial<T, SZ, modPol> a(i);
                    BasicBinPolynomial<T, SZ, modPol> b(j);
                    temp[a.value][b.value] = (a / b).getVal();
                }
            }
            divTable = temp;
            return temp;
        }

        using BasicBinPolynomial<T, SZ, modPol>::BasicBinPolynomial;

        TableBinPolynomial operator * (const TableBinPolynomial& other) const {
            return mulTable[this->getVal()][other.getVal()];
        }

        TableBinPolynomial& operator *= (const TableBinPolynomial& other) {
            this->val() = *this * other;
            return *this;
        }

        TableBinPolynomial operator / (const TableBinPolynomial& other) const {
            return divTable[this->getVal()][other.getVal()];
        }

        TableBinPolynomial& operator /= (const TableBinPolynomial& other) {
            this->val() = *this / other;
            return *this;
        }
    };
}
