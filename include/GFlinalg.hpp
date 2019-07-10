#pragma once
#include <any>
#include <array>
#include <stdexcept>
#include <iostream>
#include <memory>


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
        constexpr static size_t order = 1 << SZ;

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
        static BasicBinPolynomial polMulOld(const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            BasicBinPolynomial res(0);
            for (size_t i = 0; i < order; ++i) {
                if ((b.value >> i) & 1) {
                    res.value ^= a.val() << i;
                }
            }
            res.reduce();
            return res;
        }
        // Internal multiplication version 2
        static BasicBinPolynomial polMul(const BasicBinPolynomial &a, const BasicBinPolynomial &b) {
            BasicBinPolynomial res{};
            auto av = a.value;
            auto bv = b.value;
            while (bv > 0) {
                if (bv & 1)
                    res.val() ^= av;
                bv >>= 1;
                av <<= 1;
                if (av & BasicBinPolynomial::gfOrder())
                    av ^= BasicBinPolynomial::modpol;
            }
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
            if (b.value == 0)
                throw std::out_of_range("Division by zero");
            auto invB = pow(b, order - 2);
            return polMul(a, invB);
        }

    public:
        explicit BasicBinPolynomial() : value(0) {}
        explicit BasicBinPolynomial(const T& val, bool doReduce = true) : value(val) {
            if (doReduce) reduce();
        }

        /*
        Construct polynomial from a container (coefficients are passed in left to right)
            Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
            Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
        */
        template<typename Iter>
        explicit constexpr BasicBinPolynomial(Iter first, Iter last) : value(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
            }
            reduce();
        }

        // Returns copy of the polynomial in the contained form
        T val() const noexcept{ return value; }
        /*
        Returns polynomial in the contained form
            Note: using this operator to alter contained value violates internal
                  invariant. To resolve, use "reduce()" method.

        */
        T& val() noexcept { return value; }
        // Primitive modulus polynomial (modPol)
        static constexpr T modpol = modPol;
        // For GF(2^n) returns n
        static size_t gfDegree() { return SZ; }
        // For GF(2^n) returns 2^n
        static size_t  gfOrder() { return order; }
        // Ñalculates the degree of the polynomial
        size_t degree(size_t startPos = 1) {
            return order - leadElemPos(value, startPos);
        }
        // Reduces polynomial by modulus polynomial (modPol)
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
        explicit operator T1() { return static_cast<T1>(value); }

        friend const T& operator + (const BasicBinPolynomial& a) {
            return a.value;
        }

        friend BasicBinPolynomial operator + (const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            return BasicBinPolynomial::polSum(a, b);
        }

        BasicBinPolynomial& operator += (const BasicBinPolynomial& other) {
            *this = this->polSum(*this, other);
            return *this;
        }

        friend BasicBinPolynomial operator * (const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            return BasicBinPolynomial::polMul(a, b);
        }

        BasicBinPolynomial& operator *= (const BasicBinPolynomial& other) {
            *this = this->polMul(*this, other);
            return *this;
        }

        friend BasicBinPolynomial operator / (const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            return BasicBinPolynomial::polDiv(a, b);
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
        while (pol.val()) {
            if (deg > 0) {
                out << 'x';
                if (deg > 1)
                    out << '^' << deg;
            } else {
                out << '1';
            }
            flag = true;
            pol.val() ^= (1 << pol.degree());
            deg = pol.degree(deg - 1);
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

    template <class T, T modPol>
    class FastMultContainer;

    /*
    GF element class; uses convertion to power of primitive element Through a LUT internaly
        Operation complexity:
        " + " - O(1)
        " * " - O(1)
        " / " - O(1)
    */
    template <class T, T modPol>
    class PowBinPolynomial : public BasicBinPolynomial<T, modPol> {
    private:
        using BasicBinPolynomial<T, modPol>::BasicBinPolynomial;
        using BasicBinPolynomial<T, modPol>::order;
        using BasicBinPolynomial<T, modPol>::value;
        using BasicBinPolynomial<T, modPol>::polSum;
        struct  LUTPair {
            std::array<T, (order - 1) << 1> indToPol;
            std::array<size_t, order> polToInd;
            LUTPair() : polToInd(), indToPol() {}
            LUTPair(const std::array<T, (order - 1) * 2>& alph, const std::array<size_t, order>& ind) :
                indToPol(alph), polToInd(ind) {}
            LUTPair(const std::pair<std::array<T, (order - 1) * 2>, std::array<size_t, order>>& val) :
                indToPol(val.first), polToInd(val.second) {}
        };

        const static LUTPair alphaToIndex;
        friend FastMultContainer<T, modPol>;
        
    public:
        
        explicit PowBinPolynomial(const BasicBinPolynomial<T, modPol>& pol) : BasicBinPolynomial<T, modPol>(pol) {}
        PowBinPolynomial(const FastMultContainer<T, modPol>& cont) {
            if (cont.isZero)
                *this = PowBinPolynomial(0, false);
            else
                *this = PowBinPolynomial(alphaToIndex.indToPol[cont.power], false);
        }
        /*
        Creates a pair of look-up vectors (LUTPair):
            indToPol: power of primitive element -> polynomial
            polToInd: polynomial -> power of primitive element
        */
        static constexpr LUTPair makeAlphaToIndex() {
            LUTPair temp;
            T counter = 1;
            for (size_t i = 0; i < order - 1; ++i) {
                temp.indToPol[i] = BasicBinPolynomial<T, modPol>(counter).val();
                temp.polToInd[temp.indToPol[i]] = i;
                counter <<= 1;
            }
            // This is to avoid % operations in math operators
            for (size_t i = order - 1; i < temp.indToPol.size(); ++i)
                temp.indToPol[i] = temp.indToPol[i - order + 1];

            return temp;
        }

        /*
        Returns a pair of look-up vectors (LUTPair):
            indToPol: power of primitive element -> polynomial
            polToInd: polynomial -> power of primitive element
        */
        static LUTPair getAlphaToIndex() {
            return alphaToIndex;
        }

        T val() const noexcept {
            return value;
        }

        T& val() noexcept {
            return value;
        }

        FastMultContainer<T, modPol> operator * (const PowBinPolynomial& other) const {
            if (this->value == 0 || other.value == 0)
                return FastMultContainer<T, modPol>(0, true);
            // x * y = z
            // Let a be the primitive element: x * y -> a^u * a^v = a^(u+v) -> z
            return FastMultContainer<T, modPol>(alphaToIndex.polToInd[this->value] + alphaToIndex.polToInd[other.value]);
        }

        
        FastMultContainer<T, modPol>& operator *= (const PowBinPolynomial& other) {
            this->val() = alphaToIndex.indToPol[alphaToIndex.polToInd[this->value] +
                alphaToIndex.polToInd[other.value]];
            return FastMultContainer<T, modPol>(*this);
        }

        FastMultContainer<T, modPol> operator / (const PowBinPolynomial& other) const {
            if (value == 0)
                return PowBinPolynomial(0);
            if (other.value == 0)
                throw std::out_of_range("Division by zero");
            auto temp(alphaToIndex.polToInd[this->value]);
            if (temp < alphaToIndex.polToInd[other.value])
                temp += order - 1;
            // x / y = z
            // Let a be the primitive element: x / y -> a^u / a^v = a^(u-v) -> z
            return FastMultContainer<T, modPol>(temp - alphaToIndex.polToInd[other.value]);
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
    template <class T, T modPol>
    class FastMultContainer {
    public:
        size_t power;
        bool isZero;
        FastMultContainer() : power(0), isZero(true) {}
        FastMultContainer(const BasicBinPolynomial<T, modPol>& pol) {
            power = PowBinPolynomial<T, modPol>::alphaToIndex.polToInd[pol.val()];
            isZero = false;
            if (pol.val() == 0) 
                isZero = true;
        }
        FastMultContainer(const size_t& pow_, bool Zero = false) : power(pow_), isZero(Zero) {}
        operator PowBinPolynomial<T, modPol>() {
            if (!isZero)
                return  PowBinPolynomial<T, modPol>(PowBinPolynomial<T, modPol>::alphaToIndex.indToPol[power], false);
            else
                return  PowBinPolynomial<T, modPol>(0, false);
        }
        FastMultContainer operator * (const FastMultContainer& other) {
            if (isZero)
                return FastMultContainer(0, true);
            return power + other.power;
        }
       FastMultContainer& operator *= (const FastMultContainer& other) {
           *this = *this * other;
           return *this;
       }
       FastMultContainer operator / (const FastMultContainer& other) {
            if (isZero)
                return  FastMultContainer(0, true);
            if (other.isZero)
                throw std::out_of_range("Cannot divide by zero");
            auto temp = power;
            if (power < other.power)
                temp += PowBinPolynomial<T, modPol>::gfOrder() - 1;
            return power - other.power;
        }
        FastMultContainer& operator /= (const FastMultContainer& other) {
            *this = *this / other;
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
        static const GFtable mulTable;
        static const GFtable divTable;
    public:

        explicit TableBinPolynomial(const BasicBinPolynomial<T, modPol>& pol) : BasicBinPolynomial<T, modPol>(pol) {}


        /*
        Creates multiplication table (GFtable):
            Table[pol1][pol2] = pol1 * pol2;
        */
        static constexpr GFtable makeMulTable() {
            GFtable temp{0};
            // Note: requires checking whether we have traversed through all elements
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicBinPolynomial<T, modPol> a(i);
                    BasicBinPolynomial<T, modPol> b(j);
                    temp[a.val()][b.val()] = (a * b).val();
                    temp[b.val()][a.val()] = temp[a.val()][b.val()];
                }
            }
            return temp;
        }

        /*
        Creates division table (GFtable) using multiplication table:
            Table[pol1 * pol2][pol2] = pol1;
            Table[pol1 * pol2][pol1] = pol2;

        Note: multiplication table is required to create this table
        */
        static constexpr GFtable makeInvMulTable() {
            GFtable temp{0};
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicBinPolynomial<T, modPol> a(i);
                    BasicBinPolynomial<T, modPol> b(j);
                    temp[mulTable[a.val()][b.val()]][a.val()] = b.val();
                    temp[mulTable[a.val()][b.val()]][b.val()] = a.val();
                }
            }
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
                    temp[a.value][b.value] = (a / b).val();
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
            return TableBinPolynomial(mulTable[this->val()][other.val()], false);
        }


        TableBinPolynomial& operator *= (const TableBinPolynomial& other) {
            *this = *this * other;
            return *this;
        }

        TableBinPolynomial operator / (const TableBinPolynomial& other) const {
            if (other.value == 0)
                throw std::out_of_range("Division by zero");
            return TableBinPolynomial(divTable[this->val()][other.val()], false);
        }

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


    /*
    Basic GF element class. All math operations are done directly in polynomial form
        Operation complexity:
        " + " - O(1)
        " * " - O(n^2)
        " / " - O(log^2(order) + n^2)
        Note: primitive modulus polynomial is stored in individual instances
    */
    template <class T>
    class BasicGFElem {
    protected:
        T value;
        T modPol;
        size_t SZ;
        size_t order;
        // Returns the degree of the modulus polynomial
        uint8_t modPolDegree() {
            uint8_t pos = 0;
            uint8_t containerSize = sizeof(T) << 3;
            for (uint8_t i = 1; i < containerSize + 1; ++i) {
                pos = containerSize - i;
                if ((modPol >> pos) & 1) break;
            }
            return pos;
        }
        // Returns the position of the leading 1 in the polynomial from the left
        uint8_t leadElemPos(const T& pol, uint8_t startPos = 1) const {
            uint8_t pos = 0;
            for (uint8_t i = startPos; i < this->order + 1; ++i) {
                pos = order - i;
                if ((pol >> pos) & 1) break;
            }
            return order - pos;
        }
        // Internal multiplication version 1
        BasicGFElem polMulOld(const BasicGFElem& a, const BasicGFElem& b) const {
            BasicGFElem res(0, a.modPol);
            for (size_t i = 0; i < a.order; ++i) {
                if ((b.value >> i) & 1) {
                    res.value ^= a.value << i;
                }
            }
            res.reduce();
            return res;
        }
        // Internal multiplication version 2 
        static BasicGFElem polMul(const BasicGFElem &a, const BasicGFElem &b) {
            BasicGFElem res{0,a.getMod()};
            auto av = a.value;
            auto bv = b.value;
            while (bv > 0) {
                if (bv & 1)
                    res.val() ^= av;
                bv >>= 1;
                av <<= 1;
                if (av & a.gfOrder())
                    av ^= a.getMod();
            }
            return res;
        }
        // Internal addition
        BasicGFElem polSum(const BasicGFElem& a, const BasicGFElem& b) const {
            return BasicGFElem(a.value ^ b.value, a.modPol);
        }

        // Internal Division
        BasicGFElem polDiv(const BasicGFElem& a, const BasicGFElem& b) const {
            // Get b^-1: a / b = a * b^-1
            if (b.value == 0)
                throw std::out_of_range("Division by zero");
            auto invB = pow(b, a.order - 2);
            return polMul(a, invB);
        }


    public:
        explicit constexpr BasicGFElem() : value(), modPol(), SZ(0), order(0) {}
        explicit constexpr BasicGFElem(const T& value, const T& modulus, bool doReduce = true) : value(value), modPol(modulus), SZ(modPolDegree()), order(1 << SZ) {
            if (doReduce) reduce();
        }
        template<typename Iter>
        explicit constexpr BasicGFElem(Iter first, Iter last, const T& modulus) : value(0), modPol(modulus), SZ(0), order(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
            }
            SZ = modPolDegree();
            order = 1 << SZ;
            reduce();
        }

        template<typename Iter>
        explicit constexpr BasicGFElem(Iter first, Iter last, Iter firstMod, Iter lastMod) : value(0), modPol(0), SZ(0), order(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
            }
            while (firstMod != lastMod) {
                modPol |= (static_cast<T>(*firstMod) & 1);
                modPol <<= 1;
            }
            SZ = modPolDegree();
            order = 1 << SZ;
            reduce();
        }
        // Returns copy of the polynomial in the contained form
        T val() const { return value; }
        /*
        Returns polynomial in the contained form
            Note: using this operator to alter contained value violates internal
                  invariant. To resolve, use "reduce()" method.

        */
        T& val() { return value; }
        // For GF(2^n) returns n
        size_t gfDegree() const { return SZ; }
        // For GF(2^n) returns 2^n
        size_t  gfOrder() const { return order; }
        // Returns the primitive modulus polynomial (modPol) 
        T getMod() const { return modPol; }
        // Reduces the element by modulus polynomial (modPol)
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
        /*
        Returns inverse of polynomial
           e.g. For polynomial "a" returns "a^(-1)" such that:
           a * a^(-1) == a^(-1) * a == 1
        */
        BasicGFElem getInverse() {
            return pow(*this, order - 2);
        }
        /*
        Inverts polynomial
           e.g. For polynomial "a" calculates a^(-1) such that:
           a * a^(-1) == a^(-1) * a == 1
           and sets a = a^(-1)
       */
        BasicGFElem& invert() {
            *this = pow(*this, order - 2);
            return *this;
        }

        template <class T1>
        explicit operator T1() { return static_cast<T1>(value); }

        // Calculates the degree of the polynomial
        size_t degree(size_t startPos = 1) const {
            return order - leadElemPos(value, startPos);
        }

        friend BasicGFElem operator + (const BasicGFElem& a, const BasicGFElem& b) {
            if (a.modPol != b.modPol)
                throw std::runtime_error("Cannot perform addition for elements of different fields");
            return a.polSum(a, b);
        }

        BasicGFElem& operator += (const BasicGFElem& other) {
            *this = this->polSum(*this, other);
            return *this;
        }

        friend BasicGFElem operator * (const BasicGFElem& a, const BasicGFElem& b) {
            if (a.modPol != b.modPol)
                throw std::runtime_error("Cannot perform multiplication for elements of different fields");
            return a.polMul(a, b);
        }

        BasicGFElem& operator *= (const BasicGFElem& other) {
            *this = this->polMul(*this, other);
            return *this;
        }

        friend BasicGFElem operator / (const BasicGFElem& a, const BasicGFElem& b) {
            if (a.modPol != b.modPol)
                throw std::runtime_error("Cannot perform division for elements of different fields");
            return a.polDiv(a, b);
        }

        BasicGFElem& operator /= (const BasicGFElem& other) {
            *this = this->polDiv(*this, other);
            return *this;
        }

        template <class T1>
        friend BasicGFElem<T1> pow(BasicGFElem<T1> val, size_t power);

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

    template <class T>
    BasicGFElem<T> pow(BasicGFElem<T> val, size_t power) {
        BasicGFElem<T> res{1, val.getMod()};
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

    template <class T>
    std::ostream& operator << (std::ostream& out, BasicGFElem<T> pol) {
        size_t deg = pol.degree();
        bool flag = false;
        while (pol.val()) {
            if (deg > 0) {
                out << 'x';
                if (deg > 1)
                    out << '^' << deg;
            } else {
                out << '1';
            }
            flag = true;
            pol.val() ^= (1 << pol.degree());
            deg = pol.degree(deg - 1);
            if (pol.val())
                out << '+';
        }
        if (!flag)
            out << '0';
        return out;
    }


    /*
    GF element class; uses convertion to power of primitive element through a LUT internaly
       Operation complexity:
       " + " - O(1)
       " * " - O(1)
       " / " - O(1)
       Note 1: primitive modulus polynomial is stored in individual instances
       Note 2: uses a LUTPair pointer internally; To make the LUTPair use LUTPair<T>::makeLUT(params...)
               the LUT is generated for the specific field; use LUT-s with propper parameters.
   */
    template<class T>
    class PowGFElem : public BasicGFElem<T> {
    public:
        struct LUTPair {
            std::vector<T> indToPol;
            std::vector<T> polToInd;
            LUTPair() : polToInd(), indToPol() {}
            LUTPair(size_t length1, size_t length2) :indToPol(length1), polToInd(length2) {}
            LUTPair(size_t length) :indToPol(length << 1), polToInd(length) {}
            LUTPair(const LUTPair&& lut) : polToInd(std::move(lut.polToInd)),
                indToPol(std::move(lut.indToPol)) {
            }
        };
    protected:
        using BasicGFElem<T>::polSum;
        using BasicGFElem<T>::value;
        using BasicGFElem<T>::order;
        using BasicGFElem<T>::modPol;

        LUTPair* alphaToIndex = nullptr;

    public:
        using BasicGFElem<T>::BasicGFElem;
        explicit PowGFElem(const BasicGFElem<T>& pol) : BasicGFElem<T>(pol) {}
        explicit PowGFElem(const T& val, const T& modPol, LUTPair* lut) : BasicGFElem<T>(val, modPol) {
            alphaToIndex = lut;
        }
        explicit PowGFElem(const BasicGFElem<T>& pol, LUTPair* lut) : BasicGFElem<T>(pol) {
            alphaToIndex = lut;
        }
        /*
            Creates a two-way look-up table for a specified field
            Returns a pointer to generated LUTPair
                indToPol: power of primitive element -> polynomial
                polToInd: polynomial -> power of primitive element
        */
        static LUTPair* makeLUT(size_t order, const T& modPol) {
            LUTPair* temp = new LUTPair(order - 1);
            T counter = 1;
            for (size_t i = 0; i < order - 1; ++i) {
                temp->indToPol[i] = BasicGFElem<T>{counter, modPol}.val();
                temp->polToInd[temp->indToPol[i]] = i;
                counter <<= 1;
            }
            // This is to avoid % operations in math operators
            for (size_t i = order - 1; i < temp->indToPol.size(); ++i) {
                temp->indToPol[i] = temp->indToPol[i - order + 1];
            }
            return temp;
        }

        // Access the internal look-up table
        LUTPair* lut() const {
            return alphaToIndex;
        }
        // Access the internall look-up table
        LUTPair* lut() {
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
                             alphaToIndex->polToInd[other.value]], modPol, false);
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
            return  PowGFElem(polSum(*this, other), alphaToIndex);
        }

        PowGFElem operator += (const PowGFElem& other) {
            *this = *this + other;
            return *this;
        }
        template <class T1>
        friend PowGFElem<T1> pow(const PowGFElem<T1>& val, size_t power);
    };
    template<class T>
    PowGFElem<T> pow(const PowGFElem<T>& val, size_t power) {
        return PowGFElem<T>(
            val.alphaToIndex->indToPol[
                (val.alphaToIndex->polToInd[val.value] * power) % (val.gfOrder() - 1)], val.modPol, val.alphaToIndex);
    }


    /*
    Table based GF element class. All math operations are done using multiplication table and division table
        Operation complexity:
        " + " - O(1)
        " * " - O(1)
        " / " - O(1)
    */
    template <class T>
    class TableGFElem : public BasicGFElem<T> {
    public:
        struct GFtable {
            std::vector<std::vector<T>> data;
            GFtable(size_t rowCount, size_t colCount) : data(rowCount, std::vector<T>(colCount)) {}
            size_t rowCount() { return data.size(); }
            size_t colCount() {
                if (rowCount() != 0)
                    return data[0].size();
                return 0;
            }
        };
    private:
        GFtable* mulTable = nullptr;
        GFtable* divTable = nullptr;
        using BasicGFElem<T>::order;
        using BasicGFElem<T>::value;
        using BasicGFElem<T>::modPol;
        using BasicGFElem<T>::polSum;
        using BasicGFElem<T>::BasicGFElem;
    public:

        explicit TableGFElem(const BasicGFElem<T>& pol) : BasicGFElem<T>(pol) {}
        explicit TableGFElem(const T& val, const T& modPol, GFtable* mulTable_, GFtable* divTable_) : BasicGFElem<T>(val, modPol) {
            mulTable = mulTable_;
            divTable = divTable_;
        }
        explicit TableGFElem(const BasicGFElem<T>& pol, GFtable* mulTable_, GFtable* divTable_) : BasicGFElem<T>(pol) {
            mulTable = mulTable_;
            divTable = divTable_;
        }

        /*
        Creates multiplication table (GFtable):
            Table[pol1][pol2] = pol1 * pol2;
        */
        static GFtable* makeMulTable(size_t order, const T& modPol) {
            GFtable* temp = new GFtable(order, order);
            // Note: requires checking whether we have traversed through all elements
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = i; j < order; ++j) {
                    BasicGFElem<T> a(i, modPol);
                    BasicGFElem<T> b(j, modPol);
                    temp->data[a.val()][b.val()] = (a * b).val();
                    temp->data[b.val()][a.val()] = temp->data[a.val()][b.val()];
                }
            }
            return temp;
        }

        /*
        Creates division table (GFtable) using multiplication table:
            Table[pol1 * pol2][pol2] = pol1;
            Table[pol1 * pol2][pol1] = pol2;

        Note: multiplication table is required to create this table
        */
        static GFtable* makeInvMulTable(GFtable* mulTable, const T& modPol) {
            GFtable* temp = new GFtable(mulTable->rowCount(), mulTable->colCount());
            for (size_t i = 0; i < mulTable->rowCount(); ++i) {
                for (size_t j = i; j < mulTable->colCount(); ++j) {
                    BasicGFElem<T> a(i, modPol);
                    BasicGFElem<T> b(j, modPol);
                    temp->data[mulTable->data[a.val()][b.val()]][a.val()] = b.val();
                    temp->data[mulTable->data[a.val()][b.val()]][b.val()] = a.val();
                }
            }
            return temp;
        }
        /*
        Creates division table (GFtable) :
            Table[pol1][pol2] = pol1 / pol2;
        */
        static GFtable* makeDivTable(size_t order, const T& modPol) {
            GFtable* temp = new GFtable(order, order);
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = 0; j < order; ++j) {
                    BasicGFElem<T> a(i, modPol);
                    BasicGFElem<T> b(j, modPol);
                    temp->data[a.value][b.value] = (a / b).val();
                }
            }
            return temp;
        }
        GFtable* multable() {
            return mulTable;
        }
        GFtable* multable() const {
            return mulTable;
        }
        GFtable* divtable() {
            return divTable;
        }
        GFtable* divtable() const {
            return divTable;
        }

        TableGFElem operator + (const TableGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform addition for elements of different fields");
            return TableGFElem(polSum(*this, other), mulTable, divTable);
        }

        TableGFElem& operator += (const TableGFElem& other) {
            *this = *this + other;
            return *this;
        }

        TableGFElem operator * (const TableGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform multiplication for elements of different fields");
            return TableGFElem(mulTable->data[this->val()][other.val()], modPol, mulTable, divTable);
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
            return TableGFElem(divTable->data[this->val()][other.val()], modPol, mulTable, divTable);
        }

        TableGFElem& operator /= (const TableGFElem& other) {
            *this = *this / other;
            return *this;
        }
    };

}