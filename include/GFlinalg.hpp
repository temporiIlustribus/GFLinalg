#pragma once
#include <any>
#include <array>
#include <stdexcept>
#include <iostream>
#include <memory>


namespace GFlinalg {
    //! Polynomial based GF element class
    /**
     *
     *Basic GF element class. All math operations are done directly in polynomial form.
     *  Operation complexity:
     *  " + " - O(1)
     *  " * " - O(n^2)
     *  " / " - O(log^2(order) + n^2)
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
        //! Default constructor
        explicit BasicBinPolynomial() : value(0) {}
        //! Simple constructor.
        /*!
          \param val unsigned integer type template class. Used to store the value internaly
          \param doReduce bool parameter, by default is true. Determines reduction of the stored value.
         */
        explicit BasicBinPolynomial(const T& val, bool doReduce = true) : value(val) {
            if (doReduce) reduce();
        }
        //! Iterator based constructor
        /*!
          \param first iterator to the start of the container
          \param second iterator to the end of the container

          Construct polynomial from a container (coefficients are passed in left to right)
            Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
            Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
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
          Returns polynomial in the contained form
          Note: using this operator to alter contained value violates internal
                invariant. To resolve, use "reduce()" method.

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
            return order - leadElemPos(value, startPos);
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
        Returns inverse of polynomial
            e.g. For polynomial "a" returns "a^(-1)" such that:
            a * a^(-1) == a^(-1) * a == 1
        */
        BasicBinPolynomial getInverse() {
            return pow(*this, order - 2);
        }
        /*!
        Inverts polynomial
            e.g. For polynomial "a" calculates a^(-1) such that:
            a * a^(-1) == a^(-1) * a == 1
            and sets a = a^(-1)
        */
        BasicBinPolynomial& invert() {
            *this = pow(*this, order - 2);
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
            return BasicBinPolynomial::polSum(a, b);
        }

        BasicBinPolynomial& operator += (const BasicBinPolynomial& other) {
            *this = this->polSum(*this, other);
            return *this;
        }
        //! Multiplies elements in Galois field as polynomials
        friend BasicBinPolynomial operator * (const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            return BasicBinPolynomial::polMul(a, b);
        }
        //! Multiplies elements in Galois Field as polynomials
        BasicBinPolynomial& operator *= (const BasicBinPolynomial& other) {
            *this = this->polMul(*this, other);
            return *this;
        }
        //! Divides elements in Galois field as polynomials
        friend BasicBinPolynomial operator / (const BasicBinPolynomial& a, const BasicBinPolynomial& b) {
            return BasicBinPolynomial::polDiv(a, b);
        }
        //! Divides elements in Galois field as polynomials
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
    //! Binary power in Galois field, performing operations in polynomial form
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
    //! Look-up table based GF element class
    /**
     * GF element class; uses convertion to power of primitive element Through a LUT internaly
     * Operation complexity:
     *   " + " - O(1)
     *   " * " - O(1)
     *   " / " - O(1)
     *
     */
    template <class T, T modPol>
    class PowBinPolynomial : public BasicBinPolynomial<T, modPol> {
    private:
        using BasicBinPolynomial<T, modPol>::BasicBinPolynomial;
        using BasicBinPolynomial<T, modPol>::order;
        using BasicBinPolynomial<T, modPol>::value;
        using BasicBinPolynomial<T, modPol>::polSum;
        //! Pair of arrays, used as LUTs 
        struct  LUTPair {
            std::array<T, (order - 1) << 1> indToPol; /*<Look-up table, converts powers of primitive element to polynomials */
            std::array<size_t, order> polToInd; /*<Look-up table, converts polynomials to powers of primitive element*/
            LUTPair() : polToInd(), indToPol() {}
            LUTPair(const std::array<T, (order - 1) * 2>& alph, const std::array<size_t, order>& ind) :
                indToPol(alph), polToInd(ind) {
            }
            LUTPair(const std::pair<std::array<T, (order - 1) * 2>, std::array<size_t, order>>& val) :
                indToPol(val.first), polToInd(val.second) {
            }
        };

        const static LUTPair alphaToIndex;  /*<Internal pair of look-up arrays (LUTPair)*/
        friend FastMultContainer<T, modPol>;

    public:
        //! Default constructor
        explicit PowBinPolynomial(const BasicBinPolynomial<T, modPol>& pol) : BasicBinPolynomial<T, modPol>(pol) {}
        //! Copy constructor; used for accelecated sequencial multiplication
        PowBinPolynomial(const FastMultContainer<T, modPol>& cont) {
            if (cont.isZero)
                *this = PowBinPolynomial(0, false);
            else
                *this = PowBinPolynomial(alphaToIndex.indToPol[cont.power], false);
        }
        /*!
        Creates a pair of look-up arrays (LUTPair):
            indToPol: power of primitive element -> polynomial
            polToInd: polynomial -> power of primitive element
        */
        static constexpr LUTPair makeAlphaToIndex() {
            LUTPair temp;
            BasicBinPolynomial<T, modpol> counter{1};
            for (size_t i = 0; i < order - 1; ++i) {
                temp.indToPol[i] = counter.val();
                temp.polToInd[temp.indToPol[i]] = i;
                counter *= BasicBinPolynomial<T, modPol>(2);
            }
            // This is to avoid % operations in math operators
            for (size_t i = order - 1; i < temp.indToPol.size(); ++i)
                temp.indToPol[i] = temp.indToPol[i - order + 1];

            return temp;
        }

        /*!
        Returns a pair of look-up vectors (LUTPair):
            indToPol: power of primitive element -> polynomial
            polToInd: polynomial -> power of primitive element
        */
        static LUTPair getAlphaToIndex() {
            return alphaToIndex;
        }
        //! Returns copy of the polynomial in the contained form
        T val() const noexcept {
            return value;
        }
        /*!
          Returns polynomial in the contained form
          Note: using this operator to alter contained value violates internal
                invariant. To resolve, use "reduce()" method.

        */
        T& val() noexcept {
            return value;
        }
        //! Multiplies elements in Galois field using LUTs
        /*!
          Theory:
            a is the primitive element
            x * y = z;
            x -LUT-> a^u; y-LUT-> a^v;
            x * y = a^u * a^v = a^(u+v) -LUT-> z

        */
        FastMultContainer<T, modPol> operator * (const PowBinPolynomial& other) const {
            if (this->value == 0 || other.value == 0)
                return FastMultContainer<T, modPol>(0, true);
            return FastMultContainer<T, modPol>(alphaToIndex.polToInd[this->value] + alphaToIndex.polToInd[other.value]);
        }

        //! Multiplies elements in Galois field using LUTs
        /*!
          Theory:
            a is the primitive element
            x * y = z;
            x -LUT-> a^u; y-LUT-> a^v;
            x * y = a^u * a^v = a^(u+v) -LUT-> z

        */
        PowBinPolynomial& operator *= (const PowBinPolynomial& other) {
            this->val() = alphaToIndex.indToPol[alphaToIndex.polToInd[this->value] +
                alphaToIndex.polToInd[other.value]];
            return (*this);
        }
        //! Divides elements in Galois field using LUTs
        /*!
          Theory:
            a is the primitive element
            x / y = z;
            x -LUT-> a^u; y-LUT-> a^v;
            x / y = a^u * a^v = a^(u-v) -LUT-> z

        */
        FastMultContainer<T, modPol> operator / (const PowBinPolynomial& other) const {
            if (value == 0)
                return PowBinPolynomial(0);
            if (other.value == 0)
                throw std::out_of_range("Division by zero");
            auto temp(alphaToIndex.polToInd[this->value]);
            if (temp < alphaToIndex.polToInd[other.value])
                temp += order - 1;
            return FastMultContainer<T, modPol>(temp - alphaToIndex.polToInd[other.value]);
        }

        //! Divides elements in Galois field using LUTs
        /*!
          Theory:
            a is the primitive element
            x / y = z;
            x -LUT-> a^u; y-LUT-> a^v;
            x * y = a^u / a^v = a^(u-v) -LUT-> z

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
            //std::cout << "Converting to PowBinPol " << power << '\n';
            if (!isZero)
                return  PowBinPolynomial<T, modPol>(PowBinPolynomial<T, modPol>::alphaToIndex.indToPol[power]);
            else
                return  PowBinPolynomial<T, modPol>(0);
        }
        FastMultContainer operator * (const FastMultContainer& other) {
            if (isZero || other.isZero)
                return FastMultContainer(0, true);
            return (power + other.power) % (PowBinPolynomial<T, modPol>::order - 1);
        }
        FastMultContainer operator * (const PowBinPolynomial<T, modPol>& other) {
            if (isZero || other.val() == 0)
                return FastMultContainer(0, true);

            return (power + PowBinPolynomial<T, modPol>::alphaToIndex.polToInd[other.val()]) %
                (PowBinPolynomial<T, modPol>::gfOrder() - 1);
        }

        FastMultContainer& operator *= (const FastMultContainer& other) {
            *this = *this * other;
            return *this;
        }
        FastMultContainer& operator *= (const PowBinPolynomial<T, modPol>& other) {
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
            return temp - other.power;
        }
        FastMultContainer operator / (const PowBinPolynomial<T, modPol>& other) {
            if (isZero)
                return  FastMultContainer(0, true);
            if (other.isZero)
                throw std::out_of_range("Cannot divide by zero");
            auto otherCopy = FastMultContainer(other);
            if (otherCopy.power < power)
                return power - otherCopy.power;
            return power + PowBinPolynomial<T, modPol>::gfOrder() - 1 - otherCopy.power;
        }
        FastMultContainer& operator /= (const FastMultContainer& other) {
            *this = *this / other;
            return *this;
        }
        FastMultContainer& operator /= (const PowBinPolynomial<T, modPol>& other) {
            *this = *this / other;
            return *this;
        }
    };
    //! Table based GF element class
    /**
     * Table based GF element class. All math operations are done using multiplication table and division table
     *   Operation complexity:
     *   " + " - O(1)
     *   " * " - O(1)
     *   " / " - O(1)
     *
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
        static const GFtable mulTable; /*<Internal multiplication table */
        static const GFtable divTable; /*<Internal Division table */
    public:
        //! Simple constructor
        explicit TableBinPolynomial(const BasicBinPolynomial<T, modPol>& pol) : BasicBinPolynomial<T, modPol>(pol) {}


        /*!
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

        /*!
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
        /*!
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
        //! Multiplies elements in Galois field using multiplication table
        TableBinPolynomial operator * (const TableBinPolynomial& other) const {
            return TableBinPolynomial(mulTable[this->val()][other.val()], false);
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
            return TableBinPolynomial(divTable[this->val()][other.val()], false);
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
     *  Note: Primitive element modulus polynomial is stored in class instances.
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
        //! Default constructor
        explicit constexpr BasicGFElem() : value(), modPol(), SZ(0), order(0) {}
        //! Simple constructor.
        /*!
          \param value unsigned integer type template class. Used to store the value internaly
          \param modulus unsigned integer type template class. Modulus polynomial.
          \param doReduce bool parameter, by default is true. Determines reduction of the stored value.
         */
        explicit constexpr BasicGFElem(const T& value, const T& modulus, bool doReduce = true) : value(value), modPol(modulus), SZ(modPolDegree()), order(1 << SZ) {
            if (doReduce) reduce();
        }
        //! Iterator based constructor
        /*!
          \param first iterator to the start of the container
          \param second iterator to the end of the container
          \param modulus modulus polynomial
          Construct polynomial from a container (coefficients are passed in left to right)
            Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
            Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
        */
        template<typename Iter>
        explicit constexpr BasicGFElem(Iter first, Iter last, const T& modulus) : value(0), modPol(modulus), SZ(0), order(0) {
            while (first != last) {
                value |= (static_cast<T>(*first) & 1);
                value <<= 1;
                ++first;
            }
            SZ = modPolDegree();
            order = 1 << SZ;
            reduce();
        }
        //! Iterator based constructor
        /*!
          \param first iterator to the start of the container with GF element
          \param second iterator to the end of the container with GF element
          \param firstMod iterator to the start of the container with modulus polynomial
          \param secondMod iterator to the end of the container with modulus polynomial
          Construct polynomial from a container (coefficients are passed in left to right)
            Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
            Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
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
            SZ = modPolDegree();
            order = 1 << SZ;
            reduce();
        }
        //! Returns copy of the polynomial in the contained form
        T val() const { return value; }
        /*!
          Returns polynomial in the contained form
          Note: using this operator to alter contained value violates internal
                invariant. To resolve, use "reduce()" method.

        */
        T& val() { return value; }
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
        Returns inverse of polynomial
            e.g. For polynomial "a" returns "a^(-1)" such that:
            a * a^(-1) == a^(-1) * a == 1
        */
        BasicGFElem getInverse() {
            return pow(*this, order - 2);
        }
        /*!
        Inverts polynomial
            e.g. For polynomial "a" calculates a^(-1) such that:
            a * a^(-1) == a^(-1) * a == 1
            and sets a = a^(-1)
        */
        BasicGFElem& invert() {
            *this = pow(*this, order - 2);
            return *this;
        }
        //! Cast stored value to T1
        template <class T1>
        explicit operator T1() { return static_cast<T1>(value); }

        //! Calculates the degree of the polynomial
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
        //! Multiplies elements in Galois field as polynomials
        friend BasicGFElem operator * (const BasicGFElem& a, const BasicGFElem& b) {
            if (a.modPol != b.modPol)
                throw std::runtime_error("Cannot perform multiplication for elements of different fields");
            return a.polMul(a, b);
        }
        //! Multiplies elements in Galois field as polynomials
        BasicGFElem& operator *= (const BasicGFElem& other) {
            *this = this->polMul(*this, other);
            return *this;
        }
        //! Divides elements in Galois field as polynomials
        friend BasicGFElem operator / (const BasicGFElem& a, const BasicGFElem& b) {
            if (a.modPol != b.modPol)
                throw std::runtime_error("Cannot perform division for elements of different fields");
            return a.polDiv(a, b);
        }
        //! Divides elements in Galois field as polynomials
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

        LUTPair const* alphaToIndex = nullptr;

    public:
        using BasicGFElem<T>::BasicGFElem;
        explicit PowGFElem(const BasicGFElem<T>& pol) : BasicGFElem<T>(pol) {}
        explicit PowGFElem(const T& val, const T& modPol, LUTPair const* lut) : BasicGFElem<T>(val, modPol) {
            alphaToIndex = lut;
        }
        explicit PowGFElem(const BasicGFElem<T>& pol, LUTPair const* lut) : BasicGFElem<T>(pol) {
            alphaToIndex = lut;
        }
        /*
            Creates a two-way look-up table for a specified field
            Returns a pointer to generated LUTPair
                indToPol: power of primitive element -> polynomial
                polToInd: polynomial -> power of primitive element
        */
        static LUTPair* makeLUT(size_t order, const T& modPol) {
            LUTPair* temp = new LUTPair(2 * (order-1), order);
            BasicGFElem<T> counter{1, modPol};
            for (size_t i = 0; i < order - 1; ++i) {
                temp->indToPol[i] = counter.val();
                temp->polToInd[temp->indToPol[i]] = i;
                counter *= BasicGFElem < T>{2, modPol};
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