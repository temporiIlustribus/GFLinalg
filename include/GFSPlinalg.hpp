#pragma once

#include "GFbase.hpp"

namespace GFlinalg {
   //
   //
   // Single parameter templated versions
   //
   //


   //
   // Comments for Single parameter classes are indentical 
   // to their two template parameter analogs
   //

   //! Polynomial based single template parameter GF element class
   /**
    *
    *Basic GF element class. All math operations are done directly in polynomial form.
    *  Operation complexity:
    *  * " + " - O(1)
    *  * " * " - O(n^2)
    *  * " / " - O(log^2(order) + n^2)
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
         *    Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
         *
         *    Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
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
         *    Example 1: {1,0,1,0,0} -> x^4 + x^2 (10100)
         *
         *    Example 2: {1,1,1,0,0,1} -> x^5 + x^4 + x^3 + 1 (111001)
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
         *    e.g. For polynomial "a" returns "a^(-1)" such that:
         *
         *    a * a^(-1) == a^(-1) * a == 1
         *
         */
        BasicGFElem getInverse() {
            return pow<BasicGFElem>(*this, order - 2);
        }
        /*!
         * Inverts polynomial
         *
         *    e.g. For polynomial "a" calculates a^(-1) such that:
         *
         *    a * a^(-1) == a^(-1) * a == 1
         *
         *    and sets a = a^(-1)
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

    //! LUT based single template parameter GF element class
    /**
     * GF element class; uses convertion to power of primitive element through a LUT internaly
     *
     *   Operation complexity:
     *   * " + " - O(1)
     *   * " * " - O(1)
     *   * " / " - O(1)
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
        PowGFElem operator * (const PowGFElem& other) const {
            if (other.modPol != this->modPol)
                throw std::runtime_error("Cannot perform multiplication for elements of different fields");
            if (this->value == 0 || other.value == 0)
                return PowGFElem(0, modPol);
            return PowGFElem(alphaToIndex->indToPol[alphaToIndex->polToInd[this->value] +
                             alphaToIndex->polToInd[other.value]], modPol, alphaToIndex);
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
        PowGFElem& operator *= (const PowGFElem& other) {
            this->val() = alphaToIndex->indToPol[alphaToIndex->polToInd[this->value] +
                alphaToIndex->polToInd[other.value]];
            return *this;
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
            return PowGFElem(alphaToIndex->indToPol[temp - alphaToIndex->polToInd[other.value]], modPol, alphaToIndex);
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
            return  PowGFElem(this->val() ^ other.val(), modPol, alphaToIndex);
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

    //! Table based single template parameter GF element class
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
     *   Note 1: Primitive modulus polynomial is stored in individual instances
     *
     *   Note 2: Uses a std::vector<T> pointer internally to access the tables; To make the tables use makeMulTable(params...)
     *           the LUT is generated for the specific field; use LUT-s with propper parameters.
     *
     *   Note 3: Data access whithin the tables is done using the following method:
     *
     *       Table[a][b] -> Table[a * order +  b];
     *
     *   This documentation will imply Table[a * order +  b] when using Table[a][b] notation.
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
        explicit TableGFElem(const T& val, const T& modPol, std::vector<T> const* mulTable_, std::vector<T> const* divTable_) : BasicGFElem<T>(val, modPol) {
            mulTable = mulTable_;
            divTable = divTable_;
        }
        explicit TableGFElem(const BasicGFElem<T>& pol, std::vector<T> const* mulTable_, std::vector<T> const* divTable_) : BasicGFElem<T>(pol) {
            mulTable = mulTable_;
            divTable = divTable_;
        }

        /*!
         * Creates multiplication table (std::vector<T>):
         *
         *     Table[pol1][pol2] = pol1 * pol2;
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

        /*!
         * Creates division table ( std::vector<T>) using multiplication table:
         *
         *     Table[pol1 * pol2][pol2] = pol1;
         *     Table[pol1 * pol2][pol1] = pol2;
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
         *     Table[pol1][pol2] = pol1 / pol2;
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