#include "pch.h"
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include "GFlinalg.h"
typedef GFlinalg::PowBinPolynomial<uint8_t, 3, 11> powPol;
typedef GFlinalg::BaseBinPolynomial<uint8_t, 3, 11> basePol;
typedef GFlinalg::BasicBinPolynomial<uint8_t, 3, 11> BasicPol;
powPol::arrayPair powPol::alphaToIndex = powPol::makeAlphaToIndex();

int main() {
    powPol::printAlpha();
    size_t n;
    std::cout << "Reduction by modulus test:\n Enter test count:";
    std::cin >> n;
    {
        uint8_t temp;
        std::vector<uint8_t> tests(n);
        std::cout << "\nPowBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp = static_cast<uint32_t>(rand());
            tests[i] = temp;
            std::cout << static_cast<uint32_t>(temp) << " : ";
            powPol pol1(temp);
            std::cout << static_cast<uint32_t>(pol1.getVal()) << "\n";
        }
        std::cout << "\nBaseBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp = tests[i];
            std::cout << static_cast<uint32_t>(temp) << " : ";
            basePol pol1(temp);
            std::cout << static_cast<uint32_t>(pol1.getVal()) << "\n";
        }
        std::cout << "\nBasicBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp = tests[i];
            std::cout << static_cast<uint32_t>(temp) << " : ";
            BasicPol pol1(temp);
            std::cout << static_cast<uint32_t>(pol1.getVal()) << "\n";
        }
    }
    std::cout << "\n\nSum by modulus test:\n Enter test count:";
    std::cin >> n;
    {
        uint8_t temp1, temp2;
        std::vector<std::pair<uint8_t, uint8_t>> tests(n);
        std::cout << "\nPowBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp1 = static_cast<uint32_t>(rand());
            temp2 = static_cast<uint32_t>(rand());
            tests[i] = std::make_pair(temp1, temp2);
            std::cout << static_cast<uint32_t>(temp1) << "; " << static_cast<uint32_t>(temp2) << " : ";
            powPol pol1(temp1);
            powPol pol2(temp2);
            std::cout << static_cast<uint32_t>(pol1.getVal()) << " + " << static_cast<uint32_t>(pol2.getVal()) << " = " << static_cast<uint32_t>((pol1 + pol2).getVal()) << "\n";
        }
        std::cout << "\nBasicBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp1 = tests[i].first;
            temp2 = tests[i].second;
            std::cout << static_cast<uint32_t>(temp1) << "; " << static_cast<uint32_t>(temp2) << " : ";
            BasicPol pol1(temp1);
            BasicPol pol2(temp2);
            std::cout << static_cast<uint32_t>(pol1.getVal()) << " + " << static_cast<uint32_t>(pol2.getVal()) << " = " << static_cast<uint32_t>((pol1 + pol2).getVal()) << "\n";
        }
    }
    std::cout << "\n\nDivision test:\n Enter test count:";
    std::cin >> n;
    {
        uint8_t temp1, temp2;
        std::vector<std::pair<uint8_t, uint8_t>> tests(n);
        std::cout << "\nPowBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp1 = static_cast<uint32_t>(rand());
            temp2 = static_cast<uint32_t>(rand());
            tests[i] = std::make_pair(temp1, temp2);
            std::cout << static_cast<uint32_t>(temp1) << "; " << static_cast<uint32_t>(temp2) << " : ";
            powPol pol1(temp1);
            powPol pol2(temp2);
            if (pol2.getVal() != 0)
                std::cout << static_cast<uint32_t>(pol1.getVal()) << " / " << static_cast<uint32_t>(pol2.getVal()) << " = " << static_cast<uint32_t>((pol1 / pol2).getVal()) << "\n";
        }
        std::cout << "\nBasicBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp1 = tests[i].first;
            temp2 = tests[i].second;
            std::cout << static_cast<uint32_t>(temp1) << "; " << static_cast<uint32_t>(temp2) << " : ";
            BasicPol pol1(temp1);
            BasicPol pol2(temp2);
            if (pol2.getVal() != 0)
                std::cout << static_cast<uint32_t>(pol1.getVal()) << " / " << static_cast<uint32_t>(pol2.getVal()) << " = " << static_cast<uint32_t>((pol1 / pol2).getVal()) << "\n";
        }
    }
    std::cout << "\n\nMultiplication by modulus test:\n Enter test count:";
    std::cin >> n;
    {
        uint8_t temp1, temp2;
        std::vector<std::pair<uint8_t, uint8_t>> tests(n);
        std::cout << "\nPowBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp1 = static_cast<uint32_t>(rand());
            temp2 = static_cast<uint32_t>(rand());
            tests[i] = std::make_pair(temp1, temp2);
            std::cout << static_cast<uint32_t>(temp1) << "; " << static_cast<uint32_t>(temp2) << " : ";
            powPol pol1(temp1);
            powPol pol2(temp2);
            std::cout << static_cast<uint32_t>(pol1.getVal()) << " * " << static_cast<uint32_t>(pol2.getVal()) << " = " << static_cast<uint32_t>((pol1 * pol2).getVal()) << "\n";
        }
        std::cout << "\nBasicBinPolynomial: \n";
        for (size_t i = 0; i < n; ++i) {
            temp1 = tests[i].first;
            temp2 = tests[i].second;
            std::cout << static_cast<uint32_t>(temp1) << "; " << static_cast<uint32_t>(temp2) << " : ";
            BasicPol pol1(temp1);
            BasicPol pol2(temp2);
            std::cout << static_cast<uint32_t>(pol1.getVal()) << " * " << static_cast<uint32_t>(pol2.getVal()) << " = " << static_cast<uint32_t>((pol1 * pol2).getVal()) << "\n";
        }
    }
}
