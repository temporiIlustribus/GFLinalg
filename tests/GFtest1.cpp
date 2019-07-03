#include <iostream>
#include "GFlinalg.h"
#include <random>
#include <string>

uint8_t GFlinalg::binPolyniomial<uint8_t>::modPol = 11;

int main()
{
    size_t n;
    std::cout << "Reduction by modulus test:\n Enter test count:";
    std::cin >> n;
    {
        uint8_t temp;
        for (size_t i = 0; i < n; ++i) {
            temp = static_cast<uint32_t>(rand());
            std::cout << static_cast<uint32_t>(temp) << " : ";
            GFlinalg::binPolyniomial<uint8_t> pol1(temp);
            pol1.reduce();
            std::cout << static_cast<uint32_t>(pol1.getVal()) << "\n";
        }
    }
    std::cout << "\n\nSum by modulus test:\n Enter test count:";
    std::cin >> n;
    {
        uint8_t temp1, temp2;
        for (size_t i = 0; i < n; ++i) {
            temp1 = static_cast<uint32_t>(rand());
            temp2 = static_cast<uint32_t>(rand());
            std::cout << static_cast<uint32_t>(temp1) <<"; "<< static_cast<uint32_t>(temp2) << " : ";
            GFlinalg::binPolyniomial<uint8_t> pol1(temp1);
            pol1.reduce();
            GFlinalg::binPolyniomial<uint8_t> pol2(temp2);
            pol2.reduce();
            std::cout << static_cast<uint32_t>(pol1.getVal()) << " + " << static_cast<uint32_t>(pol2.getVal()) << " = " << static_cast<uint32_t>((pol1 + pol2).getVal()) <<"\n";
        }
    }
    std::cout << "\n\nMultiplication by modulus test:\n Enter test count:";
    std::cin >> n;
    {
        uint8_t temp1, temp2;
        for (size_t i = 0; i < n; ++i) {
            temp1 = static_cast<uint32_t>(rand());
            temp2 = static_cast<uint32_t>(rand());
            std::cout << static_cast<uint32_t>(temp1) << "; " << static_cast<uint32_t>(temp2) << " : ";
            GFlinalg::binPolyniomial<uint8_t> pol1(temp1);
            pol1.reduce();
            GFlinalg::binPolyniomial<uint8_t> pol2(temp2);
            pol2.reduce();
            std::cout << static_cast<uint32_t>(pol1.getVal()) << " * " << static_cast<uint32_t>(pol2.getVal()) << " = " << static_cast<uint32_t>((pol1 * pol2).getVal()) << "\n";
        }
    }
}
