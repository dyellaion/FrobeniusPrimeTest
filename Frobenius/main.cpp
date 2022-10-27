#include <stdint.h>
#include <cmath>

#include "Frobenius.h"
#include <boost/multiprecision/cpp_int.hpp>


int main()
{
    mp::uint1024_t input;
    int iterations;
    std::cout << "Input a number: ";
    std::cin >> input;
    std::cout << "Input iterations: ";
    std::cin >> iterations;

    bool rc = CheckIsPrime(input, iterations);

    return 0;
}