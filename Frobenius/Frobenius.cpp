#include <stdint.h>
#include <cmath>
#include <random>
#include <numeric>
#include <iostream>
#include <chrono>
#include <fstream>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random.hpp>
#include <boost/math/common_factor_rt.hpp>
#include <boost/integer/extended_euclidean.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "Frobenius.h"



bool CheckIsSquare(mp::uint1024_t number)
{
	if (number == 0 || number == 1)
	{
		return true;
	}

	mp::uint1024_t root = mp::sqrt(number);

	if (root * root == number)
	{
		return true;
	}
	return false;
}

void ABSelect(mp::uint1024_t number, mp::uint1024_t& a, mp::uint1024_t& b, mp::uint1024_t& delta)
{
	boost::random::random_device prng{};
	// распределение в промежутке [1, a]
	boost::random::uniform_int_distribution<mp::uint1024_t> distA(1, a - 1);
	// распределение в промежутке [1, b]
	boost::random::uniform_int_distribution<mp::uint1024_t> distB(1, b - 1);

	// генерация двух рандомных чисел
	a = distA(prng);
	b = distB(prng);
	// подсчет дельты
	delta = a * a - 4 * b;

	// должно удовлетворять 2 условиям:
	// 1) delta - не является полным квадратом
	// 2) НОД(number, 2ab*delta) = 1
	while (!(!CheckIsSquare(delta) && boost::math::gcd(2 * delta * a * b, number) == 1))
	{
		a = distA(prng);
		b = distB(prng);
		delta = a * a - 4 * b;
	}
}

bool CheckIsPrime(mp::uint1024_t number, int iterations)
{
	std::vector<long long> time;
	bool flag = false;
	for (int i = 0; i < iterations; ++i)
	{
		auto begin = std::chrono::high_resolution_clock::now();
		bool rc = FrobeniusAlgorithm(number);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
		time.push_back(duration);
		if (rc)
		{
			flag = true;
		}
	}
	if (flag)
	{
		std::cout << "Number is prime\n";
	}
	else
	{
		std::cout << "Number is not prime\n";
	}
	long long resTime = 0;
	for (int i = 0; i < time.size(); ++i)
	{
		resTime += time[i];
	}
	resTime /= time.size();
	std::ofstream fout("out.dat", std::ios_base::app);
	fout << number << ' ' << resTime << '\n';
	return flag;
}

mp::uint1024_t ExtendedGCD(mp::uint1024_t a, mp::uint1024_t b, mp::int1024_t& x, mp::int1024_t& y)
{
	if (a == 0) {
		x = 0; y = 1;
		return b;
	}
	mp::int1024_t x1 = 0, y1 = 0;
	mp::uint1024_t d = ExtendedGCD(b % a, a, x1, y1);
	x = y1 - (mp::int1024_t)(b / a) * x1;
	y = x1;
	return d;
}

mp::uint1024_t ModulePow(mp::uint1024_t number, mp::uint1024_t degree, mp::uint1024_t module)
{
	mp::uint1024_t res = 1;

	if ((degree >> 24) == 0)
	{
		for (int i = 0; i < degree; ++i)
		{
			res *= (number % module);
			res %= module;
		}
	}
	else
	{
		while (degree != 0)
		{
			if ((degree & 1) != 0)
			{
				res *= number;
				res %= module;
			}

			number *= number;
			number %= module;
			degree >>= 1;
		}
	}
	return res;
}

void CalculateSequence(mp::uint1024_t w1, mp::uint1024_t m, mp::uint1024_t number, mp::uint1024_t& Wm, mp::uint1024_t& Wm1)
{
	Wm = 2;
	Wm1 = w1;

	int log = 0;
	while (m)
	{
		log++;
		m >>= 1;
	}
	log -= 2;
	if (log < 0)
	{
		log = 0;
	}
	mp::uint1024_t mask = 1 << log;
	while (mask <= m)
	{
		mask <<= 1;
	}
	mask >>= 1;
	while (mask)
	{
		if (mask & m != 0)
		{
			Wm = (Wm * Wm1 - w1) % number;
			Wm1 = (Wm1 * Wm1 - 2) % number;
		}
		else
		{
			Wm = (Wm * Wm - 2) % number;
			Wm1 = (Wm * Wm1 - w1) % number;
		}
		mask >>= 1;
	}
}

int CountJacobi(mp::uint1024_t number, mp::uint1024_t delta)
{
	int res = 1;
	while (number)
	{
		// number % 2 == 0
		while (!(number & 1))
		{
			number >>= 1;
			// delta % 8 == 3 || d % 8 == 5
			if ((delta & 7) == 3 || (delta & 7) == 5)
			{
				res = -res;
			}
		}
		std::swap(number, delta);
		// delta % 4 == 3
		if ((delta & 3) == 3 && (delta & 3) == (number & 3))
		{
			res = -res;
		}
		number %= delta;
	}
	 
}


bool FrobeniusAlgorithm(mp::uint1024_t number)
{
	// первичная проверка чисел
	if (number < 3 || number % 2 == 0)
	{
		return false;
	}

	if (CheckIsSquare(number))
	{
		return false;
	}
	

	// выбор чисел А и В
	mp::uint1024_t a = 0, b = 0, delta = 0;
	ABSelect(number, a, b, delta);
	// Алгоритм Евклида и поиск коэффициентов
	mp::int1024_t x = 0, y = 0;
	mp::uint1024_t gcd = ExtendedGCD(a, b, x, y);
	// первый элемент последовательности Wj
	mp::int1024_t w1 = ((mp::int1024_t)(a % number) * (mp::int1024_t)(a % number) * (x % number)) % number - 2;
	if (w1 < 0)
		w1 += number;
	// подсчет символа Якоби для (delta/number)
	mp::uint1024_t m = (number - CountJacobi(delta, number)) >> 1;

	mp::uint1024_t wm = 0, wm1 = 0;
	// подсчет m и (m+1) символа последовательности
	CalculateSequence((mp::uint1024_t)w1, m, number, wm, wm1);


	if ((w1 * wm) != (2 * wm1 % number))
	{
		return false;
	}
	// возведение b в степень (number-1)/2
	b = ModulePow(b, (number - 1) >> 1, number);
	return ((b * wm) % number) == 2;
}

