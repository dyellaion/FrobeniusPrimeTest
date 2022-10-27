#pragma once
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;
/**
* Основная функция, проверка числа на простоту, реализуя тест Фробениуса
* \param[in] number Число для проверки
* \param[in] iterations Количество итераций проверки
* \return true, если число простое, false - иначе
*/
bool CheckIsPrime(mp::uint1024_t number, int iterations);

/**
* Функция, реализуюая алгоритм Фробениуса для проверки числа на простоту
* \param[in] number Число для проверки
* \return true, если число простое, false - иначе
*/
bool FrobeniusAlgorithm(mp::uint1024_t number);

/**
* Функция, генерирующая 2 случайных числа, удовлетворяющих условиям теста Люка
* \param[in] number Число для проверки на простоту
* \param[out] a Первое сгенерированное число
* \param[out] b Второе сгенерированное число
* \param[out] delta = a ^ 2 - 4 * b
*/
void ABSelect(mp::uint1024_t number, mp::uint1024_t& a, mp::uint1024_t& b, mp::uint1024_t& delta);

/**
* Функция, генерирующая рекуррентную последовательность последовательность Wj
* \param[in] w1 Первый элемент последовательности
* \param[in] m Индекс последовательности 
* \paran[in] number Число для проверки на простоту (по нему будет браться модуль)
* \param[out] Wm m-й элемент последовательности Wj
* \param[out] Wm1 (m+1)-й элемент последовательности Wj
*/
void CalculateSequence(mp::uint1024_t w1, mp::uint1024_t m, mp::uint1024_t number, mp::uint1024_t& Wm, mp::uint1024_t& Wm1);

/**
* Функция, проверяющая, является ли число квадратом
* \param[in] Проверяемое число
* \return true, если число является квадратом, false - иначе
*/
bool CheckIsSquare(mp::uint1024_t number);

/**
* Функция для возведения в степень числа по модулю
* \param[in] Число, возводимое в степень
* \param[in] Степень, в которую необходимо возвести число
* \param[in] Модуль, по которому будет производиться возведение
* \return Результат возведения в степень по модулю
*/
mp::uint1024_t ModulePow(mp::uint1024_t number, mp::uint1024_t degree, mp::uint1024_t module);

/**
* Функция для вычисления символа Якоби
* \param[in] number Нижний аргумент
* \param[in] delta Верхний аргумент
* \return Значение символа Якоби
*/
int CountJacobi(mp::uint1024_t number, mp::uint1024_t delta);

/**
* Расширенный алгоритм Евклида
* \param[in] a Первое число, для которого считается НОД
* \param[in] b Второе число, для которого считается НОД
* \param[out] x Множитель при числе а (ах+by=d)
* \param[out] y Множитель при числе b (ax+by=d)
* \return НОД чисел а и b
*/
mp::uint1024_t ExtendedGCD(mp::uint1024_t a, mp::uint1024_t b, mp::int1024_t& x, mp::int1024_t& y);