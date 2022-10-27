#pragma once
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;
/**
* �������� �������, �������� ����� �� ��������, �������� ���� ����������
* \param[in] number ����� ��� ��������
* \param[in] iterations ���������� �������� ��������
* \return true, ���� ����� �������, false - �����
*/
bool CheckIsPrime(mp::uint1024_t number, int iterations);

/**
* �������, ���������� �������� ���������� ��� �������� ����� �� ��������
* \param[in] number ����� ��� ��������
* \return true, ���� ����� �������, false - �����
*/
bool FrobeniusAlgorithm(mp::uint1024_t number);

/**
* �������, ������������ 2 ��������� �����, ��������������� �������� ����� ����
* \param[in] number ����� ��� �������� �� ��������
* \param[out] a ������ ��������������� �����
* \param[out] b ������ ��������������� �����
* \param[out] delta = a ^ 2 - 4 * b
*/
void ABSelect(mp::uint1024_t number, mp::uint1024_t& a, mp::uint1024_t& b, mp::uint1024_t& delta);

/**
* �������, ������������ ������������ ������������������ ������������������ Wj
* \param[in] w1 ������ ������� ������������������
* \param[in] m ������ ������������������ 
* \paran[in] number ����� ��� �������� �� �������� (�� ���� ����� ������� ������)
* \param[out] Wm m-� ������� ������������������ Wj
* \param[out] Wm1 (m+1)-� ������� ������������������ Wj
*/
void CalculateSequence(mp::uint1024_t w1, mp::uint1024_t m, mp::uint1024_t number, mp::uint1024_t& Wm, mp::uint1024_t& Wm1);

/**
* �������, �����������, �������� �� ����� ���������
* \param[in] ����������� �����
* \return true, ���� ����� �������� ���������, false - �����
*/
bool CheckIsSquare(mp::uint1024_t number);

/**
* ������� ��� ���������� � ������� ����� �� ������
* \param[in] �����, ���������� � �������
* \param[in] �������, � ������� ���������� �������� �����
* \param[in] ������, �� �������� ����� ������������� ����������
* \return ��������� ���������� � ������� �� ������
*/
mp::uint1024_t ModulePow(mp::uint1024_t number, mp::uint1024_t degree, mp::uint1024_t module);

/**
* ������� ��� ���������� ������� �����
* \param[in] number ������ ��������
* \param[in] delta ������� ��������
* \return �������� ������� �����
*/
int CountJacobi(mp::uint1024_t number, mp::uint1024_t delta);

/**
* ����������� �������� �������
* \param[in] a ������ �����, ��� �������� ��������� ���
* \param[in] b ������ �����, ��� �������� ��������� ���
* \param[out] x ��������� ��� ����� � (��+by=d)
* \param[out] y ��������� ��� ����� b (ax+by=d)
* \return ��� ����� � � b
*/
mp::uint1024_t ExtendedGCD(mp::uint1024_t a, mp::uint1024_t b, mp::int1024_t& x, mp::int1024_t& y);