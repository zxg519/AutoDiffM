/*
   一个支持多变量自动微分框架的原理代码，still too young too simple,sometimes naiive.
   by Zhang X.G.
   School of Ins.
   Southeast University
*/
#include<cstring>
#include <complex>
#include <vector>
using namespace std;

using DIFFS = vector<double>;

DIFFS operator*(const DIFFS& v1, const double d2)
{
	DIFFS ret(v1);
	for (auto& data : ret)
	{
		data *= d2;
	}
	return ret;
}
DIFFS operator/(const DIFFS& v1, const double d2)
{
	DIFFS ret(v1);
	for (auto& data : ret)
	{
		data /= d2;
	}
	return ret;
}

DIFFS operator*(const double d2, const DIFFS& v1)
{
	DIFFS ret(v1);
	for (auto& data : ret)
	{
		data *= d2;
	}
	return ret;
}
DIFFS operator+(const DIFFS& v1, const DIFFS& v2)
{
	DIFFS ret(v1);
	DIFFS::iterator it1 = ret.begin();
	DIFFS::const_iterator it2 = v2.begin();
	for (; it1 != ret.end() && it2 != v2.end(); ++it1, ++it2)
	{
		*it1 += *it2;
	}
	return ret;
}
DIFFS operator-(const DIFFS& v1, const DIFFS& v2)
{
	DIFFS ret(v1);
	DIFFS::iterator it1 = ret.begin();
	DIFFS::const_iterator it2 = v2.begin();
	for (; it1 != ret.end() && it2 != v2.end(); ++it1, ++it2)
	{
		*it1 -= *it2;
	}
	return ret;
}


class dual
{
public:
	dual(const double value = 1, const int variable_order = 0, int variable_dim = 1) :_value(value), _diffs(variable_dim, (double)0.0)
	{
		_diffs[variable_order] = 1;
	}
	dual(const double value, const DIFFS& diffs) :_value(value), _diffs(diffs)
	{
	}
	friend dual operator + (const dual& d1, const dual&d2)
	{
		return dual(d1._value + d2._value, d1._diffs + d2._diffs);
	}
	friend dual operator - (const dual& d1, const dual&d2)
	{
		return dual(d1._value + d2._value, d1._diffs - d2._diffs);
	}
	friend dual operator * (const dual& d1, const dual&d2)
	{
		return dual(d1._value * d2._value, d2._value * d1._diffs + d1._value *d2._diffs);
	}

	friend dual operator * (const double d1, const dual&d2)
	{
		return dual(d1 * d2._value, d1 *d2._diffs);
	}
	friend dual operator * (const dual& d2, const double d1)
	{
		return dual(d1 * d2._value, d1 *d2._diffs);
	}


	friend dual operator / (const dual& d1, const dual&d2)
	{
		return dual(d1._value / d2._value,
			d1._diffs / d2._value - (d1._value / (d2._value*d2._value))*d2._diffs);
	}
	friend dual operator / (const double d1, const dual&d2)
	{
		return dual(d1 / d2._value,
			-(d1 / (d2._value*d2._value))*d2._diffs);
	}
	friend dual operator / (const dual&d1, const double d2)
	{
		return dual(d1._value / d2,
			d1._diffs / d2);
	}


	friend dual sin(const dual& d)
	{
		return dual(sin(d._value), cos(d._value)*d._diffs);
	}
	friend dual cos(const dual& d)
	{
		return dual(cos(d._value), -sin(d._value)*d._diffs);
	}
	friend dual exp(const dual& d)
	{
		return dual(exp(d._value), exp(d._value)*d._diffs);
	}
	friend dual pow(const dual&d, const int N)
	{
		double value = 1;
		for (int i = 0; i < N - 1; ++i)
			value *= d._value;
		
		return dual(value*d._value, N*value*d._diffs);
	}
	friend dual operator ^ (const dual&d, const int N)
	{
		return pow(d, N);
	}
	/*
	friend dual log_e(const dual& d)
	{
	// to do
	}
	friend dual log_10(const dual& d)
	{
	// to do
	}
	*/
public:
	double get_value() const
	{
		return _value;
	}
	double get_diff(const int variable_pos) const   // 得到第variable_pos个变量对应的偏微分值
	{
		return _diffs[variable_pos];
	}
protected:
	double               _value;  // 函数值
	DIFFS       _diffs;  // 偏微分列表 存储x1,x2,x3,...,xn的微分值，计算微分和计算函数值同步。
};

int main()
{
	dual x1{ 2, 0, 3 };    // 第0个变量，所以偏微分向量中第0个单元设置为1，其余单元必须设置为0
	dual x2{ 3, 1, 3 };    // 第1个变量，所以偏微分向量中第1个单元设置为1，其余单元必须设置为0
	dual x3{ 4, 2, 3 };    // 第2个变量，所以偏微分向量中第2个单元设置为1，其余单元必须设置为0

	//dual y = 2*x1*x2*x3 + sin(x3)*cos(x3)*pow(x1, 2)/3 + 1/(x2 ^ 4);
	dual y = x2*x3 + 1 / (x1 ^ 5);  //^优先级最低，一定要加括号！
	cout << "f(x1,x2,x3)=" << y.get_value() << endl;
	cout << "dy/dx1=" << y.get_diff(0) << endl;
	cout << "dy/dx2=" << y.get_diff(1) << endl;
	cout << "dy/dx3=" << y.get_diff(2) << endl;
}
