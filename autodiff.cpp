/*
   一个支持多变量自动微分框架的原理代码， still too young too simple, sometimes naiive.
   by Zhang X.G.
   School of Ins.
   Southeast University
*/
#include<iostream>
#include<cstring>
#include <complex>
#include <vector>
#include <map>
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
		return dual(d1._value - d2._value, d1._diffs - d2._diffs);
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
	friend dual pow(const dual&x, const int N1)
	{

		double value = 1;
		if (N1 == 0)
		{
			return dual(1, vector<double>(x.diff_size(), 0));
		}
		else if (N1 > 0)
		{
			for (int i = 0; i < N1 - 1; ++i)
				value *= x._value;

			return dual(value*x._value, N1*value*x._diffs);
		}
		else //if (N < 0)
		{
			int N = -N1;

			for (int i = 0; i < N; ++i)
				value *= x._value;

			return dual(value, -N / (value*x._value) * x._diffs);
		}
	}
	friend dual pow(const dual&x1, const dual& x2) // x1^x2
	{
		return dual(exp(x2._value*log(x1._value)),
			exp(x2._value*log(x1._value))*(log(x1._value)*x2._diffs + x2._value/x1._value*x1._diffs)
			);
	}

	friend dual sqrt_N(const dual&x, const int N)
	{
		//x ^ (1 / N) = e ^ (1 / N*log(x) )
		return dual(exp(log(x._value) / N),
			exp(log(x._value) / N) / N / x._value * x._diffs
			        );
	}
	//==-----------------------------------------------------------------------------------------
	friend dual operator ^ (const dual&x1, const dual& x2) // x1^x2
	{
		return pow(x1, x2);
	}
	friend dual operator ^ (const dual&x, const int N)
	{
		return pow(x, N);
	}
	friend dual operator ^ (const double down, const dual& x)
	{
		// down^x= e ^ (x ln down)
		// to do list
		return dual(exp(x._value*log(down)), 
			        exp(x._value*log(down))*log(down)*x._diffs);
	}
	friend dual operator ^ (const dual& x, const double exponent)  //x^p
	{
		// down^x= e ^ (x ln down)
		// to do list
		return dual(exp(exponent*log(x._value)),
			exp(exponent*log(x._value))*exponent/(x._value)*x._diffs);
	}
	
	friend dual log_e(const dual& d)
	{
		return dual(log(d._value),  
		            1.0/d._value * d._diffs);
	}
	friend dual log_10(const dual& x)
	{
		// log_10(x) = ln(x) / ln(10)
		return dual(log(x._value)/log(10.0),
			1.0 / x._value /log(10.0)* x._diffs);

	}
	
public:
	double get_value() const
	{
		return _value;
	}
	const int diff_size() const
	{
		return _diffs.size();
	}
	
	double get_diff(const int variable_pos) const   // 得到第variable_pos个变量对应的偏微分值
	{
		return _diffs[variable_pos];
	}
	double get_diff(const dual& variable) const     // 求相对于某个自变量的偏微分值
	{
		int pos = 0;
		for (auto& data : variable._diffs)
		{
			if (data == 0.0)
			{
				pos++;
			}
			else if (data == 1.0)
			{
				break;
			}
			else
			{
				cout << "ERROR state in line:"<<__LINE__<< endl;
				exit(0);
			}		
		}
		return _diffs[pos];
	}
protected:
	double      _value;  // 函数值
	DIFFS       _diffs;  // 偏微分列表 存储x1,x2,x3,...,xn的微分值，计算微分和计算函数值同步。
};

// To do, regiteration mechanism, automatica set the diff vectors;
//   and automatic setting the patial diff position in the vector for each variable
/*
class argument_register
{
public:
	void begin_regist()
	{
	}
	void regist(const string& argument_name, 
	            const double argument_value) //注册变量名字及计算偏微分时该变量的取值
	{
	}
	void end_regist()
	{
	}

public:
	dual& get_argument(const string& name)
	{

	}
	dual& operator[](const string& name)
	{
	
	}
	int  diff_pos(const string& name)
	{

	}
protected:
	vector<dual> duals;
	map<const string, dual&> duals_mapping;
};
*/

int main()
{
	dual x1{ 2, 0, 3 };    // 第0个变量，所以偏微分向量中第0个单元设置为1，其余单元必须设置为0
	dual x2{ 3, 1, 3 };    // 第1个变量，所以偏微分向量中第1个单元设置为1，其余单元必须设置为0
	dual x3{ 4, 2, 3 };    // 第2个变量，所以偏微分向量中第2个单元设置为1，其余单元必须设置为0

	//dual y = 2*x1*x2*x3 + sin(x3)*cos(x3)*pow(x1, 2)/3 + 1/(x2 ^ 4);
	dual y = x2*x3 + 1 / (x1 ^ 5);  //^优先级低，一定要加括号！
	cout << "f(x1,x2,x3)=" << y.get_value() << endl;
	cout << "dy/dx1=" << y.get_diff(0) << endl;
	cout << "dy/dx2=" << y.get_diff(1) << endl;
	cout << "dy/dx3=" << y.get_diff(2) << endl;

	cout << "--------------------------------------" << endl;
	dual x5{ 2, 0, 1}; //单变量微分
	dual y1 = sqrt_N(x5,2);
	cout << y1.get_value() << endl;
	cout << "dy/dx1="<<y1.get_diff(0) << endl;


	cout << "--------------------------------------" << endl;
	dual x8{ 2, 0, 1}; //单变量微分
	dual y3 = sqrt_N(x5, 2);
	cout << y3.get_value() << endl;
	cout << "dy/dx1=" << y3.get_diff(0) << endl;

	cout << "--------------------------------------" << endl;
	dual x9{ 2, 0, 2 }; //双变量微分
	dual x10{ 2, 1, 2 }; //双变量微分
	dual y4 = x9^x10;
	cout << y4.get_value() << endl;
	cout << "dy/dx1=" << y4.get_diff(0) << endl;
	cout << "dy/dx2=" << y4.get_diff(1) << endl;

	/*
	    // to do list

		argument_register reg;
		reg.begin_regist();
			reg.register("x1",2);  // x1, partial diff on x1 = 2
			reg.register("x2",3);  // x2, partial diff on x2 = 3
		reg.end_register();

		auto& x1 = reg["x1"];
		auto& x2 = reg["x2"]

		auto y = x1*x2 + sin(x2)*log(x1);

		cout<<"dy/dx1 = "<<y.get_diff(x1)<<endl;
		cout<<"dy/dx1 = "<<y.get_diff(x2)<<endl;
	
	*/

}


