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
#include <fstream>
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

// 不鼓励使用 +=, -=
DIFFS& operator+=(DIFFS& v1, const DIFFS& v2)
{
	DIFFS::iterator it1 = v1.begin();
	DIFFS::const_iterator it2 = v2.begin();
	for (; it1 != v1.end() && it2 != v2.end(); ++it1, ++it2)
	{
		*it1 += *it2;
	}
	return v1;
}
DIFFS& operator-=(DIFFS& v1, const DIFFS& v2)
{	
	DIFFS::iterator it1 = v1.begin();
	DIFFS::const_iterator it2 = v2.begin();
	for (; it1 != v1.end() && it2 != v2.end(); ++it1, ++it2)
	{
		*it1 -= *it2;
	}
	return v1;
}



DIFFS operator+(const DIFFS& v1, const DIFFS& v2)
{
	DIFFS ret(v1);
	ret += v2;
	return ret;
}
DIFFS operator-(const DIFFS& v1, const DIFFS& v2)
{
	DIFFS ret(v1);
	ret -= v2;
	return ret;
}

class argument_register;
class dual
{
public:
	/*
	dual(const dual& d1)
	{
		//double      _value;  // 函数值
		//DIFFS       _diffs;  // 偏微分列表 存储x1,x2,x3,...,xn的微分值，计算微分和计算函数值同步。
		_value = d1._value;
	}*/
public:
	dual(const double value = 1, const int variable_order = 0, int variable_dim = 1) :_value(value), _diffs(variable_dim, (double)0.0)
	{
		if (variable_order >= 0 && variable_order<variable_dim)
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
			exp(x2._value*log(x1._value))*(log(x1._value)*x2._diffs + x2._value / x1._value*x1._diffs)
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
			exp(exponent*log(x._value))*exponent / (x._value)*x._diffs);
	}

	friend dual log_e(const dual& d)
	{
		return dual(log(d._value),
			1.0 / d._value * d._diffs);
	}
	friend dual log_10(const dual& x)
	{
		// log_10(x) = ln(x) / ln(10)
		return dual(log(x._value) / log(10.0),
			1.0 / x._value / log(10.0)* x._diffs);

	}
public:
	// disable this function !!!
	dual& operator +=(const dual& d1)
	{
		this->_value += d1._value;
		this->_diffs += d1._diffs;
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
				cout << "ERROR state in line:" << __LINE__ << endl;
				exit(0);
			}
		}
		return _diffs[pos];
	}
	vector<double> get_diff(const vector<dual>& duals)
	{
		vector<double> ret;
		for (auto& d : duals)
		{
			ret.push_back(get_diff(d));
		}

		return ret;
	}
public:
	friend class argument_register;
protected:
	// this interface is specially designed for class argument_register!!!
	DIFFS &    _get_diffs()
	{
		return _diffs;
	}
protected:
	double      _value;  // 函数值
	DIFFS       _diffs;  // 偏微分列表 存储x1,x2,x3,...,xn的微分值，计算微分和计算函数值同步。
};

// To do, regiteration mechanism, automatica set the diff vectors;
//   and automatic setting the patial diff position in the vector for each variable

class argument_register
{
public:
	void begin_regist()
	{
		// clear everything !
		_map.clear();
		_duals.clear();	
		for (const auto& p : _map_array)
			delete p.second; // 删除动态分配的 vector<int>*			
		
		_map_array.clear();
	}
	
	bool regist(const string& argument_name,
				const double argument_value=0) //注册变量名字及计算偏微分时该变量的取值
	{
		auto iter = _map.find(argument_name);
		if (iter != _map.end()) 
		{
			cout << "fail to register variable " << argument_name <<", it is already there!"<< endl;
			return false;
		}
		DIFFS dfs;
		_duals.push_back(dual(argument_value, dfs));
		_map.insert({ argument_name, _duals.size() - 1});
		return true;
	}
	bool regist(const string& argument_name,
	         	const vector<double>& init_datas)
	{
		auto iter = _map_array.find(argument_name);
		if (iter != _map_array.end())
		{
			cout << "fail to register variable " << argument_name << ", it is already there!" << endl;
			return false;
		}

		// 我们允许单变量和数组同名，因为分别通过不同接口得到dual变量的
		vector<int>* pv = new vector<int>;
		DIFFS dfs;
		for (auto& data : init_datas)
		{
			_duals.push_back(dual(data, dfs));
			int order = _duals.size() - 1;
			pv->push_back(order);
		}
		_map_array.insert( {argument_name, pv} );
		return true;
	}
	
	void end_regist()
	{
		// 把每个变量的标记位自动确定下来
		int size = _duals.size();
		int i = 0;
		for (auto& data : _duals)
		{
			DIFFS& dfs = data._get_diffs();
			dfs.resize(size);
			for (auto& data1 : dfs)
				data1 = 0;

			// set variable pos in the diff table
			dfs[i] = 1;
			++i;
		}
	}
public:
	// 返回积分计算的0
	dual zero()
	{
		return dual(0, -1, _duals.size());		
	}
public:
	// 单变量接口
	dual& get_argument(const string& name)
	{
		return _duals[_map[name]];
	}
	dual& operator[](const string& name)
	{
		return _duals[_map[name]];
	}	
	int  diff_pos(const string& name)
	{
		return _map[name];
	}
	// 数组接口
	vector<dual> get_argument_array(const string& name)
	{
		auto iter = _map_array.find(name);
		if (iter == _map_array.end())
		{
			cout << "unregistered array " << name << "[], cannot find it!" << endl;
			return vector<dual>();
		}

		vector<dual> v;
		for (auto& pos : *_map_array[name])
			v.push_back(_duals[pos]);

		return v;
	}
	inline vector<dual> operator()(const string& name)
	{
		return get_argument_array(name);
	}
	
protected:
	vector<dual> _duals;
	map<const string, int> _map;  //记录在数组中的下表
	map<const string, vector<int>*> _map_array;
};




/*
	test array
*/
void test_1()
{
	// sample code for calling the multi-dimensional auto-diff feature

	// 1. register all variables with its name and value
	argument_register reg;
	reg.begin_regist();
		reg.regist("x1", 2);  // x1, partial diff on x1 = 2
		reg.regist("x2", 3);  // x2, partial diff on x2 = 3
		reg.regist("x3", 0);
		reg.regist("x4", 1);
		reg.regist("x1", 1); // already registed, it will trig a failure!!!
	reg.end_regist();

	// 2. get the corresponding variable using its name
	auto& x1 = reg["x1"];
	auto& x2 = reg["x2"];
	auto& x3 = reg["x3"];
	auto& x4 = reg["x4"];

	// 3. computing y and its parital dirivative values on x1,x2,x3,x4
	auto y1 = x1*x2*x1*x2 + sin(x2)*log_e(x1) + 2 * x3 + x4 + (x1^x4) + log_10(x1);

	// 4. output results
	cout << "dy/dx1 = " << y1.get_diff(x1) << endl;
	cout << "dy/dx1 = " << y1.get_diff(x2) << endl;
	cout << "dy/dx3 = " << y1.get_diff(x3) << endl;
	cout << "dy/dx4 = " << y1.get_diff(x4) << endl;
	/*
	   output:
		   dy/dx1 = 37.2877
		   dy/dx1 = 23.3138
		   dy/dx3 = 2
		   dy/dx4 = 2.38629
	
	*/
}


// test DR ...
void test_array1()
{
	cout << "---------------------------------------------------------" << endl;
	cout << "test DR and its derivatives over theta and L" << endl;
	// a test of DR algorithm to find which will greatly affect the system performance.
	const double PI = 3.1415926525;
	argument_register reg1;
	reg1.begin_regist();
		reg1.regist("theta1", 6 * PI / 180.00);
		reg1.regist("theta2", 6 * PI / 180.00);
		reg1.regist("theta3", 6 * PI / 180.00);
		reg1.regist("theta4", 6 * PI / 180.00);
		reg1.regist("theta5", 6 * PI / 180.00);
	
		reg1.regist("L1", 20);
		reg1.regist("L2", 20);
		reg1.regist("L3", 20);
		reg1.regist("L4", 20);
		reg1.regist("L5", 20);
	reg1.end_regist();

	dual theta[] =
	{
		reg1["theta0"],
		reg1["theta1"],
		reg1["theta2"],
		reg1["theta3"],
		reg1["theta4"]
	};
	dual L[] =
	{
		reg1["L0"],
		reg1["L1"],
		reg1["L2"],
		reg1["L3"],
		reg1["L4"]
	};

	dual x[] = { reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero() };  // 6 data
	dual y[] = { reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero() };
	dual alpha[] = { reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero(), reg1.zero() };

	int size = sizeof(theta) / sizeof(theta[0]);

	// DR algorithms using difference-angle from gyro and distance(L) from odometers
	for (int i = 0; i < size; ++i)
	{
		alpha[i + 1] = alpha[i] + theta[i];
		x[i + 1] = x[i] + L[i] * cos(alpha[i + 1]);
		y[i + 1] = y[i] + L[i] * sin(alpha[i + 1]);
	}
	// out alpha, x, y
	ofstream os("d:\\111.csv");
	for (int i = 0; i < size; ++i)
	{

		cout << x[i + 1].get_value() << ",";
		cout << y[i + 1].get_value() << ",";
		cout << alpha[i + 1].get_value() << endl;

		os << x[i + 1].get_value() << ",";
		os << y[i + 1].get_value() << ",";
		os << alpha[i + 1].get_value() << endl;
	}
	os.clear();
	os.close();
	//output derivative value of x0 to L0,L1,L2,L3,L4, and theta1,theta2,theta3,theta4,theta5
	cout << "--------------------------------------------" << endl;
	for (int i = 0; i < size; ++i)
	{
		cout << "dx/d_L[" << i << "] = " << x[size].get_diff(L[i]) << endl;
	}
	for (int i = 0; i < size; ++i)
	{
		cout << "dx/d_theta[" << i << "] = " << x[size].get_diff(theta[i]) << endl;
	}
	cout << "--------------------------------------------" << endl;
	for (int i = 0; i < size; ++i)
	{
		cout << "dy/d_L[" << i << "] = " << y[size].get_diff(L[i]) << endl;
	}
	for (int i = 0; i < size; ++i)
	{
		cout << "dy/d_theta[" << i << "] = " << y[size].get_diff(theta[i]) << endl;
	}
}

/*
test real array
*/
void test_array2()
{
	cout << "-------------------test real array --------------------" << endl;
	argument_register reg;
	reg.begin_regist();
		reg.regist("x", 1.5);
		reg.regist("hello", { 1, 2, 3, 4, 1 });
	reg.end_regist();

	auto x = reg["x"];
	vector<dual> hellos = reg("hello");   // use [] to get singla variables, and () for array

	auto y = hellos[0] * hellos[1] + sin(hellos[3]) + (x ^ 2);

	auto diffs = y.get_diff(hellos);
	int i = 0;
	for (auto& diff : diffs)
	{
		cout << "dy/d_hello[" << i << "] = " << diff << endl;
		++i;
	}
	cout << "dy/dx=" << y.get_diff(x) << endl;
	/*
	output:
	------------------- test2 --------------------
	dy/d_hello[0] = 2
	dy/d_hello[1] = 1
	dy/d_hello[2] = 0
	dy/d_hello[3] = -0.653644
	dy/d_hello[4] = 0
	dy/dx=3

	*/
}

int main()
{
	// test with registering single-variables
	cout << "--------3. test with registering single variable -- -----------" << endl;
	test_1();

	// test with registering single-variable, and combining array 
	cout << "--------3. test with registering single variable and combining array-- -----------" << endl;
	test_array1();

	// test with registering array
	cout << "--------3. test with registering array-- -----------" << endl;
	test_array2();
	
	cout << "---------- end -----------" << endl;
}
