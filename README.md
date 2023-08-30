# AutoDiffM
automatic differential framework in C++ supporting vectors, 多变量微分自动计算框架(存粹为了展示基本原理）

# 功能
支持多变量正向自动微分计算，极端原始的原理版本，反向版本to be done.

```cpp
int main()
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
	auto y = x1*x2*x1*x2 + sin(x2)*log_e(x1) + 2 * x3 + x4 + (x1^x4) + log_10(x1);

	// 4. output results
	cout << "dy/dx1 = " << y.get_diff(x1) << endl;
	cout << "dy/dx1 = " << y.get_diff(x2) << endl;
	cout << "dy/dx3 = " << y.get_diff(x3) << endl;
	cout << "dy/dx4 = " << y.get_diff(x4) << endl;
}
```

例子2
```cpp
// test DR ...
void test()
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
···

```c++
输出结果类似如下（我改过程序参数的，不见得一样了）：
fail to register variable x1, it is already there!
dy/dx1 = 37.2877
dy/dx1 = 23.3138
dy/dx3 = 2
dy/dx4 = 2.38629
--------------------------------------------
dx0/dL0 = -18.4861
dx0/dL1 = 0.997564
dx0/dL2 = 0.994522
dx0/dL3 = 0.990268
dx0/dL4 = 0.984808
dx0/dtheta0 = -18.4861
dx0/dtheta1 = -18.4861
dx0/dtheta2 = -8.34699
dx0/dtheta3 = -6.25643
dx0/dtheta4 = -3.47296
--------------------------------------------
dy0/dL0 = 158.756
dy0/dL1 = 0.0697565
dy0/dL2 = 0.104528
dy0/dL3 = 0.139173
dy0/dL4 = 0.173648
dy0/dtheta0 = 158.756
dy0/dtheta1 = 158.756
dy0/dtheta2 = 59.392
dy0/dtheta3 = 39.5015
dy0/dtheta4 = 19.6962
```
```c++
// example 3, supporting registering an array
void test_array()
{
	cout << "---------------------------------------------------------" << endl;
	cout << "test DR and its derivatives over theta and L" << endl;
	// a test of DR algorithm to find which will greatly affect the system performance.
	const double PI = 3.1415926525;
	argument_register reg1;
	reg1.begin_regist();
                // register two arrays, with name theta and L
		reg1.regist("theta", {1.0, 1.0,    2,   3,   5, -0.5});   // array theta
		reg1.regist("L",     {10,   10, 10.5,  1.5,  0,    0});   // array L
	reg1.end_regist();
        vector<dual> theta = reg1["theta"];
        vector<dual> L = reg1["L"];

        // now use the arrays to compute ....
}
```


# 联系
别联系了，还很原始，我有功夫改进吧。 
