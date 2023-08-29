# AutoDiffM
automatic differential framework in C++ supporting vectors, 多变量微分自动计算框架(存粹为了展示基本原理）

# 功能
支持多变量正向自动微分计算，极端原始的原理版本，反向版本to be done.


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


# 联系
别联系了，还很原始，我有功夫改进吧。 
