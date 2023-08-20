# AutoDiffM
automatic differential framework in C++ supporting vectors, 多变量微分自动计算框架

# 功能
支持多变量自动微分计算，极端原始的原理版本，期望改进后最终实现如下调用形式功能

		argument_register reg;
		reg.begin_regist();
			reg.register("x1",2);  // x1, partial diff on x1 = 2
			reg.register("x2",3);  // x2, partial diff on x2 = 3
		reg.end_register();

		auto& x1 = reg["x1"];
		auto& x2 = reg["x2"];

		auto y = x1*x2 + sin(x2)*log(x1);

		cout<<"dy/dx1 = "<<y.get_diff(x1)<<endl;
		cout<<"dy/dx2 = "<<y.get_diff(x2)<<endl;


# 联系
别联系了，还很原始，我有功夫改进吧。 
