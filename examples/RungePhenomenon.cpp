/************************************************************************
 * 
 * Runge 现象
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=拉格朗日多项式#Runge
 * 
 ***********************************************************************/


#include <iostream>
#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>


double func(double x)
{
    return 1 / (1 + x * x);
}

void sample(int n, std::vector<double> & x_nodes, std::vector<double> & y_nodes)
{
    double lbound = -5.0;
    double hbound = 5.0;

    double step = (hbound - lbound) / n;
    x_nodes.resize(n+1);
    y_nodes.resize(n+1);

    for (int i = 0; i <= n; i++) {
        x_nodes[i] = lbound + i * step;
        y_nodes[i] = func(x_nodes[i]);
    }
}


int main(int argc, char * argv[])
{
    std::vector<double> x_nodes;
    std::vector<double> y_nodes;

    sample(5, x_nodes, y_nodes);
    auto p5 = xiaotu::LagrangePolynomial<xiaotu::Polynomial<double>>(x_nodes, y_nodes);
    std::cout << "p5 : " << p5 << std::endl;

    sample(7, x_nodes, y_nodes);
    auto p7 = xiaotu::LagrangePolynomial<xiaotu::Polynomial<double>>(x_nodes, y_nodes);
    std::cout << "p7 : " << p7 << std::endl;
    
    sample(9, x_nodes, y_nodes);
    auto p9 = xiaotu::LagrangePolynomial<xiaotu::Polynomial<double>>(x_nodes, y_nodes);
    std::cout << "p9 : " << p9 << std::endl;

    sample(15, x_nodes, y_nodes);
    auto p15 = xiaotu::LagrangePolynomial<xiaotu::Polynomial<double>>(x_nodes, y_nodes);
    std::cout << "p15: " << p15 << std::endl;

    sample(17, x_nodes, y_nodes);
    auto p17 = xiaotu::LagrangePolynomial<xiaotu::Polynomial<double>>(x_nodes, y_nodes);
    std::cout << "p17: " << p17 << std::endl;

    std::cout << "------------------------------------" << std::endl;
    double lbound = -5.0;
    double hbound = 5.0;

    std::cout << "x,y,p5,p7,p9,p15,p17" << std::endl;
    for (double x = lbound; x <= hbound; x += 0.1) {
        std::cout << x << ",";
        std::cout << func(x) << ",";
        std::cout << p5(x)   << ",";
        std::cout << p7(x)   << ",";
        std::cout << p9(x)   << ",";
        std::cout << p15(x)  << ",";
        std::cout << p17(x)  << std::endl;
    }

    return 0;
}
