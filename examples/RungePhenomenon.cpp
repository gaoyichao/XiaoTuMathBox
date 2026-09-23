/************************************************************************
 * 
 * Runge 现象
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=拉格朗日多项式#Runge
 * 
 ***********************************************************************/

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <matplot/matplot.h>

#include <iostream>


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
    auto p5 = xiaotu::BarycentricLagrange(x_nodes, y_nodes);
    std::vector<double> p5_y;

    sample(7, x_nodes, y_nodes);
    auto p7 = xiaotu::BarycentricLagrange(x_nodes, y_nodes);
    std::vector<double> p7_y;
    
    sample(9, x_nodes, y_nodes);
    auto p9 = xiaotu::BarycentricLagrange(x_nodes, y_nodes);
    std::vector<double> p9_y;

    sample(15, x_nodes, y_nodes);
    auto p15 = xiaotu::BarycentricLagrange(x_nodes, y_nodes);
    std::vector<double> p15_y;

    sample(17, x_nodes, y_nodes);
    auto p17 = xiaotu::BarycentricLagrange(x_nodes, y_nodes);
    std::vector<double> p17_y;

    double lbound = -5.0;
    double hbound = 5.0;

    std::vector<double> func_x;
    std::vector<double> func_y;
    for (double x = lbound; x <= hbound; x += 0.1) {
        func_x.push_back(x);
        func_y.push_back(func(x));
        p5_y.push_back(p5(x));
        p7_y.push_back(p7(x));
        p9_y.push_back(p9(x));
        p15_y.push_back(p15(x));
        p17_y.push_back(p17(x));
    }

    matplot::figure(true); 
    matplot::gcf()->size(1000, 600);
    
    matplot::plot(func_x, func_y, "-ok")->line_width(2).marker_size(8).display_name("y");
    matplot::hold(matplot::on);
    matplot::plot(func_x, p5_y, "--")->line_width(2).display_name("p5");
    matplot::plot(func_x, p7_y, "--")->line_width(2).display_name("p7");
    matplot::plot(func_x, p9_y, "--")->line_width(2).display_name("p9");
    matplot::plot(func_x, p15_y, "--s")->line_width(2).marker_size(8).marker_face_color("none").display_name("p15");
    matplot::plot(func_x, p17_y, "--^")->line_width(2).marker_size(8).marker_face_color("none").display_name("p17");

    matplot::title("龙格现象");
    matplot::xlabel("X Axis");
    matplot::ylabel("Values");
    matplot::legend();
    matplot::grid(true);
    matplot::save("龙格现象.png");

    return 0;
}
