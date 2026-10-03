#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <matplot/matplot.h>

#include <iostream>


void NatureSpline()
{
    std::vector<double> x_nodes = {0.0, 1.0, 2.0, 3.0,  4.0,  5.0};
    std::vector<double> y_nodes = {0.0, 0.0, 0.0, 0.0, 10.0, 10.0};
    // std::vector<double> d_nodes = {2.0, 4.0, 8.0};

    auto spline = xiaotu::CubicSpline(x_nodes, y_nodes);

    std::vector<double> curve_x;
    std::vector<double> curve_y;
    double step = 0.1;
    for (double x = x_nodes.front(); x <= x_nodes.back(); x += step) {
        curve_x.push_back(x);
        curve_y.push_back(spline(x));
    }

    using namespace matplot;
    auto f = figure(true);

    plot(curve_x, curve_y, "-b")->line_width(2).display_name("自然三次样条(nature cubic spline)");
    hold(on);

    scatter(x_nodes, y_nodes, 10)->marker_face(true).marker_color({1.0, 0.0, 0.0}).marker_face_color({1.0, 0.3, 0.3});
    
    title("自然三次样条(nature cubic spline)");
    xlabel("X Axis");
    ylabel("Y Axis");
    grid(on);
    legend();
    hold(off);
    f->save("自然三次样条.png");
}


void ClampedSpline()
{
    std::vector<double> x_nodes = {0.0, 1.0, 2.0, 3.0,  4.0,  5.0};
    std::vector<double> y_nodes = {0.0, 0.0, 0.0, 0.0, 10.0, 10.0};
    // std::vector<double> d_nodes = {2.0, 4.0, 8.0};

    auto spline = xiaotu::CubicSpline(x_nodes, y_nodes, 0.0, 0.0);
    auto spline1 = xiaotu::CubicSpline(x_nodes, y_nodes, 1.0, 0.0);
    auto spline2 = xiaotu::CubicSpline(x_nodes, y_nodes, 10.0, 0.0);

    std::vector<double> curve_x;
    std::vector<double> curve_y;
    std::vector<double> curve_y1;
    std::vector<double> curve_y2;
    double step = 0.1;
    for (double x = x_nodes.front(); x <= x_nodes.back(); x += step) {
        curve_x.push_back(x);
        curve_y.push_back(spline(x));
        curve_y1.push_back(spline1(x));
        curve_y2.push_back(spline2(x));
    }

    using namespace matplot;
    auto f = figure(true);

    plot(curve_x, curve_y, "-b")->line_width(2).display_name("df_x0 = 0.0, df_xn = 0.0");
    hold(on);
    plot(curve_x, curve_y1, "-r")->line_width(2).display_name("df_x0 = 1.0, df_xn = 0.0");
    plot(curve_x, curve_y2, "-g")->line_width(2).display_name("df_x0 = 10.0, df_xn = 0.0");

    scatter(x_nodes, y_nodes, 10)->marker_face(true).marker_color({1.0, 0.0, 0.0}).marker_face_color({1.0, 0.3, 0.3});
    
    title("固定三次样条(clamped cubic spline)");
    xlabel("X Axis");
    ylabel("Y Axis");
    grid(on);
    legend()->location(legend::general_alignment::topleft);
    hold(off);
    f->save("固定三次样条.png");
}


int main() {
    NatureSpline();
    ClampedSpline();

    return 0;
}


