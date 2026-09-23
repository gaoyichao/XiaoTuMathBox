#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <matplot/matplot.h>

#include <iostream>


int main() {
    std::vector<double> x_nodes = {1.0, 2.0, 4.0};
    std::vector<double> y_nodes = {1.0, 4.0, 16.0};
    std::vector<double> d_nodes = {2.0, 4.0, 8.0};

    auto npch = xiaotu::NaivePieceCubicHermite(x_nodes, y_nodes, d_nodes);

    std::vector<double> curve_x;
    std::vector<double> curve_y;
    double step = 0.1;
    for (double x = 1.0; x <= 4.0; x += step) {
        curve_x.push_back(x);
        curve_y.push_back(npch(x));
    }

    using namespace matplot;

    plot(curve_x, curve_y, "-b")->line_width(2).display_name("Hermite 插值");
    hold(on);

    scatter(x_nodes, y_nodes, 10)->marker_face(true).marker_color({1.0, 0.0, 0.0}).marker_face_color({1.0, 0.3, 0.3});
    
    title("朴素分段三次 Hermite 插值(Naive PCHI");
    xlabel("X Axis");
    ylabel("Y Axis");
    grid(on);
    legend();
    save("朴素分段三次Hermite插值.png");

    return 0;
}

