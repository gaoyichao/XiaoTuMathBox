#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <matplot/matplot.h>

#include <iostream>


void NaivePieceCubicHermite()
{
    std::vector<double> x_nodes = {1.0, 2.0, 4.0};
    std::vector<double> y_nodes = {1.0, 4.0, 16.0};
    std::vector<double> d_nodes = {2.0, 4.0, 8.0};

    auto npch = xiaotu::PieceCubicHermite(x_nodes, y_nodes, d_nodes);

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
    hold(off);
}

void FritschCarlsonPieceCubicHermite()
{
    std::vector<double> x_nodes = {0.0, 1.0, 2.0, 3.0,  4.0,  5.0};
    std::vector<double> y_nodes = {0.0, 0.0, 0.0, 0.0, 10.0, 10.0};

    auto d_fc = xiaotu::FritschCarlsonDerivatives<double>(x_nodes, y_nodes);
    auto fc_pch = xiaotu::PieceCubicHermite(x_nodes, y_nodes, d_fc);

    auto d_cr = xiaotu::CatmullRomDerivatives<double>(x_nodes, y_nodes);
    auto cr_pch = xiaotu::PieceCubicHermite(x_nodes, y_nodes, d_cr);


    std::vector<double> curve_x;
    std::vector<double> curve_fc;
    std::vector<double> curve_cr;
    double step = 0.05;
    for (double x = x_nodes.front(); x <= x_nodes.back(); x += step) {
        curve_x.push_back(x);
        curve_fc.push_back(fc_pch(x));
        curve_cr.push_back(cr_pch(x));
    }

    using namespace matplot;

    plot(curve_x, curve_fc, "-b")->line_width(2).display_name("Fritsch & Carlson 插值");
    hold(on);
    plot(curve_x, curve_cr, "-r")->line_width(2).display_name("Catmull & Rom 插值");

    scatter(x_nodes, y_nodes, 10)->marker_face(true).marker_color({1.0, 0.0, 0.0}).marker_face_color({1.0, 0.3, 0.3});
    
    title("分段三次 Hermite 插值");
    xlabel("X Axis");
    ylabel("Y Axis");
    grid(on);
    legend()->location(legend::general_alignment::topleft);
    save("fc_cr分段三次Hermite插值.png");
}

int main() {
    NaivePieceCubicHermite();
    FritschCarlsonPieceCubicHermite();
    return 0;
}

