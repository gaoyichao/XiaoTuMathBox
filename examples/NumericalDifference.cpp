#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <matplot/matplot.h>

void Plot(std::vector<double> const & x,
          std::vector<double> const & y,
          std::string const & t)
{
    using namespace matplot;
    auto f = figure(true);
    auto ax = f->current_axes();

    ax->plot(x, y, "-ok")->line_width(2).marker_size(8);
    ax->title(t);
    ax->xlabel("X Axis");
    ax->ylabel("Y Axis");
    ax->grid(on);
    f->save(t + ".png");
}

int main() {
    auto fn = [](double x) {
        return std::pow(x, 9);
    };

    std::vector<double> h_list;
    std::vector<double> e_list;
    std::vector<double> e3_list;
    std::vector<double> e3e_list;
    std::vector<double> e5_list;
    std::vector<double> e5e_list;

    double factor = 0.5;
    double h = 0.1;
    for (int i = 0; i < 10; i++) {
        auto f1 = xiaotu::ForwardDifference<double>(fn, 1.0, h);
        auto f3 = xiaotu::MidThreePoints<double>(fn, 1.0, h);
        auto f3_e = xiaotu::EndThreePoints<double>(fn, 1.0, h);
        auto f5 = xiaotu::MidFivePoints<double>(fn, 1.0, h);
        auto f5_e = xiaotu::EndFivePoints<double>(fn, 1.0, h);

        h_list.push_back(-std::log(h));
        e_list.push_back(xiaotu::RelativeError(9.0, f1));
        e3_list.push_back(xiaotu::RelativeError(9.0, f3));
        e3e_list.push_back(xiaotu::RelativeError(9.0, f3_e));
        e5_list.push_back(xiaotu::RelativeError(9.0, f5));
        e5e_list.push_back(xiaotu::RelativeError(9.0, f5_e));

        h *= factor;
    }

    Plot(h_list, e_list, "前向微分");
    Plot(h_list, e3_list, "三点中心公式");
    Plot(h_list, e3e_list, "三点端点公式");
    Plot(h_list, e5_list, "五点中心公式");
    Plot(h_list, e5e_list, "五点端点公式");


    return 0;
}
