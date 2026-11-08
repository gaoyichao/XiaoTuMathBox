#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <matplot/matplot.h>

#include <iostream>


void Plot(std::vector<double> const & t,
          std::vector<double> const & y,
          std::string title, std::string y_title)
{
    using namespace matplot;
    auto f = figure(true);

    plot(t, y, "-b")->line_width(2).display_name(title);
    hold(on);

    matplot::title(title);
    xlabel("t");
    ylabel(y_title);
    grid(on);
    legend()->location(legend::general_alignment::topleft);
    hold(off);
    f->save(title + ".png");
}


int main() {

    xiaotu::Matrix<double, 2, 2> A = {
        0, 1,
        0, 0,
    };
    xiaotu::Matrix<double, 2, 1> B = {
        0,
        0.5 
    };
    xiaotu::Matrix<double, 2, 1> x0 = {
        0,
        30,
    };
    xiaotu::LinearSystem<double> sys(x0, A, B);

    std::vector<double> t_list;
    std::vector<double> v_list;
    std::vector<double> s_list;

    double step = 0.1;
    double max_t = 30;
    double t = 0;
    while (t < max_t) {
        auto const & curr_x = sys.State();

        t_list.push_back(t);
        s_list.push_back(curr_x(0));
        v_list.push_back(curr_x(1));

        xiaotu::Matrix<double, 1, 1> u = { -curr_x(1) };
        auto next_x = sys.EulerStep(step, &u);
        sys.UpdateState(next_x);
        t += step;
    }

    Plot(t_list, v_list, "速度曲线", "v(t) m/s");
    Plot(t_list, s_list, "位移曲线", "s(t) m");

    std::cout << sys.State() << std::endl;
}

