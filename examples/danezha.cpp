#include <matplot/matplot.h>
#include <cmath>
#include <vector>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>



void plot_trapezoidal_rule()
{
    using namespace matplot;

    // 1. 准备数据
    auto f = [](double x) { return sin(x) + 2.0 + 0.1 * std::pow(x - 2, 3); };
    
    double a = 1.0;
    double b = 4.0;
    
    std::vector<double> x_curve = linspace(0.5, 4.5, 100);
    std::vector<double> y_curve = transform(x_curve, f);

    // 2. 梯形顶点坐标（用于半透明阴影填充）
    std::vector<double> x_fill = {a, b, b, a};
    std::vector<double> y_fill = {0, 0, f(b), f(a)};

    // 开始绘图
    auto h = figure();
    hold(on);

    // 3. 完全隐藏系统默认的边框、四周线、以及自带的数字刻度
    gca()->visible(false);

    // 4. 绘制填充区域
    auto fl = fill(x_fill, y_fill);
    fl->color({0.75, 0.82, 0.93, 0.6}); // 淡蓝色半透明

    // 5. 绘制函数曲线 y = f(x)
    auto p1 = plot(x_curve, y_curve);
    p1->line_width(2.5).color("blue");

    // 6. 绘制近似直线 y = P_1(x)
    auto p2 = plot(std::vector<double>{a, b}, std::vector<double>{f(a), f(b)});
    p2->line_width(2).color("black");

    // 7. 绘制两侧的垂直边界实线（对应原图的细黑实线）
    auto l1 = plot(std::vector<double>{a, a}, std::vector<double>{0, f(a)});
    l1->line_width(1).color("black");
    auto l2 = plot(std::vector<double>{b, b}, std::vector<double>{0, f(b)});
    l2->line_width(1).color("black");

    // 8. 🔥 绘制带有“负数尾巴”的 x 轴和 y 轴
    // x 轴：从 -0.5 开始，穿过原点，画到 4.8
    auto axis_x = arrow(-0.5, 0.0, 4.8, 0.0); 
    axis_x->color("black").line_width(1.5);
    text(4.8, -0.2, "x")->font_size(14); // x 轴标签

    // y 轴：从 -0.5 开始，穿过原点，画 to 3.8
    auto axis_y = arrow(0.0, -0.5, 0.0, 3.8); 
    axis_y->color("black").line_width(1.5);
    text(-0.2, 3.8, "y")->font_size(14); // y 轴标签

    // 9. 添加圆点标志（交点）
    plot(std::vector<double>{a, b}, std::vector<double>{f(a), f(b)}, "ko")->marker_face_color("black");

    // 10. 文字标注微调
    text(2.3, f(2.3) + 0.4, "y = f(x)")->font_size(12).color("blue");
    text(3.4, f(3.4) - 0.4, "y = P_1(x)")->font_size(12).color("black");
    
    // X轴对应点位置标注 (向下微调一点，避开横轴)
    text(a - 0.2, -0.25, "a = x_0")->font_size(12);
    text(b - 0.2, -0.25, "x_1 = b")->font_size(12);

    // 11. 锁定视口范围
    xlim({-0.8, 5.2});
    ylim({-0.8, 4.2});

    save("trapezoidal_rule.svg");
}

void plot_simpson_rule()
{
    using namespace matplot;

    // 1. 准备数据
    auto f = [](double x) { return sin(x) + 2.0 + 0.1 * std::pow(x - 2, 3); };
    
    double a = 1.0;
    double b = 4.0;
    double mid = (a + b) / 2.0; // 辛普森法则的中点
    
    std::vector<double> x_curve = linspace(0.5, 4.5, 100);
    std::vector<double> y_curve = transform(x_curve, f);

    // 2. 构造穿过 (a, f(a)), (mid, f(mid)), (b, f(b)) 的二次插值多项式 P_2(x)
    auto P2 = [=](double x) {
        double l0 = ((x - mid) * (x - b)) / ((a - mid) * (a - b));
        double l1 = ((x - a) * (x - b)) / ((mid - a) * (mid - b));
        double l2 = ((x - a) * (x - mid)) / ((b - a) * (b - mid));
        return f(a) * l0 + f(mid) * l1 + f(b) * l2;
    };

    // 3. 准备抛物线的数据（用于绘制近似曲线和阴影填充区域）
    std::vector<double> x_para = linspace(a, b, 50);
    std::vector<double> y_para = transform(x_para, P2);

    // 构造填充区域的封闭多边形顶点坐标
    std::vector<double> x_fill = x_para;
    x_fill.push_back(b);
    x_fill.push_back(a);
    
    std::vector<double> y_fill = y_para;
    y_fill.push_back(0.0);
    y_fill.push_back(0.0);

    // 开始绘图
    auto h = figure();
    hold(on);

    // 4. 完全隐藏系统默认的边框、四周线、以及自带的数字刻度
    gca()->visible(false);

    // 5. 绘制抛物线下方的填充区域
    auto fl = fill(x_fill, y_fill);
    fl->color({0.75, 0.82, 0.93, 0.6}); // 淡蓝色半透明

    // 6. 绘制函数曲线 y = f(x)
    auto p1 = plot(x_curve, y_curve);
    p1->line_width(2.5).color("blue");

    // 7. 绘制近似二次抛物线 y = P_2(x)
    auto p2 = plot(x_para, y_para);
    p2->line_width(2).color("black");

    // 8. 绘制两侧的垂直边界实线
    auto l1 = plot(std::vector<double>{a, a}, std::vector<double>{0, f(a)});
    l1->line_width(1).color("black");
    auto l2 = plot(std::vector<double>{b, b}, std::vector<double>{0, f(b)});
    l2->line_width(1).color("black");
    
    // 9. 绘制中间节点的垂直细虚线
    auto l3 = plot(std::vector<double>{mid, mid}, std::vector<double>{0, f(mid)}, "--");
    l3->line_width(1).color("gray");

    // 10. 🔥 绘制带有“负数尾巴”的 x 轴和 y 轴
    auto axis_x = arrow(-0.5, 0.0, 4.8, 0.0); 
    axis_x->color("black").line_width(1.5);
    text(4.8, -0.2, "x")->font_size(14); // x 轴标签

    auto axis_y = arrow(0.0, -0.5, 0.0, 3.8); 
    axis_y->color("black").line_width(1.5);
    text(-0.2, 3.8, "y")->font_size(14); // y 轴标签

    // 11. 添加圆点标志（3个插值节点交点）
    plot(std::vector<double>{a, mid, b}, std::vector<double>{f(a), f(mid), f(b)}, "ko")->marker_face_color("black");

    // 12. 文字标注微调
    text(2.1, f(2.1) + 0.4, "y = f(x)")->font_size(12).color("blue");
    text(3.4, P2(3.4) - 0.4, "y = P_2(x)")->font_size(12).color("black");
    
    // X轴对应点位置标注 (向下微调一点，避开横轴)
    text(a - 0.25, -0.25, "a = x_0")->font_size(12);
    text(mid - 0.15, -0.25, "x_1")->font_size(12);
    text(b - 0.2, -0.25, "x_2 = b")->font_size(12);

    // 13. 锁定视口范围
    xlim({-0.8, 5.2});
    ylim({-0.8, 4.2});

    save("simpson_rule.svg");
}

// 辅助函数：拉格朗日插值计算
double lagrange_interpolate(double x, const std::vector<double>& x_nodes, const std::vector<double>& y_nodes) {
    double total = 0.0;
    size_t n = x_nodes.size();
    for (size_t i = 0; i < n; ++i) {
        double term = y_nodes[i];
        for (size_t j = 0; j < n; ++j) {
            if (j != i) {
                term *= (x - x_nodes[j]) / (x_nodes[i] - x_nodes[j]);
            }
        }
        total += term;
    }
    return total;
}

void plot_closed_newton_cotes()
{
    using namespace matplot;
    // 1. 设定等距节点（严格在 [a, b] 范围内）
    double x0 = 1.2; // a
    double x1 = 2.4;
    double x2 = 3.6;
    double x3 = 4.8;
    double x4 = 6.0; // b
    std::vector<double> x_nodes = {x0, x1, x2, x3, x4};
    std::vector<double> y_nodes = {2.1, 2.5, 2.2, 2.7, 3.1};

    // 2. 模拟原函数 f(x) [蓝色]
    auto f = [&](double x) {
        double base = lagrange_interpolate(x, x_nodes, y_nodes);
        double perturb = 0.25 * std::sin(M_PI * (x - x0) / (x1 - x0));
        return base + perturb;
    };

    // 3. 🚨 关键修改：采样点范围严格限定在 [x0, x4] 之间，绝不向外延伸
    std::vector<double> x_curve = linspace(x0, x4, 200);
    std::vector<double> y_f;
    std::vector<double> y_Pn;
    for (double x : x_curve) {
        y_f.push_back(f(x));
        y_Pn.push_back(lagrange_interpolate(x, x_nodes, y_nodes));
    }

    // 4. 构造积分阴影填充区域（同样严格在 [x0, x4] 内闭合）
    std::vector<double> x_fill = x_curve;
    std::vector<double> y_fill = y_Pn;
    x_fill.push_back(x4); y_fill.push_back(0.0);
    x_fill.push_back(x0); y_fill.push_back(0.0);

    // 开始绘图
    auto h = figure();
    hold(on);

    // 5. 隐藏默认背景与边框刻度
    gca()->visible(false);

    // 6. 绘制半透明积分阴影面积
    auto fl = fill(x_fill, y_fill);
    fl->color({0.68, 0.76, 0.86, 0.7}); 

    // 7. 绘制原函数 y = f(x) [蓝色曲线]
    auto p_f = plot(x_curve, y_f);
    p_f->line_width(2.2).color("blue");

    // 8. 绘制插值多项式 y = P_n(x) [黑色曲线]
    auto p_P = plot(x_curve, y_Pn);
    p_P->line_width(2.2).color("black");

    // 9. 绘制等距节点的垂直边界线与内部中间线
    plot(std::vector<double>{x0, x0}, std::vector<double>{0, y_nodes[0]})->line_width(1.2).color("black");
    plot(std::vector<double>{x1, x1}, std::vector<double>{0, y_nodes[1]})->line_width(1.0).color("black");
    plot(std::vector<double>{x2, x2}, std::vector<double>{0, y_nodes[2]})->line_width(1.0).color("black");
    plot(std::vector<double>{x3, x3}, std::vector<double>{0, y_nodes[3]})->line_width(1.0).color("black");
    plot(std::vector<double>{x4, x4}, std::vector<double>{0, y_nodes[4]})->line_width(1.2).color("black");

    // 10. 绘制带有紧凑负数尾巴的坐标轴
    arrow(-0.3, 0.0, 6.8, 0.0)->color("black").line_width(1.5); // X 轴
    arrow(0.0, -0.3, 0.0, 4.2)->color("black").line_width(1.5); // Y 轴

    // 11. 在曲线上打上离散的数据节点
    plot(x_nodes, y_nodes, "ko")->marker_face_color("black").marker_size(5.5);

    // 12. 纯文本标签标注
    text(6.8, -0.2, "x")->font_size(14);
    text(-0.2, 4.2, "y")->font_size(14);
    
    // 曲线标签位置
    text(2.6, f(2.6) + 0.25, "y = f(x)")->font_size(12).color("blue");
    text(5.0, lagrange_interpolate(5.0, x_nodes, y_nodes) + 0.25, "y = P_n(x)")->font_size(12).color("black");
    
    // 节点标注
    text(x0 - 0.3, -0.25, "x_0 = a")->font_size(11);
    text(x1 - 0.1, -0.25, "x_1")->font_size(11);
    text(x2 - 0.1, -0.25, "x_2")->font_size(11);
    text(x3 - 0.3, -0.25, "x_{n-1}")->font_size(11);
    text(x4 - 0.1, -0.25, "x_n = b")->font_size(11);

    // 13. 锁定视口范围并导出
    xlim({-0.5, 7.2});
    ylim({-0.5, 4.5});
    
    save("closed_newton_cotes.svg");
}

void plot_open_newton_cotes()
{
    using namespace matplot;

    // 1. 设定严格等距的 Open 型区间和内部节点
    double a = 1.0;  // x_{-1}
    double x0 = 2.0; 
    double x1 = 3.0; 
    double x2 = 4.0; 
    double x3 = 5.0; // 对应图中的 x_n
    double b = 6.0;  // x_{n+1}
    
    std::vector<double> x_nodes = {x0, x1, x2, x3};
    std::vector<double> y_nodes = {2.5, 2.6, 2.3, 3.1}; // 节点高度

    // 2. 构造多项式 P_n(x) 和原函数 f(x)
    // 巧妙利用正弦交替扰动，使原函数与多项式在节点精准相交，且区间内完美呈现一上一下的波浪
    auto P_n = [&](double x) { return lagrange_interpolate(x, x_nodes, y_nodes); };
    auto f = [&](double x) { return P_n(x) - 0.25 * std::sin(M_PI * (x - x0)); };

    // 3. 准备绘图数据
    // 蓝色原函数严格限制在 [a, b]
    std::vector<double> x_blue = linspace(a, b, 200);
    std::vector<double> y_blue;
    for (double x : x_blue) y_blue.push_back(f(x));

    // 黑色多项式向外延伸，左侧超出 a，右侧超出 b
    std::vector<double> x_black = linspace(a - 0.4, b + 0.2, 200);
    std::vector<double> y_black;
    for (double x : x_black) y_black.push_back(P_n(x));

    // 4. 构造阴影区域：严格位于 [a, b] 内、蓝色原函数 f(x) 下方的面积
    std::vector<double> x_fill = x_blue;
    std::vector<double> y_fill = y_blue;
    x_fill.push_back(b); y_fill.push_back(0.0);
    x_fill.push_back(a); y_fill.push_back(0.0);

    // 开始绘图
    auto h = figure();
    hold(on);

    // 5. 隐藏默认网格与边框刻度
    gca()->visible(false);

    // 6. 绘制半透明钢蓝色阴影面积 (严格限制在 a 到 b，封顶在蓝色线上)
    auto fl = fill(x_fill, y_fill);
    fl->color({0.68, 0.76, 0.86, 0.7}); 

    // 7. 绘制原函数 y = f(x) [蓝色曲线，刚好在 a 和 b 截断]
    auto p_f = plot(x_blue, y_blue);
    p_f->line_width(2.2).color("blue");

    // 8. 绘制插值多项式 y = P_n(x) [黑色曲线，两头自然延伸出去]
    auto p_P = plot(x_black, y_black);
    p_P->line_width(2.2).color("black");

    // 9. 绘制左右两侧 [a, b] 的垂直边界实线 (封顶到蓝色原函数)
    plot(std::vector<double>{a, a}, std::vector<double>{0, f(a)})->line_width(1.2).color("black");
    plot(std::vector<double>{b, b}, std::vector<double>{0, f(b)})->line_width(1.2).color("black");
    
    // 10. 绘制内部各个等距节点的 X 轴微小刻度线 (不向上连线)
    plot(std::vector<double>{x0, x0}, std::vector<double>{0, y_nodes[0]})->line_width(1).color("black");
    plot(std::vector<double>{x1, x1}, std::vector<double>{0, y_nodes[1]})->line_width(1).color("black");
    plot(std::vector<double>{x2, x2}, std::vector<double>{0, y_nodes[2]})->line_width(1).color("black");
    plot(std::vector<double>{x3, x3}, std::vector<double>{0, y_nodes[3]})->line_width(1).color("black");

    // 11. 绘制带有负数尾巴的高清坐标轴
    arrow(-0.4, 0.0, 7.0, 0.0)->color("black").line_width(1.5); // X 轴
    arrow(0.0, -0.4, 0.0, 4.2)->color("black").line_width(1.5); // Y 轴

    // 12. 在内部交点处打上黑色实心圆点
    plot(x_nodes, y_nodes, "ko")->marker_face_color("black").marker_size(5.5);

    // 13. 纯文本标签标注 (规避 LaTeX 渲染困扰)
    text(7.0, -0.2, "x")->font_size(14);
    text(-0.2, 4.2, "y")->font_size(14);
    
    // 曲线标签 (微调位置防止遮挡)
    text(3.3, f(3.3) + 0.15, "y = f(x)")->font_size(12).color("blue");
    text(4.6, P_n(4.6) + 0.25, "y = P_n(x)")->font_size(12).color("black");
    
    // 严格对照原图的 X 轴刻度文本
    text(a - 0.2, -0.25, "x_{-1} = a")->font_size(11);
    text(x0 - 0.1, -0.25, "x_0")->font_size(11);
    text(x1 - 0.1, -0.25, "x_1")->font_size(11);
    text(x2 - 0.1, -0.25, "x_2")->font_size(11);
    text(x3 - 0.1, -0.25, "x_n")->font_size(11);
    text(b - 0.2, -0.25, "x_{n+1} = b")->font_size(11);

    // 14. 锁定视口范围并导出
    xlim({-0.6, 7.4});
    ylim({-0.6, 4.6});
    
    save("open_newton_cotes.svg");
}


double func(double x)
{
    return std::exp(x);
}

void sample(double a, double b, int n,
            std::vector<double> & x_nodes, std::vector<double> & y_nodes)
{
    double step = (b - a) / n;
    x_nodes.resize(n+1);
    y_nodes.resize(n+1);

    for (int i = 0; i <= n; i++) {
        x_nodes[i] = a + i * step;
        y_nodes[i] = func(x_nodes[i]);
    }
}

int main()
{
    // plot_trapezoidal_rule();
    // plot_simpson_rule();
    // plot_closed_newton_cotes();
    // plot_open_newton_cotes();


    {
        std::vector<double> x_nodes;
        std::vector<double> y_nodes;

        double a = 0.0;
        double b = 2.0;
        sample(a, b, 10, x_nodes, y_nodes);
        xiaotu::LagrangePolynomialIntegral lpint(a, b, x_nodes);
        std::cout << lpint.Integral(y_nodes) << std::endl;

        std::cout << func(b) - func(a)  << std::endl;
        std::cout << std::log(func(1)) << std::endl;
        // auto da = xiaotu::LagrangeBasis(a, 0);
        // XTLog(std::cout) << "da = " << da << std::endl;

    }

    return 0;
}
