#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

namespace xiaotu {

    void _douniwan_()
    {
        Matrix<double, 3, 3> A = {
            9, -3, 1,
            1,  1, 1,
            4,  2, 1
        };
        Matrix<double, 3, 1> b = { 20, 0, 10 };

        GaussJordanEliminate(A, &b);
        XTLog(std::cout) << "A = " << A << std::endl;
        XTLog(std::cout) << "b = " << b << std::endl;

        xiaotu::Polynomial<double> a({1.0, 2.0, 3.0});
        std::cout << a << std::endl;

    }

}

