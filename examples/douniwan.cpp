#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>


namespace xiaotu {

    void _douniwan_();

}

int main(int argc, char * argv[])
{
    xiaotu::_douniwan_();

    xiaotu::Polynomial<double> a({1.0, 2.0, 3.0});
    std::cout << a << std::endl;

    return 0;
}
