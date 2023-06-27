#ifndef XTMB_LA_CONSTANTS_H
#define XTMB_LA_CONSTANTS_H

namespace xiaotu {
namespace math {

    //! @brief 参考 eigen 的实现，用于动态内存的模板定义
    const int Dynamic = -1;

    //! @brief 存储方式
    enum EStorageOptions {
        //! 列优先存储
        eColMajor = 0,
        //! 行优先存储
        eRowMajor = 1
    };

}
}

#endif
