#ifndef XTMB_LA_CONSTANTS_H
#define XTMB_LA_CONSTANTS_H

namespace xiaotu::math {

    //! @brief 存储方式
    enum EStorageOptions {
        //! 列优先存储
        eColMajor = 0,
        //! 行优先存储
        eRowMajor = 1
    };

#define SMALL_VALUE 1e-6

}

#endif


