#ifndef XTMB_LA_CONSTANTS_H
#define XTMB_LA_CONSTANTS_H

namespace xiaotu::math {

    //! @brief 存储方式
    enum EAlignType : uint32_t {
        //! 列优先存储, 无存储空间
        eColMajor = 0x00,
        //! 行优先存储
        eRowMajor = 0x01
    };

    enum EStoreType : uint32_t {
        //! std::vector<double> data
        eStoreVector = 0x00,
        //! double data[n]
        eStoreArray = 0x01
    };

#define SMALL_VALUE 1e-6

}

#endif


