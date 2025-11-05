#ifndef XTMB_LA_CONSTANTS_H
#define XTMB_LA_CONSTANTS_H

namespace xiaotu {

    //! @brief 存储方式
    enum EAlignType : uint32_t {
        //! 列优先存储
        eColMajor = 0x00,
        //! 行优先存储
        eRowMajor = 0x01
    };

    enum EStoreType : uint32_t {
        eStoreNone = 0x00,
        //! std::vector<double> data
        eStoreVector = 0x01,
        //! double data[n]
        eStoreArray = 0x02,
        //! 可动态扩展 std::vector
        eStoreDyna = 0x03,
        //! 只是个代理,不管理内存,用于 SubMatrix
        eStoreProxy = 0x04
    };

    #define SMALL_VALUE 1e-12

}

#endif


