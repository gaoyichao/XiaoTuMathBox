#ifndef XTMB_LA_DECLARATIONS_H
#define XTMB_LA_DECLARATIONS_H


namespace xiaotu::math {

    //! @brief 矩阵相关的类型萃取器
    template<typename T>
    struct Traits;

    //! @brief 矩阵视图
    template <typename T, int numRows, int numCols, EStorageOptions option = EStorageOptions::eColMajor>
    class MatrixView;

    //! @brief 稠密矩阵
    template <typename T, int numRows, int numCols, EStorageOptions option = EStorageOptions::eColMajor>
    class Matrix;
}

#endif
