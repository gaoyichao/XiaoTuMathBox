#ifndef XTMB_LA_DECLARATIONS_H
#define XTMB_LA_DECLARATIONS_H


namespace xiaotu::math {

    //! @brief 矩阵视图
    template <typename T, int numRows, int numCols, EStorageOptions option = EStorageOptions::eColMajor>
    class MatrixView;

}

#endif
