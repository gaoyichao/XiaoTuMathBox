#ifndef XTMB_LA_ORTHOGONOLIZATION_H
#define XTMB_LA_ORTHOGONOLIZATION_H

#include <stdint.h>
#include <cassert>
#include <vector>
#include <iostream>


namespace xiaotu::math {


    //! @brief Gram-Schmidt 标准正交化
    //!
    //! @param [in] col_vectors 按照列向量的形式排列的基向量
    //! @param [out] ortho 标准正交基
    template <typename MatViewIn, typename MatViewOut>
    void GramSchmidt(MatViewIn const & col_vectors, MatViewOut & ortho)
    {
        int num = col_vectors.Cols();
        int dim = col_vectors.Rows();

        for (int i = 0; i < num; i++) {
            
        }
    }

}


#endif

