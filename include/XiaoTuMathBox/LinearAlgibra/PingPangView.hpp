#ifndef XTMB_LA_MATRIX_PINGPANG_VIEW_H
#define XTMB_LA_MATRIX_PINGPANG_VIEW_H

#include <cassert>
#include <iostream>
#include <initializer_list>
#include <vector>


namespace xiaotu {

    /**
     * @brief 乒乓队列, 主要用于 QR 迭代、SVD 分解
     */
    template <typename Mat, bool IsMatrix>
    class PingPangView
    {
        public:
            typedef MatrixSubView<Mat> SubView;

            /**
             * @brief 默认构造函数
             */
            PingPangView()
            {}

            PingPangView(SubView && view)
            {
                mPingPtr->push_back(view);
            }

            /**
             * @brief 初始化
             */
            PingPangView & Init(SubView && view)
            {
                Clear();
                mPingPtr->push_back(view);
                return *this;
            }

            /**
             * @brief 清空
             */
            void Clear()
            {
                mPing.clear();
                mPang.clear();
                mPingPtr = &mPing;
                mPangPtr = &mPang;
            }

            /**
             * @brief 总是处理 Ping 队列
             */
            std::vector<SubView> & Ping()
            {
                return *mPingPtr;
            }

            std::vector<SubView> * PingPtr()
            {
                return mPingPtr;
            }

            /**
             * @brief 总是在 Pang 队列中作缓存
             */
            std::vector<SubView> & Pang()
            {
                return *mPangPtr;
            }

            std::vector<SubView> * PangPtr()
            {
                return mPangPtr;
            }

            /**
             * @brief 交换乒乓队列
             */
            PingPangView & Swap()
            {
                std::swap(mPingPtr, mPangPtr);
                mPangPtr->clear();
                return *this;
            }
        private:
            std::vector<SubView> mPing;
            std::vector<SubView> mPang;
            std::vector<SubView> * mPingPtr{&mPing};
            std::vector<SubView> * mPangPtr{&mPang};


    };



}

#endif

