/*********************************************************************************************
 * 
 * 一个集成了数组和链表的数据结构
 * 
 * 主要是看《计算几何》的时候，说是用链表来表示多边形更方便一些。
 * 我想对多边形顶点的随机访问，也是比较常见的操作。
 * 
 * 所以想结合两种数据结构，用链表来描述顶点的拓扑关系，用数组提供随机访问的能力。
 * 在写 Polygon 的过程中，觉得这个数据结构可能有一定的通用性。
 * 所以专门放置到 XiaoTuDataBox 中收集起来。
 * 
 * https://gaoyichao.com/Xiaotu/?book=%E8%AE%A1%E7%AE%97%E5%87%A0%E4%BD%95_in_C&title=Chapter1.4
 *
 ***************************************************************** GAO YiChao 2025.1104 *****/
#ifndef __XTDB_VECTOR_LIST_H__
#define __XTDB_VECTOR_LIST_H__

#include <XiaoTuDataBox/Utils.hpp>

#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>

namespace xiaotu {

    template <typename DataType>
    class VectorList {
        public:
            /**
             * @brief 数组链表的一个节点
             */
            class Node {
                friend class VectorList;

                public:
                    /**
                     * @brief 拷贝 d
                     * 
                     * @param [in] d 目标数据
                     * @param [in] vl 构造本对象的容器，vl->CreateNode
                     * @param [in] idx 成功构造后的索引
                     */
                    Node(DataType const & d, VectorList *vl, size_t idx)
                        : mData(d), mContainer(vl),
                          mIndex(idx), mPrev(idx), mNext(idx)
                    {}

                    /**
                     * @brief 移动 d
                     * 
                     * @param [in] d 目标数据
                     * @param [in] vl 构造本对象的容器，vl->CreateNode
                     * @param [in] idx 成功构造后的索引
                     */
                    Node(DataType && d, VectorList * vl, size_t idx)
                        : mData(std::move(d)), mContainer(vl),
                          mIndex(idx), mPrev(idx), mNext(idx)
                    {}

                    /**
                     * @brief 以指针或者迭代器的形式解引用
                     */
                    DataType & operator * () { return mData; }
                    DataType const & operator * () const { return mData; }

                    /**
                     * @brief 以指针或者迭代器的形式解引用
                     */
                    DataType * operator -> () { return &mData; }
                    DataType const * operator -> () const { return &mData; }

                    /**
                     * @brief 获取当前节点的后继
                     */
                    Node & Next() { return (*mContainer)[mNext]; }
                    Node const & Next() const { return (*mContainer)[mNext]; }

                    /**
                     * @brief 获取当前节点的前驱
                     */
                    Node & Prev() { return (*mContainer)[mPrev]; }
                    Node const & Prev() const { return (*mContainer)[mPrev]; }

                    /**
                     * @brief 将 _new 插到 prev 和 next 之间
                     */
                    static void Insert(Node & _new, Node & prev, Node & next)
                    {
                        next.mPrev = _new.mIndex;
                        _new.mNext = next.mIndex;
                        _new.mPrev = prev.mIndex;
                        prev.mNext = _new.mIndex;
                    }

                    void InsertNext(Node & _new)
                    {
                        Insert(_new, *this, Next());
                    }

                    void InsertPrev(Node & _new)
                    {
                        Insert(_new, Prev(), *this);
                    }

                    /**
                     * @brief 调试用的字符串
                     */
                    std::string ToString() const
                    {
                        std::stringstream ss;
                        ss << mPrev << "[" << mIndex << "]" << mNext;
                        return ss.str();
                    }

                private:
                    //! 容器内容
                    DataType mData;
                    //! 管理 this 的容器对象
                    VectorList * mContainer{nullptr};
                    //! 在 mContainer->mNodes 中的索引
                    size_t mIndex{0};
                    //! 前驱节点的索引
                    size_t mPrev{0};
                    //! 后继节点的索引
                    size_t mNext{0};
            };


            /**
             * @brief 默认构造一个空的容器
             */
            VectorList()
            {}

            /**
             * @brief 拷贝构造一个节点
             */
            Node & CreateNode(DataType const & d)
            {
                mNodes.emplace_back(d, this, mNodes.size());
                return mNodes.back();
            }

            /**
             * @brief 移动构造一个节点
             */
            Node & CreateNode(DataType && d)
            {
                mNodes.emplace_back(std::move(d), this, mNodes.size());
                return mNodes.back();
            }

            /**
             * @brief 在 head 之前增加一个节点
             * 
             * @param [in] d 目标数据
             * @param [in] head 应当是 mVertices 中的一个成员，若为 nullptr 则是 mVertices[0]
             * @return 新增的节点对象
             */
            Node & Add(DataType const & d, Node * head = nullptr)
            {
                auto & node = CreateNode(d);
                return __InsertPrev__(node, head);
            }

            Node & Add(DataType && d, Node * head = nullptr)
            {
                auto & node = CreateNode(std::move(d));
                return __InsertPrev__(node, head);
            }

            /**
             * @brief 接口 Add 的辅助函数，将 n 插到 head 的前面
             */
            Node & __InsertPrev__(Node & n, Node * head = nullptr)
            {
                if (1 == mNodes.size())
                    return n;
                if (nullptr == head)
                    head = &mNodes[0];
                head->InsertPrev(n);
                return n;
            }

            Node & operator [] (int idx) { return mNodes[idx]; }
            Node const & operator [] (int idx) const { return mNodes[idx]; }

            std::vector<Node> & Nodes() { return mNodes; }
            std::vector<Node> const & Nodes() const { return mNodes; }

        private:
            //! @brief 节点列表
            std::vector<Node> mNodes;
    };


}


#endif

