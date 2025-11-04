#include <XiaoTuDataBox/VectorList.hpp>

#include <iostream>
#include <gtest/gtest.h>


TEST(VectorList, Node_double)
{
    using namespace xiaotu::data;

    VectorList<double> list;
    {
        double tmp{1.2};
        auto & node = list.Add(tmp);
        EXPECT_EQ(1, list.Nodes().size());
        std::cout << std::endl;
        std::cout << node.ToString() << "  " << *node << std::endl;
    }

    {
        VectorList<double>::Node & node = list.Add(3.14);
        EXPECT_EQ(2, list.Nodes().size());
        std::cout << std::endl;
        std::cout << node.ToString() << std::endl;
        std::cout << node.Next().ToString() << std::endl;
        std::cout << node.Prev().ToString() << std::endl;
    }

    {
        VectorList<double>::Node & node = list.Add(1.414);
        EXPECT_EQ(3, list.Nodes().size());
        std::cout << std::endl;
        std::cout << node.ToString() << std::endl;
        std::cout << node.Next().ToString() << std::endl;
        std::cout << node.Prev().ToString() << std::endl;
    }

}

struct TestStruct {
    double x;
    double y;
};


TEST(VectorList, Node_struct)
{
    using namespace xiaotu::data;

    VectorList<TestStruct> list;
    {
        TestStruct tmp{3.14, 1.5};
        auto & node = list.Add(tmp);
        EXPECT_EQ(1, list.Nodes().size());
        std::cout << node.ToString() << "  " << node->x << "," << node->y << std::endl;
    }

    {
        VectorList<TestStruct>::Node & node = list.Add({1.414, 3.15});
        EXPECT_EQ(2, list.Nodes().size());
        std::cout << node.ToString() << std::endl;
        std::cout << node.Next().ToString() << std::endl;
        std::cout << node.Prev().ToString() << std::endl;
    }

    {
        VectorList<TestStruct>::Node & node = list.Add({1.414, 3.15});
        EXPECT_EQ(3, list.Nodes().size());
        std::cout << node.ToString() << std::endl;
        std::cout << node.Next().ToString() << std::endl;
        std::cout << node.Prev().ToString() << std::endl;
    }
}

