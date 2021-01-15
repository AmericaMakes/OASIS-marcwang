#include "gtest/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(blaTest, test1)
{
    //arrange
    //act
    //assert
    EXPECT_EQ(0, 0);
}