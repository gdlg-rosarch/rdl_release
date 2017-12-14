#include <gtest/gtest.h>
#include "UnitTestUtils.hpp"

using namespace RobotDynamics;

class ReferenceFrameTest : public ::testing::Test
{
	protected:

		virtual void SetUp()
		{
			std::srand( time(NULL) );

			allFrames.push_back(root1.get());
			allFrames.push_back(frame1.get());
			allFrames.push_back(frame2.get());
			allFrames.push_back(frame3.get());

			allFrames.push_back(root2.get());
			allFrames.push_back(frame4.get());
			allFrames.push_back(frame5.get());
			allFrames.push_back(frame6.get());
			allFrames.push_back(frame7.get());
			allFrames.push_back(frame8.get());

			frames1.push_back(root1.get());
			frames1.push_back(frame1.get());
			frames1.push_back(frame2.get());
			frames1.push_back(frame3.get());

			frames2.push_back(root2.get());
			frames2.push_back(frame4.get());
			frames2.push_back(frame5.get());
			frames2.push_back(frame6.get());
			frames2.push_back(frame7.get());
			frames2.push_back(frame8.get());

			for(auto const &frame : allFrames)
			{
				frame->update();
			}
		}
		virtual void TearDown()
		{
			allFrames.clear();
			frames1.clear();
			frames2.clear();
		}

		std::shared_ptr<ReferenceFrame> root1 = ReferenceFrame::createARootFrame("root1");
		std::shared_ptr<ReferenceFrame> root2 = ReferenceFrame::createARootFrame("root2");

		std::shared_ptr<unit_test_utils::RandomUnchangingFrame> frame1 = unit_test_utils::RandomUnchangingFrame::create("frame1", root1.get(),1);
		std::shared_ptr<unit_test_utils::RandomUnchangingFrame> frame2 = unit_test_utils::RandomUnchangingFrame::create("frame2", frame1.get(),2);
		std::shared_ptr<unit_test_utils::RandomUnchangingFrame> frame3 = unit_test_utils::RandomUnchangingFrame::create("frame3", frame2.get(),3);

		std::shared_ptr<unit_test_utils::RandomlyChangingFrame> frame4 = unit_test_utils::RandomlyChangingFrame::create("frame4", root2.get(),1);
		std::shared_ptr<unit_test_utils::RandomlyChangingFrame> frame5 = unit_test_utils::RandomlyChangingFrame::create("frame5", frame4.get(),2);
		std::shared_ptr<unit_test_utils::RandomUnchangingFrame> frame6 = unit_test_utils::RandomUnchangingFrame::create("frame6", root2.get(),1);
		std::shared_ptr<unit_test_utils::RandomlyChangingFrame> frame7 = unit_test_utils::RandomlyChangingFrame::create("frame7", frame6.get(),2);
		std::shared_ptr<unit_test_utils::RandomlyChangingFrame> frame8 = unit_test_utils::RandomlyChangingFrame::create("frame8", frame7.get(),3);

		std::vector<ReferenceFrame*> allFrames;
		std::vector<ReferenceFrame*> frames1;
		std::vector<ReferenceFrame*> frames2;

		int nTests = 20;

	private:

};

TEST_F(ReferenceFrameTest, constructors)
{
	ReferenceFrame testFrame(*frame1.get());

	EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(testFrame.getTransformToRoot().toMatrix(),frame1->getTransformToRoot().toMatrix(),unit_test_utils::TEST_PREC));
	EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(testFrame.getInverseTransformToRoot().toMatrix(),frame1->getInverseTransformToRoot().toMatrix(),unit_test_utils::TEST_PREC));
	EXPECT_STREQ(testFrame.getName().c_str(),frame1->getName().c_str());
	EXPECT_EQ(testFrame.getIsBodyFrame(),frame1->getIsBodyFrame());
	EXPECT_EQ(testFrame.getIsWorldFrame(),frame1->getIsWorldFrame());
}

TEST_F(ReferenceFrameTest, testRootsHaveNullParent)
{
	EXPECT_TRUE(root1->getParentFrame() == nullptr);
	EXPECT_TRUE(root2->getParentFrame() == nullptr);

	try
	{
		std::shared_ptr<ReferenceFrame> frame(new ReferenceFrame("blah", nullptr, SpatialTransform(), false,1));
	}
	catch(ReferenceFrameException &e)
	{
		EXPECT_STREQ(e.what(),"You are not allowed to create a frame with parentFrame=nullptr. Only a root frame and the world frame may have parentFrame=nullptr");
	}
}

TEST_F(ReferenceFrameTest, testVectorOfFramesWithRoot)
{
    ReferenceFrame* worldFrame = ReferenceFrame::getWorldFrame().get();


    EXPECT_EQ(worldFrame->getFramesStartingWithRootEndingWithThis().size(),1);
}

TEST_F(ReferenceFrameTest, testWorldFramePointerStuff)
{
	const std::shared_ptr<ReferenceFrame> worldFrame1 = ReferenceFrame::getWorldFrame();
	const std::shared_ptr<ReferenceFrame> worldFrame2 = ReferenceFrame::getWorldFrame();

    EXPECT_EQ(worldFrame1.get(),worldFrame2.get());
}

TEST_F(ReferenceFrameTest, testRootFramesArentTheSame)
{
	EXPECT_FALSE(root1 == root2);
}

TEST_F(ReferenceFrameTest, testGetRootFrame)
{
	EXPECT_TRUE(frame2->getRootFrame() == root1.get());
	EXPECT_TRUE(frame7->getRootFrame() == frame5->getRootFrame());

	frame7.get()->verifyFramesHaveSameRoot(frame6.get());

	try
	{
		frame7->verifyFramesHaveSameRoot(frame1.get());
	}
	catch ( RobotDynamics::ReferenceFrameException &e)
	{
		EXPECT_STREQ("Frames frame1 and frame7 have mismatched roots!",e.what());
	}
}

TEST_F(ReferenceFrameTest, testCheckReferenceFramesMatch)
{
	try
	{
		frame2->checkReferenceFramesMatch(nullptr);
	}
	catch(ReferenceFrameException &e)
	{
		EXPECT_STREQ(e.what(),"Reference frame is nullptr!");
	}

	try
	{
		frame2->checkReferenceFramesMatch(frame7.get());
	}
	catch(ReferenceFrameException &e)
	{
		EXPECT_STREQ(e.what(),"Reference frames do not match!");
	}
}

TEST_F(ReferenceFrameTest, testGetTransformBetweenFrames)
{

	for (int i = 0; i < nTests; i++)
	{
		unit_test_utils::updateAllFrames(allFrames);

		ReferenceFrame* tmpFrame1 = unit_test_utils::getARandomFrame(frames1);
		ReferenceFrame* tmpFrame2 = unit_test_utils::getARandomFrame(frames1);

		RobotDynamics::Math::SpatialTransform transform1 = tmpFrame1->getTransformToDesiredFrame(tmpFrame2);
		RobotDynamics::Math::SpatialTransform transform2 = tmpFrame2->getTransformToDesiredFrame(tmpFrame1);

		RobotDynamics::Math::SpatialTransform shouldBeIdentity = transform1 * transform2;
		RobotDynamics::Math::SpatialTransform identityTransform;

		EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(shouldBeIdentity.toMatrix(), identityTransform.toMatrix(), unit_test_utils::TEST_PREC));
	}

	for (int i = 0; i < nTests; i++)
	{
		unit_test_utils::updateAllFrames(allFrames);

		ReferenceFrame* tmpFrame1 = unit_test_utils::getARandomFrame(frames2);
		ReferenceFrame* tmpFrame2 = unit_test_utils::getARandomFrame(frames2);

		RobotDynamics::Math::SpatialTransform transform1 = tmpFrame1->getTransformToDesiredFrame(tmpFrame2);
		RobotDynamics::Math::SpatialTransform transform2 = tmpFrame2->getTransformToDesiredFrame(tmpFrame1);

		RobotDynamics::Math::SpatialTransform shouldBeIdentity = transform1 * transform2;
		RobotDynamics::Math::SpatialTransform identityTransform;

		EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(shouldBeIdentity.toMatrix(), identityTransform.toMatrix(), unit_test_utils::TEST_PREC));
	}
}

TEST_F(ReferenceFrameTest, testGetTransformToParent)
{
	for (int i = 1; i < allFrames.size(); i++)
	{
		ReferenceFrame* tmpFrame2 = allFrames[i];

		ReferenceFrame* parentFrame = tmpFrame2->getParentFrame();

		if (parentFrame != nullptr)
		{
			Eigen::Matrix4d m1,m2;
			RobotDynamics::Math::SpatialTransform t1 = tmpFrame2->getTransformToParent();
			RobotDynamics::Math::SpatialTransform t2 = tmpFrame2->getTransformToDesiredFrame(parentFrame);
			EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(t1.toMatrix(),t2.toMatrix(),unit_test_utils::TEST_PREC));
		}
	}
}

TEST_F(ReferenceFrameTest, testGetTransformToRoot)
{
	for (int j = 0; j < nTests; j++)
	{
		unit_test_utils::updateAllFrames(allFrames);

		for (int i = 0; i < allFrames.size(); i++)
		{
			ReferenceFrame* frame = allFrames[i];
			RobotDynamics::Math::SpatialTransform transformToRoot = unit_test_utils::getTransformToRootByClimbingTree(frame);

			EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(transformToRoot.toMatrix(),frame->getTransformToRoot().toMatrix(),unit_test_utils::TEST_PREC));
		}
	}
}

TEST_F(ReferenceFrameTest, testGetTransformToSelf)
{
	for (int i = 0; i < nTests; i++)
	{
		unit_test_utils::updateAllFrames(allFrames);

		for (int j = 0; j < allFrames.size(); j++)
		{
			ReferenceFrame* tmpFrame = allFrames[j];
			RobotDynamics::Math::SpatialTransform shouldBeIdentity = tmpFrame->getTransformToDesiredFrame(tmpFrame);

            RobotDynamics::Math::SpatialMatrix identityTransform = RobotDynamics::Math::SpatialMatrixIdentity;
            EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(shouldBeIdentity.toMatrix(),identityTransform,unit_test_utils::TEST_PREC));

		}
	}
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	return RUN_ALL_TESTS();
}