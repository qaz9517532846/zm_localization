#include <zm_localization/zm_localization.h>

int main(int argc, char** argv)
{
  ros::init(argc, argv, "zm_localization");

  zmLocalization zm_localization(ros::this_node::getName());

  ros::spin();
  return 0;
}