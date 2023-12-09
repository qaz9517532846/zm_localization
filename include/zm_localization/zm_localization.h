#ifndef ZM_LOCALIZATION_H_
#define ZM_LOCALIZATION_H_

#include <ros/ros.h>

#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <geometry_msgs/PoseArray.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/OccupancyGrid.h>
#include <sensor_msgs/LaserScan.h>

#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#include <mutex>
#include <thread>
#include <angles/angles.h>
#include <random>

#include <zm_localization/zm_solver.h>
#include <zm_localization/zm_gridMap.h>

typedef enum
{
    LOCALIZATION_0D,
    LOCALIZATION_1D,
    LOCALIZATION_2D,
    LOCALIZATION_3D
} LOCALIZATION_MODE;

class zmLocalization
{
    public:
        zmLocalization(std::string name);
        ~zmLocalization();

    private:
        void MapCB(const nav_msgs::OccupancyGrid::ConstPtr& map);
        void ScanCB(const sensor_msgs::LaserScan::ConstPtr &scan, std::string topic);
        void broadcast();
        void InitialPoseCB(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& pose);
        void LocalizationUpdate(const ros::TimerEvent& event);
        std::vector<ScanPoint_t> ConvertScan(const sensor_msgs::LaserScan::ConstPtr& scan, const Eigen::Matrix4d& odomToBase);
        void PublishPointCloud2Msg(std::vector<ScanPoint_t> points);
        void UpdateMap();
        void UpdateLoop();
        LOCALIZATION_MODE Localization_Mode();

        ros::NodeHandle nh_;
        ros::Subscriber initialPoseSub;
        std::vector<ros::Subscriber> scanSub;
        ros::Subscriber mapSub;

        ros::Publisher mapTilePub;
        ros::Publisher localPosePub;
        ros::Publisher localPose2Pub;
        ros::Publisher poseArrayPub;
        ros::Publisher scanMergerPub;


        ros::Timer updateTimer_;

        std::shared_ptr<GridMap<float>> mMap;	  // map tile
        nav_msgs::OccupancyGrid::ConstPtr mWorld; // whole map

        std::mt19937 mGenerator;
        std::mutex nodeMutex_;
        std::thread mapUpdateThread_;
	    tf::TransformListener tf_;
	    tf::TransformBroadcaster tfPub_;

        std::map<std::string, sensor_msgs::LaserScan::ConstPtr> scanBuffer;

        Eigen::Matrix4d mWorldToMap;
        Eigen::Matrix4d mGridToMap;

        //ros parameter
        bool broadcastTf_;
        std::string baseFrame_;
        std::string odomFrame_;
        std::string mapFrame_;
        int scanNum_;
        int mapSize_;
        int mapDownScale_;
        double mapUpdateRate_;
        double locUpdateRate_;
        int numSmooth_;
        double minScore_;
        double solverGain_;
        double solverDamping_;
        int solverIterations_;
        int sampleRate_;
        int minPoints_;
        int updateCnt = 0;
        double updateGain;
        double confidenceGain;
        double odometryStdXY_;
        double odometryStdYaw_;
        double minSampleStdXY_;
        double minSampleStdYaw_;
        double maxSampleStdXY_;
        double maxSampleStdYaw_;
        double constrainThreshold_;
        double constrainThresholdYaw_;
        double transformTimeout_;

        ros::Time offsetTime;
        double offsetX = 0;					// current x offset between odom and map
	    double offsetY = 0;					// current y offset between odom and map
	    double offsetYaw = 0;				// current yaw offset between odom and map
        double sampleStdXY = 0;				// current sample spread in xy
        double sampleStdYaw = 0;			// current sample spread in yaw
        Eigen::Matrix<double, 4, 1> lastOdomPose;
        Solver mSolver;
};
#endif