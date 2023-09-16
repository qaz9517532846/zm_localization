#include <zm_localization/zm_localization.h>

zmLocalization::zmLocalization(std::string name) : nh_("~")
{
    nh_.param<bool>("broadcast_tf", broadcastTf_, true);
	nh_.param<std::string>("base_frame", baseFrame_, "base_link");
	nh_.param<std::string>("odom_frame", odomFrame_, "odom");
	nh_.param<std::string>("map_frame", mapFrame_, "map");
	nh_.param<int>("scan_num", scanNum_, 1);
	nh_.param<int>("map_size", mapSize_, 1000);
	nh_.param<int>("map_downscale", mapDownScale_, 0);
	nh_.param<double>("map_update_rate", mapUpdateRate_, 0.5);
	nh_.param<double>("loc_update_rate", locUpdateRate_, 10.);
	nh_.param<int>("num_smooth", numSmooth_, 5);
	nh_.param<double>("min_score", minScore_, 0.2);
	nh_.param<double>("solver_gain", mSolver.gain, 0.1);
	nh_.param<double>("solver_damping", mSolver.damping, 1000.);
	nh_.param<int>("solver_iterations", solverIterations_, 20);
	nh_.param<int>("sample_rate", sampleRate_, 10);
	nh_.param<int>("min_points", minPoints_, 20);
	nh_.param<double>("update_gain", updateGain, 0.5);
	nh_.param<double>("confidence_gain", confidenceGain, 0.01);
	nh_.param<double>("odometry_std_xy", odometryStdXY_, 0.01);
	nh_.param<double>("odometry_std_yaw", odometryStdYaw_, 0.01);
	nh_.param<double>("min_sample_std_xy", minSampleStdXY_, 0.025);
	nh_.param<double>("min_sample_std_yaw", minSampleStdYaw_, 0.025);
	nh_.param<double>("max_sample_std_xy", maxSampleStdXY_, 0.5);
	nh_.param<double>("max_sample_std_yaw", maxSampleStdYaw_, 0.5);
	nh_.param<double>("constrain_threshold", constrainThreshold_, 0.1);
	nh_.param<double>("constrain_threshold_yaw", constrainThresholdYaw_, 0.2);
	nh_.param<double>("transform_timeout", transformTimeout_, 0.2);

	scanSub.resize(scanNum_);
	for(int i = 0; i < scanNum_; i++)
	{
		std::string ScanTopicName;
		ScanTopicName = "/scan_" + std::to_string(i);
		scanSub[i] = nh_.subscribe<sensor_msgs::LaserScan>(ScanTopicName.c_str(), 10, boost::bind(&zmLocalization::ScanCB, this, _1, ScanTopicName.c_str()));
	}

    initialPoseSub = nh_.subscribe("/initialpose", 1, &zmLocalization::InitialPoseCB, this);
	mapSub = nh_.subscribe("/map", 1, &zmLocalization::MapCB, this);

	mapTilePub = nh_.advertise<nav_msgs::OccupancyGrid>("/map_tile", 1);
	localPosePub = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("/amcl_pose", 10);
	localPose2Pub = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("/map_pose", 10);
	poseArrayPub = nh_.advertise<geometry_msgs::PoseArray>("/particlecloud", 10);

	updateTimer_ = nh_.createTimer(ros::Rate(locUpdateRate_), &zmLocalization::LocalizationUpdate, this);
	mapUpdateThread_ = std::thread(&zmLocalization::UpdateLoop, this);
}

zmLocalization::~zmLocalization()
{
	if(mapUpdateThread_.joinable()) mapUpdateThread_.join();
}

void zmLocalization::InitialPoseCB(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& pose)
{
	std::lock_guard<std::mutex> lock(nodeMutex_);

	if(pose->header.frame_id != mapFrame_)
	{
		ROS_WARN_STREAM("ZM Localization Node: Invalid pose estimate frame: " << pose->header.frame_id);
		return;
	}

	tf::Transform mapPose;
	tf::poseMsgToTF(pose->pose.pose, mapPose);

	ROS_INFO_STREAM("ZM Localization Node: Got new map pose estimate: x=" << mapPose.getOrigin()[0]
					<< " m, y=" <<  mapPose.getOrigin()[1] << " m, yaw=" << tf::getYaw(mapPose.getRotation()) << " rad");

	tf::StampedTransform baseToOdom;
	try
	{
		tf_.lookupTransform(odomFrame_, baseFrame_, ros::Time(), baseToOdom);
	}
	catch(const std::exception& ex)
	{
		ROS_WARN_STREAM("ZM Localization Node: lookupTransform(baseFrame, odomFrame) failed: " << ex.what());
		return;
	}

	const Eigen::Matrix4d L = ConvertTransform25(baseToOdom);
	const Eigen::Matrix<double, 3, 1> newOffset = (ConvertTransform25(mapPose) * L.inverse()).block<3, 1>(0, 3);

	// set new offset based on given position
	offsetX = newOffset(0, 0);
	offsetY = newOffset(1, 0);
	offsetYaw = newOffset(2, 0);

	// reset particle spread to maximum
	sampleStdXY = maxSampleStdXY_;
	sampleStdYaw = maxSampleStdYaw_;

	broadcast();

	UpdateMap();
}

void zmLocalization::broadcast()
{
	if(broadcastTf_)
	{
		// compose and publish transform for tf package
		geometry_msgs::TransformStamped pose;
		// compose header
		pose.header.stamp = offsetTime;
		pose.header.frame_id = mapFrame_;
		pose.child_frame_id = odomFrame_;
		// compose data container
		pose.transform.translation.x = offsetX;
		pose.transform.translation.y = offsetY;
		pose.transform.translation.z = 0;
		tf::quaternionTFToMsg(tf::createQuaternionFromYaw(offsetYaw), pose.transform.rotation);
		tfPub_.sendTransform(pose);
	}
}

void zmLocalization::MapCB(const nav_msgs::OccupancyGrid::ConstPtr& map)
{
	std::lock_guard<std::mutex> lock(nodeMutex_);
	
	ROS_INFO_STREAM("ZM Localization Node: Got new map with dimensions " << map->info.width << " x " << map->info.height
					<< " and cell size " << map->info.resolution);

	tf::Transform tmp;
	tf::poseMsgToTF(map->info.origin, tmp);
	mWorldToMap = ConvertTransform25(tmp);
	mWorld = map;
	sampleStdXY = maxSampleStdXY_;
	sampleStdYaw = maxSampleStdYaw_;
}

void zmLocalization::ScanCB(const sensor_msgs::LaserScan::ConstPtr &scan, std::string topic)
{
	std::lock_guard<std::mutex> lock(nodeMutex_);

	if(!mMap) return;

	scanBuffer[scan->header.frame_id] = scan;
}

std::vector<ScanPoint_t> zmLocalization::ConvertScan(const sensor_msgs::LaserScan::ConstPtr& scan, const Eigen::Matrix4d& odomToBase)
{
	std::vector<ScanPoint_t> points;
	tf::StampedTransform sensorToBase, baseToOdom;

	try
	{
		tf_.lookupTransform(baseFrame_, scan->header.frame_id, ros::Time(0), sensorToBase);
	}
	catch(const std::exception& ex)
	{
		ROS_WARN_STREAM("ZMLocalizationNode: lookupTransform(scan->header.frame_id, m_base_frame) failed: " << ex.what());
		return points;
	}

	try
	{
		tf_.waitForTransform(odomFrame_, baseFrame_, scan->header.stamp, ros::Duration(transformTimeout_));
		tf_.lookupTransform(odomFrame_, baseFrame_, scan->header.stamp, baseToOdom);
	}
	catch(const std::exception& ex)
	{
		ROS_WARN_STREAM("ZMLocalizationNode: lookupTransform(baseFrame_, odomFrame_) failed: " << ex.what());
		return points;
	}

	const Eigen::Matrix4d S = ConvertTransform3(sensorToBase);
	const Eigen::Matrix4d L = ConvertTransform25(baseToOdom);

	// precompute transformation matrix from sensor to requested base
	const Eigen::Matrix4d T = odomToBase * L * S;

	for(size_t i = 0; i < scan->ranges.size(); ++i)
	{
		if(scan->ranges[i] <= scan->range_min || scan->ranges[i] >= scan->range_max)
		{
			continue;	// no actual measurement
		}

		// transform sensor points into base coordinate system
		const Eigen::Matrix<double, 4, 1> scanRange{scan->ranges[i], 0, 0, 1};
		const Eigen::Matrix<double, 4, 1> scanPos = T * Rotate3Z<double>(scan->angle_min + i * scan->angle_increment) * scanRange;
		ScanPoint_t point;
		point.x = scanPos(0, 0);
		point.y = scanPos(1, 0);
		points.emplace_back(point);
	}
	
	return points;
}

void zmLocalization::LocalizationUpdate(const ros::TimerEvent& event)
{
	std::lock_guard<std::mutex> lock(nodeMutex_);

	if(!mMap || scanBuffer.empty()) return;

	tf::StampedTransform baseToOdom;
	try
	{
		tf_.lookupTransform(odomFrame_, baseFrame_, ros::Time(0), baseToOdom);
	}
	catch(const std::exception& ex)
	{
		ROS_WARN_STREAM("ZMLocalizationNode: lookupTransform(baseFrame_, odomFrame_) failed: " << ex.what());
		return;
	}

	const Eigen::Matrix4d L = ConvertTransform25(baseToOdom);
	const Eigen::Matrix4d T = Translate25(offsetX, offsetY) * Rotate25Z(offsetYaw); // odom to map

	const Eigen::Matrix<double, 4, 1> odomPose = L.block<4, 1>(0, 3);
	const double distMoved = sqrt(pow(odomPose(0, 0) - lastOdomPose(0, 0), 2) + pow(odomPose(1, 0) - lastOdomPose(1, 0), 2));
	const double radRotated = fabs(angles::normalize_angle(odomPose(3, 0) - lastOdomPose(3, 0)));

	std::vector<ScanPoint_t> points;

	// convert all scans to current base frame
	for(const auto& scan : scanBuffer)
	{
		auto scanPoints = ConvertScan(scan.second, L.inverse());
		points.insert(points.end(), scanPoints.begin(), scanPoints.end());
	}
	
	// check for number of points
	if(points.size() < minPoints_)
	{
		ROS_WARN_STREAM("ZMLocalizationNode: Number of points too low: " << points.size());
		return;
	}

	auto poseArray = boost::make_shared<geometry_msgs::PoseArray>();
	poseArray->header.stamp = baseToOdom.stamp_;
	poseArray->header.frame_id = mapFrame_;

	// calc predicted grid pose based on odometry
	const Eigen::Matrix<double, 4, 1> gridPose = (mGridToMap.inverse() * T * L).block<4, 1>(0, 3);

	// setup distributions
	std::normal_distribution<double> distX(gridPose(0, 0), sampleStdXY);
	std::normal_distribution<double> distY(gridPose(1, 0), sampleStdXY);
	std::normal_distribution<double> distYaw(gridPose(2, 0), sampleStdYaw);

	// solve odometry prediction first
	mSolver.poseX = gridPose(0, 0);
	mSolver.poseY = gridPose(1, 0);
	mSolver.poseYaw = gridPose(2, 0);

	for(int iter = 0; iter < solverIterations_; ++iter)
	{
		mSolver.Solve<float>(*mMap, points);
	}

	double bestX = mSolver.poseX;
	double bestY = mSolver.poseY;
	double bestYaw = mSolver.poseYaw;
	double bestScore = mSolver.rNorm;
	//ROS_INFO("Best X = %f, y = %f, Yaw = %f, bestScore = %f\n", bestX, bestY, bestYaw, bestScore);
	
	std::vector<Eigen::Matrix<double, 4, 1>> seeds(sampleRate_);
	std::vector<Eigen::Matrix<double, 4, 1>> samples(sampleRate_);
	std::vector<double> sampleErrors(sampleRate_);

	for(int i = 0; i < sampleRate_; ++i)
	{
		// generate new sample
		mSolver.poseX = distX(mGenerator);
		mSolver.poseY = distY(mGenerator);
		mSolver.poseYaw = distYaw(mGenerator);

		seeds[i] = (Eigen::Matrix<double, 4, 1>() << mSolver.poseX, mSolver.poseY, mSolver.poseYaw, 1).finished();

		// solve sample
		for(int iter = 0; iter < solverIterations_; ++iter)
		{
			mSolver.Solve<float>(*mMap, points);
		}

		// save sample
		const auto sample = (Eigen::Matrix<double, 4, 1>() << mSolver.poseX, mSolver.poseY, mSolver.poseYaw, 1).finished();
		samples[i] = sample;
		sampleErrors[i] = mSolver.rNorm;

		// check if sample is better
		if(mSolver.rNorm > bestScore)
		{
			bestX = mSolver.poseX;
			bestY = mSolver.poseY;
			bestYaw = mSolver.poseYaw;
			bestScore = mSolver.rNorm;
		}

		// add to visualization
		const Eigen::Matrix<double, 4, 1> mapPose = mGridToMap * sample;
		tf::Pose pose;
		pose.setOrigin(tf::Vector3(mapPose(0, 0), mapPose(1, 0), 0));
		pose.setRotation(tf::createQuaternionFromYaw(mapPose(2, 0)));
		geometry_msgs::Pose tmp;
		tf::poseTFToMsg(pose, tmp);
		poseArray->poses.push_back(tmp);
	}
	
	//compute covariances
	double meanScore = 0;
	Eigen::Matrix<double, 4, 1> meanXYW;
	Eigen::Matrix<double, 4, 1> seedMeanXYW;
	const double varError = Compute_Variance(sampleErrors, meanScore);
	const Eigen::Matrix<double, 4, 4> varXYW = Compute_Variance(samples, meanXYW);
	const Eigen::Matrix<double, 3, 3> gradVarXYW = ComputeVirtualScanCovarianceXYW(mMap, points, (Eigen::Matrix<double, 3, 1>() << bestX, bestY, bestYaw).finished());

	// compute gradient characteristic
	std::array<Eigen::Matrix<double, 2, 1>, 2> gradEigenVectors;
	Eigen::Matrix<double, 2, 2> mat = gradVarXYW.block<2, 2>(0, 0);
	const Eigen::Matrix<double, 2, 1> gradEigenValues = Compute_Eigenvectors_2(mat, gradEigenVectors);
	const Eigen::Matrix<double, 3, 1> gradStdUVW = (Eigen::Matrix<double, 3, 1>() << sqrt(gradEigenValues(0, 0)),
																					 sqrt(gradEigenValues(1, 0)),
																					 sqrt(gradVarXYW(2, 2))).finished();
	
	// decide if we have 3D, 2D, 1D or 0D localization
	LOCALIZATION_MODE mode = LOCALIZATION_0D;
	if(bestScore > minScore_)
	{
		if(gradStdUVW(0, 0) > constrainThreshold_)
		{
			if(gradStdUVW(1, 0) > constrainThreshold_)
				mode = LOCALIZATION_3D;	// 2D position + rotation
			else if(gradStdUVW(2, 0) > constrainThresholdYaw_)
				mode = LOCALIZATION_2D;	// 1D position + rotation
			else
				mode = LOCALIZATION_1D;	// 1D position only	
		}
	}

	double newGridX = bestX;
	double newGridY = bestY;
	double newGridYaw = bestYaw;
	switch(mode)
	{
		case LOCALIZATION_0D:
		break;
		case LOCALIZATION_1D:
		case LOCALIZATION_2D:
			newGridYaw = gridPose(2, 0);	// keep old orientation
		case LOCALIZATION_3D:
			// constrain update to the good direction (ie. in direction of the eigen vector with the smaller sigma)
			const auto delta = Eigen::Matrix<double, 2, 1>{bestX, bestY} - Eigen::Matrix<double, 2, 1>{gridPose(0, 0), gridPose(1, 0)};
			const auto dist = gradEigenVectors[0].dot(delta);
			newGridX = gridPose(0, 0) + dist * gradEigenVectors[0](0, 0);
			newGridY = gridPose(1, 0) + dist * gradEigenVectors[0](1, 0);
	}

	
	if(mode > LOCALIZATION_0D)
	{
		// use best sample for update
		Eigen::Matrix<double, 4, 4> gridPoseNew = Translate25(newGridX, newGridY) * Rotate25Z(newGridYaw);

		// compute new odom to map offset from new grid pose
		const Eigen::Matrix<double, 4, 1> newOffset = (mGridToMap * gridPoseNew * L.inverse()).block<4, 1>(0, 3);

		// apply new offset with an exponential low pass filter
		offsetX += (newOffset(0, 0) - offsetX) * updateGain;
		offsetY += (newOffset(1, 0) - offsetY) * updateGain;
		offsetYaw += angles::shortest_angular_distance(offsetYaw, newOffset(2, 0)) * updateGain;
	}

	offsetTime = baseToOdom.stamp_;

	// update particle spread depending on mode
	if(mode >= LOCALIZATION_3D)
		sampleStdXY *= (1 - confidenceGain);
	else
		sampleStdXY += distMoved * odometryStdXY_;

	if(mode >= LOCALIZATION_2D)
		sampleStdYaw *= (1 - confidenceGain);
	else
		sampleStdYaw += radRotated * odometryStdYaw_;
		
	// limit particle spread
	sampleStdXY = fmin(fmax(sampleStdXY, minSampleStdXY_), maxSampleStdXY_);
	sampleStdYaw = fmin(fmax(sampleStdYaw, minSampleStdYaw_), maxSampleStdYaw_);

	// publish new transform
	broadcast();

	
	const Eigen::Matrix<double, 3, 1> newMapPose = (Translate25(offsetX, offsetY) * Rotate25Z(offsetYaw) * L).block<3, 1>(0, 3);

	// publish localization pose
	auto localPose = boost::make_shared<geometry_msgs::PoseWithCovarianceStamped>();
	localPose->header.stamp = offsetTime;
	localPose->header.frame_id = mapFrame_;
	localPose->pose.pose.position.x = newMapPose(0, 0);
	localPose->pose.pose.position.y = newMapPose(1, 0);
	localPose->pose.pose.position.z = 0;
	tf::quaternionTFToMsg(tf::createQuaternionFromYaw(newMapPose(2, 0)), localPose->pose.pose.orientation);
	for(int j = 0; j < 3; j++)
	{
		for(int i = 0; i < 3; ++i)
		{
			const int i_ = (i == 2 ? 5 : i);
			const int j_ = (j == 2 ? 5 : j);
			localPose->pose.covariance[j_ * 6 + i_] = varXYW(i, j);
		}
	}

	localPosePub.publish(localPose);
	localPose2Pub.publish(localPose);

	// publish visualization
	poseArrayPub.publish(poseArray);

	// keep last odom pose
	lastOdomPose = odomPose;

	if(updateCnt++ % 10 == 0)
	{
		updateCnt = 0;
		ROS_INFO_STREAM("ZMLocalizationNode: score = " << float(bestScore) << ", gradUVW = [" << float(gradStdUVW(0, 0)) << ", " << float(gradStdUVW(1, 0))
					<< ", " << float(gradStdUVW(2, 0)) << "], stdXY=" << float(sampleStdXY) << " m, stdYaw=" << float(sampleStdYaw)
					<< " rad, mode=" << mode << "D, " << scanBuffer.size() << " scans");
	}

	// clear scan buffer
	scanBuffer.clear();

}

void zmLocalization::UpdateMap()
{
	Eigen::Matrix4d worldToMap;
	Eigen::Matrix<double, 3, 1> worldPose;
	tf::StampedTransform baseToOdom_;
	nav_msgs::OccupancyGrid::ConstPtr world;

	{
		std::lock_guard<std::mutex> lock(nodeMutex_);
	
		if(!mWorld) return;

		try
		{
			tf_.lookupTransform(odomFrame_, baseFrame_, ros::Time(0), baseToOdom_);
		}
		catch(const std::exception& ex)
		{
			ROS_WARN_STREAM("ZM Localization Node: lookupTransform(baseFrame, odomFrame) failed: " << ex.what());
			return;
		}

		const Eigen::Matrix4d L = ConvertTransform25(baseToOdom_);
		const Eigen::Matrix4d T = Translate25(offsetX, offsetY) * Rotate25Z(offsetYaw);		// odom to map
		worldPose = (worldToMap.inverse() * T * L).block<3, 1>(0, 3);

		world = mWorld;
		worldToMap = mWorldToMap;
	}

	// compute tile origin in pixel coords
	const double worldScale = world->info.resolution;
	const int tileX = int(worldPose(0, 0) / worldScale) - mapSize_ / 2;
	const int tileY = int(worldPose(1, 0) / worldScale) - mapSize_ / 2;

	auto map = std::make_shared<GridMap<float>>(mapSize_, mapSize_, worldScale);

	// extract tile and convert to our format (occupancy between 0 and 1)
	for(int y = 0; y < map->sizeY(); ++y)
	{
		for(int x = 0; x < map->sizeX(); ++x)
		{
			const int x_ = std::min(std::max(tileX + x, 0), int(world->info.width) - 1);
			const int y_ = std::min(std::max(tileY + y, 0), int(world->info.height) - 1);
			const auto cell = world->data[y_ * world->info.width + x_];

			if(cell >= 0) (*map)(x, y) = fminf(cell / 100.f, 1.f);
			else (*map)(x, y) = 0;
		}
	}

	// optionally downscale map
	for(int i = 0; i < mapDownScale_; ++i)
		map = map->DownScale();

	// smooth map
	for(int i = 0; i < numSmooth_; ++i)
		map->Smooth33_1();

	// update map
	{
		std::lock_guard<std::mutex> lock(nodeMutex_);
		mMap = map;
		mGridToMap = worldToMap * Translate25<double>(tileX * worldScale, tileY * worldScale);
	}

	const Eigen::Matrix<double, 3, 1> tileOrigin = mGridToMap.block<3, 1>(0, 3);
	const Eigen::Matrix<double, 4, 1> mapCenter = (Eigen::Matrix<double, 4, 1>() << map->Scale() * map->sizeX() / 2,
																					map->Scale() * map->sizeY() / 2,
																					0, 1).finished();
	const Eigen::Matrix<double, 3, 1> tileCenter = (mGridToMap * mapCenter).block<3, 1>(0, 0);

	// publish new map tile for visualization
	auto rosGrid = boost::make_shared<nav_msgs::OccupancyGrid>();
	rosGrid->info.resolution = map->Scale();
	rosGrid->info.width = map->sizeX();
	rosGrid->info.height = map->sizeY();
	rosGrid->info.origin.position.x = tileOrigin(0, 0);
	rosGrid->info.origin.position.y = tileOrigin(1, 0);
	tf::quaternionTFToMsg(tf::createQuaternionFromYaw(tileOrigin(2, 0)), rosGrid->info.origin.orientation);
	rosGrid->data.resize(map->NumCells());
	for(int y = 0; y < map->sizeY(); ++y)
	{
		for(int x = 0; x < map->sizeX(); ++x)
		{
			rosGrid->data[y * map->sizeX() + x] = (*map)(x, y) * 100.f;
		}
	}
	
	mapTilePub.publish(rosGrid);
	ROS_DEBUG_STREAM("ZM Localization Node: Got new grid at offset (" << tileX << ", " << tileY << ") [iworld], "
					 "center = (" << tileCenter(0, 0) << ", " << tileCenter(1, 0) << ") [map]");
}

void zmLocalization::UpdateLoop()
{
	ros::Rate rate(mapUpdateRate_);
	while(ros::ok())
	{
		try
		{
			UpdateMap();	// get a new map tile periodically
		}
		catch(const std::exception& ex)
		{
			ROS_WARN_STREAM("ZM Localization Node: UpdateMap() failed: " << ex.what());
		}
			
		rate.sleep();
	}
}