#ifndef ZM_SOLVER_H_
#define ZM_SOLVER_H_

#include <zm_localization/zm_gridMap.h>
#include <zm_localization/zm_math.h>
#include <vector>
#include <eigen3/Eigen/Dense>

struct ScanPoint_t
{
	float x = 0;		// [m]
	float y = 0;		// [m]
};

struct ScanPointEx_t : public ScanPoint_t
{
	float w = 1;		// weight [1]
	int layer = 0;		// layer index
};


class Solver {
public:
	double poseX = 0;					// initial guess / output grid pose
	double poseY = 0;					// initial guess / output grid pose
	double poseYaw = 0;				// initial guess / output grid pose

	double gain = 0.1;					// how fast to converge (0 to 1)
	double damping = 1;					// numerical hessian damping
	double rNorm = 0;					// current error norm

	Eigen::Matrix<double, 3, 1> G;				// gradient vector
	Eigen::Matrix<double, 3, 3> H;				// Hessian matrix

	template<typename T>
	void Solve(const GridMap<T>& grid, const std::vector<ScanPoint_t> points)
	{
		Reset();

		// compute transformation matrix first
		const Eigen::Matrix<double, 3, 3> P = Transform2(poseX, poseY, poseYaw);

		for(const auto& point : points)
		{
			Eigen::Matrix<double, 3, 1> pointMatrix{point.x, point.y, 1};
			// transform sensor point to grid coordinates
			const auto q = P * pointMatrix;

			const float gridX = grid.WorldToGrid(q(0, 0));
			const float gridY = grid.WorldToGrid(q(1, 0));

			// compute error based on grid
			const float rI = grid.BilinearLookup(gridX, gridY);
			rNorm += rI * rI;

			// compute error gradient based on grid
			float dx, dy;
			grid.CalcGradient(gridX, gridY, dx, dy);

			Integrate(point.x, point.y, rI, dx, dy);
		}

		// we want average rNorm
		rNorm = sqrt(rNorm / points.size());

		Update();
	}

	template<typename T>
	void Solve(const MultiGridMap<T>& multiGrid,const std::vector<ScanPointEx_t>& points)
	{
		Reset();

		// compute transformation matrix first
		const Eigen::Matrix<double, 3, 3> P = Transform2(poseX, poseY, poseYaw);

		for(const auto& point : points)
		{
			Eigen::Matrix<double, 3, 1> pointMatrix{point.x, point.y, 1};
			auto& grid = multiGrid.layers[point.layer];

			// transform sensor point to grid coordinates
			const auto q = P * pointMatrix;
			const float gridX = grid.worldToGrid(q(0, 0));
			const float gridY = grid.worldToGrid(q(1, 0));

			// compute error based on grid
			const float rI = grid.BilinearLookup(gridX, gridY);
			rNorm += rI * rI;

			// compute error gradient based on grid
			float dx, dy;
			grid.CalGradient(gridX, gridY, dx, dy);

			Integrate(point.x, point.y, rI, dx, dy);
		}

		// we want average rNorm
		rNorm = sqrt(rNorm / points.size());

		Update();
	}

protected:
	void Reset()
	{
		Eigen::Matrix<double, 3, 1> G;
		Eigen::Matrix<double, 3, 3> H;
		rNorm = 0;
	}

	void Integrate(const float pX, const float pY, const float rI, const float dx, const float dy)
	{
		const float jX = dx * 1.f;
		const float jY = dy * 1.f;
		const float jYaw = dx * (-sinf(poseYaw) * pX - cosf(poseYaw) * pY) + dy * (cosf(poseYaw) * pX - sinf(poseYaw) * pY);

		// direct gradient vector summation
		G[0] += jX * rI;
		G[1] += jY * rI;
		G[2] += jYaw * rI;

		// direct Hessian matrix summation
		H(0, 0) += jX * jX;
		H(1, 1) += jY * jY;
		H(2, 2) += jYaw * jYaw;

		H(0, 1) += jX * jY;
		H(1, 0) += jX * jY;

		H(0, 2) += jX * jYaw;
		H(2, 0) += jX * jYaw;

		H(1, 2) += jY * jYaw;
		H(2, 1) += jY * jYaw;
	}

	void Update()
	{
		// add Hessian damping
		H(0, 0) += damping;
		H(1, 1) += damping;
		H(2, 2) += damping;

		// solve Gauss-Newton step
		const auto X = H.inverse() * G;

		// apply new solution with a gain (optimize max. rNorm)
		poseX += gain * X[0];
		poseY += gain * X[1];
		poseYaw += gain * X[2];
	}

};


/*
 * Computes a "virtual" covariance matrix based on second order gradients.
 * A higher covariance means a larger gradient, so the meaning of "covariance" is inverted here.
 * A higher gradient is better for localization accuracy.
 */
inline Eigen::Matrix<double, 3, 3> ComputeVirtualScanCovarianceXYW(std::shared_ptr<const GridMap<float>> grid, const std::vector<ScanPoint_t>& points, const Eigen::Matrix<double, 3, 1>& pose)
{
	Eigen::Matrix<double, 3, 3> varXYW;
	const Eigen::Matrix<double, 3, 3> P = Transform2(pose);	// pre-compute transformation matrix

	for(const auto& point : points)
	{
		// transform sensor point to grid coordinates
        const Eigen::Matrix<double, 3, 1> pointMatrix{point.x, point.y, 1};
		const auto q = P * pointMatrix;
		const float grid_x = grid->WorldToGrid(q[0]);
		const float grid_y = grid->WorldToGrid(q[1]);

		float ddx, ddy;
		grid->CalcGradient2(grid_x, grid_y, ddx, ddy);

		const float dir = pose[2] + atan2f(point.y, point.x);
		const float length = hypotf(point.x, point.y);
		const float ddyaw = (sinf(dir) * ddx + cosf(dir) * ddy) * length;

		varXYW(0, 0) += ddx * ddx;
		varXYW(1, 0) += ddx * ddy;
		varXYW(0, 1) += ddy * ddx;
		varXYW(1, 1) += ddy * ddy;
		varXYW(2, 2) += ddyaw * ddyaw;
	}
	varXYW *= 1. / points.size();
	return varXYW;
}



#endif /* INCLUDE_SOLVER_H_ */