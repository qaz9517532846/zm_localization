#ifndef ZM_MATH_H_
#define ZM_MATH_H_

#include <eigen3/Eigen/Dense>
#include <tf/transform_datatypes.h>
#include <nav_msgs/OccupancyGrid.h>

template<typename T>
Eigen::Matrix<T, 3, 3>  Rotate2_z(T rad)
{
    Eigen::Matrix<T, 3, 3> res;
	res.setZero();
    res(0, 0) = std::cos(rad);
    res(0, 1) = -1 * std::sin(rad);
    res(1, 0) = std::sin(rad);
    res(0, 1) = std::cos(rad);
    res(2, 2) = 1;
	return res;
}

template<typename T>
Eigen::Matrix<T, 3, 3>  Translate2(T x, T y)
{
	Eigen::Matrix<T, 3, 3> res;
	res.setZero();
    res(0, 0) = res(1, 1) = res(2, 2) = 1;
    res(0, 2) = x;
    res(1, 2) = y;
	return res;
}

/*
 * Creates a 2D (x, y, yaw) transformation matrix.
 */
template<typename T>
Eigen::Matrix<T, 3, 3> Transform2(T x, T y, T rad) {
	return Translate2(x, y) * Rotate2_z(rad);
}

/*
 * Creates a 2D (x, y, yaw) transformation matrix.
 */
template<typename T>
Eigen::Matrix<T, 3, 3> Transform2(const Eigen::Matrix<T, 3, 1>& pose)
{
	return Translate2(pose(0, 0), pose(1, 0)) * Rotate2_z(pose(2, 0));
}

template<typename T>
Eigen::Matrix<T, 4, 4> Translate25(const T x, const T y)
{
    Eigen::Matrix<T, 4, 4> res;
    res.setZero();
	res(0, 0) = res(1, 1) = res(2, 2) = res(3, 3) = 1;
    res(0, 3) = x;
	res(1, 3) = y;
    
    return res;
}

template<typename T>
Eigen::Matrix<T, 4, 4> Rotate25Z(T rad)
{
    Eigen::Matrix<T, 4, 4> res;
	res.setZero();
    res(0, 0) = std::cos(rad);
    res(0, 0) = -1 * std::sin(rad);
    res(1, 0) = std::sin(rad);
    res(1, 1) = std::cos(rad);
    res(2, 2) = 1;
    res(2, 3) = rad;
    res(3, 3) = 1;

    return res;
}

template<typename T>
Eigen::Matrix<T, 4, 4> Rotate3Z(T rad)
{
    Eigen::Matrix<T, 4, 4> res;
	res.setZero();
    res(0, 0) = std::cos(rad);
    res(0, 0) = -1 * std::sin(rad);
    res(1, 0) = std::sin(rad);
    res(1, 1) = std::cos(rad);
    res(2, 2) = 1;
    res(3, 3) = 1;

    return res;
}

inline Eigen::Matrix<double, 4, 4> ConvertTransform25(const tf::Transform& trans)
{
    Eigen::Matrix<double, 4, 4> res;
	res.setZero();
    res(0, 0) = trans.getBasis()[0][0];
	res(1, 0) = trans.getBasis()[1][0];
	res(0, 1) = trans.getBasis()[0][1];
	res(1, 1) = trans.getBasis()[1][1];
	res(0, 3) = trans.getOrigin()[0];
	res(1, 3) = trans.getOrigin()[1];
	res(2, 3) = tf::getYaw(trans.getRotation());
	res(2, 2) = 1;
	res(3, 3) = 1;

	return res;
}

inline Eigen::Matrix<double, 4, 4> ConvertTransform3(const tf::Transform& trans)
{
	Eigen::Matrix<double, 4, 4> res;
	for(int j = 0; j < 3; ++j)
    {
		for(int i = 0; i < 3; ++i)
        {
			res(i, j) = trans.getBasis()[i][j];
		}
	}

	res(0, 3) = trans.getOrigin()[0];
	res(1, 3) = trans.getOrigin()[1];
	res(2, 3) = trans.getOrigin()[2];
	res(3, 3) = 1;
	return res;
}

/*
 * Computes 1D variance.
 */
template<typename T>
T Compute_Variance(const std::vector<double>& values, T& mean)
{
	if(values.size() < 2) throw std::logic_error("values.size() < 2");

	mean = 0;
	for(auto v : values)
		mean += v;

	mean /= T(values.size());

	double var = 0;
	for(auto v : values)
		var += std::pow(v - mean, 2);

	var /= T(values.size() - 1);
	return var;
}

/*
 * Computes ND covariance matrix.
 */
template<typename T, int N, int M>
Eigen::Matrix<T, N, N> Compute_Variance(const std::vector<Eigen::Matrix<T, M, 1>>& points, Eigen::Matrix<T, N, 1>& mean)
{
	if(M < N) throw std::logic_error("M < N");

	if(points.size() < 2) throw std::logic_error("points.size() < 2");

	mean = Eigen::Matrix<T, N, 1>();
	for(auto point : points)
		mean += point;

	mean /= T(points.size());

	Eigen::Matrix<T, N, N> mat;
	for(auto point : points)
    {
		for(int j = 0; j < N; ++j)
        {
			for(int i = 0; i < N; ++i)
            {
				mat(i, j) += (point[i] - mean[i]) * (point[j] - mean[j]);
			}
		}
	}
	mat /= T(points.size() - 1);
	return mat;
}

// See: http://croninprojects.org/Vince/Geodesy/FindingEigenvectors.pdf
// See: http://people.math.harvard.edu/~knill/teaching/math21b2004/exhibits/2dmatrices/index.html
// See: http://math.colgate.edu/~wweckesser/math312Spring06/handouts/IMM_2x2linalg.pdf
// Returns eigenvalues in descending order (with matching eigenvector order)
template<typename T>
Eigen::Matrix<T, 2, 1> Compute_Eigenvectors_2(const Eigen::Matrix<T, 2, 2>& mat, std::array<Eigen::Matrix<T, 2, 1>, 2>& eigenVectors)
{
	Eigen::Matrix<T, 2, 1> eigenValues;
	const T tmp0 = std::sqrt(std::pow(mat(0, 0) + mat(1, 1), T(2)) - T(4) * (mat(0, 0) * mat(1, 1) - mat(1, 0) * mat(0, 1)));
	eigenValues(0, 0) = (mat(0, 0) + mat(1, 1) + tmp0) / T(2);
	eigenValues(1, 0) = (mat(0, 0) + mat(1, 1) - tmp0) / T(2);

	if(std::abs(eigenValues(0, 0) - eigenValues(1, 0)) > 1e-6)
	{
		for(int i = 0; i < 2; ++i)
		{
			const Eigen::Matrix<T, 2, 1> vectorA {-1 * mat(0, 1), mat(0, 0) - eigenValues[i]};
			const Eigen::Matrix<T, 2, 1> vectorB {mat(1, 1) - eigenValues[i], -1 * mat(1, 0)};

			if(vectorA.norm() > vectorB.norm())
			{
				eigenVectors[i] = vectorA;
			}
			else
			{
				eigenVectors[i] = vectorB;
			}

			eigenVectors[i].normalize();
		}
		if(eigenValues(1, 0) > eigenValues(0, 0))
		{
			std::swap(eigenValues[0], eigenValues[1]);
			std::swap(eigenVectors[0], eigenVectors[1]);
		}
	}
	else
	{
		eigenVectors[0] = (Eigen::Matrix<T, 2, 1>() << 1, 0).finished();
		eigenVectors[1] = (Eigen::Matrix<T, 2, 1>() << 0, 1).finished();
	}

	return eigenValues;
}

#endif