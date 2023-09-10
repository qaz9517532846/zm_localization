#ifndef ZM_GRID_MAP_H_
#define ZM_GRID_MAP_H_

#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdexcept>
#include <memory>
#include <vector>

template<typename T>
class GridMap
{
public:
    /*
	 * @param sizeX Size of map in pixels
	 * @param sizeY Size of map in pixels
	 * @param scale Size of one pixel in meters
	 */
	GridMap(int sizeX, int sizeY, float scale) :   mSizeX(sizeX),
                                                   mSizeY(sizeY),
			                                       mScale(scale),
			                                       mInvScale(1 / scale)
	{
		mMap = new T[size_t(sizeX) * sizeY];
	}

	/*
	 * Deep copy constructor.
	 */
	GridMap(const GridMap& other) : GridMap(other.mSizeX, other.mSizeY, other.mScale)
	{
		*this = other;
	}

	~GridMap()
	{
		delete [] mMap;
		mMap = 0;
	}

	/*
	 * Deep assignment operator.
	 */
	GridMap& operator=(const GridMap& other)
	{
		if(mSizeX != other.mSizeX || mSizeY != other.mSizeY)
        {
			throw std::logic_error("grid size mismatch");
		}

		mScale = other.mScale;
		mInvScale = other.mInvScale;
		::memcpy(mMap, other.mMap, NumCells() * sizeof(T));
		return *this;
	}

	int sizeX() const
    {
		return mSizeX;
	}

	int sizeY() const
    {
		return mSizeY;
	}

	float Scale() const
    {
		return mScale;
	}

	float InvScale() const
    {
		return mInvScale;
	}

	size_t NumCells() const
    {
		return size_t(mSizeX) * mSizeY;
	}

	/*
	 * Sets all pixels to given value.
	 */
	void Clear(const T& value) const
    {
		for(size_t i = 0; i < NumCells(); ++i) {
			mMap[i] = value;
		}
	}

	/*
	 * Access a given cell by index.
	 */
	T& operator[](size_t index)
	{
		return mMap[index];
	}

	const T& operator[](size_t index) const
	{
		return mMap[index];
	}

	/*
	 * Access a given cell by coordinate.
	 */
	T& operator()(int x, int y)
	{
		return mMap[size_t(y) * mSizeX + x];
	}

	const T& operator()(int x, int y) const
	{
		return mMap[size_t(y) * mSizeX + x];
	}

	/*
	 * Converts a coordinate in meters to a grid coordinate in pixels.
	 */
	float WorldToGrid(float meters) const
	{
		return meters * mInvScale - 0.5f;
	}

	/*
	 * Bilinear interpolation at given pixel position.
	 * A coordinate of (0, 0) gives the exact value of the first pixel.
	 */
	float BilinearLookup(float x, float y) const
	{
		const float a = x - floorf(x);
		const float b = y - floorf(y);

		return BilinearLookupEx(x, y, a, b);
	}

	/*
	 * Same as bilinear_lookup() but with pre-computed offsets a and b.
	 */
	float BilinearLookupEx(int x, int y, float a, float b) const
	{
		const int x0 = std::min(std::max(x, 0), mSizeX - 1);
		const int x1 = std::min(x0 + 1, mSizeX - 1);
		const int y0 = std::min(std::max(y, 0), mSizeY - 1);
		const int y1 = std::min(y0 + 1, mSizeY - 1);

		return		(*this)(x0, y0) * ((1.f - a) * (1.f - b))
				+	(*this)(x1, y0) * (a * (1.f - b))
				+	(*this)(x0, y1) * ((1.f - a) * b)
				+	(*this)(x1, y1) * (a * b);
	}

	/*
	 * Bilinear summation at a given pixel position.
	 * A coordinate of (0, 0) will sum at the first pixel exclusively.
	 */
	void BilinearSummation(float x, float y, const T& value)
	{
		const float a = x - floorf(x);
		const float b = y - floorf(y);

		BilinearSummationEx(x, y, a, b, value);
	}

	/*
	 * Same as bilinear_summation() but with pre-computed offsets a and b.
	 */
	void BilinearSummationEx(int x, int y, float a, float b, const T& value)
	{
		const int x0 = std::min(std::max(x, 0), mSizeX - 1);
		const int x1 = std::min(x0 + 1, mSizeX - 1);
		const int y0 = std::min(std::max(y, 0), mSizeY - 1);
		const int y1 = std::min(y0 + 1, mSizeY - 1);

		(*this)(x0, y0) += value * ((1.f - a) * (1.f - b));
		(*this)(x1, y0) += value * (a * (1.f - b));
		(*this)(x0, y1) += value * ((1.f - a) * b);
		(*this)(x1, y1) += value * (a * b);
	}

	/*
	 * Computes gauss-filtered and bilinear-interpolated first-order x and y gradient
	 * at given pixel position.
	 */
	void CalcGradient(float x, float y, float& dx, float& dy) const
	{
		static const float coeff_33_dxy[3][3] =
        {
				{-0.139505, 0, 0.139505},
				{-0.220989, 0, 0.220989},
				{-0.139505, 0, 0.139505}
		};

		const int x0 = x;
		const int y0 = y;

		const float a = x - floorf(x);
		const float b = y - floorf(y);

		dx = 0;
		dy = 0;

		for(int j = -1; j <= 1; ++j)
        {
			for(int i = -1; i <= 1; ++i)
            {
				const float value = BilinearLookupEx(x0 + i, y0 + j, a, b);
				dx += coeff_33_dxy[j+1][i+1] * value;
				dy += coeff_33_dxy[i+1][j+1] * value;
			}
		}

		dx /= 2 * mScale;
		dy /= 2 * mScale;
	}

	/*
	 * Computes gauss-filtered and bilinear-interpolated second-order x and y gradient
	 * at given pixel position.
	 */
	void CalcGradient2(float x, float y, float& ddx, float& ddy) const
	{
		static const float coeff_33_ddxy[3][3] =
        {
				{0.069752, -0.139504, 0.069752},
				{0.110494, -0.220989, 0.110494},
				{0.069752, -0.139504, 0.069752}
		};

		const int x0 = x;
		const int y0 = y;

		const float a = x - floorf(x);
		const float b = y - floorf(y);

		ddx = 0;
		ddy = 0;

		for(int j = -1; j <= 1; ++j)
        {
			for(int i = -1; i <= 1; ++i)
            {
				const float value = BilinearLookupEx(x0 + i, y0 + j, a, b);
				ddx += coeff_33_ddxy[j+1][i+1] * value;
				ddy += coeff_33_ddxy[i+1][j+1] * value;
			}
		}

		ddx /= 2 * mScale;
		ddy /= 2 * mScale;
	}

	/*
	 * Applies one smoothing iteration using a 3x3 gaussian kernel with sigma 1.
	 */
	void Smooth33_1()
	{
		static const float coeff_33_1[3][3] =
        {
				{0.077847, 0.123317, 0.077847},
				{0.123317, 0.195346, 0.123317},
				{0.077847, 0.123317, 0.077847}
		};

		GridMap<T> tmp(mSizeX, mSizeY, mScale);

		for(int y = 0; y < mSizeY; ++y)
        {
			for(int x = 0; x < mSizeX; ++x)
            {
				float sum = 0;
				for(int j = -1; j <= 1; ++j)
                {
					const int y_ = std::min(std::max(y + j, 0), mSizeY - 1);
					for(int i = -1; i <= 1; ++i)
                    {
						const int x_ = std::min(std::max(x + i, 0), mSizeX - 1);
						sum += coeff_33_1[j+1][i+1] * (*this)(x_, y_);
					}
				}
				tmp(x, y) = sum;
			}
		}
		*this = tmp;
	}

	/*
	 * Applies one smoothing iteration using a 5x5 gaussian kernel with sigma 2.
	 */
	void Smooth55_2()
	{
		static const float coeff_55_2[5][5] =
        {
				{0.0232468, 0.033824, 0.0383276, 0.033824, 0.0232468},
				{0.033824, 0.0492136, 0.0557663, 0.0492136, 0.033824},
				{0.0383276, 0.0557663, 0.0631915, 0.0557663, 0.0383276},
				{0.033824, 0.0492136, 0.0557663, 0.0492136, 0.033824},
				{0.0232468, 0.033824, 0.0383276, 0.033824, 0.0232468}
		};

		GridMap<T> tmp(mSizeX, mSizeY, mScale);

		for(int y = 0; y < mSizeY; ++y)
        {
			for(int x = 0; x < mSizeX; ++x)
            {
				float sum = 0;
				for(int j = -2; j <= 2; ++j)
                {
					const int y_ = std::min(std::max(y + j, 0), mSizeY - 1);
					for(int i = -2; i <= 2; ++i)
                    {
						const int x_ = std::min(std::max(x + i, 0), mSizeX - 1);
						sum += coeff_55_2[j+2][i+2] * (*this)(x_, y_);
					}
				}
				tmp(x, y) = sum;
			}
		}
		*this = tmp;
	}

	/*
	 * Returns a 2x downscaled map using a 4x4 gaussian filter with sigma 1.
	 */
	std::shared_ptr<GridMap<T>> DownScale()
	{
		static const float coeff_44_1[4][4]
        {
			{0.0180824, 0.049153, 0.049153, 0.0180824},
			{0.049153, 0.133612, 0.133612, 0.049153},
			{0.049153, 0.133612, 0.133612, 0.049153},
			{0.0180824, 0.049153, 0.049153, 0.0180824}
		};

		auto res = std::make_shared<GridMap<T>>(mSizeX / 2, mSizeY / 2, mScale * 2);

		for(int y = 0; y < mSizeY / 2; ++y)
        {
			for(int x = 0; x < mSizeX / 2; ++x)
            {
				float sum = 0;
				for(int j = -1; j <= 2; ++j)
                {
					const int y_ = std::min(std::max(y * 2 + j, 0), mSizeY - 1);
					for(int i = -1; i <= 2; ++i)
                    {
						const int x_ = std::min(std::max(x * 2 + i, 0), mSizeX - 1);
						sum += coeff_44_1[j+1][i+1] * (*this)(x_, y_);
					}
				}
				(*res)(x, y) = sum;
			}
		}
		return res;
	}

private:
	int mSizeX = 0;
	int mSizeY = 0;

	float mScale = 0;
	float mInvScale = 0;

	T* mMap = 0;

};


template<typename T>
class MultiGridMap
{
public:
	std::vector<std::shared_ptr<GridMap<T>>> layers;

	/*
	 * @param sizeXm Size of map in meters
	 * @param sizeYm Size of map in meters
	 * @param scale Size of one pixel in meters
	 */
	MultiGridMap(float sizeXm, float sizeYm, float scale, int numLayers)
	{
		layers.resize(numLayers);
		const int baseSizeX = sizeXm / (scale * (1 << (numLayers - 1))) + 0.5f;
		const int baseSizeY = sizeYm / (scale * (1 << (numLayers - 1))) + 0.5f;

		for(int i = 0; i < numLayers; ++i)
        {
			layers[i] = std::make_shared<GridMap<T>>(baseSizeX * (1 << (numLayers - i - 1)),
													 baseSizeY * (1 << (numLayers - i - 1)),
													 scale * (1 << i));
		}
	}


};


#endif