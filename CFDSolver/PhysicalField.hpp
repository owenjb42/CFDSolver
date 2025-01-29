#pragma once
#include <vector>

class PhysicalField
{
public:
	PhysicalField(int nx, int ny) : nx(nx), ny(ny) { values.resize(nx * ny); }

	inline const double& operator() (int i, int j) const { return values[i + nx * j]; }
	inline double& operator() (int i, int j) { return values[i + nx * j]; }

	std::vector<double>::iterator begin() { return values.begin(); }
	std::vector<double>::iterator end() { return values.end(); }

	void operator+= (const PhysicalField& other) { for (size_t i = 0; i < other.values.size(); ++i) { values[i] += other.values[i]; } }

	friend PhysicalField operator-(const PhysicalField& lhs, const PhysicalField& rhs)
	{
		PhysicalField tmp(lhs.nx, rhs.ny);
		for (int i = 0; i < tmp.values.size(); ++i) { tmp.values[i] = lhs.values[i] - rhs.values[i]; }
		return tmp;
	}

	friend PhysicalField operator/(const PhysicalField& lhs, const PhysicalField& rhs)
	{
		PhysicalField tmp(lhs.nx, rhs.ny);
		for (int i = 0; i < tmp.values.size(); ++i) { tmp.values[i] = lhs.values[i] / rhs.values[i]; }
		return tmp;
	}

	void reset()
	{
		for (auto& value : values) value = 0.0;
	}

	int nx, ny;
	std::vector<double> values;
};