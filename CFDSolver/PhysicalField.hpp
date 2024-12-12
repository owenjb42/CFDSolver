#pragma once
#include <vector>

class PhysicalField
{
public:
	PhysicalField(int nx, int ny) : nx(nx), ny(ny) { values.resize(nx * ny); }

	const double& operator() (int i, int j) const { return values[i + nx * j]; }
	double& operator() (int i, int j) { return values[i + nx * j]; }

	void operator+= (const PhysicalField& other) { for (size_t i = 0; i < other.values.size(); ++i) { values[i] += other.values[i]; } }

	int nx, ny;
	std::vector<double> values;
};