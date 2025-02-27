#pragma once
#include <vector>

template <typename T>
class Field
{
public:
	Field() = default;
	Field(int nx, int ny) : nx(nx), ny(ny) { values.resize(nx * ny); }

	inline const T& operator() (int i, int j) const { return values[i + nx * j]; }
	inline T& operator() (int i, int j) { return values[i + nx * j]; }

	std::vector<T>::iterator begin() { return values.begin(); }
	std::vector<T>::iterator end() { return values.end(); }

	std::vector<T>::const_iterator begin() const { return values.begin(); }
	std::vector<T>::const_iterator end() const { return values.end(); }

	Field& operator= (const Field& other) = default;

	void operator+= (const Field& other) { for (size_t i = 0; i < other.values.size(); ++i) { values[i] += other.values[i]; } }
	void operator-= (const Field& other) { for (size_t i = 0; i < other.values.size(); ++i) { values[i] -= other.values[i]; } }

	friend Field operator-(const Field& lhs, const Field& rhs)
	{
		Field tmp(lhs.nx, rhs.ny);
		for (int i = 0; i < tmp.values.size(); ++i) { tmp.values[i] = lhs.values[i] - rhs.values[i]; }
		return tmp;
	}

	friend Field operator+(const Field& lhs, const Field& rhs)
	{
		Field tmp(lhs.nx, rhs.ny);
		for (int i = 0; i < tmp.values.size(); ++i) { tmp.values[i] = lhs.values[i] + rhs.values[i]; }
		return tmp;
	}

	friend Field operator/(const Field& lhs, const Field& rhs)
	{
		Field tmp(lhs.nx, rhs.ny);
		for (int i = 0; i < tmp.values.size(); ++i) { tmp.values[i] = lhs.values[i] / (rhs.values[i] + 1e-20); }
		return tmp;
	}

	friend Field operator*(const Field& lhs, const Field& rhs)
	{
		Field tmp(lhs.nx, rhs.ny);
		for (int i = 0; i < tmp.values.size(); ++i) { tmp.values[i] = lhs.values[i] * rhs.values[i]; }
		return tmp;
	}

	void reset()
	{
		for (auto& value : values) value = 0.0;
	}

private:
	int nx, ny;
	std::vector<T> values;
};

using PhysicalField = Field<double>;