#pragma once
#include <vector>

enum class Flag : uint32_t {
    Open = 1 << 0,  // 0001
    Unused0 = 1 << 1,  // 0010
    Unused1 = 1 << 2,  // 0100
    Unused2 = 1 << 3   // 1000
};

template<typename T>
constexpr auto operator|(T lhs, T rhs) -> T {
    return static_cast<T>(static_cast<uint32_t>(lhs) | static_cast<uint32_t>(rhs));
}
template<typename T>
constexpr auto operator&(T lhs, T rhs) -> bool {
    return (static_cast<uint32_t>(lhs) & static_cast<uint32_t>(rhs)) != 0;
}

class FlagField 
{
public:
    FlagField(int nx, int ny) : nx(nx), ny(ny) { values.resize(nx * ny); }

    void setFlag(int i, int j, Flag flag)
    {
        values[j * nx + i] |= static_cast<uint32_t>(flag);
    }

    void setFlag(Flag flag)
    {
        for (auto& elm : values)
            elm |= static_cast<uint32_t>(flag);
    }

    void clearFlag(int i, int j, Flag flag)
    {
        values[j * nx + i] &= ~static_cast<uint32_t>(flag);
    }

    bool operator()(int i, int j, Flag flag) const
    {
        return (values[j * nx + i] & static_cast<uint32_t>(flag)) != 0;
    }
private:
    std::vector<uint32_t> values;  // Each element stores flag bits efficiently
    size_t nx, ny;
};



//class FlagField
//{
//public:
//    //PhysicalField(int nx, int ny) : nx(nx), ny(ny) { values.resize(nx * ny); }
//
//    inline const int& operator() (int i, int j) const { return values[i + nx * j]; }
//    inline int& operator() (int i, int j) { return values[i + nx * j]; }
//
//    std::vector<int>::iterator begin() { return values.begin(); }
//    std::vector<int>::iterator end() { return values.end(); }
//
//    void reset()
//    {
//        for (auto& value : values) value = 0.0;
//    }
//
//private:
//    int nx, ny;
//    std::vector<int> values;
//};
//
//enum class Flags
//{
//    H
//};
