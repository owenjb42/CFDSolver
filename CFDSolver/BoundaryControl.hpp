#pragma once
#include <vector>

class SolverStaggeredIMEXTemp;

struct BoundaryCondition
{
    BoundaryCondition(bool direction, bool component) : direction(direction), component(component) {}

    bool direction; // 0 : +, 1 : -
    bool component; // 0 : u, 1 : v
    std::vector<std::pair<int, int>> boundary_faces;
};

struct FixedInletBoundaryCondition : public BoundaryCondition
{
    FixedInletBoundaryCondition(bool direction, bool component) : BoundaryCondition(direction, component) {}

    void ApplyForVelocity(SolverStaggeredIMEXTemp& solver);
    void ApplyForTemperature(SolverStaggeredIMEXTemp& solver);
    void CorrectBoundaryCellVelocities(SolverStaggeredIMEXTemp& solver);

    double inlet_temperature{};
    double inlet_velocity{};
};

struct FixedOutletBoundaryCondition : public BoundaryCondition
{
    FixedOutletBoundaryCondition(bool direction, bool component) : BoundaryCondition(direction, component) {}

    void ApplyForVelocity(SolverStaggeredIMEXTemp& solver);
    void CorrectBoundaryCellVelocities(SolverStaggeredIMEXTemp& solver);

    double outlet_velocity{};
};

struct OpenBoundaryCondition : public BoundaryCondition
{
    OpenBoundaryCondition(bool direction, bool component) : BoundaryCondition(direction, component) {}

    void ApplyForPressure(SolverStaggeredIMEXTemp& solver);
    void ApplyForTemperature(SolverStaggeredIMEXTemp& solver);
    void CorrectBoundaryCellVelocities(SolverStaggeredIMEXTemp& solver);

    double temperature{};
    double pressure{};
};

struct FrictionBoundaryCondition : public BoundaryCondition
{

};