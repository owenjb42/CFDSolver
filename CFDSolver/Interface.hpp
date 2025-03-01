#pragma once

#include "raygui.h"
#include "raylib.h"
#include <vector>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <random>

#include "FairMutex.hpp"
#include "Helper.hpp"

#include "PhysicalField.hpp"

struct Face
{
    bool operator==(const Face& other) const = default;
    bool component_dir /*0: u, 1: v*/; int i; int j;
};
struct Boundary : public Face
{
    bool operator==(const Boundary& other) const { return Face::operator==(other); }
    enum Type{FixedInflow = 1, FixedOutflow = 2, Open = 3};
    Type type{ Open };
    double optional_velocity{ 0.0 }, optional_temp{ 0.0 }, optional_pressure{ 0.0 };
    bool boundary_dir{ 0 };// 0: +, 1: -
};
namespace std
{
    template <> struct hash<Face> { size_t operator()(const Face& f) const { return (hash<bool>()(f.component_dir) ^ (hash<int>()(f.i) << 1) ^ (hash<int>()(f.j) << 2)); } };
    template <> struct hash<Boundary> { size_t operator()(const Boundary& f) const { return (hash<bool>()(f.component_dir) ^ (hash<int>()(f.i) << 1) ^ (hash<int>()(f.j) << 2)); } };
}

class SolverStaggeredIMEXTemp;

class Interface 
{
public:
    Interface()
    {
        SetTraceLogLevel(LOG_ERROR);
        SetConfigFlags(FLAG_WINDOW_RESIZABLE);
        InitWindow(1000, 1000, "CFD Simulation Visualization");
        SetWindowMinSize(1000, 1000);
        SetTargetFPS(60);
        ResetView();
    }

    ~Interface()
    {
        CloseWindow();
    }

    void SetData(SolverStaggeredIMEXTemp& solver);
    void GetDataFromBuffer();

    void RenderModel()
    {
        in_model_creator_mode = true;
        is_solving = false;

        while (!WindowShouldClose() && in_model_creator_mode)
        {
            HandleInput();

            BeginDrawing();
            ClearBackground(RAYWHITE);

            DrawGrid();
            HandlePannel();

            EndDrawing();
        }

        GenerateSpeedPoints();
    }

    void RenderResults()
    {
        in_render_mode = true;
        is_solving = false;

        while (!WindowShouldClose() && in_render_mode)
        {
            GetDataFromBuffer();

            // setup only needed if there is new data 
            if (recalculate_auxiliary_data) { SetAuxData(); }

            HandleInput();

            BeginDrawing();
            ClearBackground(RAYWHITE);

            if (enable_flags[2]) DrawField();
            DrawGrid();
            if (enable_flags[1]) DrawVelocityVectors();
            if (enable_flags[0]) DrawStreamlines();
            HandlePannel();

            EndDrawing();
        }
    }

    // UI
    float sidebar_width = 400.0f;

    bool in_model_creator_mode{ false };
    bool in_render_mode{ false };
    bool is_solving{ false };

    // Model creator UI
    char fixed_in_options[2][32] = { "Velocity", "Temperature" };
    bool fixed_in_edit_mode[2] = { false, false };
    char fixed_in_inputs[2][32] = { "1", "20" };
    float fixed_in_vel{ 1.0f }, fixed_in_temp{ 20.0f };
    bool fixed_in_dir{ false };

    char open_options[2][32] = { "Pressure", "Temperature" };
    bool open_edit_mode[2] = { false, false };
    char open_inputs[2][32] = { "0", "20" };
    float open_pressure{ 0.0f }, open_temp{ 20.0f };

    char fixed_out_options[32] = { "Velocity" };
    bool fixed_out_edit_mode = { false };
    char fixed_out_inputs[32] = { "1" };
    float fixed_out_vel{ 1.0f };
    bool fixed_out_dir{ false };

    // Selection UI
    float key_hold_time[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

    Vector2 drag_start{};
    bool dragging{ false };

    std::vector<Face> selected_faces;
    std::unordered_set<Face> blocked_faces;
    std::unordered_set<Boundary> boundary_faces;

    // Results UI
    bool enable_flags[3] = { true, true, true };

    char stremline_options[4][32] = { "num", "steps", "thickness", "reverse" };
    bool stremline_option_editmode[3] = { false, false, false };
    char stremline_option_inputs[3][32] = { "100", "100", "2" };
    int streamline_num{ 100 }, streamline_steps{ 100 };
    float streamline_thickness{ 2.0f };
    bool reverse_streamlines{ false };

    char v_arrow_options[32] = { "scale" };
    bool v_arrow_editmode = { false };
    char v_arrow_inputs[32] = { "1" };
    float v_arrow_scale{ 1.0f };

    char field_options[2][32] = { "Presure", "Temperature" };
    bool plot_p_or_t{ false };

    // Result infomation
    int nx{ 5 }, ny{ 5 };
    float dx{ 0.005f }, dy{ 0.005f };
    PhysicalField u_buffer, v_buffer, p_buffer, t_buffer, divergence_buffer, u_face_buffer, v_face_buffer;
    PhysicalField u, v, p, t, divergence, u_face, v_face;
    FairMutex mutex;
    float cellWidth{ 0 }, cellHeight{ 0 };
    float minPressure{ 0.0 }, maxPressure{ 0.0 }, minTemp{ 0.0 }, maxTemp{ 0.0 }, minVelocity{ 0.0 }, maxVelocity{ 0.0 };
      
    // Residuals
    char residuals[3][32] = { "Pressure", "Divergence", "Temperature" };
    float residual_values_buffer[3] = { 0.0, 0.0, 0.0 };
    char residual_values[3][32] = {"", "", ""};

    // Render settings
    float offsetX{ 0 }, offsetY{ 0 }; // offset in the screen co-ords
    float zoom{ 1.0f };
    float zoom_sensitivity{ 0.1f };
    float minZoomLevel{ 0.01f };
    float velocityScale{ 0.01f };

    bool needs_update{ true };
    bool recalculate_auxiliary_data{ false };

    void HandleInput();
    void HandlePannel();

    Vector2 ScreenToGridPos(Vector2 mouse_pos)
    {
        return Vector2{ (mouse_pos.x - offsetX) / zoom, (mouse_pos.y - offsetY) / zoom };
    }

    Vector2 GridToScreenPos(Vector2 grid_pos)
    {
        return Vector2{ grid_pos.x * zoom + offsetX, grid_pos.y * zoom + offsetY };
    }

    float GetViewScreenWidth()
    {
        return GetScreenWidth() - sidebar_width;
    }

    bool IsMouseInView()
    {
        auto pos = GetMousePosition();
        return pos.x > 0.0f && pos.x < GetViewScreenWidth() && pos.y > 0.0f && pos.y < GetScreenHeight();
    }

    void ResetView()
    {
        offsetX = offsetY = 0.0f;
        double max_scale = std::max(nx * dx / GetViewScreenWidth(), ny * dy / GetScreenHeight());
        zoom = 1 / max_scale;
    }

    void SanatiseFaceSets()
    {
        // Remove faces which no longer conform to the current grid size
        for (auto it = boundary_faces.begin(); it != boundary_faces.end(); )
        {
            bool is_pos_flow = (it->type == Boundary::FixedInflow && it->boundary_dir == 0);
            if (it->component_dir == 0)
            {
                if (it->i > (is_pos_flow ? nx - 1 : nx) || it->j > ny - 1)
                    it = boundary_faces.erase(it);
                else ++it;
            }
            else
            {
                if (it->i > nx - 1 || it->j > (is_pos_flow ? ny - 1 : ny))
                    it = boundary_faces.erase(it);
                else ++it;
            }
        }

        for (auto it = blocked_faces.begin(); it != blocked_faces.end(); ) 
        {
            if (it->component_dir == 0)
            {
                if (it->i > nx || it->j > ny - 1)
                    it = blocked_faces.erase(it);
                else ++it;
            }
            else
            {
                if (it->i > nx - 1 || it->j > ny)
                    it = blocked_faces.erase(it);
                else ++it;
            }
        }
    }

    void SetAuxData()
    {
        minPressure = *std::ranges::min_element(p.begin(), p.end());
        maxPressure = *std::ranges::max_element(p.begin(), p.end());

        maxVelocity = std::numeric_limits<float>::min();
        minVelocity = std::numeric_limits<float>::max();

        minTemp = *std::ranges::min_element(t.begin(), t.end());
        maxTemp = *std::ranges::max_element(t.begin(), t.end());

        double minDivergence = *std::ranges::min_element(divergence.begin(), divergence.end());
        double maxDivergence = *std::ranges::max_element(divergence.begin(), divergence.end());

        double minMaxDivergence = std::max(std::abs(minDivergence), std::abs(maxDivergence));

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                float uVal = static_cast<float>(u(i, j));
                float vVal = static_cast<float>(v(i, j));

                // Calculate the magnitude and apply logarithmic scaling
                float magnitude = sqrt(uVal * uVal + vVal * vVal);
                if (magnitude > maxVelocity) maxVelocity = magnitude;
                if (magnitude < minVelocity) minVelocity = magnitude;
            }
        }

        GenerateStreamlines();

        recalculate_auxiliary_data = false;
    }

    void DrawGrid()
    {
        // Draw Grid Lines
        for (int i = 0; i <= nx; ++i)
        {
            float xPos = i * cellWidth + offsetX;
            DrawLineV({ xPos, offsetY }, { xPos, offsetY + ny * cellHeight }, GRAY);
        }

        for (int j = 0; j <= ny; ++j)
        {
            float yPos = j * cellHeight + offsetY;
            DrawLineV({ offsetX, yPos }, { offsetX + nx * cellWidth, yPos }, GRAY);
        }

        // Draw Boundary & Blocked Faces
        float bold_width = std::max(zoom / 4000.0f, 2.0f);
        float ofset = 0.5f * bold_width;

        // Border outline
        DrawLineEx({ offsetX, offsetY - ofset },                   { offsetX, offsetY + ny * cellHeight + ofset }, bold_width, BLACK);
        DrawLineEx({ offsetX + nx * cellWidth, offsetY - ofset },  { offsetX + nx * cellWidth, offsetY + ny * cellHeight + ofset }, bold_width, BLACK);
        DrawLineEx({ offsetX - ofset, offsetY },                   { offsetX + nx * cellWidth + ofset, offsetY }, bold_width, BLACK);
        DrawLineEx({ offsetX - ofset, offsetY + ny * cellHeight }, { offsetX + nx * cellWidth + ofset, offsetY + ny * cellHeight }, bold_width, BLACK);

        for (const auto& face : blocked_faces)
        {
            float xPos = face.i * cellWidth + offsetX;
            float yPos = face.j * cellHeight + offsetY;
            if (face.component_dir == 0)
                DrawLineEx({ xPos, yPos - ofset }, { xPos, yPos + cellHeight + ofset }, bold_width, BLACK);
            else
                DrawLineEx({ xPos - ofset, yPos }, { xPos + cellHeight + ofset, yPos }, bold_width, BLACK);
        }

        for (const auto& boundary : boundary_faces)
        {
            float xPos = boundary.i * cellWidth + offsetX;
            float yPos = boundary.j * cellHeight + offsetY;

            DrawBoundary(boundary.component_dir, boundary.boundary_dir, xPos, yPos, cellWidth, cellHeight, bold_width, boundary.type);
        }

        for (const auto& face : selected_faces)
        {
            float xPos = face.i * cellWidth + offsetX;
            float yPos = face.j * cellHeight + offsetY;
            if (face.component_dir == 0)
                DrawLineEx({ xPos, yPos - ofset }, { xPos, yPos + cellHeight + ofset }, bold_width, RED);
            else
                DrawLineEx({ xPos - ofset, yPos }, { xPos + cellHeight + ofset, yPos }, bold_width, RED);
        }
    }

    void DrawVelocityVectors()
    {
        auto DrawArrowHead = [](Vector2 centre, Vector2 direction, float scale, Color color)
        {
            float angle = atan2(direction.y, direction.x);
            float arrowHeadAngle = 0.8f;

            Vector2 arrowPoint1 = Vector2{ centre.x + direction.x * scale, centre.y + direction.y * scale };
            Vector2 arrowPoint2 = Vector2{ centre.x - scale * cos(angle + arrowHeadAngle), centre.y - scale * sin(angle + arrowHeadAngle) };
            Vector2 arrowPoint3 = Vector2{ centre.x - scale * cos(angle - arrowHeadAngle), centre.y - scale * sin(angle - arrowHeadAngle) };
            DrawTriangle(arrowPoint1, arrowPoint2, arrowPoint3, color);
            DrawTriangleLines(arrowPoint1, arrowPoint2, arrowPoint3, BLACK);
        };

        float scale = velocityScale * zoom * v_arrow_scale;

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                float uVal = static_cast<float>(u(i, j));
                float vVal = static_cast<float>(v(i, j));

                // Calculate the magnitude and apply logarithmic scaling
                float magnitude = sqrt(uVal * uVal + vVal * vVal);
                uVal /= magnitude; vVal /= magnitude;
                float logMagnitude = log10(magnitude + 1.0f);

                // Calculate the starting position of the vector (center of the cell)
                float xPos = i * cellWidth + cellWidth / 2 + offsetX;
                float yPos = j * cellHeight + cellHeight / 2 + offsetY;

                DrawArrowHead(Vector2{ (float)xPos, (float)yPos }, Vector2{ (float)uVal, (float)vVal }, logMagnitude * scale, MapToColor(magnitude, minVelocity, maxVelocity));
            }
        }
    }

    void DrawField()
    {
        const PhysicalField& f = plot_p_or_t == 0 ? p : t;
        const float minValue = plot_p_or_t == 0 ? minPressure : minTemp;
        const float maxValue = plot_p_or_t == 0 ? maxPressure : maxTemp;

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                Color color = MapToColor(f(i, j), minValue, maxValue);
                float xPos = i * cellWidth + offsetX;
                float yPos = j * cellHeight + offsetY;
                DrawRectangleV({ xPos, yPos }, { cellWidth, cellHeight }, color);
            }
        }
    }

    struct Streamline
    {
        bool is_reverse{ false };
        std::vector<Vector2> points;
    };
    std::vector<Streamline> streamlines;
    std::vector<std::pair<float, float>> seeds;

    void GenerateSpeedPoints()
    {
        seeds.clear();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dis(0.0f, 1.0f);
        for (int seed = 0; seed < streamline_num; ++seed)
        {
            // Generate random seed points within the field
            float x = dis(gen) * (float)nx;
            float y = dis(gen) * (float)ny;
            seeds.push_back({ x, y });
        }
    }

    void GenerateStreamlines()
    {
        auto sampleVelocity = [&](float x, float y) -> Vector2
        {
            if (x <= 0.5f || y <= 0.5f || x >= (float)(nx) - 0.5f || y >= (float)(ny) - 0.5f)
                return { 0.0f, 0.0f }; // Out of bounds

            x -= 0.5f;
            y -= 0.5f;

            int i = static_cast<int>(x);
            int j = static_cast<int>(y);

            // Bilinear interpolation
            float tx = x - i;
            float ty = y - j;

            float vx = (1 - tx) * (1 - ty) * u(i, j) +
                tx * (1 - ty) * u(i + 1, j) +
                (1 - tx) * ty * u(i, j + 1) +
                tx * ty * u(i + 1, j + 1);

            float vy = (1 - tx) * (1 - ty) * v(i, j) +
                tx * (1 - ty) * v(i + 1, j) +
                (1 - tx) * ty * v(i, j + 1) +
                tx * ty * v(i + 1, j + 1);

            return { vx, vy };
        };

        streamlines.clear();

        float stepSize = 0.1f;

        for (auto [x, y] : seeds)
        {
            float sign = reverse_streamlines ? -1.0f : 1.0f;

            Streamline streamline;

            for (int step = 0; step < streamline_steps; ++step)
            {
                Vector2 velocity = sampleVelocity(x, y);

                float magnitude = std::pow(std::pow(velocity.x, 2) + std::pow(velocity.y, 2), 0.5);

                if (magnitude < 1e-10)
                    break;

                // Normalize velocity and move
                velocity = { velocity.x * (stepSize / magnitude), velocity.y * (stepSize / magnitude) };

                float nextX = x + sign * velocity.x;
                float nextY = y + sign * velocity.y;

                // Store point in streamline
                streamline.points.push_back(Vector2{ x * dx, y * dy });

                x = nextX;
                y = nextY;

                // Stop if out of bounds
                if (x <= 0 || y <= 0 || x >= nx || y >= ny)
                    break;
            }
            streamline.is_reverse = reverse_streamlines;
            streamlines.push_back(streamline);
        }
    }

    void DrawStreamlines()
    {
        for (const auto& streamline : streamlines)
        {
            for (size_t i = 1; i < streamline.points.size(); ++i)
            {
                auto screenPoint = GridToScreenPos(streamline.points[i]);
                auto otherScreenPoint = GridToScreenPos(streamline.points[i - 1]);

                DrawLineEx(
                    { screenPoint.x, screenPoint.y },
                    { otherScreenPoint.x, otherScreenPoint.y },
                    streamline_thickness,
                    BLACK
                );

                //if ((i + 5) % 10 == 0)
                //{
                //    if (streamline.is_reverse) std::swap(screenPoint, otherScreenPoint);
                //
                //    float arrowAngle = PI / 6;
                //    float angle = atan2(screenPoint.y - otherScreenPoint.y, screenPoint.x - otherScreenPoint.x);
                //    float scale = zoom * 0.0005;
                //
                //    Vector2 arrowPoint2 = Vector2{ screenPoint.x - scale * cos(angle + arrowAngle), screenPoint.y - scale * sin(angle + arrowAngle) };
                //    Vector2 arrowPoint3 = Vector2{ screenPoint.x - scale * cos(angle - arrowAngle), screenPoint.y - scale * sin(angle - arrowAngle) };
                //
                //    DrawLine(
                //        screenPoint.x, screenPoint.y,
                //        arrowPoint2.x, arrowPoint2.y,
                //        BLACK
                //    );
                //
                //    DrawLine(
                //        screenPoint.x, screenPoint.y,
                //        arrowPoint3.x, arrowPoint3.y,
                //        BLACK
                //    );
                //}
            }
        }
    }
};