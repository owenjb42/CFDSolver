#pragma once

#include "raylib.h"
#include <vector>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include "FairMutex.hpp"

#include "PhysicalField.hpp"

struct Face
{
    bool operator==(const Face& other) const = default;
    bool dir; int i; int j;
};
namespace std
{
    template <> struct hash<Face> { size_t operator()(const Face& f) const { return (hash<bool>()(f.dir) ^ (hash<int>()(f.i) << 1) ^ (hash<int>()(f.j) << 2)); } };
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
        SetTargetFPS(60);
        ResetView();
    }

    ~Interface()
    {
        CloseWindow();
    }

    void SetData(SolverStaggeredIMEXTemp& solver);

    void RenderModel()
    {
        in_model_creator_mode = true;
        while (!WindowShouldClose() && in_model_creator_mode)
        {
            HandleInput();
            BeginDrawing();
            ClearBackground(RAYWHITE);
            DrawGrid();
            EndDrawing();
        }
    }

    void RenderResults()
    {
        double minTemp{ 0.0 };
        double maxTemp{ 0.0 };

        while (!WindowShouldClose())
        {
            std::lock_guard lk(mutex);

            // setup only needed if there is new data 
            if (needs_update)
            {
                minPressure = *std::ranges::min_element(p.begin(), p.end());
                maxPressure = *std::ranges::max_element(p.begin(), p.end());

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

                //GenerateStreamlines(200, 200, 0.1);

                needs_update = false;
            }

            HandleInput();

            BeginDrawing();
            ClearBackground(RAYWHITE);

            DrawField(t, minTemp, maxTemp);
            DrawGrid();

            DrawVelocityVectors();
            //DrawStreamlines();

            EndDrawing();
        }
    }

    FairMutex mutex;

    // Model creator
    bool in_model_creator_mode{ false };
    float key_hold_time[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

    Vector2 drag_start{};
    bool dragging{ false };

    std::vector<Face> selected_faces;
    std::unordered_set<Face> blocked_faces;

    // Result infomation
    int nx{ 5 }, ny{ 5 };
    float dx{ 0.005f }, dy{ 0.005f };
    PhysicalField u, v, p, t, divergence, u_face, v_face;

    float cellWidth{ 0 }, cellHeight{ 0 };

    bool needs_update{ true };

    // Render settings
    float offsetX{ 0 }, offsetY{ 0 }; // offset in the screen co-ords
    float zoom{ 1.0f };
    float zoom_sensitivity{ 0.1f };
    float minZoomLevel{ 0.01f };
    float velocityScale{ 0.01f };

    float minPressure{ 0.0 }, maxPressure{ 0.0 };
    float minVelocity{ std::numeric_limits<float>::max() }, maxVelocity{ std::numeric_limits<float>::lowest() };

    void HandleInput()
    {
        auto mouse_pos = GetMousePosition();

        // Zoom in/out using mouse wheel
        auto new_zoom = zoom * (1.0f + GetMouseWheelMove() * zoom_sensitivity);

        // Keep the mouse position pointing at the same position
        auto old_mouse_grid = ScreenToGridPos(mouse_pos);
        offsetX *= new_zoom / zoom;
        offsetY *= new_zoom / zoom;
        zoom = new_zoom;
        auto new_mouse_grid = ScreenToGridPos(mouse_pos);
        offsetX -= (old_mouse_grid.x - new_mouse_grid.x) * zoom;
        offsetY -= (old_mouse_grid.y - new_mouse_grid.y) * zoom;

        cellWidth = zoom * dx;
        cellHeight = zoom * dy;

        // Pan using right mouse button drag
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT))
        {
            Vector2 mouseDelta = GetMouseDelta();
            offsetX += mouseDelta.x;
            offsetY += mouseDelta.y;
        }

        if (IsKeyPressed(KEY_R))
            ResetView();

        if (in_model_creator_mode)
        {
            if (IsKeyPressed(KEY_S))
                in_model_creator_mode = false;

            key_hold_time[0] = IsKeyDown(KEY_LEFT) ? key_hold_time[0] + GetFrameTime() : 0.0f;
            key_hold_time[1] = IsKeyDown(KEY_RIGHT) ? key_hold_time[1] + GetFrameTime() : 0.0f;
            key_hold_time[2] = IsKeyDown(KEY_UP) ? key_hold_time[2] + GetFrameTime() : 0.0f;
            key_hold_time[3] = IsKeyDown(KEY_DOWN) ? key_hold_time[3] + GetFrameTime() : 0.0f;

            if (IsKeyPressed(KEY_LEFT) || key_hold_time[0] > 0.3f)
                nx = std::max(3, nx - 1);
            if (IsKeyPressed(KEY_RIGHT) || key_hold_time[1] > 0.3f)
                nx = std::min(200, nx + 1);
            if (IsKeyPressed(KEY_UP) || key_hold_time[2] > 0.3f)
                ny = std::max(3, ny - 1);
            if (IsKeyPressed(KEY_DOWN) || key_hold_time[3] > 0.3f)
                ny = std::min(200, ny + 1);

            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
            {
                if (!dragging)
                {
                    drag_start = mouse_pos;
                    dragging = true;
                }
                else 
                {
                    DrawRectangle(mouse_pos.x, mouse_pos.y, drag_start.x - mouse_pos.x, drag_start.y - mouse_pos.y, SKYBLUE);
                    DrawRectangleLines(mouse_pos.x, mouse_pos.y, drag_start.x - mouse_pos.x, drag_start.y - mouse_pos.y, BLUE);

                    auto grid_pos_1 = ScreenToGridPos(drag_start);
                    auto grid_pos_2 = ScreenToGridPos(mouse_pos);

                    int i1 = std::ceil(std::min(grid_pos_1.x, grid_pos_2.x) / dx);
                    int i2 = std::floor(std::max(grid_pos_1.x, grid_pos_2.x) / dx);
                    int j1 = std::ceil(std::min(grid_pos_1.y, grid_pos_2.y) / dy);
                    int j2 = std::floor(std::max(grid_pos_1.y, grid_pos_2.y) / dy);

                    selected_faces.clear();
                    if (i1 < i2 || j1 < j2)
                    {
                        for (int i = std::max(i1, 0); i <= std::min(i2, nx); ++i)
                            for (int j = std::max(j1, 0); j < std::min(j2, ny); ++j)
                                selected_faces.push_back({ 0, i, j });
                        for (int i = std::max(i1, 0); i < std::min(i2, nx); ++i)
                            for (int j = std::max(j1, 0); j <= std::min(j2, ny); ++j)
                                selected_faces.push_back({ 1, i, j });
                    }
                }
            }
            else
            {
                dragging = false;
            }

            // toggle blocked face
            if (IsKeyPressed(KEY_B))
            {
                for (const auto& face : selected_faces) 
                {
                    if (IsKeyDown(KEY_LEFT_SHIFT)) blocked_faces.erase(face);
                    else blocked_faces.insert(face);
                }
            }
        }
        else
        {
            // Scale Velocity
            if (IsKeyDown(KEY_A))
                velocityScale *= 1.0f + zoom_sensitivity;
            if (IsKeyDown(KEY_S))
                velocityScale *= 1.0f / (1.0f + zoom_sensitivity);
        }
    }

    Vector2 ScreenToGridPos(Vector2 mouse_pos)
    {
        return Vector2{ (mouse_pos.x - offsetX) / zoom, (mouse_pos.y - offsetY) / zoom };
    }

    Vector2 GridToScreenPos(Vector2 grid_pos)
    {
        return Vector2{ grid_pos.x * zoom + offsetX, grid_pos.y * zoom + offsetY };
    }

    void ResetView()
    {
        offsetX = offsetY = 0.0f;
        double max_scale = std::max(nx * dx / GetScreenWidth(), ny * dy / GetScreenHeight());
        zoom = 1 / max_scale;
    }

    void DrawGrid()
    {
        // Draw grid lines
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

        float bold_width = zoom / 4000.0f;
        float ofset = 0.5f * bold_width;

        for (const auto& face : blocked_faces)
        {
            float xPos = face.i * cellWidth + offsetX;
            float yPos = face.j * cellHeight + offsetY;
            if (face.dir == 0)
                DrawLineEx({ xPos, yPos - ofset }, { xPos, yPos + cellHeight + ofset }, bold_width, BLACK);
            else
                DrawLineEx({ xPos - ofset, yPos }, { xPos + cellHeight + ofset, yPos }, bold_width, BLACK);
        }

        for (const auto& face : selected_faces)
        {
            float xPos = face.i * cellWidth + offsetX;
            float yPos = face.j * cellHeight + offsetY;
            if (face.dir == 0)
                DrawLineEx({ xPos, yPos - ofset }, { xPos, yPos + cellHeight + ofset }, bold_width, RED);
            else
                DrawLineEx({ xPos - ofset, yPos }, { xPos + cellHeight + ofset, yPos }, bold_width, RED);
        }
    }

    Color MapToColor(float value, float minValue, float maxValue)
    {
        if (std::isnan(value) || value < minValue || value > maxValue || minValue == maxValue)
            return WHITE;

        float normalized = (value - minValue) / (maxValue - minValue);

        float r = std::clamp(1.5 - std::abs(4.0 * normalized - 3.0), 0.0, 1.0);
        float g = std::clamp(1.5 - std::abs(4.0 * normalized - 2.0), 0.0, 1.0);
        float b = std::clamp(1.5 - std::abs(4.0 * normalized - 1.0), 0.0, 1.0);

        return Color{ (unsigned char)(255 * r), (unsigned char)(255 * g), (unsigned char)(255 * b), 255 };
    }

    void DrawVelocityVectors()
    {
        float scale = velocityScale * zoom;

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

    void DrawArrowHead(Vector2 centre, Vector2 direction, float scale, Color color)
    {
        float angle = atan2(direction.y, direction.x);
        float arrowHeadAngle = 0.8f;

        Vector2 arrowPoint1 = Vector2{ centre.x + direction.x * scale, centre.y + direction.y * scale };
        Vector2 arrowPoint2 = Vector2{ centre.x - scale * cos(angle + arrowHeadAngle), centre.y - scale * sin(angle + arrowHeadAngle) };
        Vector2 arrowPoint3 = Vector2{ centre.x - scale * cos(angle - arrowHeadAngle), centre.y - scale * sin(angle - arrowHeadAngle) };
        DrawTriangle(arrowPoint1, arrowPoint2, arrowPoint3, color);
        DrawTriangleLines(arrowPoint1, arrowPoint2, arrowPoint3, BLACK);
    }

    void DrawField(const PhysicalField& f, float minValue, float maxValue)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                float value = f(i, j);

                // Map the value to a color (heatmap)
                Color color = MapToColor(value, minValue, maxValue);

                // Draw a rectangle for each cell representing the value
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

    void GenerateStreamlines(int numSeeds, int steps, float stepSize)
    {
        auto sampleVelocity = [&](float x, float y) -> Vector2
        {
            if (x <= 0.5f || y <= 0.5f || x >= nx - 0.5f || y >= ny - 0.5f)
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

        for (int seed = 0; seed < numSeeds; ++seed)
        {
            // Generate random seed points within the field
            float x = GetRandomValue(0, nx);
            float y = GetRandomValue(0, ny);

            bool reverse = GetRandomValue(0, 1) > 0.5f;
            float sign = reverse ? -1.0f : 1.0f;

            Streamline streamline;

            for (int step = 0; step < steps; ++step)
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
            streamline.is_reverse = reverse;
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

                DrawLine(
                    screenPoint.x, screenPoint.y,
                    otherScreenPoint.x, otherScreenPoint.y,
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