#pragma once

#include "raylib.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "FairMutex.hpp"

#include "PhysicalField.hpp"

class CFDVisualizer {
public:
    CFDVisualizer(int nx, int ny, float dx, float dy, const PhysicalField& u, const PhysicalField& v, const PhysicalField& p, FairMutex& cell_data_mutex)
        : nx(nx), ny(ny), dx(dx), dy(dy), u(u), v(v), p(p), mutex(cell_data_mutex)
    {
        SetTraceLogLevel(LOG_ERROR);
        SetConfigFlags(FLAG_WINDOW_RESIZABLE);
        InitWindow(1000, 1000, "CFD Simulation Visualization");
        SetTargetFPS(60);

        minPressure = *std::ranges::min_element(p.values);
        maxPressure = *std::ranges::max_element(p.values);

        double max_scale = std::max(nx * dx / 1000, ny * dy / 1000);
        zoom = 1 / max_scale;
        velocityScale = 1000.0f * max_scale;
    }

    ~CFDVisualizer() 
    {
        CloseWindow();
    }

    void Render() 
    {
        GenerateStreamlines(100, 5000, 0.1);

        while (true)
        {
            std::lock_guard lk(mutex);

            HandleInput();

            BeginDrawing();
            ClearBackground(RAYWHITE);

            DrawField(p, minPressure, maxPressure);
            DrawGrid();
            //DrawVelocityVectors();
            DrawStreamlines();

            EndDrawing();
        }
    }

private:
    // Result infomation
    int nx, ny;
    float dx, dy;
    FairMutex& mutex;
    const PhysicalField& u, & v, & p;
    int cellWidth{ 0 }, cellHeight{ 0 };

    // Render settings
    float offsetX{ 0 }, offsetY{ 0 }; // offset in the screen co-ords
    float zoom{ 1.0f };
    float zoom_sensitivity{ 0.1f };
    float minZoomLevel{ 0.01f };
    float velocityScale{ 300.0f };

    float minPressure{ 0.0 }, maxPressure{ 0.0 };

    void HandleInput() 
    {
        auto mouse_pos = GetMousePosition();

        // Zoom in/out using mouse wheel
        auto new_zoom = zoom * (1.0f + GetMouseWheelMove() * zoom_sensitivity);

        // Scale Velocity
        if (IsKeyDown(KEY_A))
            velocityScale *= 1.0f + zoom_sensitivity;
        else if (IsKeyDown(KEY_S))
            velocityScale *= 1.0f - zoom_sensitivity;

        // Keep the mouse position pointing at the same position
        auto old_mouse_grid = ScreenToGridPos(mouse_pos);
        offsetX *= new_zoom / zoom;
        offsetY *= new_zoom / zoom;
        zoom = new_zoom;
        auto new_mouse_grid = ScreenToGridPos(mouse_pos);
        offsetX -= (old_mouse_grid.x - new_mouse_grid.x) * zoom;
        offsetY -= (old_mouse_grid.y - new_mouse_grid.y) * zoom;

        // Pan using right mouse button drag
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) 
        {
            Vector2 mouseDelta = GetMouseDelta();
            offsetX += mouseDelta.x;
            offsetY += mouseDelta.y;
        }
    }

    Vector2 ScreenToGridPos(Vector2 mouse_pos)
    {
        return Vector2{ (mouse_pos.x + offsetX) / zoom, (mouse_pos.y + offsetY) / zoom };
    }

    Vector2 GridToScreenPos(Vector2 grid_pos)
    {
        return Vector2{ grid_pos.x * zoom + offsetX, grid_pos.y * zoom + offsetY };
    }

    void DrawGrid() 
    {
        cellWidth = zoom * dx;
        cellHeight = zoom * dy;

        // Draw grid lines
        for (int i = 0; i <= nx; ++i) 
        {
            float xPos = i * cellWidth + offsetX;
            DrawLine(xPos, offsetY, xPos, offsetY + ny * cellHeight, GRAY);
        }

        for (int j = 0; j <= ny; ++j) 
        {
            float yPos = j * cellHeight + offsetY;
            DrawLine(offsetX, yPos, offsetX + nx * cellWidth, yPos, GRAY);
        }
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
                if (magnitude > 0) { magnitude = log10(magnitude + 1.0f); }

                // Calculate the starting position of the vector (center of the cell)
                float xPos = i * cellWidth + cellWidth / 2 + offsetX;
                float yPos = j * cellHeight + cellHeight / 2 + offsetY;
                
                DrawArrowHead(Vector2{ (float)xPos, (float)yPos }, Vector2{ (float)uVal, (float)vVal }, magnitude * scale, BLACK);
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
    }

    void DrawField(const PhysicalField& f, double minValue, double maxValue)
    {
        for (int j = 0; j < ny; ++j) 
        {
            for (int i = 0; i < nx; ++i) 
            {
                double value = f(i, j);

                // Map the value to a color (heatmap)
                Color color = MapToColor(value, minValue, maxValue);

                // Draw a rectangle for each cell representing the value
                float xPos = i * cellWidth + offsetX;
                float yPos = j * cellHeight + offsetY;
                DrawRectangle(xPos, yPos, cellWidth, cellHeight, color);
            }
        }
    }

    Color MapToColor(double value, double minValue, double maxValue)
    {
        // Map the pressure value to a color (blue = low, red = high)
        float normValue = static_cast<float>(value);
        normValue = normValue < minValue ? minValue : (normValue > maxValue ? maxValue : normValue);
        normValue = (normValue - minValue) / (maxValue - minValue);
        return Color{ (unsigned char)(255 * (1.0f - normValue)),
                     (unsigned char)(255 * normValue),
                     0, 255 };
    }

    struct Streamline
    {
        std::vector<Vector2> points;
    };

    std::vector<Streamline> streamlines;

    void GenerateStreamlines(int numSeeds, int steps, float stepSize)
    {
        auto sampleVelocity = [&](float x, float y) -> Vector2 
        {
            int i = static_cast<int>(x);
            int j = static_cast<int>(y);

            if (i < 0 || j < 0 || i >= nx - 1 || j >= ny - 1)
                return { 0.0f, 0.0f }; // Out of bounds

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

        for (int seed = 0; seed < numSeeds; ++seed) 
        {
            // Generate random seed points within the field
            float x = GetRandomValue(0, nx - 1);
            float y = GetRandomValue(0, ny - 1);

            Streamline streamline;

            for (int step = 0; step < steps; ++step) 
            {
                Vector2 velocity = sampleVelocity(x, y);

                float magnitude = std::pow(std::pow(velocity.x, 2) + std::pow(velocity.y, 2), 0.5);

                if (magnitude < 1e-10)
                    break;

                // Normalize velocity and move
                velocity = { velocity.x * (stepSize / magnitude), velocity.y * (stepSize / magnitude) };

                float nextX = x + velocity.x;
                float nextY = y + velocity.y;

                // Store point in streamline
                streamline.points.push_back(Vector2{ x*dx, y*dy });

                x = nextX;
                y = nextY;

                // Stop if out of bounds
                if (x <= 0 || y <= 0 || x >= nx || y >= ny)
                    break;
            }

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
            }
        }
    }
};