#pragma once

#include "raylib.h"
#include <vector>
#include <iostream>
#include <algorithm>

#include "PhysicalField.hpp"

class CFDVisualizer {
public:
    CFDVisualizer(int nx, int ny, float dx, float dy, const PhysicalField& u, const PhysicalField& v, const PhysicalField& p)
        : nx(nx), ny(ny), dx(dx), dy(dy), u(u), v(v), p(p) 
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

    void Render() {
        while (true)
        {
            HandleInput();

            BeginDrawing();
            ClearBackground(RAYWHITE);

            DrawField(p, minPressure, maxPressure);
            DrawGrid();
            DrawVelocityVectors();

            EndDrawing();
        }
    }

private:
    // Result infomation
    int nx, ny;
    double dx, dy;
    const PhysicalField& u, v, p;
    int cellWidth{ 0 }, cellHeight{ 0 };

    // Render settings
    int offsetX{ 0 }, offsetY{ 0 }; // offset in the screen co-ords
    float zoom{ 1.0f };
    float zoom_sensitivity{ 0.1f };
    float minZoomLevel{ 0.01f };
    float velocityScale{ 300.0f };

    float minPressure{ 0.0 }, maxPressure{ 0.0 };

    void HandleInput() 
    {
        auto mouse_pos = GetMousePosition();

        // Zoom in/out using mouse wheel
        //auto new_zoom_level = zoom_level * (1.0f + GetMouseWheelMove() * zoom_sensitivity);
        //if (new_zoom_level < min_zoom_level) new_zoom_level = min_zoom_level;

        float new_zoom{ zoom };
        if (IsKeyPressed(KEY_Z))
            new_zoom *= 1.0f + zoom_sensitivity;  // Zoom in
        else if (IsKeyPressed(KEY_X))
            new_zoom *= 1.0f - zoom_sensitivity;  // Zoom out

        // Scale Velocity
        if (IsKeyPressed(KEY_A))
            velocityScale *= 1.0f + zoom_sensitivity;
        else if (IsKeyPressed(KEY_S))
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
            offsetX += (int)mouseDelta.x;
            offsetY += (int)mouseDelta.y;
        }
    }

    Vector2 ScreenToGridPos(Vector2 mouse_pos)
    {
        return Vector2{ (mouse_pos.x + offsetX) / zoom, (mouse_pos.y + offsetY) / zoom };
    }

    Vector2 GridToScreenPos(Vector2 grid_pos)
    {
        return Vector2{ grid_pos.x * zoom - offsetX, grid_pos.y * zoom - offsetY };
    }

    void DrawGrid() 
    {
        cellWidth = zoom * dx;
        cellHeight = zoom * dy;

        // Draw grid lines
        for (int i = 0; i <= nx; ++i) 
        {
            int xPos = i * cellWidth + offsetX;
            DrawLine(xPos, offsetY, xPos, offsetY + ny * cellHeight, GRAY);
        }

        for (int j = 0; j <= ny; ++j) 
        {
            int yPos = j * cellHeight + offsetY;
            DrawLine(offsetX, yPos, offsetX + nx * cellWidth, yPos, GRAY);
        }
    }

    void DrawVelocityVectors() 
    {
        float scale = velocityScale * zoom;

        for (int j = 1; j < ny - 1; ++j) 
        {
            for (int i = 1; i < nx - 1; ++i) 
            {
                float uVal = static_cast<float>(u(i, j));
                float vVal = static_cast<float>(v(i, j));

                // Calculate the magnitude and apply logarithmic scaling
                float magnitude = sqrt(uVal * uVal + vVal * vVal);
                if (magnitude > 0) { magnitude = log10(magnitude + 1.0f); }

                // Calculate the starting position of the vector (center of the cell)
                int xPos = i * cellWidth + cellWidth / 2 + offsetX;
                int yPos = j * cellHeight + cellHeight / 2 + offsetY;

                float len = magnitude * scale;

                // Draw the arrow using two lines
                DrawArrow(Vector2{ (float)xPos, (float)yPos }, Vector2{ xPos + uVal * len, yPos + vVal * len }, len, BLACK);
            }
        }
    }

    void DrawArrow(Vector2 start, Vector2 end, float scale, Color color)
    {
        float angle = atan2(end.y - start.y, end.x - start.x);
        float arrowheadSize = scale / 2.0f;
        float lineThickness = scale / 10.0f;
        float arrowHeadAngle = 0.4f;

        Vector2 arrowPoint1 = Vector2{ end.x - arrowheadSize * cos(angle - arrowHeadAngle), end.y - arrowheadSize * sin(angle - arrowHeadAngle) };
        Vector2 arrowPoint2 = Vector2{ end.x - arrowheadSize * cos(angle + arrowHeadAngle), end.y - arrowheadSize * sin(angle + arrowHeadAngle) };
        Vector2 lineEnd = Vector2{ end.x - arrowheadSize * cos(arrowHeadAngle) * cos(angle), end.y - arrowheadSize * cos(arrowHeadAngle) * sin(angle) };

        // Draw the main arrow line
        DrawLineEx(start, lineEnd, lineThickness, color);

        // Draw the arrowhead
        DrawTriangle(end, arrowPoint2, arrowPoint1, color);
    }

    void DrawField(const PhysicalField& f, double minValue, double maxValue)
    {
        for (int j = 1; j < ny - 1; ++j) 
        {
            for (int i = 1; i < nx - 1; ++i) 
            {
                double value = f(i, j);

                // Map the value to a color (heatmap)
                Color color = MapToColor(value, minValue, maxValue);

                // Draw a rectangle for each cell representing the value
                int xPos = i * cellWidth + offsetX;
                int yPos = j * cellHeight + offsetY;
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
};