#include "Helper.hpp"
#include <algorithm>
#include <cmath>

void DrawArrow(float x1, float y1, float x2, float y2, float size, bool double_headed, Color color)
{
    DrawLineEx(Vector2(x1, y1), Vector2(x2, y2), 2, color);
    auto DrawHead = [&size, &color](float x1, float y1, float x2, float y2)
        {
            float angle = atan2f(y2 - y1, x2 - x1);
            float ax = x1 + 0.75f * (x2 - x1);
            float ay = y1 + 0.75f * (y2 - y1);

            Vector2 p1 = { ax - size * sinf(angle), ay - size * cosf(angle) };
            Vector2 p2 = { ax + size * sinf(angle), ay + size * cosf(angle) };

            DrawTriangleSafe(Vector2(x2, y2), p2, p1, color);
        };
    DrawHead(x1, y1, x2, y2);
    if (double_headed) DrawHead(x2, y2, x1, y1);
}

void DrawDashedRectangle(float x, float y, float width, float height, int dash, int gap, float bold_width, Color color)
{
    if (width < 0.0f) { x += width; width = -width; }
    if (height < 0.0f) { y += height; height = -height; }

    for (int i = 0; i < width; i += dash + gap) 
    {
        DrawLineEx({ x + i, y }, { x + i + dash, y }, bold_width, color);                  // Top
        DrawLineEx({ x + i, y + height }, { x + i + dash, y + height }, bold_width, color); // Bottom
    }
    for (int i = 0; i < height; i += dash + gap) 
    {
        DrawLineEx({ x, y + i }, { x, y + i + dash }, bold_width, color);                  // Left
        DrawLineEx({ x + width, y + i }, { x + width, y + i + dash }, bold_width, color);   // Right
    }
}

void DrawBoundary(bool component_dir, bool boundary_dir, float xPos, float yPos, float cellWidth, float cellHeight, float bold_width, int type)
{
    float offset = 0.0;// 0.5f * bold_width;

    auto colour = BLACK;
    if (type == 1) // Fixed inflow
        colour = VIOLET;
    else if (type == 2) // Fixed outflow
        colour = LIME;
    else if (type == 3) // Open
        colour = GOLD;

    int dash_size = std::max((int)(0.01f * std::min(cellWidth, cellHeight)), 2);


    if (component_dir == 0)
    {
        DrawDashedRectangle(xPos, yPos - offset, (boundary_dir == 0 ? 1 : -1) * cellWidth / 2, cellHeight + 2 * offset, dash_size, dash_size, bold_width, colour);
        if (type == 3) DrawDashedRectangle(xPos, yPos - offset, (boundary_dir == 0 ? -1 : 1) * cellWidth / 2, cellHeight + 2 * offset, dash_size, dash_size, bold_width, colour);
        DrawLineEx({ xPos, yPos - offset }, { xPos, yPos + cellHeight + offset }, bold_width, colour);
        if ((type == 1 && boundary_dir == 0) || (type == 2 && boundary_dir == 1))
            DrawArrow(xPos - cellWidth / 4, yPos + cellHeight / 2, xPos + cellWidth / 4, yPos + cellHeight / 2, cellHeight * 0.1, type == 3, colour);
        else
            DrawArrow(xPos + cellWidth / 4, yPos + cellHeight / 2, xPos - cellWidth / 4, yPos + cellHeight / 2, cellHeight * 0.1, type == 3, colour);
    }
    else
    {
        DrawDashedRectangle(xPos - offset, yPos, cellWidth + 2 * offset, (boundary_dir == 0 ? 1 : -1) * cellHeight / 2, dash_size, dash_size, bold_width, colour);
        if (type == 3) DrawDashedRectangle(xPos - offset, yPos, cellWidth + 2 * offset, (boundary_dir == 0 ? -1 : 1) * cellHeight / 2, dash_size, dash_size, bold_width, colour);
        DrawLineEx({ xPos - offset, yPos }, { xPos + cellHeight + offset, yPos }, bold_width, colour);
        if ((type == 1 && boundary_dir == 0) || (type == 2 && boundary_dir == 1))
            DrawArrow(xPos + cellWidth / 2, yPos - cellHeight / 4, xPos + cellWidth / 2, yPos + cellHeight / 4, cellWidth * 0.1, type == 3, colour);
        else
            DrawArrow(xPos + cellWidth / 2, yPos + cellHeight / 4, xPos + cellWidth / 2, yPos - cellHeight / 4, cellWidth * 0.1, type == 3, colour);
    }
}

void DrawRectangleSafe(int x, int y, int width, int height, Color color)
{
    if (width < 0)
    {
        x += width;
        width = -width;
    }
    if (height < 0)
    {
        y += height;
        height = -height;
    }
    DrawRectangle(x, y, width, height, color);
}

void DrawTriangleSafe(Vector2 p1, Vector2 p2, Vector2 p3, Color color)
{
    if (((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)) > 0)
    {
        Vector2 temp = p2;
        p2 = p3;
        p3 = temp;
    }
    DrawTriangle(p1, p2, p3, color);
}

bool IsValidFloat(const char* text, bool pos) 
{
    bool hasDecimal = false;
    bool hasDigit = false;
    for (int i = 0; text[i] != '\0'; i++) 
    {
        if (text[i] == '.' && !hasDecimal) 
        {
            hasDecimal = true;
        }
        else if (isdigit(text[i])) 
        {
            hasDigit = true;
        }
        else if (i == 0 && text[i] == '-') 
        {
            if (pos)
                return false;
            else
                continue;
        }
        else 
        {
            return false;
        }
    }
    return hasDigit;
}

bool IsValidInt(const char* text, bool pos) 
{
    bool hasDigit = false;
    for (int i = 0; text[i] != '\0'; i++)
    {
        if (isdigit(text[i]))
        {
            hasDigit = true;
        }
        else if (i == 0 && text[i] == '-')
        {
            if (pos)
                return false;
            else
                continue;
        }
        else
        {
            return false;
        }
    }
    return hasDigit;
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

void DrawLegend(const char* text, float x1, float y1, float xsize, float ysize, float minValue, float maxValue)
{
    DrawText(text, x1 + xsize / 2 - MeasureText(text, 20) / 2, y1 - 30, 20, BLACK);

    float y2 = y1 + ysize;
    float x2 = x1 + xsize;

    int steps = 100;
    float stepSize = (y2 - y1) / steps;

    for (int i = 0; i < steps; i++) 
    {
        float yPos = y1 + i * stepSize;
        float value = maxValue + ((minValue - maxValue) * (float)i / steps);
        Color color = MapToColor(value, minValue, maxValue);
        DrawRectangle(x1, yPos, x2 - x1, stepSize + 1, color);
    }

    // Draw text labels for values
    int n = 4;
    for (int i = 0; i < n; i++) 
    {
        float t = (float)i / (n - 1);
        float pointY = y1 + t * (y2 - y1);
        DrawText(TextFormat("%.1e", maxValue - t * (maxValue - minValue)), x2 + 5, pointY - 10, 20, BLACK);
    }
}