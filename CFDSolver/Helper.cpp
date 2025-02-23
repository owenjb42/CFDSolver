#include "Helper.hpp"

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

bool IsValidFloat(const char* text, bool pos) {
    bool hasDecimal = false;
    bool hasDigit = false;
    for (int i = 0; text[i] != '\0'; i++) {
        if (text[i] == '.' && !hasDecimal) {
            hasDecimal = true;
        }
        else if (isdigit(text[i])) {
            hasDigit = true;
        }
        else if (i == 0 && text[i] == '-') {
            if (pos)
                return false;
            else
                continue;
        }
        else {
            return false;
        }
    }
    return hasDigit;
}