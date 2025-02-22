#include "Helper.hpp"

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