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

bool IsValidFloat(const char* text) {
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
            continue;
        }
        else {
            return false;
        }
    }
    return hasDigit;
}