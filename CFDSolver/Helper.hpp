#pragma once
#include "raylib.h"

void DrawRectangleSafe(int x, int y, int width, int height, Color color);

#include "raylib.h"
#include <string>

struct Button 
{
    Rectangle rect;
    Color color;
    std::string label;
    bool isHovered;
    bool isPressed;

    void Draw() 
    {
        DrawRectangleRec(rect, isHovered ? LIGHTGRAY : color);
        DrawRectangleLines(rect.x, rect.y, rect.width, rect.height, BLACK);
        DrawText(label.c_str(), rect.x + 10, rect.y + 10, 20, BLACK);
    }

    void Update(Vector2 mousePos) 
    {
        isHovered = CheckCollisionPointRec(mousePos, rect);
        isPressed = isHovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON);
    }
};

struct InputBox 
{
    Rectangle rect;
    std::string text;
    bool isActive;

    void Draw() 
    {
        DrawRectangleRec(rect, WHITE);
        DrawRectangleLines(rect.x, rect.y, rect.width, rect.height, BLACK);
        DrawText(text.c_str(), rect.x + 5, rect.y + 5, 20, BLACK);
    }

    void Update(Vector2 mousePos) 
    {
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) 
        {
            isActive = CheckCollisionPointRec(mousePos, rect);
        }

        if (isActive) 
        {
            int key = GetCharPressed();
            while (key > 0) {
                if ((key >= '0' && key <= '9') || key == '.' || key == '-') 
                {
                    text += (char)key;
                }
                key = GetCharPressed();
            }
            if (IsKeyPressed(KEY_BACKSPACE) && !text.empty()) 
            {
                text.pop_back();
            }
        }
    }
};