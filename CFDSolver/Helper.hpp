#pragma once

#include "raylib.h"
#include <string>

void DrawArrow(float x1, float y1, float x2, float y2, float size, bool double_headed, Color color);

void DrawBoundary(bool component_dir, bool boundary_dir, float xPos, float yPos, float cellWidth, float cellHeight, float bold_width, int type);

void DrawRectangleSafe(int x, int y, int width, int height, Color color);

void DrawTriangleSafe(Vector2 p1, Vector2 p2, Vector2 p3, Color color);

bool IsValidFloat(const char* text, bool pos = false);

bool IsValidInt(const char* text, bool pos = false);

Color MapToColor(float value, float minValue, float maxValue);