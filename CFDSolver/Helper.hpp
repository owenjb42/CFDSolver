#pragma once

#include "raylib.h"
#include <string>

void DrawRectangleSafe(int x, int y, int width, int height, Color color);

bool IsValidFloat(const char* text);