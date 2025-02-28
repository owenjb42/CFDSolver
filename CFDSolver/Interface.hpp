#pragma once

#include "raygui.h"
#include "raylib.h"
#include <vector>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <random>

#include "FairMutex.hpp"
#include "Helper.hpp"

#include "PhysicalField.hpp"

struct Face
{
    bool operator==(const Face& other) const = default;
    bool component_dir /*0: u, 1: v*/; int i; int j;
};
struct Boundary : public Face
{
    bool operator==(const Boundary& other) const { return Face::operator==(other); }
    enum Type{FixedInflow = 1, FixedOutflow = 2, Open = 3};
    Type type{ Open };
    double optional_velocity{ 0.0 }, optional_temp{ 0.0 };
    bool boundary_dir{ 0 };// 0: +, 1: -
};
namespace std
{
    template <> struct hash<Face> { size_t operator()(const Face& f) const { return (hash<bool>()(f.component_dir) ^ (hash<int>()(f.i) << 1) ^ (hash<int>()(f.j) << 2)); } };
    template <> struct hash<Boundary> { size_t operator()(const Boundary& f) const { return (hash<bool>()(f.component_dir) ^ (hash<int>()(f.i) << 1) ^ (hash<int>()(f.j) << 2)); } };
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
        SetWindowMinSize(1000, 1000);
        SetTargetFPS(60);
        ResetView();
    }

    ~Interface()
    {
        CloseWindow();
    }

    void SetData(SolverStaggeredIMEXTemp& solver);
    void GetDataFromBuffer();

    void RenderModel()
    {
        in_model_creator_mode = true;
        is_solving = false;

        while (!WindowShouldClose() && in_model_creator_mode)
        {
            HandleInput();

            BeginDrawing();
            ClearBackground(RAYWHITE);

            DrawGrid();
            HandlePannel();

            EndDrawing();
        }

        GenerateSpeedPoints();
    }

    void RenderResults()
    {
        in_render_mode = true;
        is_solving = false;

        while (!WindowShouldClose() && in_render_mode)
        {
            GetDataFromBuffer();

            // setup only needed if there is new data 
            if (recalculate_auxiliary_data) { SetAuxData(); }

            HandleInput();

            BeginDrawing();
            ClearBackground(RAYWHITE);

            if (enable_flags[2]) DrawField();
            DrawGrid();
            if (enable_flags[1]) DrawVelocityVectors();
            if (enable_flags[0]) DrawStreamlines();
            HandlePannel();

            EndDrawing();
        }
    }

    // UI
    float sidebar_width = 400.0f;

    bool in_model_creator_mode{ false };
    bool in_render_mode{ false };
    bool is_solving{ false };

    // Model creator UI
    float key_hold_time[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    bool text_box_edit_mode[3] = { false, false, false };

    char velocity_input[32] = "1";
    char velocity_input2[32] = "1";
    char temp_input[32] = "20";
    bool direction_toggle[2] = { false, false };

    Vector2 drag_start{};
    bool dragging{ false };

    std::vector<Face> selected_faces;
    std::unordered_set<Face> blocked_faces;
    std::unordered_set<Boundary> boundary_faces;

    // Result infomation
    int nx{ 5 }, ny{ 5 };
    float dx{ 0.005f }, dy{ 0.005f };
    PhysicalField u_buffer, v_buffer, p_buffer, t_buffer, divergence_buffer, u_face_buffer, v_face_buffer;
    PhysicalField u, v, p, t, divergence, u_face, v_face;
    FairMutex mutex;
    float cellWidth{ 0 }, cellHeight{ 0 };
    float minPressure{ 0.0 }, maxPressure{ 0.0 }, minTemp{ 0.0 }, maxTemp{ 0.0 }, minVelocity{ 0.0 }, maxVelocity{ 0.0 };

    // Results UI
    bool enable_flags[3] = { true, true, true };

    char stremline_options[4][20] = { "num", "steps", "thickness", "reverse" };
    bool stremline_option_editmode[3] = { false, false, false };
    char stremline_option_inputs[3][32] = { "100", "100", "2" };
    int streamline_num{ 100 }, streamline_steps{ 100 };
    float streamline_thickness{ 2.0f };
    bool reverse_streamlines{ false };

    char v_arrow_options[20] = { "scale" };
    bool v_arrow_editmode = { false };
    char v_arrow_inputs[32] = { "1" };
    float v_arrow_scale{ 1.0f };

    char field_options[2][20] = { "Presure", "Temperature" };
    bool plot_p_or_t{ false };
    //bool plot_temperature_field{ true };
    //char field_options[21] = "Pressure;Temperature";
    //int field_combo_option;
      
    // Render settings
    float offsetX{ 0 }, offsetY{ 0 }; // offset in the screen co-ords
    float zoom{ 1.0f };
    float zoom_sensitivity{ 0.1f };
    float minZoomLevel{ 0.01f };
    float velocityScale{ 0.01f };

    bool needs_update{ true };
    bool recalculate_auxiliary_data{ false };

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
            if (IsKeyPressed(KEY_ENTER))
            {
                in_model_creator_mode = false;
                selected_faces.clear();
            }

            bool grid_shrank = false;

            key_hold_time[0] = IsKeyDown(KEY_LEFT) ? key_hold_time[0] + GetFrameTime() : 0.0f;
            key_hold_time[1] = IsKeyDown(KEY_RIGHT) ? key_hold_time[1] + GetFrameTime() : 0.0f;
            key_hold_time[2] = IsKeyDown(KEY_UP) ? key_hold_time[2] + GetFrameTime() : 0.0f;
            key_hold_time[3] = IsKeyDown(KEY_DOWN) ? key_hold_time[3] + GetFrameTime() : 0.0f;

            if (IsKeyPressed(KEY_LEFT) || key_hold_time[0] > 0.3f)
            {
                nx = std::max(3, nx - 1);
                grid_shrank = true;
            }
            if (IsKeyPressed(KEY_RIGHT) || key_hold_time[1] > 0.3f)
            {
                nx = std::min(200, nx + 1);
            }
            if (IsKeyPressed(KEY_UP) || key_hold_time[2] > 0.3f)
            {
                ny = std::max(3, ny - 1);
                grid_shrank = true;
            }
            if (IsKeyPressed(KEY_DOWN) || key_hold_time[3] > 0.3f)
            {
                ny = std::min(200, ny + 1);
            }

            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
            {
                if (!dragging)
                {
                    if (IsMouseInView())
                    {
                        drag_start = mouse_pos;
                        dragging = true;
                    }
                }
                else
                {
                    auto drag_end = Vector2(std::clamp(mouse_pos.x, 0.0f, GetViewScreenWidth()), std::clamp(mouse_pos.y, 0.0f, (float)GetScreenHeight()));

                    DrawRectangleSafe(drag_end.x, drag_end.y, drag_start.x - drag_end.x, drag_start.y - drag_end.y, SKYBLUE);
                    DrawRectangleLines(drag_end.x, drag_end.y, drag_start.x - drag_end.x, drag_start.y - drag_end.y, BLUE);

                    auto grid_pos_1 = ScreenToGridPos(drag_start);
                    auto grid_pos_2 = ScreenToGridPos(drag_end);

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

            if (grid_shrank)
                SanatiseFaceSets();
        }
        else
        {
            if (IsKeyPressed(KEY_ENTER))
                in_render_mode = false;

            // Scale Velocity
            if (IsKeyDown(KEY_UP))
                velocityScale *= 1.0f + zoom_sensitivity;
            if (IsKeyDown(KEY_DOWN))
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

    float GetViewScreenWidth()
    {
        return GetScreenWidth() - sidebar_width;
    }

    bool IsMouseInView()
    {
        auto pos = GetMousePosition();
        return pos.x > 0.0f && pos.x < GetViewScreenWidth() && pos.y > 0.0f && pos.y < GetScreenHeight();
    }

    void ResetView()
    {
        offsetX = offsetY = 0.0f;
        double max_scale = std::max(nx * dx / GetViewScreenWidth(), ny * dy / GetScreenHeight());
        zoom = 1 / max_scale;
    }

    void SanatiseFaceSets()
    {
        // Remove faces which no longer conform to the current grid size
        for (auto it = boundary_faces.begin(); it != boundary_faces.end(); )
        {
            bool is_pos_flow = (it->type == Boundary::FixedInflow && it->boundary_dir == 0);
            if (it->component_dir == 0)
            {
                if (it->i > (is_pos_flow ? nx - 1 : nx) || it->j > ny - 1)
                    it = boundary_faces.erase(it);
                else ++it;
            }
            else
            {
                if (it->i > nx - 1 || it->j > (is_pos_flow ? ny - 1 : ny))
                    it = boundary_faces.erase(it);
                else ++it;
            }
        }

        for (auto it = blocked_faces.begin(); it != blocked_faces.end(); ) 
        {
            if (it->component_dir == 0)
            {
                if (it->i > nx || it->j > ny - 1)
                    it = blocked_faces.erase(it);
                else ++it;
            }
            else
            {
                if (it->i > nx - 1 || it->j > ny)
                    it = blocked_faces.erase(it);
                else ++it;
            }
        }
    }

    void HandlePannel()
    {
        int x1 = GetViewScreenWidth(), x2 = GetScreenWidth(), y1 = 0.0f, y2 = GetScreenHeight();
        int width = x2 - x1;
        Rectangle panel = { x1, y1, x2 - x1, y2 - y1 };
        DrawRectangleRec(panel, LIGHTGRAY);
        DrawRectangleLinesEx(panel, 2, DARKGRAY);

        float current_y = y1 + 50;
        int text_size = 20;
        float local_size{ 0.0f };

        if (in_model_creator_mode)
        {
            // Blocked and Clear buttons
            int xoffset = 30, xsize = 150, ysize = 60;
            GuiSetStyle(DEFAULT, TEXT_SIZE, 20);
            if (GuiButton(Rectangle(x1 + xoffset, current_y, xsize, ysize), "Blocked")) 
            {
                for (const auto& face : selected_faces)
                    blocked_faces.insert(face);
            }
            if (GuiButton(Rectangle(x1 + width - xsize - xoffset, current_y, xsize, ysize), "Clear")) 
            {
                for (const auto& face : selected_faces)
                    blocked_faces.erase(face);

                for (const auto& face : selected_faces)
                    boundary_faces.erase(Boundary(face));
            }
            current_y += ysize + 50;

            DrawRectangleLinesEx(Rectangle(x1 + xoffset / 2, current_y, x2 - x1 - xoffset, 200), 2, BLACK);
            DrawText("Fixed Inflow", x1 + xoffset, y1 + 10 + current_y, 20, BLACK);
            DrawText("Velocity", x1 + xoffset + 140, y1 + 60 + current_y, 20, BLACK);
            DrawText("Temperature", x1 + xoffset + 140, y1 + 100 + current_y, 20, BLACK);
            if (GuiTextBox(Rectangle(x1 + xoffset, y1 + 55 + current_y, 120, 30), velocity_input, 32, text_box_edit_mode[0]))
            {
                text_box_edit_mode[0] = !text_box_edit_mode[0];
                if (!IsValidFloat(velocity_input, true)) strcpy(velocity_input, "");
            }
            if (GuiTextBox(Rectangle(x1 + xoffset, y1 + 95 + current_y, 120, 30), temp_input, 32, text_box_edit_mode[1]))
            {
                text_box_edit_mode[1] = !text_box_edit_mode[1];
                if (!IsValidFloat(temp_input, true)) strcpy(temp_input, "");
            }
            if (GuiButton(Rectangle(x1 + xoffset, y1 + 150 + current_y, 150, 40), direction_toggle[0] == false ? "+ve direction" : "-ve direction"))
            { 
                direction_toggle[0] = !direction_toggle[0];
            }
            if (GuiButton(Rectangle(x2 - 100 - xoffset, y1 + 150 + current_y, 100, 40), "ADD")) 
            {
                for (const auto& face : selected_faces)
                    boundary_faces.erase(Boundary(face));

                char* endptr;
                float temp = strtof(temp_input, &endptr);
                float vel = strtof(velocity_input, &endptr);

                for (const auto& face : selected_faces)
                {
                    // Check for invalid boundary selection and skip in this case
                    bool is_pos = direction_toggle[0] == 0;
                    if (face.component_dir == 0)
                    {
                        if (face.i > (is_pos ? nx - 1 : nx) || face.i < (is_pos ? 0 : 1)) continue;
                    }
                    else
                    {
                        if (face.j > (is_pos ? ny - 1 : ny) || face.j < (is_pos ? 0 : 1)) continue;
                    }

                    auto new_boundary = Boundary(face);
                    new_boundary.optional_temp = temp;
                    new_boundary.optional_velocity = vel;
                    new_boundary.boundary_dir = direction_toggle[0];
                    new_boundary.type = Boundary::FixedInflow;
                    boundary_faces.insert(new_boundary);
                }
            }
            current_y += 200 + 50;

            DrawRectangleLinesEx(Rectangle(x1 + xoffset / 2, current_y, x2 - x1 - xoffset, 100), 2, BLACK);
            DrawText("Open", x1 + xoffset, y1 + 10 + current_y, 20, BLACK);
            if (GuiButton(Rectangle(x2 - 100 - xoffset, y1 + 50 + current_y, 100, 40), "ADD"))
            {
                for (const auto& face : selected_faces)
                    boundary_faces.erase(Boundary(face));

                for (const auto& face : selected_faces)
                {
                    auto new_boundary = Boundary(face);
                    new_boundary.type = Boundary::Open;
                    boundary_faces.insert(new_boundary);
                }
            }
            current_y += 100 + 50;

            DrawRectangleLinesEx(Rectangle(x1 + xoffset / 2, current_y, x2 - x1 - xoffset, 150), 2, BLACK);
            DrawText("Fixed Outflow", x1 + xoffset, y1 + 10 + current_y, 20, BLACK);
            DrawText("Velocity", x1 + xoffset + 140, y1 + 60 + current_y, 20, BLACK);
            if (GuiTextBox(Rectangle(x1 + xoffset, y1 + 55 + current_y, 120, 30), velocity_input2, 32, text_box_edit_mode[2]))
            {
                text_box_edit_mode[2] = !text_box_edit_mode[2];
                if (!IsValidFloat(velocity_input2, true)) strcpy(velocity_input2, "");
            }
            if (GuiButton(Rectangle(x1 + xoffset, y1 + 100 + current_y, 150, 40), direction_toggle[1] == false ? "+ve direction" : "-ve direction"))
            {
                direction_toggle[1] = !direction_toggle[1];
            }
            if (GuiButton(Rectangle(x2 - 100 - xoffset, y1 + 100 + current_y, 100, 40), "ADD"))
            {
                for (const auto& face : selected_faces)
                    boundary_faces.erase(Boundary(face));

                for (const auto& face : selected_faces)
                {
                    char* endptr;
                    float vel = strtof(velocity_input2, &endptr);

                    // Check for invalid boundary selection and skip in this case
                    bool is_pos = direction_toggle[1] == 0;
                    if (face.component_dir == 0)
                    {
                        if (face.i > (is_pos ? nx - 1 : nx) || face.i < (is_pos ? 0 : 1)) continue;
                    }
                    else
                    {
                        if (face.j > (is_pos ? ny - 1 : ny) || face.j < (is_pos ? 0 : 1)) continue;
                    }

                    auto new_boundary = Boundary(face);
                    new_boundary.boundary_dir = direction_toggle[1];
                    new_boundary.optional_velocity = vel;
                    new_boundary.type = Boundary::FixedOutflow;
                    boundary_faces.insert(new_boundary);
                }
            }
        }
        else if (in_render_mode)
        {
            float x_total{ (float)width }, x_offset_outer{ 10.0f }, x_offset_inner{ 10.0f }, x_offset_block{ 200.0f }, x_rec_size = 60.0f, text_offset{ 5.0f };
            float size = 30.0f;
            float y_offset_outer{ 30.0f }, y_offset_inner{ 15.0f };
            float bold{ 2.0f };

            // Streamlines
            local_size = 5 * y_offset_inner + 4 * size;
            DrawRectangleLinesEx(Rectangle(x1 + x_offset_outer, current_y, x_total - 2 * x_offset_outer, local_size), bold, BLACK);
            GuiCheckBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + local_size / 2 - size / 2, size, size), "Streamlines", &enable_flags[0]);
            if (GuiTextBox(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 1 + size * 0, x_rec_size, size), stremline_option_inputs[0], 32, stremline_option_editmode[0]))
            {
                stremline_option_editmode[0] = !stremline_option_editmode[0];
                if (!IsValidInt(stremline_option_inputs[0], true)) strcpy(stremline_option_inputs[0], "100");
                int new_streamline_num = std::stoi(stremline_option_inputs[0]);
                if (streamline_num != new_streamline_num) 
                {
                    streamline_num = new_streamline_num;
                    GenerateSpeedPoints();
                }
            }
            DrawText(stremline_options[0], x1 + x_offset_block + x_rec_size + text_offset, current_y + y_offset_inner * 1 + size * 0, text_size, BLACK);
            if (GuiTextBox(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 2 + size * 1, x_rec_size, size), stremline_option_inputs[1], 32, stremline_option_editmode[1]))
            {
                stremline_option_editmode[1] = !stremline_option_editmode[1];
                if (!IsValidInt(stremline_option_inputs[1], true)) strcpy(stremline_option_inputs[1], "100");
                streamline_steps = std::stoi(stremline_option_inputs[1]);
            }
            DrawText(stremline_options[1], x1 + x_offset_block + x_rec_size + text_offset, current_y + y_offset_inner * 2 + size * 1, text_size, BLACK);
            if (GuiTextBox(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 3 + size * 2, x_rec_size, size), stremline_option_inputs[2], 32, stremline_option_editmode[2]))
            {
                stremline_option_editmode[2] = !stremline_option_editmode[2];
                if (!IsValidFloat(stremline_option_inputs[2], true)) strcpy(stremline_option_inputs[2], "2");
                streamline_thickness = std::stof(stremline_option_inputs[2]);
            }
            DrawText(stremline_options[2], x1 + x_offset_block + x_rec_size + text_offset, current_y + y_offset_inner * 3 + size * 2, text_size, BLACK);
            GuiCheckBox(Rectangle(x1 + x_offset_block + (x_rec_size - size) / 2, current_y + y_offset_inner * 4 + size * 3, size, size), stremline_options[3], &reverse_streamlines);
            current_y += local_size + y_offset_outer;

            // Cell velocity plot
            local_size = 2 * y_offset_inner + size;
            DrawRectangleLinesEx(Rectangle(x1 + x_offset_outer, current_y, x_total - 2 * x_offset_outer, local_size), bold, BLACK);
            GuiCheckBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + local_size / 2 - size / 2, size, size), "Cell Velocity", &enable_flags[1]);
            if (GuiTextBox(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 1 + size * 0, x_rec_size, size), v_arrow_inputs, 32, v_arrow_editmode))
            {
                v_arrow_editmode = !v_arrow_editmode;
                if (!IsValidFloat(v_arrow_inputs, true)) strcpy(v_arrow_inputs, "1");
                v_arrow_scale = std::stof(v_arrow_inputs);
            }
            current_y += local_size + y_offset_outer;

            // Field result plot
            local_size = 2 * y_offset_inner + 1 * size;
            DrawRectangleLinesEx(Rectangle(x1 + x_offset_outer, current_y, x_total - 2 * x_offset_outer, local_size), bold, BLACK);
            GuiCheckBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + local_size / 2 - size / 2, size, size), "Field Plot", &enable_flags[2]);
            //if (GuiCheckBox(Rectangle(x1 + x_offset_block + (x_rec_size - size) / 2, current_y + y_offset_inner * 1 + size * 0, size, size), field_options[0], &plot_pressure_field)) { plot_temperature_field = false; }
            //if (GuiCheckBox(Rectangle(x1 + x_offset_block + (x_rec_size - size) / 2, current_y + y_offset_inner * 2 + size * 1, size, size), field_options[1], &plot_temperature_field)) { plot_pressure_field = false; }
            GuiToggle(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 1 + size * 0, x_rec_size * 2.5, size), field_options[plot_p_or_t], &plot_p_or_t);
            current_y += local_size + y_offset_outer;

            // TODO ledged plot
        }

        int yoffset = 30, ysize = 50, xsize = 100;
        if (GuiButton(Rectangle((width - xsize) / 2 + x1, y2 - yoffset - ysize, xsize, ysize), in_model_creator_mode ? "Start" : (is_solving ? "Stop" : "Rest")))
        {
            in_render_mode = in_model_creator_mode = false;
            selected_faces.clear();
        }
    }

    void SetAuxData()
    {
        minPressure = *std::ranges::min_element(p.begin(), p.end());
        maxPressure = *std::ranges::max_element(p.begin(), p.end());

        maxVelocity = std::numeric_limits<float>::min();
        minVelocity = std::numeric_limits<float>::max();

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

        GenerateStreamlines();

        recalculate_auxiliary_data = false;
    }

    void DrawGrid()
    {
        // Draw Grid Lines
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

        // Draw Boundary & Blocked Faces
        float bold_width = std::max(zoom / 4000.0f, 2.0f);
        float ofset = 0.5f * bold_width;

        // Border outline
        DrawLineEx({ offsetX, offsetY - ofset },                   { offsetX, offsetY + ny * cellHeight + ofset }, bold_width, BLACK);
        DrawLineEx({ offsetX + nx * cellWidth, offsetY - ofset },  { offsetX + nx * cellWidth, offsetY + ny * cellHeight + ofset }, bold_width, BLACK);
        DrawLineEx({ offsetX - ofset, offsetY },                   { offsetX + nx * cellWidth + ofset, offsetY }, bold_width, BLACK);
        DrawLineEx({ offsetX - ofset, offsetY + ny * cellHeight }, { offsetX + nx * cellWidth + ofset, offsetY + ny * cellHeight }, bold_width, BLACK);

        for (const auto& face : blocked_faces)
        {
            float xPos = face.i * cellWidth + offsetX;
            float yPos = face.j * cellHeight + offsetY;
            if (face.component_dir == 0)
                DrawLineEx({ xPos, yPos - ofset }, { xPos, yPos + cellHeight + ofset }, bold_width, BLACK);
            else
                DrawLineEx({ xPos - ofset, yPos }, { xPos + cellHeight + ofset, yPos }, bold_width, BLACK);
        }

        for (const auto& boundary : boundary_faces)
        {
            float xPos = boundary.i * cellWidth + offsetX;
            float yPos = boundary.j * cellHeight + offsetY;

            DrawBoundary(boundary.component_dir, boundary.boundary_dir, xPos, yPos, cellWidth, cellHeight, bold_width, boundary.type);
        }

        for (const auto& face : selected_faces)
        {
            float xPos = face.i * cellWidth + offsetX;
            float yPos = face.j * cellHeight + offsetY;
            if (face.component_dir == 0)
                DrawLineEx({ xPos, yPos - ofset }, { xPos, yPos + cellHeight + ofset }, bold_width, RED);
            else
                DrawLineEx({ xPos - ofset, yPos }, { xPos + cellHeight + ofset, yPos }, bold_width, RED);
        }
    }

    void DrawVelocityVectors()
    {
        auto DrawArrowHead = [](Vector2 centre, Vector2 direction, float scale, Color color)
        {
            float angle = atan2(direction.y, direction.x);
            float arrowHeadAngle = 0.8f;

            Vector2 arrowPoint1 = Vector2{ centre.x + direction.x * scale, centre.y + direction.y * scale };
            Vector2 arrowPoint2 = Vector2{ centre.x - scale * cos(angle + arrowHeadAngle), centre.y - scale * sin(angle + arrowHeadAngle) };
            Vector2 arrowPoint3 = Vector2{ centre.x - scale * cos(angle - arrowHeadAngle), centre.y - scale * sin(angle - arrowHeadAngle) };
            DrawTriangle(arrowPoint1, arrowPoint2, arrowPoint3, color);
            DrawTriangleLines(arrowPoint1, arrowPoint2, arrowPoint3, BLACK);
        };

        float scale = velocityScale * zoom * v_arrow_scale;

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

    void DrawField()
    {
        const PhysicalField& f = plot_p_or_t == 0 ? p : t;
        const float minValue = plot_p_or_t == 0 ? minPressure : minTemp;
        const float maxValue = plot_p_or_t == 0 ? maxPressure : maxTemp;

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                Color color = MapToColor(f(i, j), minValue, maxValue);
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
    std::vector<std::pair<float, float>> seeds;

    void GenerateSpeedPoints()
    {
        seeds.clear();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dis(0.0f, 1.0f);
        for (int seed = 0; seed < streamline_num; ++seed)
        {
            // Generate random seed points within the field
            float x = dis(gen) * (float)nx;
            float y = dis(gen) * (float)ny;
            seeds.push_back({ x, y });
        }
    }

    void GenerateStreamlines()
    {
        auto sampleVelocity = [&](float x, float y) -> Vector2
        {
            if (x <= 0.5f || y <= 0.5f || x >= (float)(nx) - 0.5f || y >= (float)(ny) - 0.5f)
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

        float stepSize = 0.1f;

        for (auto [x, y] : seeds)
        {
            float sign = reverse_streamlines ? -1.0f : 1.0f;

            Streamline streamline;

            for (int step = 0; step < streamline_steps; ++step)
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
            streamline.is_reverse = reverse_streamlines;
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

                DrawLineEx(
                    { screenPoint.x, screenPoint.y },
                    { otherScreenPoint.x, otherScreenPoint.y },
                    streamline_thickness,
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