#include "Interface.hpp"
#include "Solver_Staggered_IMEX_Temp.hpp"

void Interface::SetData(SolverStaggeredIMEXTemp& solver)
{
    std::lock_guard lk(mutex);

    u_buffer = solver.cell_u;
    v_buffer = solver.cell_v;
    p_buffer = solver.p;
    t_buffer = solver.t;
    divergence_buffer = solver.divergence;
    u_face_buffer = solver.u;
    v_face_buffer = solver.v;

    residual_values_buffer[0] = solver.maxPressureResidual;
    residual_values_buffer[1] = solver.maxDivergence;
    residual_values_buffer[2] = solver.maxTemperatureResidual;

    needs_update = true;
    is_solving = true;
}

void Interface::GetDataFromBuffer()
{
    std::lock_guard lk(mutex);
    if (needs_update)
    {
        u = u_buffer;
        v = v_buffer;
        p = p_buffer;
        t = t_buffer;
        divergence = divergence_buffer;
        u_face = u_face_buffer;
        v_face = v_face_buffer;

        std::sprintf(residual_values[0], "%.1e", residual_values_buffer[0]);
        std::sprintf(residual_values[1], "%.1e", residual_values_buffer[1]);
        std::sprintf(residual_values[2], "%.1e", residual_values_buffer[2]);

        needs_update = false;
        recalculate_auxiliary_data = true;
    }
}

void Interface::HandlePannel()
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

        // Boundary Conditions
        float x_total{ (float)width }, x_offset_outer{ 10.0f }, x_offset_inner{ 15.0f }, x_offset_block{ 200.0f }, x_rec_size = 120.0f, text_offset{ 15.0f };
        float size = 30.0f;
        float y_offset_outer{ 30.0f }, y_offset_inner{ 10.0f }, y_offset_inner_small{ 10.0f };
        float bold{ 2.0f };
        float text_size{ 20.0f };

        // Fixed Inflow
        local_size = 5 * y_offset_inner + 4 * size;
        DrawRectangleLinesEx(Rectangle(x1 + x_offset_outer, current_y, x_total - 2 * x_offset_outer, local_size), bold, BLACK);
        DrawText("Fixed Inflow", x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner_small, text_size, BLACK);
        if (GuiTextBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner * 2 + size * 1, x_rec_size, size), fixed_in_inputs[0], 32, fixed_in_edit_mode[0]))
        {
            fixed_in_edit_mode[0] = !fixed_in_edit_mode[0];
            if (fixed_in_edit_mode[0]) strcpy(fixed_in_inputs[0], "");
            else if (IsValidFloat(fixed_in_inputs[0], true)) fixed_in_vel = std::stof(fixed_in_inputs[0]);
            else std::sprintf(fixed_in_inputs[0], "%.f", fixed_in_vel);
        }
        DrawText(fixed_in_options[0], x1 + x_offset_outer + x_offset_inner + x_rec_size + text_offset, current_y + y_offset_inner * 2 + size * 1, text_size, BLACK);
        if (GuiTextBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner * 3 + size * 2, x_rec_size, size), fixed_in_inputs[1], 32, fixed_in_edit_mode[1]))
        {
            fixed_in_edit_mode[1] = !fixed_in_edit_mode[1];
            if (fixed_in_edit_mode[1])  strcpy(fixed_in_inputs[1], "");
            else if (IsValidFloat(fixed_in_inputs[1], false)) fixed_in_temp = std::stof(fixed_in_inputs[1]);
            else std::sprintf(fixed_in_inputs[1], "%.f", fixed_in_temp);
        }
        DrawText(fixed_in_options[1], x1 + x_offset_outer + x_offset_inner + x_rec_size + text_offset, current_y + y_offset_inner * 3 + size * 2, text_size, BLACK);
        if (GuiButton(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner * 4 + size * 3, x_rec_size * 1.4, size), fixed_in_dir ? "-ve direction" : "+ve direction"))
        {
            fixed_in_dir = !fixed_in_dir;
        }
        if (GuiButton(Rectangle(x2 - x_offset_outer - x_offset_inner - x_rec_size * 0.7, current_y + y_offset_inner * 4 + size * 3, x_rec_size * 0.7, size), "ADD"))
        {
            for (const auto& face : selected_faces)
            {
                // Check for invalid boundary selection and skip in this case
                if ((face.component_dir == 0) && (face.i > (fixed_in_dir ? nx : nx - 1) || face.i < (fixed_in_dir ? 1 : 0)) ||
                    (face.component_dir == 1) && (face.j > (fixed_in_dir ? ny : ny - 1) || face.j < (fixed_in_dir ? 1 : 0)))
                    continue;

                boundary_faces.erase(Boundary(face));

                auto new_boundary = Boundary(face);
                new_boundary.optional_temp = fixed_in_temp;
                new_boundary.optional_velocity = fixed_in_vel;
                new_boundary.boundary_dir = fixed_in_dir;
                new_boundary.type = Boundary::FixedInflow;
                boundary_faces.insert(new_boundary);
            }
        }
        current_y += local_size + y_offset_outer;

        // Open
        local_size = 5 * y_offset_inner + 4 * size;
        DrawRectangleLinesEx(Rectangle(x1 + x_offset_outer, current_y, x_total - 2 * x_offset_outer, local_size), bold, BLACK);
        DrawText("Open", x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner_small, text_size, BLACK);
        if (GuiTextBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner * 2 + size * 1, x_rec_size, size), open_inputs[0], 32, open_edit_mode[0]))
        {
            open_edit_mode[0] = !open_edit_mode[0];
            if (open_edit_mode[0]) strcpy(open_inputs[0], "");
            else if (IsValidFloat(open_inputs[0], false)) open_pressure = std::stof(open_inputs[0]);
            else std::sprintf(open_inputs[0], "%.f", open_pressure);
        }
        DrawText(open_options[0], x1 + x_offset_outer + x_offset_inner + x_rec_size + text_offset, current_y + y_offset_inner * 2 + size * 1, text_size, BLACK);
        if (GuiTextBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner * 3 + size * 2, x_rec_size, size), open_inputs[1], 32, open_edit_mode[1]))
        {
            open_edit_mode[1] = !open_edit_mode[1];
            if (open_edit_mode[1])  strcpy(open_inputs[1], "");
            else if (IsValidFloat(open_inputs[1], false)) open_temp = std::stof(open_inputs[1]);
            else std::sprintf(open_inputs[1], "%.f", open_temp);
        }
        DrawText(open_options[1], x1 + x_offset_outer + x_offset_inner + x_rec_size + text_offset, current_y + y_offset_inner * 3 + size * 2, text_size, BLACK);
        if (GuiButton(Rectangle(x2 - x_offset_outer - x_offset_inner - x_rec_size * 0.7, current_y + y_offset_inner * 4 + size * 3, x_rec_size * 0.7, size), "ADD"))
        {
            for (const auto& face : selected_faces)
            {
                boundary_faces.erase(Boundary(face));

                auto new_boundary = Boundary(face);
                new_boundary.optional_temp = open_temp;
                new_boundary.optional_pressure = open_pressure;
                new_boundary.type = Boundary::Open;
                boundary_faces.insert(new_boundary);
            }
        }
        current_y += local_size + y_offset_outer;

        // Fixed Outflow
        local_size = 4 * y_offset_inner + 3 * size;
        DrawRectangleLinesEx(Rectangle(x1 + x_offset_outer, current_y, x_total - 2 * x_offset_outer, local_size), bold, BLACK);
        DrawText("Fixed Outflow", x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner_small, text_size, BLACK);
        if (GuiTextBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner * 2 + size * 1, x_rec_size, size), fixed_out_inputs, 32, fixed_out_edit_mode))
        {
            fixed_out_edit_mode = !fixed_out_edit_mode;
            if (fixed_out_edit_mode) strcpy(fixed_out_inputs, "");
            else if (IsValidFloat(fixed_out_inputs, true)) fixed_out_vel = std::stof(fixed_out_inputs);
            else std::sprintf(fixed_out_inputs, "%.f", fixed_out_vel);
        }
        DrawText(fixed_out_options, x1 + x_offset_outer + x_offset_inner + x_rec_size + text_offset, current_y + y_offset_inner * 2 + size * 1, text_size, BLACK);
        if (GuiButton(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + y_offset_inner * 3 + size * 2, x_rec_size * 1.4, size), fixed_out_dir ? "-ve direction" : "+ve direction"))
        {
            fixed_out_dir = !fixed_out_dir;
        }
        if (GuiButton(Rectangle(x2 - x_offset_outer - x_offset_inner - x_rec_size * 0.7, current_y + y_offset_inner * 3 + size * 2, x_rec_size * 0.7, size), "ADD"))
        {
            for (const auto& face : selected_faces)
            {
                // Check for invalid boundary selection and skip in this case
                if ((face.component_dir == 0) && (face.i > (fixed_out_dir ? nx : nx - 1) || face.i < (fixed_out_dir ? 1 : 0)) ||
                    (face.component_dir == 1) && (face.j > (fixed_out_dir ? ny : ny - 1) || face.j < (fixed_out_dir ? 1 : 0)))
                    continue;

                boundary_faces.erase(Boundary(face));

                auto new_boundary = Boundary(face);
                new_boundary.optional_velocity = fixed_out_vel;
                new_boundary.boundary_dir = fixed_out_dir;
                new_boundary.type = Boundary::FixedOutflow;
                boundary_faces.insert(new_boundary);
            }
        }
        current_y += local_size + y_offset_outer;
    }
    else if (in_render_mode)
    {
        float x_total{ (float)width }, x_offset_outer{ 10.0f }, x_offset_inner{ 10.0f }, x_offset_block{ 200.0f }, x_rec_size = 60.0f, text_offset{ 5.0f };
        float size = 30.0f;
        float y_offset_outer{ 30.0f }, y_offset_inner{ 15.0f }, y_legend_size{ 300.0f };
        float bold{ 2.0f };

        // Streamlines
        local_size = 5 * y_offset_inner + 4 * size;
        DrawRectangleLinesEx(Rectangle(x1 + x_offset_outer, current_y, x_total - 2 * x_offset_outer, local_size), bold, BLACK);
        GuiCheckBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + local_size / 2 - size / 2, size, size), "Streamlines", &enable_flags[0]);
        if (GuiTextBox(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 1 + size * 0, x_rec_size, size), stremline_option_inputs[0], 32, stremline_option_editmode[0]))
        {
            stremline_option_editmode[0] = !stremline_option_editmode[0];
            if (stremline_option_editmode[0]) strcpy(stremline_option_inputs[0], "");
            else if (!IsValidInt(stremline_option_inputs[0], true)) std::sprintf(stremline_option_inputs[0], "%.d", streamline_num);
            else
            {
                int new_streamline_num = std::stoi(stremline_option_inputs[0]);
                if (streamline_num != new_streamline_num)
                {
                    streamline_num = new_streamline_num;
                    GenerateSpeedPoints();
                }
            }
        }
        DrawText(stremline_options[0], x1 + x_offset_block + x_rec_size + text_offset, current_y + y_offset_inner * 1 + size * 0, text_size, BLACK);
        if (GuiTextBox(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 2 + size * 1, x_rec_size, size), stremline_option_inputs[1], 32, stremline_option_editmode[1]))
        {
            stremline_option_editmode[1] = !stremline_option_editmode[1];
            if (stremline_option_editmode[1]) strcpy(stremline_option_inputs[1], "");
            else if (IsValidInt(stremline_option_inputs[1], true)) streamline_steps = std::stoi(stremline_option_inputs[1]);
            else std::sprintf(stremline_option_inputs[1], "%.d", streamline_steps);
        }
        DrawText(stremline_options[1], x1 + x_offset_block + x_rec_size + text_offset, current_y + y_offset_inner * 2 + size * 1, text_size, BLACK);
        if (GuiTextBox(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 3 + size * 2, x_rec_size, size), stremline_option_inputs[2], 32, stremline_option_editmode[2]))
        {
            stremline_option_editmode[2] = !stremline_option_editmode[2];
            if (stremline_option_editmode[2]) strcpy(stremline_option_inputs[2], "");
            else if (IsValidFloat(stremline_option_inputs[2], true)) streamline_thickness = std::stoi(stremline_option_inputs[2]);
            else std::sprintf(stremline_option_inputs[2], "%.f", streamline_thickness);
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
            if (v_arrow_editmode) strcpy(v_arrow_inputs, "");
            else if (IsValidFloat(v_arrow_inputs, true)) v_arrow_scale = std::stof(v_arrow_inputs);
            else std::sprintf(v_arrow_inputs, "%.f", v_arrow_scale);
        }
        current_y += local_size + y_offset_outer;

        // Field result plot
        local_size = 2 * y_offset_inner + 1 * size;
        DrawRectangleLinesEx(Rectangle(x1 + x_offset_outer, current_y, x_total - 2 * x_offset_outer, local_size), bold, BLACK);
        GuiCheckBox(Rectangle(x1 + x_offset_outer + x_offset_inner, current_y + local_size / 2 - size / 2, size, size), "Field Plot", &enable_flags[2]);
        GuiToggle(Rectangle(x1 + x_offset_block, current_y + y_offset_inner * 1 + size * 0, x_rec_size * 2.5, size), field_options[plot_p_or_t], &plot_p_or_t);
        current_y += local_size + y_offset_outer * 2;

        // Ledged Plot
        if (enable_flags[1] && enable_flags[2])
        {
            DrawLegend(plot_p_or_t ? "Temperature" : "Presssure", x1 + x_offset_outer + x_offset_inner * 4, current_y, x_rec_size * 1.5, y_legend_size, plot_p_or_t ? minTemp : minPressure, plot_p_or_t ? maxTemp : maxPressure);
            DrawLegend("Velocity", x2 - x_offset_outer - x_offset_inner * 5 - x_rec_size * 1.5, current_y, x_rec_size * 1.5, y_legend_size, minVelocity, maxVelocity);
        }
        else if (enable_flags[2])
        {
            DrawLegend(plot_p_or_t ? "Temperature" : "Presssure", x1 + x_total / 2 - x_rec_size * 0.75, current_y, x_rec_size * 1.5, y_legend_size, plot_p_or_t ? minTemp : minPressure, plot_p_or_t ? maxTemp : maxPressure);
        }
        else if (enable_flags[1])
        {
            DrawLegend("Velocity", x1 + x_total / 2 - x_rec_size * 0.75, current_y, x_rec_size * 1.5, y_legend_size, minVelocity, maxVelocity);
        }
        current_y += y_legend_size + y_offset_outer * 2;

        // Residuals
        const char* title = "Residuals";
        float title_x = x1 + (x_total - MeasureText(title, text_size)) / 2;
        float title_y = current_y - text_size - 10;
        DrawText(title, title_x, title_y, text_size, BLACK);
        float title_width = MeasureText(title, text_size);
        DrawLine(title_x, title_y + text_size + 2, title_x + title_width, title_y + text_size + 2, BLACK);
        current_y += y_offset_outer;

        float box_x_size = (x_total - x_offset_outer * 4) / 3;
        for (int i = 0; i < 3; ++i) 
        {
            float box_x = x1 + x_offset_outer * (i + 1) + box_x_size * i;
            DrawRectangleLinesEx(Rectangle(box_x, current_y, box_x_size, size), bold, BLACK);

            float text_x = box_x + box_x_size / 2 - MeasureText(residuals[i], text_size) / 2;
            float text_y = current_y - text_size - 5;
            DrawText(residuals[i], text_x, text_y, text_size, BLACK);

            float value_x = box_x + box_x_size / 2 - MeasureText(residual_values[i], text_size) / 2;
            float value_y = current_y + (size - text_size) / 2;
            DrawText(residual_values[i], value_x, value_y, text_size, BLACK);
        }
    }

    int yoffset = 30, ysize = 50, xsize = 100;
    if (GuiButton(Rectangle((width - xsize) / 2 + x1, y2 - yoffset - ysize, xsize, ysize), in_model_creator_mode ? "Start" : (is_solving ? "Stop" : "Rest")))
    {
        in_render_mode = in_model_creator_mode = false;
        selected_faces.clear();
    }
}

void Interface::HandleInput()
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