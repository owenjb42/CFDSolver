#pragma once

#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <utility>

#include "PhysicalField.hpp"

// Function to generate an SVG with a grid and arrows based on the vector field
void generateVectorFieldWithGrid(const PhysicalField& u, const PhysicalField& v, const std::string& filename = "pic.svg", int cell_size = 50, double scale = 0.003) {
    if (u.nx != v.nx || u.ny != v.ny) {
        std::cerr << "Error: Field dimensions do not match!" << std::endl;
        return;
    }

    int nx = u.nx, ny = u.ny;
    int canvas_width = nx * cell_size + 100;
    int canvas_height = ny * cell_size + 100;

    std::ofstream svgFile(filename);
    if (!svgFile.is_open()) {
        std::cerr << "Error: Could not open file for writing!" << std::endl;
        return;
    }

    // Begin SVG file
    svgFile << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << canvas_width
        << "\" height=\"" << canvas_height << "\" viewBox=\"0 0 "
        << canvas_width << " " << canvas_height << "\">\n";
    svgFile << "<rect width=\"" << canvas_width << "\" height=\"" << canvas_height << "\" fill=\"white\"/>\n";

    // Draw the grid
    svgFile << "<g>\n";
    for (int i = 0; i <= nx; ++i) {
        int x = 50 + i * cell_size;
        std::string color = (i == 0 || i == nx) ? "red" : "black";
        svgFile << "<line x1=\"" << x << "\" y1=\"50\" x2=\"" << x << "\" y2=\"" << canvas_height - 50
            << "\" stroke=\"" << color << "\" stroke-width=\"1\"/>\n";
    }
    for (int j = 0; j <= ny; ++j) {
        int y = 50 + j * cell_size;
        std::string color = (j == 0 || j == ny) ? "red" : "black";
        svgFile << "<line x1=\"50\" y1=\"" << y << "\" x2=\"" << canvas_width - 50 << "\" y2=\"" << y
            << "\" stroke=\"" << color << "\" stroke-width=\"1\"/>\n";
    }
    svgFile << "</g>\n";

    // Draw the arrows with scaled arrowheads and line thickness
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double ux = u(i, j);
            double vy = -v(i, j);
            double magnitude = std::sqrt(ux * ux + vy * vy);

            // Skip drawing if the velocity is zero
            if (magnitude < 1e-6) continue;

            // Center of the cell
            double x1 = 50 + i * cell_size + cell_size / 2;
            double y1 = 50 + j * cell_size + cell_size / 2;

            // Arrow endpoint, scaled by vector magnitude
            double x2 = x1 + scale * ux * cell_size;
            double y2 = y1 - scale * vy * cell_size;

            // Scale stroke width and arrowhead size based on magnitude
            double stroke_width = magnitude * 4.0;
            double arrowhead_size = magnitude * 10.0;

            // Inline marker definition for dynamic scaling
            svgFile << "<defs>\n"
                << "<marker id=\"arrow" << j << "_" << i
                << "\" viewBox=\"0 0 10 10\" refX=\"5\" refY=\"5\" orient=\"auto\" markerWidth=\""
                << arrowhead_size << "\" markerHeight=\"" << arrowhead_size << "\">\n"
                << "<polygon points=\"0 0, 10 5, 0 10\" fill=\"black\"/>\n"
                << "</marker>\n"
                << "</defs>\n";

            svgFile << "<line x1=\"" << x1 << "\" y1=\"" << y1
                << "\" x2=\"" << x2 << "\" y2=\"" << y2
                << "\" stroke-width=\"" << stroke_width
                << "\" marker-end=\"url(#arrow" << j << "_" << i << ")\"/>\n";
        }
    }
    svgFile << "</svg>\n";

    svgFile.close();
}