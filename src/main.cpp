/**************************************************************************
 *   This file is part of SME.                                            *
 *                                                                        *
 *   Author: Ivo Filot <ivo@ivofilot.nl>                                  *
 *                                                                        *
 *   SME is free software:                                                *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   SME is distributed in the hope that it will be useful,               *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <cairo/cairo.h>

// assume that all celestial bodies are in the same plane

// gravitational constant
static const double G = 6.673E-11;

// starting positions (aphelion) in m
static const std::vector<glm::dvec2> r0 = {
glm::dvec2(0.0, 0.0),                           // sun
glm::dvec2(149597870700 * 1.02, 0.0),           // earth
glm::dvec2(149597870700 * 1.02 + 405400000, 0.0)     // moon
};

// initial velocities in m/s
static const std::vector<glm::dvec2> v0 = {
glm::dvec2(0.0, 0.0),                     // sun
glm::dvec2(0.0, 29300.0),                 // earth
glm::dvec2(0.0, 29300.0 + 970.0)          // moon
};

// masses
static const std::vector<double> m =  {
1.989e30,       // sun
5.9722e24,      // earth
7.342e22        // moon
};

/**
 * @brief      calculate acceleration
 *
 * @param      a     force variable (output)
 * @param[in]  r     position
 * @param[in]  m     masses
 */
void calculate_acc(std::vector<glm::dvec2>& a,
                   const std::vector<glm::dvec2>& r,
                   const std::vector<double>& m) {

    for(unsigned int i=0; i<a.size(); i++) {
        glm::dvec2 acc(0.0, 0.0);

        for(unsigned int j=0; j<a.size(); j++) {
            if(i == j) {    // ignore self-interaction
                continue;
            }

            acc += G * m[j] * glm::normalize(r[j] - r[i]) / glm::distance2(r[i],r[j]);
        }

        a[i] = acc;
    }
}

/**
 * @brief      calculate kinetic energy
 *
 * @param      ekin  kinetic energy
 * @param[in]  v     velocities
 * @param[in]  v     velocities in previous timestep
 * @param[in]  m     masses
 */
void calculate_kin(std::vector<double>& ekin,
                   const std::vector<glm::dvec2>& v,
                   const std::vector<glm::dvec2>& vprev,
                   const std::vector<double>& m) {

    for(unsigned int i=0; i<v.size(); i++) {
        ekin[i] = 0.5 * m[i] * glm::length2(0.5 * (v[i] + vprev[i]));
    }
}

/**
 * @brief      create graph from simulation data
 *
 * @param[in]  pos       vector with position data
 * @param[in]  filename  filename
 */
void output_cairo_graph(const std::vector<glm::dvec2>& pos,
                        const std::string& filename) {
    static const double canvas_size = 1000;

    // create canvas
    cairo_surface_t *surface;
    cairo_t *cr;
    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, (unsigned int)canvas_size, (unsigned int)canvas_size);
    cr = cairo_create (surface);

    // set background to white
    cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
    cairo_paint(cr);

    // print circle in the center for reference
    cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, canvas_size / 2.0, canvas_size / 2.0);
    cairo_arc(cr, canvas_size / 2.0, canvas_size / 2.0, 10.0, 0, 2.0*M_PI);
    cairo_fill(cr);

    // establish maximum and minimum values
    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    for(unsigned int i=0; i<pos.size(); i++) {
        xmin = std::min(xmin, pos[i].x);
        xmax = std::max(xmax, pos[i].x);

        ymin = std::min(ymin, pos[i].y);
        ymax = std::max(ymax, pos[i].y);
    }
    const double range = std::max(xmax - xmin, ymax - ymin) * 1.1;

    // print data points as pixels on canvas
    const double r = 0.5f;                      // two pixel radius

    for(unsigned int i=0; i<pos.size(); i++) {
        double x = canvas_size / 2.0 + (pos[i].x / range) * canvas_size;
        double y = canvas_size / 2.0 + (pos[i].y / range) * canvas_size;

        // set color for the data points based on the index in the
        // vector range
        double col = 1.0 - (double)i / (double)pos.size();
        cairo_set_source_rgb (cr, 0.0, 0.0, col);   // draw in black

        cairo_move_to(cr, x, y);
        cairo_arc(cr, x, y, r, 0, 2.0*M_PI);
        cairo_fill(cr);
    }

    // Clean-up and create picture file
    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, filename.c_str());
    cairo_surface_destroy (surface);
}

/**
 * @brief      output data to data file
 *
 * @param[in]  t         vector of time steps
 * @param[in]  r         vector of positions (or velocities)
 * @param[in]  filename  filename
 */
void output_to_datafile(const std::vector<double>& t,
                        const std::vector<glm::dvec2>& r,
                        const std::string& filename) {
    std::ofstream out(filename);

    for(unsigned int i=0; i<t.size(); i++) {
        out << t[i] << "\t" << r[i].x << "\t" << r[i].y << std::endl;
    }

    out.close();
}

// step size in seconds
static const double tstep = 3600.0;

int main() {

    // initialize working variables
    std::vector<glm::dvec2> a(m.size());
    auto r = r0;
    auto v = v0;
    auto vprev = v;
    std::vector<double> ekin(m.size());
    std::vector<double> epot(m.size());

    // construct vectors to store data
    std::vector<double> data_t;
    std::vector<glm::dvec2> data_r_earth;
    std::vector<glm::dvec2> data_r_moon;
    std::vector<glm::dvec2> data_v;
    std::vector<double> data_kin;

    // loop over span of 10 years
    for(double t=0; t < 3600*24*365.25 * 10; t += tstep) {

        // calculate acceleration
        calculate_acc(a, r, m);

        // update velocity
        for(unsigned int j=0; j<v.size(); j++) {
            vprev[j] = v[j];
            v[j] += a[j] * tstep;
        }

        // update positions
        for(unsigned int j=0; j<v.size(); j++) {
            r[j] += v[j] * tstep;
        }

        // calculate kinetic energy
        calculate_kin(ekin, v, vprev, m);

        // collect data of interest; change indices here
        // to grab what you need!
        data_t.push_back(t);
        data_r_earth.push_back(r[1]);
        data_r_moon.push_back(r[2] - r[1]);
        data_v.push_back(v[1]);
        data_kin.push_back(ekin[2]);
    }

    // output to data files (can be opened in Origin or so)
    output_to_datafile(data_t, data_r_earth, "earth.dat");
    output_to_datafile(data_t, data_r_moon, "moon.dat");

    // output to graphs
    output_cairo_graph(data_r_earth, "earth.png");
    output_cairo_graph(data_r_moon, "moon.png");
}
