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

#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

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

// step size in seconds
static const double tstep = 360.0;

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
    std::vector<glm::dvec2> data_r;
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
        data_r.push_back(r[2] - r[1]);
        data_v.push_back(v[1]);
        data_kin.push_back(ekin[2]);
    }

    // output data (this is of course not the best way, but it works)
    for(unsigned int i=0; i<data_t.size(); i++) {
        std::cout << data_t[i] << "\t" << data_r[i].x << std::endl;
    }
}
