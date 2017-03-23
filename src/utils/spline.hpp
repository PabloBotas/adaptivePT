/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */


#ifndef __TK_SPLINE_Hpp__
#define __TK_SPLINE_Hpp__

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <valarray>


// band matrix solver
class band_matrix
{
private:
    std::vector< std::vector<float> > m_upper;  // upper band
    std::vector< std::vector<float> > m_lower;  // lower band
public:
    band_matrix() {};                             // constructor
    band_matrix(int dim, int n_u, int n_l);       // constructor
    ~band_matrix() {};                            // destructor
    void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
    int dim() const;                             // matrix dimension
    int num_upper() const
    {
        return m_upper.size()-1;
    }
    int num_lower() const
    {
        return m_lower.size()-1;
    }
    // access operator
    float& operator () (int i, int j);        // write
    float  operator () (int i, int j) const;  // read
    // we can store an additional diogonal (in m_lower)
    float& saved_diag(int i);
    float  saved_diag(int i) const;
    void    lu_decompose();
    std::vector<float> r_solve(const std::vector<float>& b) const;
    std::vector<float> l_solve(const std::vector<float>& b) const;
    std::vector<float> lu_solve(const std::vector<float>& b,
                                 bool is_lu_decomposed=false);

};


// spline interpolation
class Spline_t
{
public:
    enum bd_type {
        first_deriv  = 1,
        second_deriv = 2
    };

private:
    std::vector<float> m_x,m_y;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
    std::vector<float> m_a, m_b, m_c;      // spline coefficients
    float  m_b0, m_c0;                     // for left extrapolation
    bd_type m_left, m_right;
    float  m_left_value, m_right_value;
    bool    m_force_linear_extrapolation;

public:
    Spline_t();

    // optional, but if called it has to come before set_points()
    void set_boundary(bd_type left, float left_value,
                      bd_type right, float right_value,
                      bool force_linear_extrapolation=false);
    void set_points(const std::vector<float>& x,
                    const std::vector<float>& y, bool cubic_spline=true);
    void set_points(const std::valarray<float>& x,
                    const std::valarray<float>& y, bool cubic_spline=true);
    float operator() (float x) const;
};


#endif
