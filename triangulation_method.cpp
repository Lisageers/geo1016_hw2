/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;


/// convert a 3 by 3 matrix of type 'Matrix<double>' to mat3
mat3 to_mat3(Matrix<double> &M) {
    mat3 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}


/// convert M of type 'matN' (N can be any positive integer) to type 'Matrix<double>'
template<typename mat>
Matrix<double> to_Matrix(const mat &M) {
    const int num_rows = M.num_rows();
    const int num_cols = M.num_columns();
    Matrix<double> result(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}

vec3 find_centroid(const std::vector<vec3> &points)
{
    float maxx, maxy;
    // initialise min like this bc otherwise init value is too small
    float minx = points[0][0];
    float miny = points[0][1];
    // go through all points
    for (int i=0;i < points.size(); i++)
    {
        // update min x if found
        if (points[i][0] < minx)
        {
            minx = points[i][0];
        }
        // update max x if found
        if (points[i][0] > maxx)
        {
            maxx = points[i][0];
        }
        // update max y if found
        if (points[i][1] < miny)
        {
            miny = points[i][1];
        }
        // update min z if found
        if (points[i][1] > maxy)
        {
            maxy = points[i][1];
        }
    }
    // calculate centers in both directions
    float x_center = ((maxx-minx)/2) + minx;
    float y_center = ((maxy-miny)/2) + miny;
    vec3 centroid = {x_center, y_center, 1};
    return centroid;
}

float compute_scaling_factor(const std::vector<vec3> &points, vec3 centroid)
{
    float scaling_factor;
    float total_distance = 0;
    // add all distances
    for (int i=0; i < points.size(); i++)
    {
        total_distance += norm(centroid-points[i]);
    }
    // calculate scaling factor sqrt(2)/mean distance
    scaling_factor = sqrt(2) / (total_distance/points.size());
    return scaling_factor;
}

const std::vector<vec3> normalising(const std::vector<vec3> &points, vec3 &centroid, float &scale)
{
    std::vector<vec3> norm_points;
    vec3 new_point;
    float x, y, z;
    // calculate qi = pi*T and put in vector
    for (int i=0; i < points.size(); i++)
    {
        x = scale*points[i][0] + centroid[0];
        y = scale*points[i][1] + centroid[1];
        z = 1;
        new_point = {x, y, z};
        norm_points.push_back(new_point);
    }
    return norm_points;
}

Matrix<double> fundamental_matrix_estimation(const std::vector<vec3> &points_0, const std::vector<vec3> &points_1)
{
    // normalisation
    Matrix<double> T_p0(3, 3, 0.0);
    Matrix<double> T_p1(3, 3, 0.0);

    // find centroids
    vec3 centroid_p0 = find_centroid(points_0);
    vec3 centroid_p1 = find_centroid(points_1);
    // find scaling factor
    float scale_p0 = compute_scaling_factor(points_0, centroid_p0); // or should i do for each dimension?
    float scale_p1 = compute_scaling_factor(points_1, centroid_p1);
    // compute normalised points
    std::vector<vec3> normalised_p0 = normalising(points_0, centroid_p0, scale_p0);
    std::vector<vec3> normalised_p1 = normalising(points_1, centroid_p1, scale_p1);

    // eight point algorithm
    // create W
    Matrix<double> W(points_0.size(), 9, 0.0);
    float v0, v1, v2, v3, v4, v5, v6, v7;
    for (int i=0; i < points_0.size(); i++)
    {
        v0 = points_0[i][0] * points_1[i][0];
        v1 = points_0[i][1] * points_1[i][0];
        v2 = points_1[i][0];
        v3 = points_0[i][0] * points_1[i][1];
        v4 = points_0[i][1] * points_1[i][1];
        v5 = points_1[i][1];
        v6 = points_0[i][0];
        v7 = points_0[i][1];
        W.set_row({v0, v1, v2, v3, v4, v5, v6, v7, 1}, i);
    }
    // calculate F
    // solve with svd
    Matrix<double> Uw(points_0.size(), points_0.size(), 0.0);   // initialized with 0s
    Matrix<double> Sw(points_0.size(), 9, 0.0);   // initialized with 0s
    Matrix<double> Vw(9, 9, 0.0);   // initialized with 0s
    // Single Value Decomposition into U, S, and V
    svd_decompose(W, Uw, Sw, Vw);
    const auto F_vector = Vw.get_column(9 - 1);

    // Reshape vector f to matrix F
    Matrix<double> F(3, 3, F_vector.data());

    std::cout << F << "F\n\n";

    // solve with svd
    Matrix<double> U(3, 3, 0.0);   // initialized with 0s
    Matrix<double> S(3, 3, 0.0);   // initialized with 0s
    Matrix<double> V(3, 3, 0.0);   // initialized with 0s
    // Single Value Decomposition into U, S, and V
    svd_decompose(F, U, S, V);
    // make rank 2 approximation
    S(2, 2) = 0;

    // calculate Fq
    F = U * S * V;
    std::cout << F << "Fq\n\n";

    // denormalise using F= T'^T * Fq * T
    // construct T
    T_p0.set_row({scale_p0, 0, centroid_p0[0]}, 0);
    T_p0.set_row({0, scale_p0, centroid_p0[1]}, 1);
    T_p0.set_row({0,0,1}, 1);
    // construct T'
    T_p1.set_row({scale_p1, 0, centroid_p1[0]}, 0);
    T_p1.set_row({0, scale_p1, centroid_p1[1]}, 1);
    T_p1.set_row({0,0,1}, 1);
    // calculate F
    F = transpose(T_p1) * F * T_p0;
    F(2, 2) = 1;
    std::cout << F << "F\n";
    return F;
}

/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'.
 */
bool Triangulation::triangulation(
        float fx, float fy,     /// input: the focal lengths (same for both cameras)
        float cx, float cy,     /// input: the principal point (same for both cameras)
        const std::vector<vec3> &points_0,    /// input: image points (in homogenous coordinates) in the 1st image.
        const std::vector<vec3> &points_1,    /// input: image points (in homogenous coordinates) in the 2nd image.
        std::vector<vec3> &points_3d,         /// output: reconstructed 3D points
        mat3 &R,   /// output: recovered rotation of 2nd camera (used for updating the viewer and visual inspection)
        vec3 &t    /// output: recovered translation of 2nd camera (used for updating the viewer and visual inspection)
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for each sub-task. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or feel free to put them in one or multiple separate files.

    /// Easy3D provides fixed-size matrix types, e.g., mat2 (2x2), mat3 (3x3), mat4 (4x4), mat34 (3x4).
    /// To use these matrices, their sizes should be known to you at the compile-time (i.e., when compiling your code).
    /// Once defined, their sizes can NOT be changed.
    /// In 'Triangulation/matrix.h', another templated 'Matrix' type is also provided. This type can have arbitrary
    /// dimensions and their sizes can be specified at run-time (i.e., when executing your program).
    /// Below are a few examples showing some of these data structures and related APIs.

    /// ----------- fixed-size matrices

    /// define a 3 by 4 matrix M (you can also define 3 by 4 matrix similarly)
    mat34 M(1.0f);  /// entries on the diagonal are initialized to be 1 and others to be 0.

    /// set the first row of M
    M.set_row(0, vec4(1,1,1,1));    /// vec4 is a 4D vector.

    /// set the second column of M
    M.set_col(1, vec4(2,2,2,2));

    /// get the 3 rows of M
    vec4 M1 = M.row(0);
    vec4 M2 = M.row(1);
    vec4 M3 = M.row(2);

    /// ----------- fixed-size vectors

    /// how to quickly initialize a std::vector
    std::vector<double> rows = {0, 1, 2, 3,
                                4, 5, 6, 7,
                                8, 9, 10, 11};
    /// get the '2'-th row of M
    const vec4 b = M.row(2);    // it assigns the requested row to a new vector b

    /// get the '1'-th column of M
    const vec3 c = M.col(1);    // it assigns the requested column to a new vector c

    /// modify the element value at row 2 and column 1 (Note the 0-based indices)
    M(2, 1) = b.x;

    /// apply transformation M on a 3D point p (p is a 3D vector)
    vec3 p(222, 444, 333);
    vec3 proj = M * vec4(p, 1.0f);  // use the homogenous coordinates. result is a 3D vector

    /// the length of a vector
    float len = p.length();
    /// the squared length of a vector
    float sqr_len = p.length2();

    /// the dot product of two vectors
    float dot_prod = dot(p, proj);

    /// the cross product of two vectors
    vec3 cross_prod = cross(p, proj);

    /// normalize this vector
    cross_prod.normalize();

    /// a 3 by 3 matrix (all entries are intentionally NOT initialized for efficiency reasons)
//    mat3 F;
    /// ... here you compute or initialize F.
    /// compute the inverse of K
//    mat3 invF = inverse(F);

    /// ----------- dynamic-size matrices

    /// define a non-fixed size matrix
    Matrix<double> W(2, 3, 0.0); // all entries initialized to 0.0.

    /// set its first row by a 3D vector (1.1, 2.2, 3.3)
    W.set_row({ 1.1, 2.2, 3.3 }, 0);   // here "{ 1.1, 2.2, 3.3 }" is of type 'std::vector<double>'

    /// get the last column of a matrix
    std::vector<double> last_column = W.get_column(W.cols() - 1);

    // TODO: delete all above demo code in the final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // check if the input is valid
    if (points_0.size() != points_1.size())
    {
        std::cout << "Input Validation FAIL: The input files are not the same size" << std::endl;
        return false;
    }
    else if (points_0.size() < 8 or points_1.size() < 8)
    {
        // Less than 8 points
        std::cout << "Input Validation FAIL: The input contains less than 8 points" << std::endl;
        return false;
    }
    else
    {
        // Everything all right
        std::cout << "Input Validation PASS." << std::endl;
    }

    // estimate the fundamental matrix F
    Matrix<double> F = fundamental_matrix_estimation(points_0, points_1);

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    //testing github push xenia again again...

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you to check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       However, there are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}