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

vec2 find_centroid(const std::vector<vec3> &points)
{
    // initialise
    float totalx = 0;
    float totaly = 0;

    // go through all points
    for (int i=0;i < points.size(); i++)
    {
        // add value of cartesian x and y to total
        totalx += points[i][0] / points[i][2];
        totaly += points[i][1] / points[i][2];

    }
    // calculate mean in both dimensions
    vec2 centroid(totalx/points.size(), totaly/points.size());
    return centroid;
}

float compute_scaling_factor(const std::vector<vec3> &points, vec2 centroid)
{
    float scaling_factor = 0;
    float total_distance = 0;
    // add all distances
    for (int i=0; i < points.size(); i++)
    {
		// use cartesian coordinates
        vec2 d2_point(points[i][0]/points[i][2], points[i][1]/points[i][2]);
		total_distance += distance(centroid, d2_point);

    }
    // calculate scaling factor sqrt(2)/mean distance
    scaling_factor = sqrt(2) / (total_distance/points.size());
    return scaling_factor;
}

Matrix<double> fundamental_matrix_estimation(const std::vector<vec3> &points_0, const std::vector<vec3> &points_1)
{
    // initialise
    Matrix<double> Transformation_p0(3, 3, 0.0);
    Matrix<double> Transformation_p1(3, 3, 0.0);
	Matrix<double> scaling_p0(3, 3, 0.0);
	Matrix<double> scaling_p1(3, 3, 0.0);

    // find centroids
    vec2 centroid_p0 = find_centroid(points_0);
    vec2 centroid_p1 = find_centroid(points_1);
    // find scaling factors
    float scale_p0 = compute_scaling_factor(points_0, centroid_p0);
    float scale_p1 = compute_scaling_factor(points_1, centroid_p1);

    // construct transformation matrix p0 and p1
    Transformation_p0.set_row({1, 0, -centroid_p0[0]}, 0);
    Transformation_p0.set_row({0, 1, -centroid_p0[1]}, 1);
    Transformation_p0.set_row({0,0,1}, 2);

    Transformation_p1.set_row({1, 0, -centroid_p1[0]}, 0);
    Transformation_p1.set_row({0, 1, -centroid_p1[1]}, 1);
    Transformation_p1.set_row({0,0,1}, 2);

	// construct scaling matrix p0 and p1
    scaling_p0.set_row({scale_p0, 0, 0}, 0);
	scaling_p0.set_row({0, scale_p0, 0}, 1);
	scaling_p0.set_row({0,0,1}, 2);

	scaling_p1.set_row({scale_p1, 0, 0}, 0);
	scaling_p1.set_row({0, scale_p1, 0}, 1);
	scaling_p1.set_row({0,0,1}, 2);

    // convert Mat to mat3
	mat3 T_p0_mat3 = to_mat3(Transformation_p0);
    mat3 T_p1_mat3 = to_mat3(Transformation_p1);
	mat3 S_p0_mat3 = to_mat3(scaling_p0);
	mat3 S_p1_mat3 = to_mat3(scaling_p1);

	// normalise
    std::vector<vec2> normalised_p0;
    std::vector<vec2> normalised_p1;
    for (int i=0; i < points_0.size(); i++)
    {
        // apply transformation
        vec3 a = S_p0_mat3 * T_p0_mat3 * points_0[i];
        vec3 b = S_p1_mat3 * T_p1_mat3 * points_1[i];
        // convert to cartesian coordinates
        normalised_p0.push_back(vec2(a.x/a.z, a.y/a.z));
        normalised_p1.push_back(vec2(b.x/b.z, b.y/b.z));
    }


    // eight point algorithm
    // create W
    Matrix<double> W(normalised_p0.size(), 9, 0.0);
    for (int i=0; i < normalised_p0.size(); i++)
    {
        // get x and y of both points
        double u0 = normalised_p0[i][0];
        double v0 = normalised_p0[i][1];
        double u1 = normalised_p1[i][0];
        double v1 = normalised_p1[i][1];
        // calculate W values
        W.set_row({u0*u1, v0*u1, u1, u0*v1, v0*v1, v1, u0, v0, 1}, i);
    }

    // calculate F with svd
    Matrix<double> Uw(normalised_p0.size(), normalised_p0.size(), 0.0);   // initialized with 0s
    Matrix<double> Sw(normalised_p0.size(), 9, 0.0);   // initialized with 0s
    Matrix<double> Vw(9, 9, 0.0);   // initialized with 0s
    // Single Value Decomposition into U, S, and V
    svd_decompose(W, Uw, Sw, Vw);
    const auto F_vector = Vw.get_column(9-1);
    // Reshape vector f to matrix F
    Matrix<double> F(3, 3, F_vector.data());

    // solve with svd
    Matrix<double> U(3, 3, 0.0);   // initialized with 0s
    Matrix<double> S(3, 3, 0.0);   // initialized with 0s
    Matrix<double> V(3, 3, 0.0);   // initialized with 0s
    // Single Value Decomposition into U, S, and V
    svd_decompose(F, U, S, V);
    // make rank 2 approximation
    S(2, 2) = 0;
    // calculate Fq
    F = U * S * transpose(V);

    // denormalise using F= T'^T * Fq * T
    // calculate F
    F = transpose(scaling_p1 * Transformation_p1) * F * (scaling_p0 * Transformation_p0);
    // make last element 1
    float last_element_F = F(2, 2);
    for (int i=0;i<3;i++)
    {
        F(i, 0) = F(i, 0) / last_element_F;
        F(i, 1) = F(i, 1) / last_element_F;
        F(i, 2) = F(i, 2) / last_element_F;
    }

    std::cout << "F:\n" << F << "\n";
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
    /// get the '2'-nd row of M
//    const vec4 b = M.row(2);    // it assigns the requested row to a new vector b

    /// get the '1'-st column of M
    const vec3 c = M.col(1);    // it assigns the requested column to a new vector c

    /// modify the element value at row 2 and column 1 (Note the 0-based indices)
//    M(2, 1) = b.x;

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

    // Essential matrix E = F * K * transpose(K). where K is the intrinsic matrix
    // Define K matrix. Requires first defining skew?
    // K = | fx skew  cx |
    //     | 0   fy   cy |
    //     | 0   0     1 |
    // For now, I have ignored skew, and replaced it with 0, as it was not working out

    //Define Extrinsic Matrix, E, as 3x3 matrix, initialised with 0's
    // TODO: E should probably be mat3 instead of matrix double? Does this matter/affect anything?
    Matrix<double> E(3, 3, 0.0);

    std::vector<double> K_array = {fx, 0, cx, 0, fy, cy, 0, 0, 1}; // fx, fy, cx, cy are the focal lengths and principal points of both cameras
    Matrix<double> K(3, 3, K_array.data());
    std::cout << "Camera Matrix K: \n" << K << std::endl;

    //Compute matrix K's inverse, invK
    Matrix<double> invK(3, 3);
    inverse(K, invK);

//    Check to see if inverse is correct
//    std::cout << "Test to see if inverse is correct, K * invK: \n" << K * invK << std::endl;

    // Calculate E using Fundamental matrix, F, and camera matrix K and its inverse, K'
    E = (F * K * transpose(K));
    std::cout << "Essential Matrix E: \n" << E << std::endl;

    // Calculating R and t based on performing an SVD, where E = U * Sigma * transpose(V)
    // where U and V are orthogonal 3x3 matrices, and Sigma is a 3x3 diagonal matrix, diag(1,1,0)
    // SVD of matrix E:
    Matrix<double> U(E.rows(),E.rows(), 0.0);
    Matrix<double> Sig(E.rows(),E.rows(), 0.0);
    Matrix<double> V(E.rows(),E.rows(), 0.0);
    svd_decompose(E, U, Sig, V);
//    std::cout << "U: \n" << U << std::endl;
//    std::cout << "Sig: \n" << Sig << std::endl;
//    std::cout << "V: \n" << V << std::endl;

    // Find the 4 candidate relative poses (based on SVD)
    // Initialise W and Z matrices for R and t calculation
    Matrix<double> W_matrix(3, 3, 0.0);
    W_matrix.set_row({0, -1, 0}, 0);
    W_matrix.set_row({1, 0, 0}, 1);
    W_matrix.set_row({0, 0, 1}, 2);
    std::cout << "W: \n" << W_matrix << std::endl;

    Matrix<double> Z_matrix(3, 3, 0.0);
    Z_matrix.set_row({0, 1, 0}, 0);
    Z_matrix.set_row({-1, 0, 0}, 1);
    Z_matrix.set_row({0, 0, 0}, 2);
    std::cout << "Z: \n" << Z_matrix << std::endl;

    // R_1 = UWV^(T) AND R_2 = UW^(T)V^(T)
    Matrix<double> R_1 (U.cols(), U.cols() , 0.0);
    Matrix<double> R_2 (U.cols(), U.cols() , 0.0);
    R_1 = (U * W_matrix * transpose(V) * determinant(U * W_matrix * transpose(V)));
    R_2 = (U * transpose(W_matrix) * transpose(V) * determinant(U * transpose(W_matrix) * transpose(V)));

    // [t]_X = UZU^(T), t = ±u_3
    std::vector<double> t_1 = U.get_column(U.cols() - 1);
    std::vector<double> t_2 = U.get_column(U.cols() - 1) * (-1);

    std::cout << "First solution Rotation Matrix, R: \n" << R_1 << std::endl;
    std::cout << "Second solution Rotation Matrix, R: \n" << R_2 << std::endl;
    std::cout << "First solution Translation vector, t: \n" << t_1 << std::endl;
    std::cout << "Second solution Translation vector, t:: \n" << t_2 << std::endl;

    // TODO: Determinants of both R's results in -1, when they should be 1.
    //  Fixed this by using determinant in R_1 and R_2 calculations, but not sure if this is the correct method.
    std::cout << "determinant R_1: \n" << determinant(R_1) << std::endl;
    std::cout << "determinant R_2: \n" << determinant(R_2) << std::endl;



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
