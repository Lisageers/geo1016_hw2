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


/// convert a 3 by 1 vector of type 'std::vector<double>' to vec3
vec3 to_vec3(std::vector<double> &V) {
    vec3 result;
    for (int i = 0; i < 3; ++i) {
        result[i] = V[i];
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

Matrix<double> fundamental_matrix_estimation(const std::vector<vec3> &points_0, const std::vector<vec3> &points_1) {
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

    std::cout << "Fundamental Matrix F:\n" << F << "\n";
    return F;
}

mat34 computeProjection(mat3 K, mat3 &R, vec3 &t) {
    // COMPUTEPROJECTION
    // Takes
    //  - the 3x3 intrinsic parameters matrix
    //  - the 3x3 extrinsic rotation matrix
    //  - the 3x1 extrinsic translation vector
    // Outputs:
    //  - Projection Matrix

    // Initialize a 3x4 Matrix
    mat34 Rt(0.0);

    // Fill Matrix with columns from R and t
    for (size_t i = 0; i < R.num_columns(); ++i) {
        Rt.set_col(i, R.col(i));
    }
    Rt.set_col(R.num_columns(), t);

    // Calculate Projection Matrix
    mat34 M(0.0);
    M = K * Rt;

    return M;
}

vec3 triangulate(vec3 pt0, vec3 pt1, mat34 M0, mat34 M1) {
    /* TRIANGULATE
     * Triangulate a pair of points (pt0 and pt1) with the help of the projection matrices (M0 and M1)
     * Output: 3d coordinates of triangulated point
     */

    // Matrix A
    mat4 A(0.0);
    A.set_row(0, vec4(pt0.x * M0.row(2) - M0.row(0)));
    A.set_row(1, vec4(pt0.y * M0.row(2) - M0.row(1)));
    A.set_row(2, vec4(pt1.x * M1.row(2) - M1.row(0)));
    A.set_row(3, vec4(pt1.y * M1.row(2) - M1.row(1)));

    // Solution based on linear method
    // Compute SVD decomposition and construct point P
    Matrix<double> U(4, 4, 0.0), S(4, 4, 0.0), V(4, 4, 0.0);
    svd_decompose(to_Matrix(A), U, S, V);

    vec3 P = vec3(V(0,3),
                  V(1,3),
                  V(2,3))
                          / V(3,3);

    return P;
}

// TODO: Uncomment once issues with mat3 are resolved
//
//vec3 triangulate_step2 (vec3 pt0,
//                       vec3 pt1,
//                       mat3 K,
//                       mat3 R,
//                       vec3 t) {
//    // Triangulate a pair of points (pt0 and pt1) with projection matrices (M0 and M1), using linear method
//
//    mat34 M0 (1.0f);
//    M0 (2,3) = 0;
//    M0 = (K * M0);
//
//    mat34 M1(1.0f);
//    for (int i = 0; i < R.num_rows(); ++i) {
//        vec3 row = R.row(i);
//        M1.set_row(i, vec4(row[0], row[1], row[2], 0.0));
//    }
//
//    M1.set_col(3, t);
//    M1 = K * M1;
//
//    Matrix<double> A(4, 4, 0.0);
//    vec4 elem00 = pt0.x * M0.row(2) - M0.row(0);
//    vec4 elem10 = pt0.y * M0.row(2) - M0.row(1);
//    vec4 elem20 = pt1.x * M1.row(2) - M1.row(0);
//    vec4 elem30 = pt1.y * M1.row(2) - M1.row(1);
//
//    for (int i = 0; i < 4; ++i) {
//        A(0,i) = elem00[i];
//        A(1,i) = elem10[i];
//        A(2,i) = elem20[i];
//        A(3,i) = elem30[i];
//    }
//
//    Matrix<double> U (A.rows(), A.rows(), 0.0);
//    Matrix<double> S (A.rows(), A.cols(), 0.0);
//    Matrix<double> V(A.cols(), A.cols(), 0.0);
//
//    svd_decompose(A, U, S, V);
//
//    Matrix<double> P(4, 1, 0.0);
//    for (int i = 0; i < 4; ++i) {
//        P(i, 0) = V (i, A.cols() - 1);
//    }
//
//    Matrix<double> normalised = (P / P(3, 0));
//    vec3 point_3d = {float (normalised(0, 0)), float (normalised(1, 0)), float (normalised(2, 0))};
//
//    return point_3d;
//}


std::tuple<mat3, vec3> best_relative_pose (mat3 R_solution1,
                                           mat3 R_solution2,
                                           vec3 t_solution1,
                                           vec3 t_solution2,
                                           mat3 matrix_K,
                                           vec3 pts0,
                                           vec3 pts1) {

    int max_pts_infront = 0;
    mat3 R_final;
    vec3 t_final;
    mat3 next_R;
    vec3 next_t;

    for (int combinations = 0; combinations < 4; ++combinations) {
        int num_pts_infront = 0;
        if (combinations == 0) {
            next_R = R_solution1;
            next_t = t_solution1;
        }
        if (combinations == 1) {
            next_R = R_solution1;
            next_t = t_solution2;
        }
        if (combinations == 2) {
            next_R = R_solution2;
            next_t = t_solution1;
        }
        if (combinations == 3) {
            next_R = R_solution2;
            next_t = t_solution2;
        }

        // TODO: Uncomment once issues with mat3 are resolved
//        for (int i = 0; i < pts0.size(); ++i) {
//            vec3 first_pt = triangulate_step2 (pts0[i], pts1[i], matrix_K, next_R, next_t);
//            vec3 second_pt = (next_R * first_pt) + next_t;
//            if (first_pt.z > 0 && second_pt.z > 0) {
//                num_pts_infront ++;
//            }
//        }
//        if (num_pts_infront > max_pts_infront){
//            max_pts_infront = num_pts_infront;
//            R_final = next_R;
//            t_final = next_t;
//        }
    }
//    return std::make_tuple(R_final, t_final);
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
    /*
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
*/
    // TODO: delete all above demo code in the final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // check if the input is valid
    if (points_0.size() != points_1.size())
    {
        std::cout << "Input Validation FAIL: The input files are not the same size" << std::endl;
        return false;
    }
    else if (points_0.size() < 8 || points_1.size() < 8)
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
    mat3 matrix_F = to_mat3(F);

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // Essential matrix E = transpose(K) * F * K. where K is the intrinsic matrix
    // Define K matrix. Requires first defining skew?
    // For now, I have ignored skew, and replaced it with 0, as it does not seem to be relevant
    // K = | fx skew  cx |
    //     | 0   fy   cy |
    //     | 0   0     1 |
    // Define array for the values to be inserted into Matrix K
    std::vector<double> K_array = {fx, 0, cx, 0, fy, cy, 0, 0, 1}; // fx, fy, cx, cy are the focal lengths and principal points of both cameras
    // Define matrix K and insert values, print matrix K
    Matrix<double> K (3, 3, K_array.data());
    mat3 matrix_K = to_mat3(K);
    std::cout << "Camera Matrix K: \n" << K << std::endl;

    //Compute matrix K's inverse, invK
    Matrix<double> invK (3, 3);
    inverse(K, invK);
    mat3 matrix_invK = to_mat3(invK);

//    Check to see if inverse is correct
//    std::cout << "Test to see if inverse is correct, K * invK: \n" << K * invK << std::endl;

    // Calculate E using Fundamental matrix, F, and camera matrix K and its transpose
    mat3 E = transpose(matrix_K) * matrix_F * matrix_K;
    std::cout << "Essential Matrix E: \n" << E << std::endl;

    // Calculating R and t based on performing an SVD, where E = U * Sigma * transpose(V)
    // where U and V are orthogonal 3x3 matrices, and Sigma is a 3x3 diagonal matrix, diag(1,1,0)
    Matrix<double> U (3, 3, 0.0);
    Matrix<double> Sig (3, 3, 0.0);
    Matrix<double> V (3, 3, 0.0);
    Matrix<double> original_E = to_Matrix(E);

    // SVD of matrix E:
    svd_decompose(original_E, U, Sig, V);
//    std::cout << "U: \n" << U << std::endl;
//    std::cout << "Sig: \n" << Sig << std::endl;
//    std::cout << "V: \n" << V << std::endl;

    // Find the 4 candidate relative poses (based on SVD)
    // Initialise W and Z matrices for R and t calculation
    Matrix<double> W (3, 3, 0.0);
    W.set_row({0, -1, 0}, 0);
    W.set_row({1, 0, 0}, 1);
    W.set_row({0, 0, 1}, 2);
    mat3 W_matrix = to_mat3(W);
    std::cout << "W matrix: \n" << W_matrix << std::endl;

    Matrix<double> Z (3, 3, 0.0);
    Z.set_row({0, 1, 0}, 0);
    Z.set_row({-1, 0, 0}, 1);
    Z.set_row({0, 0, 0}, 2);
    mat3 Z_matrix = to_mat3(Z);
    std::cout << "Z matrix: \n" << Z_matrix << std::endl;

    mat3 matrix_U = to_mat3(U);
    mat3 matrix_V = to_mat3(V);

    // R_1 = UWV^(T) AND R_2 = UW^(T)V^(T)
    mat3 R_1 = matrix_U * W_matrix * transpose(matrix_V) * (-1);
    mat3 R_2 = matrix_U * transpose(W_matrix) * transpose(matrix_V) * (-1);

    // [t]_X = UZU^(T), t = ±u_3
    vec3 t_1 = matrix_U.col(2);
    vec3 t_2 = matrix_U.col(2) * (-1);

    std::cout << "First solution Rotation Matrix, R: \n" << R_1 << std::endl;
    std::cout << "Second solution Rotation Matrix, R: \n" << R_2 << std::endl;
    std::cout << "First solution Translation vector, t: \n" << t_1 << std::endl;
    std::cout << "Second solution Translation vector, t: \n" << t_2 << std::endl;

//     Check: Determinant of R's ought to be 1
    std::cout << "determinant R_1: \n" << determinant(R_1) << std::endl;
    std::cout << "determinant R_2: \n" << determinant(R_2) << std::endl;

    // Determining correct relative pose (4 options) by finding which R and t combination has the most points in front of the camera
    // relative position between the two cameras

    // TODO: Uncomment once issues with mat3 versus matrix double are fixed
    std::tuple<mat3, vec3> correct_pose = best_relative_pose (R_1, R_2, t_1, t_2, matrix_K, points_0, points_1);
    R = std::get<0>(correct_pose);
    t = std::get<1>(correct_pose);




    // PART 3 -- DMITRI
    // TODO: Uncomment once mat3 issues fixed
//    // Convert the matrices
//    R = to_mat3(R_2);
//    t = to_vec3(t_2);
//
//    // Calculate projection matrix M for camera 0
//    // M0 = K [I | 0]
//    mat3 R0 = {1,0,0,
//               0,1,0,
//               0,0,1};
//    vec3 t0 = vec3(0.0f);
//    mat34 M0 = computeProjection(to_mat3(K), R0, t0);
//    std::cout << M0 << std::endl;
//
//    // Calculate projection matrix for camera 1
//    // M1 = K [R | t]
//    mat34 M1 = computeProjection(to_mat3(K), R, t);
//    std::cout << M1 << std::endl;
//
//    // Reconstruct 3D points by triangulating the pairs
//    for (size_t i = 0; i < points_0.size(); ++i) {
//        // For every pair:
//
//        vec3 pt3d;
//        pt3d = triangulate(points_0[i], points_1[i], M0, M1);
//
//        std::cout << "Point " << i << ": \t" << pt3d << std::endl;
//
//        points_3d.push_back(pt3d);
//    }

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
