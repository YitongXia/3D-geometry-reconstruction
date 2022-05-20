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


// TODO: check if the input is valid (always good because you never known how others will call your function).
bool if_input_valid(const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1)
{
    auto pt0_size=points_0.size();
    auto pt1_size=points_1.size();
    if(pt0_size == pt1_size && pt0_size>=8 && pt1_size>=8)
        return true;
    else return false;
}

// TODO: Estimate relative pose of two views. This can be subdivided into
//      - estimate the fundamental matrix F;
//      - compute the essential matrix E;
//      - recover rotation R and t.

//1 - calculate the centers (corresponding pixel centers in each image separately)
//2 - mean distance to center (in each image separately)
//3 - compute scaling factor using mean distance  (separately slide 14)
//4 - compute initial Fundamental matrix with norm points (slide 12-14)
//5 - rank 2 enforcement  (slide 16,17)
//6 - calculate T1, and T2 (slide 14) think about how you are integrating it to F

//7 - scale scale_invariant F   where F(2,2) = 1.

//8 - calculate E and 4 Rt settings from it (slide 20 -27)
//9 - triangulate & compute inliers (slide 27)

//10 - choose best RT setting & evaluate

std::vector<Vector2D> normalization(const std::vector<Vector2D> &points,Matrix33 &T)
{
    std::vector<Vector2D> normalized_pt;
    // calculate centers
    double tx=0;
    double ty=0;
    for (auto & i : points)
    {
        tx+=i.x();
        ty+=i.y();
    }
    // get mean value
    tx=tx/points.size();
    ty=ty/points.size();

    // mean distance to center
    double dist1=0;

    for(auto & i : points)
    {
        double transformed_1_x=i.x()-tx;
        double transformed_1_y=i.x()-ty;
        dist1+= sqrt(pow(transformed_1_x,2)+ pow(transformed_1_y,2));
    }

    double mean_dist1=dist1/points.size();

    //compute scaling factor
    double scale = sqrt(2)/mean_dist1;

    T=(scale,0,-tx * scale,
            0,scale,-ty*scale,
            0,0,1);

    for(int i=0;i<points.size();++i) {
        normalized_pt.emplace_back((points[i].x() - tx) *scale,(points[i].y() - ty)*scale);
    }

    // for validate the mean distance is sqrt(2)
    double dist=0;
    for(auto point: normalized_pt)
    {
        dist+= sqrt(pow(point.x(),2)+ pow(point.y(),2));
    }
    dist=dist/normalized_pt.size();
    std::cout<<dist<<std::endl;
    return normalized_pt;
}


// compute initial fundamental matrix
Matrix33 estimate_fundamental_F(const std::vector<Vector2D> &norm_pt_0,const std::vector<Vector2D> &norm_pt_1)
{
    Matrix W(160,9,0.0);
    for(int i=0;i<norm_pt_0.size();++i)
    {
        Vector2D pt0=norm_pt_0[i];
        Vector2D pt1=norm_pt_1[i];
        std::vector<double> w = {pt0.x() * pt1.x(),pt0.y()*pt1.x(),pt1.x(),pt0.x() * pt1.y(),pt0.y()*pt1.y(),pt1.y(),pt0.x(),pt0.y(),1};
        W.set_row(i,w);

    }

    int num_rows = W.rows();
    /// get the number of columns.
    int num_cols = W.cols();

    Matrix u(num_rows, num_rows, 0.0);   // initialized with 0s
    Matrix s(num_rows, num_cols, 0.0);   // initialized with 0s
    Matrix v(num_cols, num_cols, 0.0);   // initialized with 0s

    /// Compute the SVD decomposition of A.
    svd_decompose(W, u, s, v);

    s(2,2)=0;

    Matrix F = u * s * v;

    return F;
}

// denormalization
Matrix33 denormalization(Matrix33 &F, Matrix33 &T0, Matrix33 &T1)
{
    Matrix33 denormalization_F= transpose(T1) * F * T0;
    return denormalization_F;
}
//7 - scale scale_invariant F   where F(2,2) = 1.
Matrix33 scale_F(Matrix33 &F)
{
    F(2,2)=1;
    return F;
}

//8 - calculate E and 4 Rt settings from it (slide 20 -27)
Matrix compute_E(Matrix33 &F, double &fx, double &fy,double &cx, double &cy)
{
    Matrix33 K=(fx,0,cx,0,fy,cy,0,0,1);
    Matrix33 E= transpose(K) * F * K;
    return E;
}

/** @brief  get relative position
*	@param 	E	essential matrix
*	@param 	R	output rotation matrix
*	@param 	t	translate matrix
*/
void R_t(Matrix &E, Matrix33 &R, Vector3D &t) {

    int num_rows = E.rows();
    /// get the number of columns.
    int num_cols = E.cols();

    Matrix33 W(0,-1,0,1,0,0,0,0,1);

    Matrix33 Z(0,1,0,-1,0,0,0,0,0);

    Matrix U(num_rows, num_rows, 0.0);

    Matrix33 diag(1,0,0,0,1,0,0,0,0);

    Matrix V(num_cols, num_cols, 0.0);

    svd_decompose(E,U,diag,V);

    // R1=U WT VT or R2=U W VT

    Matrix R1= determinant(U *W *V) *U *W*V;
    Matrix R2= determinant(U* transpose(W)*V)*U* transpose(W)*V;

    Vector t1 = U.get_column(U.cols()-1);
    Vector t2 = -(U.get_column(U.cols() - 1));

    Matrix34 projection_matrix1 = R1;
    projection_matrix1.set_column(3,t1);

}



// TODO: Reconstruct 3D points. The main task is
//      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)
void reconstruct() {

}
// TODO: Don't forget to
//          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
//          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
//            which can help you check if R and t are correct).
//       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
//       viewer will be notified to visualize the 3D points and update the view).
//       There are a few cases you should return 'false' instead, for example:
//          - function not implemented yet;
//          - input not valid (e.g., not enough points, point numbers don't match);
//          - encountered failure in any step.



/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;

    Matrix33 T0;
    Matrix33 T1;
    std::vector<Vector2D> normalized_pt0= normalization(points_0,T0);
    std::vector<Vector2D> normalized_pt1= normalization(points_1,T1);
    Matrix33 F = estimate_fundamental_F(normalized_pt0,normalized_pt1);
    denormalization(F,T0,T1);
    std::cout<<"F22 is "<<F(2,2)<<std::endl;


    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the matrix E;
    //      - recover rotation R and t.

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}