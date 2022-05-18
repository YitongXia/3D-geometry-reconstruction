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

void normalization(std::vector<Vector2D> &points_0,std::vector<Vector2D> &points_1)
{
    // calculate centers
    double t0x=0;
    double t0y=0;
    double t1x=0;
    double t1y=0;
    for (int i=0;i<points_0.size();++i)
    {
        t0x+=points_0[i].x();
        t0y+=points_0[i].y();
        t1x+=points_1[i].x();
        t1y+=points_1[i].y();
    }
    // get mean value
    t0x=t0x/points_0.size();
    t0y=t0y/points_0.size();
    t1x=t1x/points_1.size();
    t1y=t1y/points_1.size();

    // mean distance to center
    double dist0=0;
    double dist1=0;

    for(int i=0;i<points_0.size();++i)
    {
        double transformed_0_x=points_0[i].x()-t0x;
        double transformed_0_y=points_0[i].y()-t0y;
        double transformed_1_x=points_1[i].x()-t1x;
        double transformed_1_y=points_1[i].x()-t1y;

        dist0+= sqrt(sqrt(pow(transformed_0_x,2)+ pow(transformed_0_y,2)));
        dist1+= sqrt(sqrt(pow(transformed_1_x,2)+ pow(transformed_1_y,2)));
    }
    double mean_dist0=dist0/points_0.size();
    double mean_dist1=dist1/points_1.size();

    //compute scaling factor
    double scale_0=mean_dist0 * (1/ sqrt(2));
    double scale_1=mean_dist1 * (1/ sqrt(2));
    Matrix33 T_0=(1/scale_0,0,-t0x/scale_0,
                  0,1/scale_0,-t0y/scale_0,
                  0,0,1);
    Matrix33 T_1=(1/scale_1,0,-t0x/scale_1,
                  0,1/scale_1,-t0y/scale_1,
                  0,0,1);
    for(int i=0;i<points_0.size();++i)
    {
        std::cout<<"the original point is: "<<points_0[i].x()<<", "<<points_0[i].y()<<std::endl;
        std::cout<<"the original point is: "<<points_1[i].x()<<", "<<points_1[i].y()<<std::endl;
        points_0[i].x()=(points_0[i].x()-t0x)/scale_0;
        points_0[i].y()=(points_0[i].y()-t0y)/scale_0;
        points_1[i].x()=(points_1[i].x()-t1x)/scale_1;
        points_1[i].y()=(points_1[i].y()-t1y)/scale_1;
        std::cout<<"the scaled point is: "<<points_0[i].x()<<", "<<points_0[i].y()<<std::endl;
        std::cout<<"the acaled point is: "<<points_1[i].x()<<", "<<points_1[i].y()<<std::endl;
    }
}

Matrix33 estimate_fundamental_F(const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1)
{
    Matrix W;
    for(int i=0;i<points_0.size();++i)
    {
        Vector2D pt0=points_0[i];
        Vector2D pt1=points_1[i];
        Vector w = (pt0.x() * pt1.x(),pt0.y()*pt1.x(),pt1.x(),pt0.x() * pt1.y(),pt0.y()*pt1.y(),pt1.y(),pt0.x(),pt0.y(),1);
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

    Vector f = v.get_column(v.cols()-1);

    Matrix33 F;
    F.set_row(0,{f[0],f[1],f[2]});
    F.set_row(1,{f[3],f[4],f[5]});
    F.set_row(2,{f[6],f[7],f[8]});
    return F;
}

Matrix compute_essential_E(const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1)
{
    Matrix E;
    for(int i=0;i<points_0.size();++i)
    {
        Vector2D pt0=points_0[i];
        Vector2D pt1=points_1[i];
        Vector e = (pt0.x()*pt1.x(),pt1.x()*pt0.y(),pt1.x(),pt1.y()*pt0.x(),pt1.y()*pt0.y(),pt1.y(),pt0.x(),pt0.y(),1);
        E.set_row(i,e);
    }
    return E;
}

void R_t(Matrix &E) {
    // E=Udiag(1,1,0)VT=SR
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
    Matrix R1 = U* transpose(W) * V;
    Matrix R2=U * W * V;

    Vector t = U.get_column(U.cols()-1);
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

    std::vector<Vector2D> pt0=points_0;
    std::vector<Vector2D> pt1=points_1;
    normalization(pt0,pt1);
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