
#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;

void print_matrix33(Matrix33 &T0)
{
    std::cout<<"Matrix:"<<std::endl;
    std::cout<<"the size is "<<T0.rows()<<" * "<<T0.cols()<<std::endl;
    std::cout<<T0(0,0) <<", "<<T0(0,1) <<", "<<T0(0,2) <<"\n";
    std::cout<<T0(1,0) <<", "<<T0(1,1) <<", "<<T0(1,2) <<"\n";
    std::cout<<T0(2,0) <<", "<<T0(2,1) <<", "<<T0(2,2) <<std::endl;
    std::cout<<"\n";
}

void print_matrix44(Matrix44 &T0)
{
    std::cout<<"Matrix:"<<std::endl;
    std::cout<<"the size is "<<T0.rows()<<" * "<<T0.cols()<<std::endl;
    std::cout<<T0(0,0) <<", "<<T0(0,1) <<", "<<T0(0,2)<<","<<T0(0,3) <<"\n";
    std::cout<<T0(1,0) <<", "<<T0(1,1) <<", "<<T0(1,2)<<","<<T0(1,3)  <<"\n";
    std::cout<<T0(2,0) <<", "<<T0(2,1) <<", "<<T0(2,2)<<","<<T0(2,3)  <<"\n";
    std::cout<<T0(3,0) <<", "<<T0(3,1) <<", "<<T0(3,2)<<","<<T0(3,3)  <<std::endl;
    std::cout<<"\n";
}

void print_matrix34(Matrix34 &T0)
{
    std::cout<<"Matrix:"<<std::endl;
    std::cout<<"the size is "<<T0.rows()<<" * "<<T0.cols()<<std::endl;
    std::cout<<T0(0,0) <<", "<<T0(0,1) <<", "<<T0(0,2)<<","<<T0(0,3) <<"\n";
    std::cout<<T0(1,0) <<", "<<T0(1,1) <<", "<<T0(1,2)<<","<<T0(1,3)  <<"\n";
    std::cout<<T0(2,0) <<", "<<T0(2,1) <<", "<<T0(2,2)<<","<<T0(2,3)   <<std::endl;
    std::cout<<"\n";
}

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

std::vector<Vector2D> normalization(const std::vector<Vector2D> &points, Matrix33 &T)
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
    double dist=0;

    for(auto & i : points)
    {
        double transformed_1_x=i.x()-tx;
        double transformed_1_y=i.y()-ty;

        dist+= sqrt(pow(transformed_1_x,2)+ pow(transformed_1_y,2));
        //dist+= pow(transformed_1_x,2)+ pow(transformed_1_y,2);
    }

    double mean_dist=dist/points.size();

    //compute scaling factor
    double scale = sqrt(2)/mean_dist;

    std::vector<double> t0={sqrt(2)/mean_dist,0,-tx * sqrt(2)/mean_dist};
    std::vector<double> t1={0,scale,-ty*scale};
    std::vector<double> t2={0,0,1};

    T.set_row(0,t0);
    T.set_row(1,t1);
    T.set_row(2,t2);

    for(int i=0;i<points.size();++i) {
        normalized_pt.emplace_back((points[i].x() - tx) *scale,(points[i].y() - ty)*scale);
    }

    // for validate the mean distance is sqrt(2)
    double dist_1=0;
    for(auto point: normalized_pt)
    {
        dist_1+= sqrt(pow(point.x(),2)+ pow(point.y(),2));
    }

    dist_1=dist_1/normalized_pt.size();

    return normalized_pt;
}


// compute initial fundamental matrix
Matrix estimate_fundamental_F(const std::vector<Vector2D> &norm_pt_0,const std::vector<Vector2D> &norm_pt_1)
{
    int num_rows = norm_pt_0.size();
    int num_cols = 9;

    Matrix W(num_rows,num_cols,0.0);
    for(int i=0;i<norm_pt_0.size();++i)
    {
        Vector2D pt0=norm_pt_0[i];
        Vector2D pt1=norm_pt_1[i];

        double u1=norm_pt_0[i].x();
        double v1=norm_pt_0[i].y();

        double u1_=norm_pt_1[i].x();
        double v1_=norm_pt_1[i].y();

        std::vector<double> w ={u1*u1_,v1 * u1_,u1_,u1*v1_,v1*v1_,v1_,u1,v1,1};
        //std::vector<double> w = {pt0.x() * pt1.x(),pt0.y()*pt1.x(),pt1.x(),pt0.x() * pt1.y(),pt0.y()*pt1.y(),pt1.y(),pt0.x(),pt0.y(),1};
        W.set_row(i,w);
    }


    Matrix U(num_rows, num_rows, 0.0);   // initialized with 0s
    Matrix S(num_rows, num_cols, 0.0);   // initialized with 0s
    Matrix V(num_cols, num_cols, 0.0);   // initialized with 0s

    svd_decompose(W, U, S, V);


    Vector f=V.get_column(V.cols()-1);
    std::vector<double> t0={f[0],f[1],f[2]};
    std::vector<double> t1={f[3],f[4],f[5]};
    std::vector<double> t2={f[6],f[7],f[8]};

    Matrix33 initial_F;
    initial_F.set_row(0,t0);
    initial_F.set_row(1,t1);
    initial_F.set_row(2,t2);

    /// Compute the SVD decomposition of A.

    Matrix u(num_rows, num_rows, 0.0);   // initialized with 0s
    Matrix s(num_rows, num_cols, 0.0);   // initialized with 0s
    Matrix v(num_cols, num_cols, 0.0);   // initialized with 0s

    svd_decompose(initial_F,u,s,v);
    s(2,2)=0;
    Matrix33 F = u * s * v.transpose();
    return F;
}

// denormalization
Matrix33 denormalization(const Matrix33 &F, const Matrix33 &T0, const Matrix33 &T1)
{
    Matrix33 denormalization_F= T1.transpose() * F * T0;
    return denormalization_F;
}

//7 - scale scale_invariant F   where F(2,2) = 1.
Matrix33 scale_F(Matrix33 &F)
{
    double scale=F(2,2);
    return F/ scale;
}

//8 - calculate E and 4 Rt settings from it (slide 20 -27)

Matrix33 intrinsic(const double &fx, const double &fy,const double &cx, const double &cy)
{
    return Matrix33(fx,0,cx,
                    0,fy,cy,
                    0,0,1);
}

Matrix33 compute_E(Matrix33 &F, Matrix33 & K)
{
    Matrix33 E= K.transpose() * F * K;
    return E;
}

Matrix34 create_M(Matrix33 &R1, Vector t, Matrix33 &K)
{
    Matrix34 Rt_matrix1(R1(0,0),R1(0,1),R1(0,2),t[0],
                        R1(1,0),R1(1,1),R1(1,2),t[1],
                        R1(2,0),R1(2,1),R1(2,2),t[2]);
    return K*Rt_matrix1 ;
}

std::vector<Vector3D> compute_3d_coord(Matrix & R1, Vector &t, Matrix33 &K, const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1)
{
    Matrix34 Rt_matrix1(R1(0,0),R1(0,1),R1(0,2),t[0],
                        R1(1,0),R1(1,1),R1(1,2),t[1],
                        R1(2,0),R1(2,1),R1(2,2),t[2]);


    Matrix34 original_Rt(1,0,0,0,
                         0,1,0,0,
                         0,0,1,0);


    Matrix M = K * original_Rt;
    Matrix34 M_ = K * Rt_matrix1;

    Vector4D M1(M(0,0), M(0,1), M(0,2),M(0,3));
    Vector4D M2(M(1,0), M(1,1), M(1,2),M(1,3));
    Vector4D M3(M(2,0), M(2,1), M(2,2),M(2,3));

    Vector4D M1_(M_(0,0),M_(0,1),M_(0,2),M_(0,3));
    Vector4D M2_(M_(1,0),M_(1,1),M_(1,2),M_(1,3));
    Vector4D M3_(M_(2,0),M_(2,1),M_(2,2),M_(2,3));

    std::vector<Vector3D> pt_3;
    for(int i=0;i<points_0.size();++i)
    {
        Matrix44 A;
        A.set_row(0,points_0[i].x() * M3-M1);
        A.set_row(1,points_0[i].y() * M3-M2);
        A.set_row(2,points_1[i].x() * M3_ - M1_);
        A.set_row(3,points_1[i].y() * M3_ - M2_);

//        A.set_row(0, points_0[i].x() * M.get_row(2) - M.get_row(0));
//        A.set_row(1,points_0[i].y() * M.get_row(2) - M.get_row(1));
//        A.set_row(2,points_1[i].x() * M_.get_row(2) - M_.get_row(0));
//        A.set_row(3,points_1[i].y() * M_.get_row(2) - M_.get_row(1));

        int num_rows = A.rows();
        int num_cols=A.cols();

        Matrix u(num_rows, num_rows, 0.0);   // initialized with 0s
        Matrix s(num_rows, num_cols, 0.0);   // initialized with 0s
        Matrix v(num_cols, num_cols, 0.0);   // initialized with 0s

        svd_decompose(A,u,s,v);

        Vector4D last_v=v.get_column(v.cols()-1);

        pt_3.emplace_back();
        pt_3.back()=last_v.cartesian();
    }
    return pt_3;
}

int point_count(std::vector<Vector3D> pt_3d)
{
    int count=0;
    for(auto &item:pt_3d)
    {
        if(item.z() >0) count++;
    }
    return count;
}

std::vector<int> positive_z(const Matrix33 &R, const Vector3D &t, const std::vector<Vector3D> &pt_3d)
{
    int count_0 = 0;
    int count_1 = 0;
    for(auto &item:pt_3d)
    {
        if(item.z() >0) count_0++;
    }

    for(auto &p:pt_3d)
    {
        Vector3D p1 = R * p + t;
        if(p1.z()>0)
            count_1++;
    }
    std::vector<int> positive;
    positive.emplace_back(count_0);
    positive.emplace_back(count_1);
    return positive;
}


void R_t(Matrix &E, Matrix33 &R, Vector3D &t, double &fx, double &fy,double &cx, double &cy,const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1, std::vector<Vector3D> &points_3d)
{

    int num_rows = E.rows();
    /// get the number of columns.
    int num_cols = E.cols();

    Matrix33 W(0,-1,0,1,0,0,0,0,1);

    Matrix33 Z(0,1,0,-1,0,0,0,0,0);

    Matrix U(num_rows, num_rows, 0.0);

    Matrix S(num_cols, num_cols, 0.0);

    Matrix V(num_cols, num_cols, 0.0);

    svd_decompose(E,U,S,V);

    // R1=U WT VT or R2=U W VT

    Matrix33 R1= determinant(U * W * V.transpose()) * U * W * V.transpose();

    Matrix33 R2= determinant(U * W.transpose() * V.transpose()) * U * W.transpose() * V.transpose();

    Vector3D t1 = U.get_column(U.cols()-1);
    Vector3D t2 = -(U.get_column(U.cols() - 1));

    Matrix33 K= intrinsic(fx,fy,cx,cy);

    // for situation1:
    std::vector<std::vector<Vector3D>> pt_3d;

    pt_3d.emplace_back(compute_3d_coord(R1,t1,K,points_0,points_1));
    pt_3d.emplace_back(compute_3d_coord(R1,t2,K,points_0,points_1));
    pt_3d.emplace_back(compute_3d_coord(R2,t1,K,points_0,points_1));
    pt_3d.emplace_back(compute_3d_coord(R2,t2,K,points_0,points_1));
    R=R1;
    t=t2;
    points_3d=pt_3d[1];

    std::vector<std::vector<int>> all_point;
    all_point.emplace_back(positive_z(R1,t1,pt_3d[0]));
    all_point.emplace_back(positive_z(R1,t2,pt_3d[1]));
    all_point.emplace_back(positive_z(R2,t1,pt_3d[2]));
    all_point.emplace_back(positive_z(R2,t2,pt_3d[3]));

    if(all_point[0][0] == points_0.size() && all_point[0][1]==points_0.size())
    {
        R=R1;
        t=t1;
        points_3d=pt_3d[0];
    }
    else if(all_point[1][0] == points_0.size() && all_point[1][1]==points_0.size())
    {
        R=R1;
        t=t2;
        points_3d=pt_3d[1];
    }
    if(all_point[2][0] == points_0.size() && all_point[2][1]==points_0.size())
    {
        R=R2;
        t=t1;
        points_3d=pt_3d[2];
    }if(all_point[3][0] == points_0.size() && all_point[3][1]==points_0.size())
    {
        R=R2;
        t=t2;
        points_3d=pt_3d[3];
    }
}

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

    Matrix33 F1 = estimate_fundamental_F(normalized_pt0,normalized_pt1);

    Matrix33 F2=denormalization(F1,T0,T1);
    Matrix33 F= scale_F(F2);

    Matrix33 K= intrinsic(fx,fy,cx,cy);
    Matrix33 E = compute_E(F,K);

    R_t(E,R,t,fx,fy,cx,cy,points_0,points_1,points_3d);

    return points_3d.size() > 0;
}
