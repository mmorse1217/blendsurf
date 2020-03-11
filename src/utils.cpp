
#include "nummatrix.hpp"
#include "numvector.hpp"
#include "vec2t.hpp"
#include "vec3t.hpp"
#include "utils.hpp"

NumVector bary_weights_1d(NumVector interpolation_nodes){
    int n = interpolation_nodes.length();
    NumVector weights(n);
    setvalue(weights, 1.0);

    // Compute w_i = 1/prod_{j \not = i} (x_i - x_j)
    for(int i = 0; i < n; i++){

        for(int j = 0; j < n; j++){
            if(j != i)
                weights(i) *= 1./(interpolation_nodes(j) - interpolation_nodes(i));
        }
    }
    return weights;

}

/*void bary2d(
        int dof,
        NumMatrix& xy_coordinates, 
        NumMatrix& function_values,
        int num_samples_1d,
        NumMatrix& refined_xy_coordinates,
        NumMatrix& refined_function_values){
    assert(refined_xy_coordinates.n() == refined_function_values.n());
    assert(xy_coordinates.n() == function_values.n());
    setvalue(refined_function_values, 0.);

     // Let s_ij = (x_i, y_j) be the sample points in xy_coordinates and
     // refined_xy_coordinates. We assume the sample points are on a 
     // tensor-product grid, i.e. x_i = y_i and (s_ik)_1 = (s_ij)_1 for any j,k


    NumVector interpolation_nodes_x(num_samples_1d);
    NumVector interpolation_nodes_y(num_samples_1d);
    for(int i =0; i < num_samples_1d; i++){
        interpolation_nodes_x(i) = xy_coordinates(0,i);
        interpolation_nodes_y(i) = xy_coordinates(1,i*num_samples_1d); 
    }


    // TODO implement 2d barycentric interpolation

    // O(n^2) precomputation of weights; we assume that the x and y
    // interpolation nodes are the same.
    NumVector barycentric_weights_x = 
        bary_weights_1d(interpolation_nodes_x);
    NumVector barycentric_weights_y = 
        bary_weights_1d(interpolation_nodes_y);

    // iterate  over coefficients of function values.... TODO bring inside loop
    // iteration over desired samples of interpolant
    for(int si =0; si < refined_xy_coordinates.n(); si++){
        Point2 eval_point(refined_xy_coordinates.clmdata(si));

        NumVector distance_to_nodes_x(num_samples_1d);
        NumVector distance_to_nodes_y(num_samples_1d);
        for(int i =0; i < num_samples_1d; i++){
            double node_x= interpolation_nodes_x(i);
            double node_y= interpolation_nodes_y(i); 
            Point2 interp_point(node_x,node_y);
            Point2 difference = eval_point - interp_point;

            // specify a minimum distance to prevent swamping 
            if(fabs(difference.x()) <= 1e-15)
                difference.x() = 1e-15;
            if(fabs(difference.y()) <= 1e-15)
                difference.y() = 1e-15;

            // invert to use dot product below
            distance_to_nodes_x(i) = 1./difference.x();
            distance_to_nodes_y(i) = 1./difference.y();
        }

        double denominator = 
            dot(barycentric_weights_x, distance_to_nodes_x)*
            dot(barycentric_weights_y, distance_to_nodes_y);

        NumVector numerator(dof);
        for(int i =0; i < num_samples_1d; i++){
            for(int j =0; j < num_samples_1d; j++){
                    double w_i = barycentric_weights_x(i);
                    double w_j = barycentric_weights_y(j);
                    
                    double x_minus_x_i_inv = distance_to_nodes_x(j);
                    double y_minus_y_j_inv = distance_to_nodes_y(i);
                    for(int d =0; d < dof; d++){
                        double f_ij = function_values(d, i*num_samples_1d + j);
                        
                        numerator(d) += w_i * w_j * 
                            x_minus_x_i_inv * y_minus_y_j_inv * f_ij;
                    }
            }
        }
        
        for(int d =0; d < dof; d++){
            refined_function_values(d,si) = numerator(d)/denominator;
        }
    }
}*/

void bary2d_2(Index2 idx,
        int dof,
        NumMatrix& xy_coordinates, 
        NumMatrix& function_values,
        //int num_samples_1d,
        //NumVector& barycentric_weights_x,
        //NumVector& barycentric_weights_y,
        double* barycentric_weights_x,
        //double* barycentric_weights_y,
        Point2& eval_point,
        Point3& interpolated_value){
    //assert(refined_xy_coordinates.n() == refined_function_values.n());
    assert(xy_coordinates.n() == function_values.n());
    //setvalue(refined_function_values, 0.);
    const int num_interp_nodes = 4;
    int num_total_samples_1d = int(sqrt(xy_coordinates.n()));

    // Let s_ij = (x_i, y_j) be the sample points in xy_coordinates and
    // refined_xy_coordinates. We assume the sample points are on a 
    // tensor-product grid, i.e. x_i = y_i and (s_ik)_1 = (s_ij)_1 for any j,k

    int base_idx = idx.x() -1;
    int base_idy = idx.y() -1;
    /*
       double interpolation_nodes_x[num_interp_nodes];
       double interpolation_nodes_y[num_interp_nodes];

       for(int i =0; i < num_interp_nodes; i++){
       interpolation_nodes_x[i] = 
       xy_coordinates(0,(base_idx+i)*num_total_samples_1d + base_idy);
       interpolation_nodes_y[i] = 
       xy_coordinates(1,(base_idx)*num_total_samples_1d + base_idy+i); 
       }*/


    // iterate  over coefficients of function values.... TODO bring inside loop
    // iteration over desired samples of interpolant
    double distance_to_nodes_x[num_interp_nodes];
    double distance_to_nodes_y[num_interp_nodes];
    for(int i =0; i < num_interp_nodes; i++){
        double node_x=
            xy_coordinates(0,(base_idx+i)*num_total_samples_1d + base_idy);

        double node_y= 
            xy_coordinates(1,(base_idx)*num_total_samples_1d + base_idy+i); 

        Point2 interp_point(node_x,node_y);
        Point2 difference = eval_point - interp_point;

        // specify a minimum distance to prevent swamping 
        if(fabs(difference.x()) < 1e-15)
            difference.x() = 1e-15;
        if(fabs(difference.y()) < 1e-15)
            difference.y() = 1e-15;

        //difference.x() = max(fabs(difference.x()),1e-15);;

        //difference.y() = max(fabs(difference.y()),1e-15);

        // invert to use dot product below
        distance_to_nodes_x[i] = 1./difference.x();
        distance_to_nodes_y[i] = 1./difference.y();
    }

    double denominator = 0;
    double weight_times_dist_x = 0;
    double weight_times_dist_y = 0;

    for(int i =0; i < num_interp_nodes; i++){
        weight_times_dist_x += barycentric_weights_x[i]*distance_to_nodes_x[i];
        weight_times_dist_y += barycentric_weights_x[i]*distance_to_nodes_y[i];
    }

    denominator = weight_times_dist_x*weight_times_dist_y;
    //double numerator[3] = {0.,0.,0.};
    Point3 numerator(0.);
    for(int i =0; i < 4; i++){

        double w_i = barycentric_weights_x[i];
        double x_minus_x_i_inv = distance_to_nodes_x[i];
        int x_index = (base_idx+i)*num_total_samples_1d;
        for(int j =0; j < 4; j++){

            double w_j = barycentric_weights_x[j];
            double y_minus_y_j_inv = distance_to_nodes_y[j];
            
            double coeff = w_i * w_j * 
                x_minus_x_i_inv * y_minus_y_j_inv;
            int func_idx = x_index+ j+base_idy;
            for(int d =0; d < 3; d++){

                double f_ij = function_values(d, func_idx);
                numerator(d) +=  coeff * f_ij;;
            }
        }
    }
    for(int d =0; d < 3; d++){

        interpolated_value(d)  =numerator(d)/denominator;
    }
}
void bary2d_subarray(
        int dof,
        int* idx, //starting indices of interpolant patch
        //Index2 width, 
        NumMatrix& xy_coordinates, 
        NumMatrix& function_values,
        NumVector& barycentric_weights_x,
        NumVector& barycentric_weights_y,
        Point2& eval_point,
        Point3& interpolant_value){
    assert(xy_coordinates.n() == function_values.n());
    interpolant_value *= 0.;
    const int num_samples_1d = 4;


    // iterate  over coefficients of function values.... 
    // iteration over desired samples of interpolant
    //for(int si =0; si < refined_xy_coordinates.n(); si++){
    int ix = idx[0];
    int iy = idx[1];
            cout << barycentric_weights_x << endl;
            cout << barycentric_weights_y<< endl;
        for(int d =0; d < dof; d++){

            double distance_to_nodes_x[num_samples_1d];
            double distance_to_nodes_y[num_samples_1d];
            for(int i =0; i < num_samples_1d; i++){
                distance_to_nodes_x[i] = 0.;
                distance_to_nodes_y[i] = 0.;
            }
            for(int i =0; i < num_samples_1d; i++){
                double node_x= xy_coordinates(0,ix+i);
                double node_y= xy_coordinates(1,(iy+i)*num_samples_1d); 
                Point2 interp_point(node_x,node_y);
                cout << interp_point << endl;
                Point2 difference = eval_point - interp_point;

                // specify a minimum distance to prevent swamping 
                if(fabs(difference.x()) <= 1e-15)
                    difference.x() = 1e-15;
                if(fabs(difference.y()) <= 1e-15)
                    difference.y() = 1e-15;

                // invert to use dot product below
                distance_to_nodes_x[i] = 1./difference.x();
                distance_to_nodes_y[i] = 1./difference.y();
            }
            double denominator = 0.;
            double x_weight = 0.;
            double y_weight = 0.;
            cout << "dist to nodes_x" << endl;
            for(int i =0; i < num_samples_1d; i++){
                cout << distance_to_nodes_x[i] << ", " ;
            }
            cout << endl;
            cout << "dist to nodes_y" << endl;
            for(int i =0; i < num_samples_1d; i++){
                cout << distance_to_nodes_y[i] << ", " ;
            }
            cout << endl;
            for(int i =0; i < num_samples_1d; i++){
                x_weight += barycentric_weights_x(i)*distance_to_nodes_x[i];
                y_weight += barycentric_weights_y(i)*distance_to_nodes_y[i];
                //dot(barycentric_weights_x, distance_to_nodes_x)*
                //dot(barycentric_weights_y, distance_to_nodes_y);

            }
            denominator = x_weight*y_weight;
            /*
            double denominator = 
                dot(barycentric_weights_x, distance_to_nodes_x)*
                dot(barycentric_weights_y, distance_to_nodes_y);
                */

            //NumVector numerator(dof);
            double numerator = 0.;
            for(int i =0; i < num_samples_1d; i++){
                for(int j =0; j < num_samples_1d; j++){
                    double w_i = barycentric_weights_x(i);
                    double w_j = barycentric_weights_y(j);

                    double x_minus_x_i_inv = distance_to_nodes_x[j];
                    double y_minus_y_j_inv = distance_to_nodes_y[i];
                    //for(int d =0; d < dof; d++){
                    double f_ij = function_values(d, (ix+i)*num_samples_1d + (iy+j));

                    numerator += w_i * w_j * 
                        x_minus_x_i_inv * y_minus_y_j_inv * f_ij;
                    //}
                }
            }
            cout << numerator << ", "<< denominator << endl;
            interpolant_value(d) = numerator/denominator;
        }
    //}
}
