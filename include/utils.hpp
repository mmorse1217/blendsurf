#ifndef __BARY2D_HPP__
#define __BARY2D_HPP__
NumVector bary_weights_1d(NumVector interpolation_nodes);
void bary2d_2(Index2 idx,
        int dof,
        NumMatrix& xy_coordinates, 
        NumMatrix& function_values,
        //int num_samples_1d,
        double* barycentric_weights_x,
        //double* barycentric_weights_y,
        Point2& eval_point,
        Point3& interpolated_value);
void bary2d_subarray(
        int dof,
        int* idx, //starting indices of interpolant patch
        //Index2 width, 
        NumMatrix& xy_coordinates, 
        NumMatrix& function_values,
        NumVector& barycentric_weights_x,
        NumVector& barycentric_weights_y,
        Point2& eval_point,
        Point3& interpolant_value);
#endif
