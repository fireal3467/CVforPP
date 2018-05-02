#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <algorithm>
#include <queue>          // std::priority_queue
#include <set>


void dlt(std::vector<double> points)
{
    //two numbers per 2d point and 2points per translation 
    int num_points = points.size()/4;  //should be 4 for this case
    std::vector<R2Point> start_points; 
    std::vector<R2Point> end_points;
    //add the first 0-7 entries as the first 4 points
    for(int i = 0; i < num_points; i++){
        start_points.push_back(R2Point(points[(2*i)], points[(2*i)+1]));
    }
    //add the next 8 - 15 entries as the next 
    for(int i = num_points; i < 2*num_points; i++){
        end_points.push_back(R2Point(points[(2*i)], points[(2*i)+1]));
    }

    // build the 5x6 matrix of equations
    int num_rows = 2*num_points;
    int num_cols = 9;
    double** linEquations = dmatrix(1,num_rows,1,num_cols);

    for(int i = 1; i <= num_rows; i++){
        int current_point = (i-1)/2;
        if(i%2 == 1){
            linEquations[i][1] = 0.0;
            linEquations[i][2] = 0.0;
            linEquations[i][3] = 0.0;
            linEquations[i][4] = -1.0 * start_points[current_point][0];
            linEquations[i][5] = -1.0 * start_points[current_point][1];
            linEquations[i][6] = -1.0;
            linEquations[i][7] = end_points[current_point][1] * start_points[current_point][0];
            linEquations[i][8] = end_points[current_point][1] * start_points[current_point][1];
            linEquations[i][9] = end_points[current_point][1];
        } else {
            linEquations[i][1] = 1.0 * start_points[current_point][0];
            linEquations[i][2] = 1.0 * start_points[current_point][1];
            linEquations[i][3] = 1.0;
            linEquations[i][4] = 0.0;
            linEquations[i][5] = 0.0;
            linEquations[i][6] = 0.0;
            linEquations[i][7] = -end_points[current_point][0] * start_points[current_point][0];
            linEquations[i][8] = -end_points[current_point][0] * start_points[current_point][1];
            linEquations[i][9] = -end_points[current_point][0];
        }
    }

    double singularValues[10]; 
    double** nullspaceMatrix = dmatrix(1,num_rows,1,num_cols);
    svdcmp(linEquations, num_rows, num_cols, singularValues, nullspaceMatrix);

    // get the result
    printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f \n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6],singularValues[7],singularValues[8],singularValues[9]);

    // find the smallest singular value:
    int smallestIndex = 1;
    for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

    printf("\n smallest Index %d \n", smallestIndex);

    return; 
}


//takes in 4 inputs, returns the homogenous matrix 
int main(int argc, char **argv)
{
	double myPoints[] = {0,0,1,0,1,1,0,1,1,0,2,0,2,1,1,1};
  	std::vector<double> v (myPoints, myPoints + sizeof(myPoints) / sizeof(double));
	dlt(v);

}