// Source file for image class

// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <algorithm>
#include <queue>          // std::priority_queue
#include <set>



////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
    // fit a 2D conic to five points
    R2Point p1(1.2,3.5);
    R2Point p2(2.1,2.2);
    R2Point p3(0.2,1.6);
    R2Point p4(0.0,0.5);
    R2Point p5(-0.2,4.2);

    // build the 5x6 matrix of equations
    double** linEquations = dmatrix(1,5,1,6);

    linEquations[1][1] = p1[0]*p1[0];
    linEquations[1][2] = p1[0]*p1[1];
    linEquations[1][3] = p1[1]*p1[1];
    linEquations[1][4] = p1[0];
    linEquations[1][5] = p1[1];
    linEquations[1][6] = 1.0;

    linEquations[2][1] = p2[0]*p2[0];
    linEquations[2][2] = p2[0]*p2[1];
    linEquations[2][3] = p2[1]*p2[1];
    linEquations[2][4] = p2[0];
    linEquations[2][5] = p2[1];
    linEquations[2][6] = 1.0;

    linEquations[3][1] = p3[0]*p3[0];
    linEquations[3][2] = p3[0]*p3[1];
    linEquations[3][3] = p3[1]*p3[1];
    linEquations[3][4] = p3[0];
    linEquations[3][5] = p3[1];
    linEquations[3][6] = 1.0;
    
    linEquations[4][1] = p4[0]*p4[0];
    linEquations[4][2] = p4[0]*p4[1];
    linEquations[4][3] = p4[1]*p4[1];
    linEquations[4][4] = p4[0];
    linEquations[4][5] = p4[1];
    linEquations[4][6] = 1.0;

    linEquations[5][1] = p5[0]*p5[0];
    linEquations[5][2] = p5[0]*p5[1];
    linEquations[5][3] = p5[1]*p5[1];
    linEquations[5][4] = p5[0];
    linEquations[5][5] = p5[1];
    linEquations[5][6] = 1.0;

    printf("\n Fitting a conic to five points:\n");
    printf("Point #1: %f,%f\n",p1[0],p1[1]);
    printf("Point #2: %f,%f\n",p2[0],p2[1]);
    printf("Point #3: %f,%f\n",p3[0],p3[1]);
    printf("Point #4: %f,%f\n",p4[0],p4[1]);
    printf("Point #5: %f,%f\n",p5[0],p5[1]);

    // compute the SVD
    double singularValues[7]; // 1..6
    double** nullspaceMatrix = dmatrix(1,6,1,6);
    svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

    // get the result
    printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

    // find the smallest singular value:
    int smallestIndex = 1;
    for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

    // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
    printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

    // make sure the solution is correct:
    printf("Equation #1 result: %f\n",  p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] + 
                                        p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] + 
                                        p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] + 
                                        p1[0]*nullspaceMatrix[4][smallestIndex] + 
                                        p1[1]*nullspaceMatrix[5][smallestIndex] + 
                                        nullspaceMatrix[6][smallestIndex]);

    printf("Equation #2 result: %f\n",  p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] + 
                                        p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] + 
                                        p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] + 
                                        p2[0]*nullspaceMatrix[4][smallestIndex] + 
                                        p2[1]*nullspaceMatrix[5][smallestIndex] + 
                                        nullspaceMatrix[6][smallestIndex]);

    printf("Equation #3 result: %f\n",  p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] + 
                                        p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] + 
                                        p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] + 
                                        p3[0]*nullspaceMatrix[4][smallestIndex] + 
                                        p3[1]*nullspaceMatrix[5][smallestIndex] + 
                                        nullspaceMatrix[6][smallestIndex]);

    printf("Equation #4 result: %f\n",  p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] + 
                                        p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] + 
                                        p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] + 
                                        p4[0]*nullspaceMatrix[4][smallestIndex] + 
                                        p4[1]*nullspaceMatrix[5][smallestIndex] + 
                                        nullspaceMatrix[6][smallestIndex]);

    printf("Equation #5 result: %f\n",  p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] + 
                                        p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] + 
                                        p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] + 
                                        p5[0]*nullspaceMatrix[4][smallestIndex] + 
                                        p5[1]*nullspaceMatrix[5][smallestIndex] + 
                                        nullspaceMatrix[6][smallestIndex]);

    R2Point test_point(0.34,-2.8);

    printf("A point off the conic: %f\n",   test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] + 
                                            test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] + 
                                            test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] + 
                                            test_point[0]*nullspaceMatrix[4][smallestIndex] + 
                                            test_point[1]*nullspaceMatrix[5][smallestIndex] + 
                                            nullspaceMatrix[6][smallestIndex]);

    return; 
}

//DLT Algorithm
//input format 
//a,b,c,d,a',b',c',d'
std::vector<double> dlt(std::vector<double> points)
{
    // printf("\n starting \n");
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

    //option to print the matrix after everything
    // for(int i = 1; i<=num_rows;i++){
    //     for(int j = 1; j<=num_cols;j++){
    //         printf("%d ",(int)linEquations[i][j]);
    //     }
    //     printf("\n\n");
    // }


    double singularValues[10]; 
    double** nullspaceMatrix = dmatrix(1,num_rows+1,1,num_cols);

    // printf("\n finished creating pre-svdcmp\n");

    svdcmp(linEquations, num_rows, num_cols, singularValues, nullspaceMatrix);

    // get the result
    // printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f \n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6],singularValues[7],singularValues[8],singularValues[9]);

    // find the smallest singular value:
    int smallestIndex = 1;
    for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

    // printf("\n smallest Index %d \n", smallestIndex);

    printf("\n NullMatrix is \n %f, %f, %f, \n %f, %f, %f, \n %f, %f, %f \n ",
        nullspaceMatrix[1][smallestIndex],
        nullspaceMatrix[2][smallestIndex],
        nullspaceMatrix[3][smallestIndex],
        nullspaceMatrix[4][smallestIndex],
        nullspaceMatrix[5][smallestIndex],
        nullspaceMatrix[6][smallestIndex],
        nullspaceMatrix[7][smallestIndex],
        nullspaceMatrix[8][smallestIndex],
        nullspaceMatrix[9][smallestIndex]);

    std::vector<double> toReturn;
    for(int i = 1; i<10;i++){
        toReturn.push_back((double)nullspaceMatrix[i][smallestIndex]);
    }


    return toReturn; 
}

// std::vector<double> inverse(std::vector<double> v){

//     float determinant = 0;
//     for(i = 0; i < 3; i++)
//         determinant = determinant + (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));
    
//     cout<<"\n\ndeterminant: "<<determinant;
    
//     cout<<"\n\nInverse of matrix is: \n";
//     for(i = 0; i < 3; i++){
//         for(j = 0; j < 3; j++)
//             cout<<((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/ determinant<<"\t";
        
//         cout<<"\n";
//     }
// }

////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelX(void)
{
    // Apply the Sobel oprator to the image in X direction
    R2Pixel *new_pixels;
    new_pixels = new R2Pixel [ npixels ];
    assert(new_pixels);
    
    for (int i = 0; i < width; i++){
        for (int j = 0; j < height; j++){
            //TODO: Make this prettier
            //bottom left
            if (i == 0 && j == 0) {
                new_pixels[i*height + j] = -2*Pixel(i+1,j) + -1*Pixel(i+1,j+1);
            }
            //bottom right
            else if (i == width-1 && j == 0){
                new_pixels[i*height + j] = 2*Pixel(i-1,j) + 1*Pixel(i-1,j+1);
            }
            //top left
            else if (i == 0 && j == height-1){
                new_pixels[i*height + j] = -2*Pixel(i+1,j) + -1*Pixel(i+1,j-1);
            }
            //top right
            else if (i == width-1 && j == height-1 ){
                new_pixels[i*height + j] = 2*Pixel(i-1,j) + 1*Pixel(i-1,j-1);
            }
            //left edge
            else if (i == 0) {
                new_pixels[i*height + j] = -1*Pixel(i+1,j-1) + -2*Pixel(i+1,j) + -1*Pixel(i+1,j+1);
            }
            //right edge
            else if (i == width-1) {
                new_pixels[i*height + j] = 1*Pixel(i-1,j-1) + 2*Pixel(i-1,j) + 1*Pixel(i-1,j+1);
            }
            //bottom edge
            else if(j == 0){
                new_pixels[i*height + j] = 2*Pixel(i-1,j) + Pixel(i-1,j+1) + -2*Pixel(i+1,j) + -1*Pixel(i+1,j+1);
            }
            //top edge
            else if (j == height-1){
                new_pixels[i*height + j] = 2*Pixel(i-1,j) + Pixel(i-1,j-1) + -2*Pixel(i+1,j) + -1*Pixel(i+1,j-1);
            }
            //Everything else
            else {
                new_pixels[i*height + j] = Pixel(i-1,j-1) + 2*Pixel(i-1,j) + Pixel(i-1,j+1) +
                    -1*Pixel(i+1,j-1) + -2*Pixel(i+1,j) + -1*Pixel(i+1,j+1);
            }
           // new_pixels[i*height + j].Clamp();
        }
    }

    pixels = new_pixels;
}


void R2Image::
SobelY(void)
{
    // Apply the Sobel oprator to the image in Y direction
    R2Pixel *new_pixels;
    new_pixels = new R2Pixel [ npixels ];
    assert(new_pixels);
    
    for (int i = 0; i < width; i++){
        for (int j = 0; j < height; j++){
            //TODO: Make this prettier
            //bottom left
            if (i == 0 && j == 0) {
                new_pixels[i*height + j] = -2*Pixel(i,j+1) + -1*Pixel(i+1,j+1);
            }
            //bottom right
            else if (i == width-1 && j == 0){
                new_pixels[i*height + j] = -2*Pixel(i,j+1) + -1*Pixel(i-1,j+1);
            }
            //top left
            else if (i == 0 && j == height-1){
                new_pixels[i*height + j] = 2*Pixel(i,j-1) + 1*Pixel(i+1,j-1);
            }
            //top right
            else if (i == width-1 && j == height-1 ){
                new_pixels[i*height + j] = 2*Pixel(i,j-1) + 1*Pixel(i-1,j-1);
            }
            //left edge
            else if (i == 0) {
                new_pixels[i*height + j] = -2*Pixel(i,j+1) + -1*Pixel(i+1,j+1) + -2*Pixel(i,j-1) + -1*Pixel(i+1,j-1);
            }
            //right edge
            else if (i == width-1) {
                new_pixels[i*height + j] = -2*Pixel(i,j+1) + -1*Pixel(i-1,j+1) + -2*Pixel(i,j-1) + -1*Pixel(i-1,j-1);
            }
            //bottom edge
            else if(j == 0){
                new_pixels[i*height + j] = -1*Pixel(i-1,j+1) + -2*Pixel(i,j+1) + -1*Pixel(i+1,j+1);
            }
            //top edge
            else if (j == height-1){
                new_pixels[i*height + j] = 1*Pixel(i-1,j-1) + 2*Pixel(i,j-1) + 1*Pixel(i,j-1);
            }
            //Everything else
            else {
                new_pixels[i*height + j] = 1*Pixel(i-1,j-1) + 2*Pixel(i,j-1) + Pixel(i+1,j-1) +
                -1*Pixel(i-1,j+1) + -2*Pixel(i,j+1) + -1*Pixel(i+1,j+1);
            }
           // new_pixels[i*height + j].Clamp();
        }
    }
    
    pixels = new_pixels;
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image
  
    double myPoints[] = {0,0,1,0,1,1,0,1,1,0,2,0,2,1,1,1};
    //double myPoints[] = {0,0,1,0,1,1,0,1,1,2,1,1,3,1,3,2};
    //double myPoints[] = {318,367,907,302,49,326,74,436,299,379,887,306,30,341,56,451};
    std::vector<double> v (myPoints, myPoints + sizeof(myPoints) / sizeof(double));
    std::vector<double> HMatrix = dlt(v);

    double x = 1;
    double y = 0;

    double transformed_x = (HMatrix[0] * x) + (HMatrix[1] * y) + (HMatrix[2] * 1);
    double transformed_y = (HMatrix[3] * x) + (HMatrix[4] * y) + (HMatrix[5] * 1);
    double transformed_z = (HMatrix[6] * x) + (HMatrix[7] * y) + (HMatrix[8] * 1);

    fprintf(stderr, "Hmatrix things %f %f %f \n", HMatrix[0], HMatrix[1], HMatrix[2]);
    fprintf(stderr, "Hmatrix things %f %f %f \n", HMatrix[3], HMatrix[4], HMatrix[5]);

    fprintf(stderr, "Hmatrix things %f %f %f \n", HMatrix[6], HMatrix[7], HMatrix[8]);

    fprintf(stderr, "transfered x,y,z = %f %f %f \n", transformed_x, transformed_y, transformed_z);

    double normalized_x = transformed_x/transformed_z;
    double normalized_y = transformed_y/transformed_z;

    fprintf(stderr, "normalized x,y = %f %f \n", normalized_x, normalized_y);


    // double difference_in_x = normalized_x - translated_features[i].x;
    // double difference_in_y = normalized_y - translated_features[i].y;


  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}

void R2Image::
ChangeSaturation(double factor)
{
    // Changes the saturation of an image
    // Find a formula that changes the saturation without affecting the image brightness
    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            double red = Pixel(i,j).Red();
            double green = Pixel(i,j).Green();
            double blue = Pixel(i,j).Blue();
            
            double gray = 0.2989*red + 0.5870*green + 0.1140*blue;
            
            double new_red = gray + (red - gray)*factor;
            double new_green = gray + (green - gray)*factor;
            double new_blue = gray + (blue - gray)*factor;
            
            Pixel(i,j).Reset(new_red,new_green,new_blue,1);
            Pixel(i,j).Clamp();
        }
    }
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
    // Gaussian blur of the image. Separable solution is preferred
    
    R2Image myTempImage(*this);
    
    int size = 6*sigma + 1;
    int midpoint = size/2;
    
    //create the kernel
    double *kernel;
    kernel = (double*)malloc(size*sizeof(double));
    
    //calculate gaussian filter
    double sum = 0;
    for (int i = 0; i < size; i++){
        int distance = abs(i - midpoint) ;
        
        double weight = ( 1 / sigma * (sqrt(2*M_PI) ) ) * exp( -0.5 * pow( distance / sigma , 2.0 ) );
        kernel[i] = weight;
        sum += weight;
    }
    
    //normalize
    for(int i = 0; i< size; i++){
        kernel[i] = kernel[i]/sum;
    }

    //first pass in x
    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            R2Pixel val;
            for(int k = 0; k < size; k++)
            {
                int xCord = i-midpoint+k;
                if (xCord >= 0 && xCord < width)
                {
                    val += kernel[k]*Pixel(xCord,j);
                }
            }
            myTempImage.Pixel(i,j) = val;
        }
    }
    (*this) = myTempImage;

    //second pass in y
    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            R2Pixel val;
            for(int k = 0; k < size; k++)
            {
                int yCord = j-midpoint+k;
                if (yCord >= 0 && yCord < height)
                {
                    val += kernel[k]*Pixel(i,yCord);
                }
                myTempImage.Pixel(i,j) = val;
            }
        }
    }
    (*this) = myTempImage;
}

struct S
{
    int x;
    int y;
    double brightness;
    
    S(int i, int j, R2Pixel pixel)
    {
        x = i;
        y = j;
        brightness = pixel.Red();
    }
    
    bool operator<(const struct S& other) const
    {
        //Your priority logic goes here
        return brightness < other.brightness;
    }
};


void R2Image::
Harris(double sigma)
{
    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
    // Output should be 50% grey at flat regions, white at corners and black/dark near edges
  
    R2Image Ix(*this);
    R2Image Iy(*this);
    R2Image Ixy(*this);
    R2Image final_img(*this);
    R2Image original_copy(*this);
    
    //TODO: Unclamp these to make them work better
    Ix.ChangeSaturation(0.0);
    Iy.ChangeSaturation(0.0);
    Ixy.ChangeSaturation(0.0);
    
    Ix.SobelX();
    Iy.SobelY();
    Ixy.SobelX();
    
    for(int i = 0; i < width ; i++){
        for(int j = 0; j < height; j++){
            //set Ixy
            Ixy.Pixel(i,j) = Ixy.Pixel(i,j) * Iy.Pixel(i,j);
            
            //set Ix
            Ix.Pixel(i,j) = Ix.Pixel(i,j) * Ix.Pixel(i,j);

            //set Iy
            Iy.Pixel(i,j) = Iy.Pixel(i,j) * Iy.Pixel(i,j);
        }
    }
    
    Ix.Blur(sigma);
    Iy.Blur(sigma);
    Ixy.Blur(sigma);
    
    std::priority_queue<S> pq;
    
    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            final_img.Pixel(i,j) = (Ix.Pixel(i,j) * Iy.Pixel(i,j)) - (Ixy.Pixel(i,j) * Ixy.Pixel(i,j)) -
                0.04 * ((Ix.Pixel(i,j) + Iy.Pixel(i,j)) * (Ix.Pixel(i,j) + Iy.Pixel(i,j)));

            final_img.Pixel(i,j).Reset(final_img.Pixel(i,j).Red() + 0.5, final_img.Pixel(i,j).Green() + 0.5, final_img.Pixel(i,j).Blue() + 0.5, 1);
            final_img.Pixel(i,j).Clamp();
            pq.push(S(i,j,final_img.Pixel(i,j)));
        }
    }
    
    short used[width][height];

    //distance to next closest feature
    int bounds = 20;
    //the feature will be a marked by a bsize x bsize box centered on the feature
    int bsize = 10;
    //to find out how many features we are marking
    int num_features = 0;
    
    while(num_features < 150){
        S next_feature = pq.top();
        pq.pop();
        int i = next_feature.x;
        int j = next_feature.y;
        
        if(used[i][j] == 1){
            continue;
        }
        
        int min_x = fmax(i-bounds,0);
        int max_x = fmin(i+bounds,width-1);
        int min_y = fmax(j-bounds,0);
        int max_y = fmin(j+bounds,height-1);
        
        //negate all of short in the range
        for(int a = min_x; a <= max_x; a++){
            for(int b = min_y; b <= max_y; b++){
                used[a][b] = 1;
            }
        }
        
        //find the bounds of the box
        min_x = fmax(i-bsize/2,0);
        max_x = fmin(i+bsize/2,width-1);
        min_y = fmax(j-bsize/2,0);
        max_y = fmin(j+bsize/2,height-1);
        
        //draw a box around it
        for(int a = min_x; a <= max_x; a++){
            for(int b = min_y; b <= max_y; b++){
                if((a == min_x || a == max_x) || (b == min_y || b == max_y))
                {
                    original_copy.Pixel(a,b).Reset(1.0,0,0,1.0);
                }
            }
        }
        num_features += 1;
    }

    fprintf(stderr, "Num Features: %d \n", num_features);
    
    (*this) = original_copy;
}

void R2Image::
Sharpen()
{
    // Sharpen an image using a linear filter. Use a kernel of your choosing.

    R2Pixel *new_pixels;
    new_pixels = new R2Pixel [ npixels ];
    assert(new_pixels);
    
    for (int i = 1; i < width-1; i++){
        for (int j = 1; j < height-1; j++){
            
            new_pixels[i*height + j] =
                0*Pixel(i-1,j+1) + -1*Pixel(i,j+1) + 0*Pixel(i+1,j+1) +
                -1*Pixel(i-1,j) + 5*Pixel(i,j) + -1*Pixel(i+1,j) +
                0*Pixel(i-1,j-1) + -1*Pixel(i,j-1) + 0*Pixel(i+1,j-1);
            
            new_pixels[i*height + j].Clamp();
        }
    }
    
    pixels = new_pixels;
}

void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
    // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
    // compute the matching translation (pixel precision is OK), and blend the translated "otherImage" 
    // into this image with a 50% opacity.
    double sigma = 2.0;
    
    R2Image Ix(*this);
    R2Image Iy(*this);
    R2Image Ixy(*this);
    R2Image final_img(*this);
    R2Image original_copy(*this);
    R2Image otherImageCopy(*otherImage);

    //TODO: Unclamp these to make them work better
    Ix.ChangeSaturation(0.0);
    Iy.ChangeSaturation(0.0);
    Ixy.ChangeSaturation(0.0);
    
    Ix.SobelX();
    Iy.SobelY();
    Ixy.SobelX();
    
    for(int i = 0; i < width ; i++){
        for(int j = 0; j < height; j++){
            //set Ixy
            Ixy.Pixel(i,j) = Ixy.Pixel(i,j) * Iy.Pixel(i,j);
            
            //set Ix
            Ix.Pixel(i,j) = Ix.Pixel(i,j) * Ix.Pixel(i,j);
            
            //set Iy
            Iy.Pixel(i,j) = Iy.Pixel(i,j) * Iy.Pixel(i,j);
        }
    }
    
    Ix.Blur(sigma);
    Iy.Blur(sigma);
    Ixy.Blur(sigma);
    
    std::priority_queue<S> pq;
    //distance to next closest feature
    int bounds = 20;
    //the feature will be a marked by a bsize x bsize box centered on the feature
    int bsize = 20;
    int bsizeRadius = bsize/2;
    //to find out how many features we are marking
    int num_features = 0;
    
    
    for(int i = bsize/2; i < width-bsize/2; i++){
        for(int j = bsize/2; j < height-bsize/2; j++){
            final_img.Pixel(i,j) = (Ix.Pixel(i,j) * Iy.Pixel(i,j)) - (Ixy.Pixel(i,j) * Ixy.Pixel(i,j)) -
            0.04 * ((Ix.Pixel(i,j) + Iy.Pixel(i,j)) * (Ix.Pixel(i,j) + Iy.Pixel(i,j)));
            
            final_img.Pixel(i,j).Reset(final_img.Pixel(i,j).Red() + 0.5, final_img.Pixel(i,j).Green() + 0.5, final_img.Pixel(i,j).Blue() + 0.5, 1);
            final_img.Pixel(i,j).Clamp();
            pq.push(S(i,j,final_img.Pixel(i,j)));
        }
    }
    
    short used[width][height];
    std::vector<S> features;
    
    //find 150 features, and store their locations.
    while(num_features < 150){
        S next_feature = pq.top();
        pq.pop();
        int i = next_feature.x;
        int j = next_feature.y;
        
        if(used[i][j] == 1){
            continue;
        }
        
        int min_x = fmax(i-bounds,0);
        int max_x = fmin(i+bounds,width-1);
        int min_y = fmax(j-bounds,0);
        int max_y = fmin(j+bounds,height-1);
        
        //negate all of short in the range
        for(int a = min_x; a <= max_x; a++){
            for(int b = min_y; b <= max_y; b++){
                used[a][b] = 1;
            }
        }
        
        //add it to the list of features to be compared
        features.push_back(next_feature);
        num_features += 1;
    }
    
    fprintf(stderr,"found 150 features\n");
    
    std::vector<S> translated_features;
    double searchSize = 0.15;
    
    //for every feature centered at pixel x,y - look for the smallest SSD
    for(int f = 0; f < (int) features.size(); f++){
        int x = features[f].x;
        int y = features[f].y;
       
        //look in a box around the pixel for the best matching bsize x bsize square.
        int search_x_max = fmin(x+(searchSize * width),width-1);
        int search_x_min = fmax(x-(searchSize * width),0);
        int search_y_max = fmin(y+(searchSize * height),height-1);
        int search_y_min = fmax(y-(searchSize * height),0);
        
        //search the box for every bsize x bsize square, and compare it to the original.
        double min_ssd = std::numeric_limits<int>::max();
        int min_ssd_x = x;
        int min_ssd_y = y;
        for(int i = search_x_min+bsizeRadius; i < search_x_max-bsizeRadius; i++){
            for(int j = search_y_min+bsizeRadius; j < search_y_max-bsizeRadius; j++){
                double SSD = 0;

                for(int a = -bsizeRadius; a < bsizeRadius; a++){
                    for(int b = -bsizeRadius; b < bsizeRadius;b++){
                        R2Pixel original_pixel = Pixel(x+a,y+b);
                        R2Pixel new_pixel = otherImageCopy.Pixel(i+a,j+b);
                        // double red_difference = original_pixel.Red()  - new_pixel.Red();
                        // double green_difference = original_pixel.Green() - new_pixel.Green();
                        // double blue_difference = original_pixel.Blue() - new_pixel.Blue();
                        // SSD += pow(red_difference,2) + pow(green_difference,2) + pow(blue_difference,2);
                        SSD += ((original_pixel.Red()  - new_pixel.Red()) * (original_pixel.Red()  - new_pixel.Red())) + 
                                ((original_pixel.Green() - new_pixel.Green()) * (original_pixel.Green() - new_pixel.Green())) +
                                ((original_pixel.Blue() - new_pixel.Blue()) * (original_pixel.Blue() - new_pixel.Blue()));
                    }
                }

                if(SSD < min_ssd){
                    min_ssd = SSD;
                    min_ssd_x = i;
                    min_ssd_y = j;
                }
            }
        }
        translated_features.push_back(S(min_ssd_x,min_ssd_y,otherImageCopy.Pixel(min_ssd_x,min_ssd_y)));
    }
        
    std::set<int> inliers;
    double threshold = 6;
    int N = 1000;

    for(int t = 0; t < N; t++){
        int track_1 = rand() % 150;
        int track_2 = rand() % 150;
        while(track_2 == track_1){
            track_2 = rand() % 150; 
        }
        int track_3 = rand() % 150;
        while(track_3 == track_1 || track_3 == track_2){
            track_3 = rand() % 150; 
        }
        int track_4 = rand() % 150;
        while(track_4 == track_1 || track_4 == track_2 || track_4 == track_3){
            track_4 = rand() % 150; 
        }

        std::vector<int> tracks;
        tracks.push_back(track_1);
        tracks.push_back(track_2);
        tracks.push_back(track_3);
        tracks.push_back(track_4);

        std::vector<double> HMatrix;
        std::vector<double> inputs;
        //add the original feature positions
        for(int i = 0;i<4;i++){
            inputs.push_back(features[tracks[i]].x);
            inputs.push_back(features[tracks[i]].y);
        }
        //add the new feature position
        for(int i = 0;i<4;i++){
            inputs.push_back(translated_features[tracks[i]].x);
            inputs.push_back(translated_features[tracks[i]].y);
        }

        HMatrix = dlt(inputs);

        std::set<int> temp_inliers;
        for(int i = 0; i < (int) features.size(); i++) {
            double transformed_x = (HMatrix[0] * features[i].x) + (HMatrix[1] * features[i].y) + (HMatrix[2] * 1);
            double transformed_y = (HMatrix[3] * features[i].x) + (HMatrix[4] * features[i].y) + (HMatrix[5] * 1);
            double transformed_z = (HMatrix[6] * features[i].x) + (HMatrix[7] * features[i].y) + (HMatrix[8] * 1);

            double normalized_x = transformed_x/transformed_z;
            double normalized_y = transformed_y/transformed_z;

            double difference_in_x = normalized_x - translated_features[i].x;
            double difference_in_y = normalized_y - translated_features[i].y;


            double distance = sqrt(pow(difference_in_x,2) + pow(difference_in_y,2));

            if(distance < threshold){
                temp_inliers.insert(i);
            }
        }

        if(temp_inliers.size() > inliers.size()){
            inliers = temp_inliers;
            fprintf(stderr,"inliers %d \n", (int)inliers.size());
        }
    }

    printf("\n num inliers final -  %d \n", (int) inliers.size());

    std::vector<double> inputs;
    std::vector<double> HMatrix;
    for(int i = 0; i < (int) features.size(); i++){
        if(inliers.find(i) != inliers.end()){
            inputs.push_back(features[i].x);
            inputs.push_back(features[i].y);
        }
    }
    for(int i = 0; i < (int) features.size(); i++){
        if(inliers.find(i) != inliers.end()){
            inputs.push_back(translated_features[i].x);
            inputs.push_back(translated_features[i].y);
        }
    }

    fprintf(stderr,"\n In Hmatrx - inputs size: %d, inliers %d \n\n\n", (int) inputs.size(), (int) inliers.size());

    HMatrix = dlt(inputs);

    printf("here finished DLT matrix\n\n");

    printf("width %d height %d\n\n", width, height);
    for(int x = 0; x < width; x++){
        for(int y = 0; y < height; y++){
            double transformed_x = (HMatrix[0] * x) + (HMatrix[1] * y) + (HMatrix[2] * 1);
            double transformed_y = (HMatrix[3] * x) + (HMatrix[4] * y) + (HMatrix[5] * 1);
            double transformed_z = (HMatrix[6] * x) + (HMatrix[7] * y) + (HMatrix[8] * 1);

            double normalized_x = transformed_x/transformed_z;
            double normalized_y = transformed_y/transformed_z;

            if(normalized_x >= width || normalized_x < 0){
                continue;
            }
            if(normalized_y >= height || normalized_y < 0){
                continue;
            }

            // otherImageCopy.Pixel(normalized_x,normalized_y) = otherImageCopy.Pixel(normalized_x,normalized_y)*.5 + Pixel(x,y)*.5;
            otherImageCopy.Pixel(x,y) = Pixel(x,y)*.5 + otherImage->Pixel(normalized_x,normalized_y)*.5;
        }
    }

    printf("finished everything \n");

    (*this) = otherImageCopy;
    
    fprintf(stderr, "fit other image using translation and blend imageB over imageA\n");

    return;
}


void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
    // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
    // compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.
    fprintf(stderr, "fit other image using a homography and blend imageB over imageA\n");
    return;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
    
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
    
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
    // Read pixel values
    int red, green, blue;
    if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
      fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
      fclose(fp);
      return 0;
    }

    // Assign pixel values
    double r = (double) red / max_value;
    double g = (double) green / max_value;
    double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


    

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width;    /* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;       /* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB;   /* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 100, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
    
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}






