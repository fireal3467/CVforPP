// Include file for image class
#ifndef R2_IMAGE_INCLUDED
#define R2_IMAGE_INCLUDED



// Constant definitions

typedef enum {
  R2_IMAGE_RED_CHANNEL,
  R2_IMAGE_GREEN_CHANNEL,
  R2_IMAGE_BLUE_CHANNEL,
  R2_IMAGE_ALPHA_CHANNEL,
  R2_IMAGE_NUM_CHANNELS
} R2ImageChannel;

typedef enum {
  R2_IMAGE_POINT_SAMPLING,
  R2_IMAGE_BILINEAR_SAMPLING,
  R2_IMAGE_GAUSSIAN_SAMPLING,
  R2_IMAGE_NUM_SAMPLING_METHODS
} R2ImageSamplingMethod;

typedef enum {
  R2_IMAGE_OVER_COMPOSITION,
  R2_IMAGE_IN_COMPOSITION,
  R2_IMAGE_OUT_COMPOSITION,
  R2_IMAGE_ATOP_COMPOSITION,
  R2_IMAGE_XOR_COMPOSITION,
} R2ImageCompositeOperation;



// Class definition

class R2Image {
 public:
  // Constructors/destructor
  R2Image(void);
  R2Image(const char *filename);
  R2Image(int width, int height);
  R2Image(int width, int height, const R2Pixel *pixels);
  R2Image(const R2Image& image);
  ~R2Image(void);

  // Image properties
  int NPixels(void) const;
  int Width(void) const;
  int Height(void) const;

  // Pixel access/update
  R2Pixel& Pixel(int x, int y);
  R2Pixel *Pixels(void);
  R2Pixel *Pixels(int row);
  R2Pixel *operator[](int row);
  const R2Pixel *operator[](int row) const;
  void SetPixel(int x, int y,  const R2Pixel& pixel);

  // Image processing
  R2Image& operator=(const R2Image& image);

  // Per-pixel operations
  void Brighten(double factor);
  void ChangeSaturation(double factor);

  // show how SVD works
  void svdTest();

  // Linear filtering operations
  void SobelX();
  void SobelY();
  void LoG();
  void Blur(double sigma);
  void Harris(double sigma);
  void Sharpen(void);

  // further operations
  void blendOtherImageTranslated(R2Image * otherImage);
  void blendOtherImageHomography(R2Image * otherImage);

  // File reading/writing
  int Read(const char *filename);
  int ReadBMP(const char *filename);
  int ReadPPM(const char *filename);
  int ReadJPEG(const char *filename);
  int Write(const char *filename) const;
  int WriteBMP(const char *filename) const;
  int WritePPM(const char *filename, int ascii = 0) const;
  int WriteJPEG(const char *filename) const;


  //added functions by Alan 
  void line(int x0, int x1, int y0, int y1, float r, float g, float b);

 private:
  // Utility functions
  void Resize(int width, int height);
  R2Pixel Sample(double u, double v,  int sampling_method);

 private:
  R2Pixel *pixels;
  int npixels;
  int width;
  int height;
};

//added functions
void R2Image::line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
    if(x0>x1)
    {
        int x=y1;
        y1=y0;
        y0=x;
        
        x=x1;
        x1=x0;
        x0=x;
    }
    int deltax = x1 - x0;
    int deltay = y1 - y0;
    float error = 0;
    float deltaerr = 0.0;
    if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
    // note that this division needs to be done in a way that preserves the fractional part
    int y = y0;
    for(int x=x0;x<=x1;x++)
    {
        Pixel(x,y).Reset(r,g,b,1.0);
        error = error + deltaerr;
        if(error>=0.5)
        {
            if(deltay>0) y = y + 1;
            else y = y - 1;
            
            error = error - 1.0;
        }
    }
    if(x0>3 && x0<width-3 && y0>3 && y0<height-3)
    {
        for(int x=x0-3;x<=x0+3;x++)
        {
            for(int y=y0-3;y<=y0+3;y++)
            {
                Pixel(x,y).Reset(r,g,b,1.0);
            }
        }
    }
}



// Inline functions

inline int R2Image::
NPixels(void) const
{
  // Return total number of pixels
  return npixels;
}



inline int R2Image::
Width(void) const
{
  // Return width
  return width;
}



inline int R2Image::
Height(void) const
{
  // Return height
  return height;
}



inline R2Pixel& R2Image::
Pixel(int x, int y)
{
  // Return pixel value at (x,y)
  // (pixels start at lower-left and go in row-major order)
  return pixels[x*height + y];
}



inline R2Pixel *R2Image::
Pixels(void)
{
  // Return pointer to pixels for whole image 
  // (pixels start at lower-left and go in row-major order)
  return pixels;
}



inline R2Pixel *R2Image::
Pixels(int x)
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}



inline R2Pixel *R2Image::
operator[](int x) 
{
  // Return pixels pointer for row at x
  return Pixels(x);
}



inline const R2Pixel *R2Image::
operator[](int x) const
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}



inline void R2Image::
SetPixel(int x, int y, const R2Pixel& pixel)
{
  // Set pixel
  pixels[x*height + y] = pixel;
}



#endif
