#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c, int sample_no) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
      if (sample_no < 0) {
          for (int i = 0; i < sample_rate; i++) {
              sample_buffer[y*(sample_rate) * width + x*(sample_rate) + i] = c;
          }
      }
      else {
          sample_buffer[y*(sample_rate) * width + x*(sample_rate) + sample_no] = c;
      }
    
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color, -1);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
      
      rasterize_line(x0, y0, x1, y1, color);
      rasterize_line(x1, y1, x2, y2, color);
      rasterize_line(x2, y2, x0, y0, color);
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    float left_bound = std::min({ x0, x1, x2 });
    float right_bound = std::max({ x0, x1, x2 });
    float t_bound = std::min({ y0, y1, y2 });
    float b_bound = std::max({ y0, y1, y2 });
    for (float i = t_bound; i < b_bound; i += 1) {
        for (float j = left_bound; j < right_bound; j += 1) {
            for (int x_sub = 0; x_sub < sqrt(sample_rate); x_sub ++) {
                for (int y_sub = 0; y_sub < sqrt(sample_rate); y_sub ++) {
                    float x = j + .5*(1/sqrt(this->sample_rate)) + x_sub*1/sqrt(this->sample_rate);
                    float y = i + .5*(1/sqrt(this->sample_rate)) + y_sub*1/sqrt(this->sample_rate);
                    if (inside(x, y, x0, y0, x1, y1, x2, y2)) {
                        fill_pixel(j, i, color, x_sub+y_sub*sqrt(sample_rate));
                    }
                }
            }
            
        }
    }
        // TODO: Task 2: Update to implement super-sampled rasterization

  }

  bool RasterizerImp::inside(float x, float y, float x0, float y0, float x1, float y1, float x2, float y2) {

      float first = l_func(x, y, x0, y0, x1, y1);
      float second = l_func(x, y, x1, y1, x2, y2);
      float third = l_func(x, y, x2, y2, x0, y0);
      return (first >= 0 && second >= 0 && third >= 0) || (first <= 0 && second <= 0 && third <= 0);

  }

  float RasterizerImp::l_func(float sample_x, float sample_y, float line_x0, float line_y0, float line_x1, float line_y1) {
      float dx = (line_x1 - line_x0);
      float dy = (line_y1 - line_y0);
      return -1 * (sample_x - line_x0) * dy + (sample_y - line_y0) * dx;
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle



  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support
      
    this->sample_rate = rate;

    this->sample_buffer.resize(width * height * rate, Color::White);
    clear_buffers();
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

      

    this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
      clear_buffers();
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
          Color my_col = Color(0,0,0);
          for (int sub = 0; sub < sample_rate; sub ++) {
                  //temp_col = sample_buffer[(int)((y+y_sub)*width*sample_rate + (x+x_sub)*sqrt(sample_rate))];
                  Color temp_col = sample_buffer[(int)(y*width*(sample_rate) + x*(sample_rate) + sub)];
                  my_col.r += temp_col.r * 1/sample_rate;
                  my_col.g += temp_col.g * 1/sample_rate;
                  my_col.b += temp_col.b * 1/sample_rate;
          }
        //Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&my_col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
