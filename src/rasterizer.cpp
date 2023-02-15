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
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
      for (int i = 0; i < get_sample_rate(); i++){
          sample_buffer[(y * width + x) * get_sample_rate() + i] = c;
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

    fill_pixel(sx, sy, color);
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
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
      float minx = floor(min(x0, min(x1, x2)));
      float maxx = ceil(max(x0, max(x1, x2)));
      float miny = floor(min(y0, min(y1, y2)));
      float maxy = ceil(max(y0, max(y1, y2)));

      double offset = 1.0 / (2 * sqrt(get_sample_rate()));
      float tick = 2.0 * offset;

      for (double x = minx+offset; x < maxx; x+=tick) {
          for (double y = miny+offset; y < maxy; y+=tick) {

              float flag1 = (-(x - x0) * (y1 - y0) + (y - y0) * (x1 - x0));
              float flag2 = (-(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1));
              float flag3 = (-(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2));

              if (flag1 >= 0 && flag2 >= 0 && flag3 >= 0) {
                  if (x < 0 || x >= width) return;
                  if (y < 0 || y >= height) return;

                  float fbx = floor(x);
                  float fby = floor(y);
                  float offset = 1.0 / (2.0 * sqrt(get_sample_rate()));
                  float pixelwidth = 1.0 * floor(sqrt(get_sample_rate()));

                  int pixelx = floor((x - fbx - offset) * pixelwidth);
                  int pixely = floor((y - fby - offset) * pixelwidth);

                  sample_buffer[(fby * width + fbx) * get_sample_rate() + (pixely * pixelwidth) + pixelx] = color;
              }
              else if (flag1 <= 0 && flag2 <= 0 && flag3 <= 0) {
                  if (x < 0 || x >= width) return;
                  if (y < 0 || y >= height) return;

                  float fbx = floor(x);
                  float fby = floor(y);
                  float offset = 1.0 / (2.0 * sqrt(get_sample_rate()));
                  float pixelwidth = 1.0 * floor(sqrt(get_sample_rate()));

                  int pixelx = floor((x - fbx - offset) * pixelwidth);
                  int pixely = floor((y - fby - offset) * pixelwidth);

                  sample_buffer[(fby * width + fbx) * get_sample_rate() + (pixely * pixelwidth) + pixelx] = color;
              }
          }
      }

    // TODO: Task 2: Update to implement super-sampled rasterization

      

  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
      float minx = floor(min(x0, min(x1, x2)));
      float maxx = ceil(max(x0, max(x1, x2)));
      float miny = floor(min(y0, min(y1, y2)));
      float maxy = ceil(max(y0, max(y1, y2)));

      double offset = 1 / (2 * sqrt(get_sample_rate()));
      float tick = 2 * offset;

      for (double x = minx + offset; x < maxx; x += tick) {
          for (double y = miny + offset; y < maxy; y += tick) {

              Color total = Color(0, 0, 0);

              double w0 = ((-(x - x1) * (y2 - y1)) + ((y - y1) * (x2 - x1))) / ((-(x0 - x1) * (y2 - y1)) + ((y0 - y1) * (x2 - x1)));
              double w1 = ((-(x - x2) * (y0 - y2)) + ((y - y2) * (x0 - x2))) / ((-(x1 - x2) * (y0 - y2)) + ((y1 - y2) * (x0 - x2)));
              double w2 = 1.0 - w0 - w1;

              total = (w0 * c0) + (w1 * c1) + (w2 * c2);


              float flag1 = (-(x - x0) * (y1 - y0) + (y - y0) * (x1 - x0));
              float flag2 = (-(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1));
              float flag3 = (-(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2));

              if (flag1 >= 0 && flag2 >= 0 && flag3 >= 0) {
                  if (x < 0 || x >= width) return;
                  if (y < 0 || y >= height) return;

                  float fbx = floor(x);
                  float fby = floor(y);
                  float offset = 1.0 / (2.0 * sqrt(get_sample_rate()));
                  double pixelwidth = 1.0 * floor(sqrt(get_sample_rate()));

                  int pixelx = floor((x - fbx - offset) * pixelwidth);
                  int pixely = floor((y - fby - offset) * pixelwidth);

                  sample_buffer[(fby * width + fbx) * get_sample_rate() + (pixely * pixelwidth) + pixelx] = total;
              }
              else if (flag1 <= 0 && flag2 <= 0 && flag3 <= 0) {
                  if (x < 0 || x >= width) return;
                  if (y < 0 || y >= height) return;

                  float fbx = floor(x);
                  float fby = floor(y);
                  float offset = 1.0 / (2.0 * sqrt(get_sample_rate()));
                  float pixelwidth = 1.0 * floor(sqrt(get_sample_rate()));

                  int pixelx = floor((x - fbx - offset) * pixelwidth);
                  int pixely = floor((y - fby - offset) * pixelwidth);

                  sample_buffer[(fby * width + fbx) * get_sample_rate() + (pixely * pixelwidth) + pixelx] = total;
              }
          }
      }


  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
      float minx = floor(min(x0, min(x1, x2)));
      float maxx = ceil(max(x0, max(x1, x2)));
      float miny = floor(min(y0, min(y1, y2)));
      float maxy = ceil(max(y0, max(y1, y2)));

      double offset = 1 / (2 * sqrt(get_sample_rate()));
      float tick = 2 * offset;
      
      int level = 0;
      SampleParams sp = SampleParams();
      sp.lsm = lsm;
      sp.psm = psm;

      for (double x = minx + offset; x < maxx; x += tick) {
          for (double y = miny + offset; y < maxy; y += tick) {

              double w0 = ((-(x - x1) * (y2 - y1)) + ((y - y1) * (x2 - x1))) / ((-(x0 - x1) * (y2 - y1)) + ((y0 - y1) * (x2 - x1)));
              double w1 = ((-(x - x2) * (y0 - y2)) + ((y - y2) * (x0 - x2))) / ((-(x1 - x2) * (y0 - y2)) + ((y1 - y2) * (x0 - x2)));
              double w2 = 1.0 - w0 - w1;

              float u = (w0 * u0) + (w1 * u1) + (w2 * u2);
              float v = (w0 * v0) + (w1 * v1) + (w2 * v2);

              sp.p_uv = Vector2D(u, v);

              double dxw0 = ((-((x+1) - x1) * (y2 - y1)) + ((y - y1) * (x2 - x1))) / ((-(x0 - x1) * (y2 - y1)) + ((y0 - y1) * (x2 - x1)));
              double dxw1 = ((-(x + 1 - x2) * (y0 - y2)) + ((y - y2) * (x0 - x2))) / ((-(x1 - x2) * (y0 - y2)) + ((y1 - y2) * (x0 - x2)));
              double dxw2 = 1.0 - dxw0 - dxw1;

              float dxu = (dxw0 * u0) + (dxw1 * u1) + (dxw2 * u2);
              float dxv = (dxw0 * v0) + (dxw1 * v1) + (dxw2 * v2);

              sp.p_dx_uv = Vector2D(dxu, dxv);

              double dyw0 = ((-(x - x1) * (y2 - y1)) + ((y + 1 - y1) * (x2 - x1))) / ((-(x0 - x1) * (y2 - y1)) + ((y0 - y1) * (x2 - x1)));
              double dyw1 = ((-(x - x2) * (y0 - y2)) + ((y + 1 - y2) * (x0 - x2))) / ((-(x1 - x2) * (y0 - y2)) + ((y1 - y2) * (x0 - x2)));
              double dyw2 = 1.0 - w0 - w1;

              float dyu = (w0 * u0) + (w1 * u1) + (w2 * u2);
              float dyv = (w0 * v0) + (w1 * v1) + (w2 * v2);

              sp.p_dy_uv = Vector2D(dyu, dyv);

              Color c = Color(0, 0, 0);

              float flag1 = (-(x - x0) * (y1 - y0) + (y - y0) * (x1 - x0));
              float flag2 = (-(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1));
              float flag3 = (-(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2));

              if (flag1 >= 0 && flag2 >= 0 && flag3 >= 0) {

                  c = tex.sample(sp);

                  if (x < 0 || x >= width) return;
                  if (y < 0 || y >= height) return;

                  float fbx = floor(x);
                  float fby = floor(y);
                  float offset = 1.0 / (2.0 * sqrt(get_sample_rate()));
                  float pixelwidth = 1.0 * floor(sqrt(get_sample_rate()));

                  int pixelx = floor((x - fbx - offset) * pixelwidth);
                  int pixely = floor((y - fby - offset) * pixelwidth);

                  sample_buffer[(fby * width + fbx) * get_sample_rate() + (pixely * pixelwidth) + pixelx] = c;
              }
              else if (flag1 <= 0 && flag2 <= 0 && flag3 <= 0) {
                  /*if (psm == P_NEAREST) {
                      c = tex.sample_nearest(Vector2D(u, v), 0);
                  }
                  else if (psm == P_LINEAR) {
                      c = tex.sample_bilinear(Vector2D(u, v), 0);
                  }*/

                  c = tex.sample(sp);

                  if (x < 0 || x >= width) return;
                  if (y < 0 || y >= height) return;

                  float fbx = floor(x);
                  float fby = floor(y);
                  float offset = 1.0 / (2.0 * sqrt(get_sample_rate()));
                  float pixelwidth = 1.0 * floor(sqrt(get_sample_rate()));

                  int pixelx = floor((x - fbx - offset) * pixelwidth);
                  int pixely = floor((y - fby - offset) * pixelwidth);

                  sample_buffer[(fby * width + fbx) * get_sample_rate() + (pixely * pixelwidth) + pixelx] = c;
              }
          }
      }

    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height * rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
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

              int pixelbegin = (y * width + x) * get_sample_rate();

              Color total = Color::Black;

              for (int i = pixelbegin; i < pixelbegin + get_sample_rate(); i++) {
                  total += sample_buffer[i];
              }

              Color average = total *(1.0 / (get_sample_rate()));
              Color col = average;// sample_buffer[(y * width + x) * get_sample_rate()];

              for (int k = 0; k < 3; ++k) {
                  this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
              }
          }
      }
  }

  Rasterizer::~Rasterizer() { }


}// CGL
