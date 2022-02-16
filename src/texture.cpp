#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
   
      //return sample_nearest(sp.p_uv, 0);
      //return sample_bilinear(sp.p_uv, 0);
      if (sp.psm == 0) {
          if (sp.lsm == 0) {
              return sample_nearest(sp.p_uv, 0);
          }
          else if (sp.lsm == 1) {
              // nearest appropriate level and pass that to nearest
              //std::cout << "NEAREST_SAMPLE, LEVEL CHANGED TO 1";
              //return sample_nearest(sp.p_uv, 0);
              return sample_nearest(sp.p_uv, round(get_level(sp)));
            //  return sample_nearest(sp.p_uv, 14);

          }
          else {
              // continuous number?? then weighted sum
              return sample_nearest(sp.p_uv, round(get_level(sp)));
          }
      }
      else {
          if (sp.lsm == 0) {
              return sample_bilinear(sp.p_uv, 0);
          }
          else if (sp.lsm == 1) {
              // nearest appropriate level and pass that to bilinear
              return sample_bilinear(sp.p_uv, round(get_level(sp)));
              //return sample_bilinear(sp.p_uv, get_level(sp));
          }
          else {
              // continuous number?? then weighted sum
              return sample_bilinear(sp.p_uv, 0);
          }
      }
// return magenta for invalid level
  //  return Color(1, 0, 1);
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
      float du_dx = (sp.p_dx_uv.x - sp.p_uv.x)*width;
      float dv_dx = (sp.p_dx_uv.y - sp.p_uv.y)*height;
      float du_dy = (sp.p_dy_uv.x - sp.p_uv.x)*width;
      float dv_dy = (sp.p_dy_uv.y - sp.p_uv.y)*height;
      float dx = sqrt(pow(du_dx, 2) + pow(dv_dx, 2));
      float dy = sqrt(pow(du_dy, 2) + pow(dv_dy, 2));
      float l = std::max(dx, dy);
      int d = round(log2(l));
      //std::cout << " d" << d;
      if (d >= kMaxMipLevels) {
          return kMaxMipLevels -1;
      }
      else if (d < 0) {
          return 0;
      } 
      return d;
  }

  Color MipLevel::get_texel(int tx, int ty, MipLevel& mip) {
      //return Color(1, 0, 1);
     // std::cout << "TX" << tx << " ";
      
      //std::cout << "TY" << ty << " ";
      if (tx < 0) {
          tx = 0;
      }
      if (ty < 0) {
          ty = 0;
      }
      if (tx >= mip.width) {
          tx = mip.width - 1;
      }
      if (ty >= mip.height) {
          ty = mip.height - 1;
      }
      return Color(&mip.texels[tx * 3 + ty * mip.width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];
    return mip.get_texel(round(uv.x * (mip.width -1)), round(uv.y * (mip.height-1)), mip);



    // return magenta for invalid level
    //return Color(1, 0, 1);
  }

    Color Texture::sample_bilinear(Vector2D uv, int level) {
        auto& mip = mipmap[level];
        float u = uv.x*mip.width - 0.5;
        float v = uv.y*mip.height - 0.5;
        int x = floor(u);
        int y = floor(v);
        double u_ratio = u - x;
        double v_ratio = v - y;
        double u_opp = 1 - u_ratio;
        double v_opp = 1 - v_ratio;
        return (mip.get_texel(x, y, mip)*u_opp + mip.get_texel(x+1, y, mip)*u_ratio)*v_opp + (mip.get_texel(x, y+1, mip)*u_opp + mip.get_texel(x+1, y+1, mip)*u_ratio)*v_ratio;
    }

  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
