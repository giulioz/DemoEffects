#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <tuple>
#include "buffer.hh"

/*
 *  2019, Giulio Zausa
 *  gcc checkers.cc -o checkers -l stdc++ -l SDL2 -std=c++17 -Ofast
 */

#define WIDTH 640
#define HEIGHT 480

#define convCoord(x, y) (y * WIDTH + x)

struct Vector3 {
  float x, y, z;

  static inline Vector3 cross(const Vector3 a, const Vector3 b) {
    Vector3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
  }
};

struct Vector4 {
  float x, y, z, w;
};

struct Point {
  float x, y;
};

Vector3 barycentric(const std::array<Point, 3> pts, const Point P) {
  Vector3 a = {pts[2].x - pts[0].x, pts[1].x - pts[0].x, pts[0].x - P.x};
  Vector3 b = {pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - P.y};
  Vector3 u = Vector3::cross(a, b);

  Vector3 dest;
  if (fabsf(u.z) < 1) {
    // degenerate triangle
    dest.x = -1;
    dest.y = 1;
    dest.z = 1;
  } else {
    dest.x = 1.f - (u.x + u.y) / u.z;
    dest.y = u.y / u.z;
    dest.z = u.x / u.z;
  }
  return dest;
}

struct Matrix4 {
  float m11, m21, m31, m41;
  float m12, m22, m32, m42;
  float m13, m23, m33, m43;
  float m14, m24, m34, m44;

  static auto identity() {
    return Matrix4{
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    };
  }

  static auto scale(const Vector3 &scale) {
    return Matrix4{
        scale.x, 0.0, 0.0,     0.0, 0.0, scale.y, 0.0, 0.0,
        0.0,     0.0, scale.z, 0.0, 0.0, 0.0,     0.0, 1.0,
    };
  }

  static auto translation(const Vector3 &t) {
    return Matrix4{
        1.0, 0.0, 0.0, t.x, 0.0, 1.0, 0.0, t.y,
        0.0, 0.0, 1.0, t.z, 0.0, 0.0, 0.0, 1.0,
    };
  }

  static auto rotationX(const float radians) {
    auto c = cos(radians);
    auto s = sin(radians);
    return Matrix4{
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, c,    s,    0.0f,
        0.0f, -s,   c,    0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
    };
  }

  static auto rotationY(const float radians) {
    auto c = cos(radians);
    auto s = sin(radians);
    return Matrix4{
        c, 0.0f, -s, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        s, 0.0f, c,  0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
    };
  }

  static auto rotationZ(const float radians) {
    auto c = cos(radians);
    auto s = sin(radians);
    return Matrix4{
        c,    s,    0.0f, 0.0f, -s,   c,    0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
    };
  }

  static auto perspective() {
    return Matrix4{
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    };
  }

  static auto transformVector3(const Vector3 &v, const Matrix4 &m,
                               float w = 1.0) {
    return Vector3{v.x * m.m11 + v.y * m.m21 + v.z * m.m31 + m.m41 * w,
                   v.x * m.m12 + v.y * m.m22 + v.z * m.m32 + m.m42 * w,
                   v.x * m.m13 + v.y * m.m23 + v.z * m.m33 + m.m43 * w};
  }

  static auto transformVector4(const Vector3 &v, const Matrix4 &m,
                               float w = 1.0) {
    return Vector4{v.x * m.m11 + v.y * m.m21 + v.z * m.m31 + m.m41 * w,
                   v.x * m.m12 + v.y * m.m22 + v.z * m.m32 + m.m42 * w,
                   v.x * m.m13 + v.y * m.m23 + v.z * m.m33 + m.m43 * w,
                   v.x * m.m14 + v.y * m.m24 + v.z * m.m34 + m.m44 * w};
  }

  Matrix4 operator*(Matrix4 value2) {
    Matrix4 m;

    m.m11 = m11 * value2.m11 + m12 * value2.m21 + m13 * value2.m31 +
            m14 * value2.m41;
    m.m12 = m11 * value2.m12 + m12 * value2.m22 + m13 * value2.m32 +
            m14 * value2.m42;
    m.m13 = m11 * value2.m13 + m12 * value2.m23 + m13 * value2.m33 +
            m14 * value2.m43;
    m.m14 = m11 * value2.m14 + m12 * value2.m24 + m13 * value2.m34 +
            m14 * value2.m44;

    m.m21 = m21 * value2.m11 + m22 * value2.m21 + m23 * value2.m31 +
            m24 * value2.m41;
    m.m22 = m21 * value2.m12 + m22 * value2.m22 + m23 * value2.m32 +
            m24 * value2.m42;
    m.m23 = m21 * value2.m13 + m22 * value2.m23 + m23 * value2.m33 +
            m24 * value2.m43;
    m.m24 = m21 * value2.m14 + m22 * value2.m24 + m23 * value2.m34 +
            m24 * value2.m44;

    m.m31 = m31 * value2.m11 + m32 * value2.m21 + m33 * value2.m31 +
            m34 * value2.m41;
    m.m32 = m31 * value2.m12 + m32 * value2.m22 + m33 * value2.m32 +
            m34 * value2.m42;
    m.m33 = m31 * value2.m13 + m32 * value2.m23 + m33 * value2.m33 +
            m34 * value2.m43;
    m.m34 = m31 * value2.m14 + m32 * value2.m24 + m33 * value2.m34 +
            m34 * value2.m44;

    m.m41 = m41 * value2.m11 + m42 * value2.m21 + m43 * value2.m31 +
            m44 * value2.m41;
    m.m42 = m41 * value2.m12 + m42 * value2.m22 + m43 * value2.m32 +
            m44 * value2.m42;
    m.m43 = m41 * value2.m13 + m42 * value2.m23 + m43 * value2.m33 +
            m44 * value2.m43;
    m.m44 = m41 * value2.m14 + m42 * value2.m24 + m43 * value2.m34 +
            m44 * value2.m44;

    return m;
  }
};

template <int width, int height>
inline auto toScreenSpace(const Vector4 &v) {
  auto aspect = (float)height / width;
  auto px = ((v.x / v.w) * aspect * width / 2) + width / 2;
  auto py = height / 2 - ((v.y / v.w) * height / 2);
  return Point{px, py};
}

template <size_t nVertices, size_t nIndices>
struct Mesh {
  std::array<Vector3, nVertices> vertices;
  std::array<int, nIndices> indices;

  Mesh(const std::array<Vector3, nVertices> &vertices,
       const std::array<int, nIndices> &indices) {
    this->vertices = vertices;
    this->indices = indices;
  }
};

template <typename arr>
auto ptsBBox(arr pts) {
  Point bboxmin = {WIDTH - 1, HEIGHT - 1};
  Point bboxmax = {0, 0};
  Point clamp = {WIDTH - 1, HEIGHT - 1};
  for (const auto &pt : pts) {
    bboxmin.x = std::max(0.0f, std::min(bboxmin.x, pt.x));
    bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pt.x));
    bboxmin.y = std::max(0.0f, std::min(bboxmin.y, pt.y));
    bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pt.y));
  }
  return std::tuple(floor(bboxmin.x), floor(bboxmin.y), ceil(bboxmax.x),
                    ceil(bboxmax.y));
}

template <size_t nVertices, size_t nIndices>
void drawMesh(uint32_t *buf, float *zbuffer, uint32_t color,
              const Mesh<nVertices, nIndices> &mesh,
              Matrix4 transform = Matrix4::identity()) {
  auto transformedV = std::array<Vector4, nVertices>{};
  std::transform(
      mesh.vertices.begin(), mesh.vertices.end(), transformedV.begin(),
      [&](auto v) { return Matrix4::transformVector4(v, transform); });

  auto drawTriangle = [&](auto a, auto b, auto c) {
    // Bail out if poligon outside of screen
    // TODO: clipping triangles
    if (transformedV[a].z < 0 || transformedV[b].z < 0 ||
        transformedV[c].z < 0) {
      return;
    }

    const auto pa = toScreenSpace<WIDTH, HEIGHT>(transformedV[a]);
    const auto pb = toScreenSpace<WIDTH, HEIGHT>(transformedV[b]);
    const auto pc = toScreenSpace<WIDTH, HEIGHT>(transformedV[c]);
    const auto points = std::array{pa, pb, pc};
    const auto [bboxminX, bboxminY, bboxmaxX, bboxmaxY] = ptsBBox(points);

    Point p = {0, 0};
    for (p.y = bboxminY; p.y <= bboxmaxY; p.y += 1.0f) {
      for (p.x = bboxminX; p.x <= bboxmaxX; p.x += 1.0f) {
        const auto bary = barycentric(points, p);
        if (bary.x >= 0 && bary.y >= 0 && bary.z >= 0) {
          const int coord = (int)convCoord(p.x, p.y);
          const bool zbufferOver = zbuffer[coord] < bary.z;
          if (zbufferOver) {
            zbuffer[coord] = bary.z;
            buf[coord] = color;
          }
        }
      }
    }
  };

  for (size_t i = 0; i < mesh.indices.size(); i += 3) {
    drawTriangle(mesh.indices[i], mesh.indices[i + 1], mesh.indices[i + 2]);
  }
}

void drawChecker(uint32_t *buf, float *zbuffer, uint32_t colorA,
                 uint32_t colorB, float w, float h,
                 Matrix4 transform = Matrix4::identity()) {
  for (float x = 0; x < w; x += 2.0) {
    for (float z = 0; z < h; z += 2.0) {
      auto mesh = Mesh(
          std::array{
              Vector3{-1.0f + x - w / 2.0f, 0.0f, 1.0f + z - h / 2.0f},
              Vector3{1.0f + x - w / 2.0f, 0.0f, 1.0f + z - h / 2.0f},
              Vector3{1.0f + x - w / 2.0f, 0.0f, -1.0f + z - h / 2.0f},
              Vector3{-1.0f + x - w / 2.0f, 0.0f, -1.0f + z - h / 2.0f},
          },
          std::array{0, 1, 3, 1, 2, 3});
      auto mode = (int)(x + z) % 4 == 0;
      drawMesh(buf, zbuffer, mode ? colorA : colorB, mesh, transform);
    }
  }
}

int main() {
  auto wnd = BufferWindow<WIDTH, HEIGHT>();
  auto zbuffer = new float[WIDTH * HEIGHT];

  auto time = 0.0;
  while (!wnd.checkExit()) {
    wnd.startDraw();
    memset(wnd.pixels, 0, WIDTH * HEIGHT * 4);

    auto transform = Matrix4::rotationY(time / 500.0) *
                     Matrix4::translation(Vector3{0.0, -4.0, 0.0}) *
                     Matrix4::scale(Vector3{0.2, 0.2, 0.2}) *
                     Matrix4::perspective();

    drawChecker(wnd.pixels, zbuffer, 0xFFFF0000, 0xFF0000FF, 100.0f, 100.0f,
                transform);

    std::fill(zbuffer, zbuffer + (WIDTH * HEIGHT), 0);
    wnd.endDraw();

    time += wnd.getDeltaTime();
  }

  return 0;
}
