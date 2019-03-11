#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <array>
#include <cmath>
#include <tuple>
#include "buffer.hh"

/*
 *  Spinning 3D Wireframe Cube
 *  gcc cube.cc -o cube -l stdc++ -l SDL2 -std=c++17 -Ofast
 */

#define WIDTH 640
#define HEIGHT 480

#define convCoord(x, y) (y * WIDTH + x)

struct Vector3 {
  float x, y, z;
};

struct Vector4 {
  float x, y, z, w;
};

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
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
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
  return std::tuple(px, py);
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

// Bresenham algorithm
inline void line(uint32_t *buf, uint32_t color, int x1, int y1, int x2,
                 int y2) {
  const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
  if (steep) {
    std::swap(x1, y1);
    std::swap(x2, y2);
  }

  if (x1 > x2) {
    std::swap(x1, x2);
    std::swap(y1, y2);
  }

  const float dx = x2 - x1;
  const float dy = fabs(y2 - y1);

  float error = dx / 2.0f;
  const int ystep = (y1 < y2) ? 1 : -1;
  int y = (int)y1;

  const int maxX = (int)x2;

  for (int x = (int)x1; x < maxX; x++) {
    if (steep) {
      if (x > 0 && x < HEIGHT && y > 0 && y < WIDTH)
        buf[convCoord(y, x)] = color;
    } else {
      if (x > 0 && x < WIDTH && y > 0 && y < HEIGHT)
        buf[convCoord(x, y)] = color;
    }

    error -= dy;
    if (error < 0) {
      y += ystep;
      error += dx;
    }
  }
}

template <size_t nVertices, size_t nIndices>
void drawMeshWireframe(uint32_t *buf, uint32_t color,
                       const Mesh<nVertices, nIndices> &mesh,
                       Matrix4 transform = Matrix4::identity()) {
  auto drawEdge = [&](auto a, auto b) {
    const auto vA = Matrix4::transformVector4(mesh.vertices[a], transform);
    const auto [px, py] = toScreenSpace<WIDTH, HEIGHT>(vA);
    const auto vB = Matrix4::transformVector4(mesh.vertices[b], transform);
    const auto [dx, dy] = toScreenSpace<WIDTH, HEIGHT>(vB);
    line(buf, color, px, py, dx, dy);
  };

  for (size_t i = 0; i < mesh.indices.size() - 1; i += 2) {
    drawEdge(mesh.indices[i], mesh.indices[i + 1]);
  }
}

int main() {
  auto wnd = BufferWindow<WIDTH, HEIGHT>();

  auto time = 0.0;
  while (!wnd.checkExit()) {
    wnd.startDraw();
    memset(wnd.pixels, 0, WIDTH * HEIGHT * 4);

    auto transform =
        Matrix4::rotationX(time / 500) * Matrix4::rotationY(time / 700) *
        Matrix4::rotationZ(time / 600) *
        Matrix4::translation(Vector3{0.0, 0.0, 3.0}) *
        Matrix4::scale(Vector3{0.2, 0.2, 0.2}) * Matrix4::perspective();

    auto cubeMesh = Mesh(
        std::array{
            Vector3{-1.0, 1.0, 1.0},
            Vector3{1.0, 1.0, 1.0},
            Vector3{1.0, -1.0, 1.0},
            Vector3{-1.0, -1.0, 1.0},
            Vector3{-1.0, 1.0, -1.0},
            Vector3{1.0, 1.0, -1.0},
            Vector3{1.0, -1.0, -1.0},
            Vector3{-1.0, -1.0, -1.0},
        },
        std::array{0, 1, 3, 2, 4, 5, 7, 6, 2, 6, 1, 5,
                   0, 4, 7, 3, 0, 3, 1, 2, 4, 7, 5, 6});
    drawMeshWireframe(wnd.pixels, 0xFFFFFFFF, cubeMesh, transform);

    wnd.endDraw();

    time += wnd.getDeltaTime();
  }

  return 0;
}
