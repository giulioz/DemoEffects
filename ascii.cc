#include "asciiBuffer.hh"
#include <cmath>
#include <ncurses.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define WIDTH 120
#define HEIGHT 50

const char *asciiLuminance = " .:-=+*#%@";
#define getLumAscii(l) (asciiLuminance[(int)floor(l * 10)])

int main() {
  auto wnd = ASCIIBuffer<WIDTH, HEIGHT>();

  float time = 0.0;
  while (1) {
    wnd.startDraw();

    for (int x = 0; x < WIDTH; x++) {
      for (int y = 0; y < HEIGHT; y++) {
        float px = (float)(x - WIDTH / 2) / WIDTH;
        float py = (float)(-y + HEIGHT / 2) / HEIGHT;
        float d1 =
            sin(time / 100.0 +
                sqrt((px + 0.2) * (px + 0.2) + (py + 0.2) * (py + 0.2)) * 16);
        float d = abs(d1);
        wnd.print(x, y, getLumAscii(d));
      }
    }

    wnd.endDraw();
    usleep(30000);
    time += 30;
  }

  return 0;
}
