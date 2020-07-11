#pragma once

#include <ncurses.h>

template <int width, int height> class ASCIIBuffer {
public:
  int maxX = 0, maxY = 0;

  ASCIIBuffer() {
    initscr();
    start_color();
    nodelay(stdscr, TRUE);
    noecho();
    curs_set(FALSE);

    getmaxyx(stdscr, maxX, maxY);
  }

  ~ASCIIBuffer() { endwin(); }

  void startDraw() { clear(); }

  void endDraw() {
    refresh();
    timeout(30);
  }

  void print(const int x, const int y, const char c) { mvaddch(y, x, c); }

  void printLine(const int y, const chtype *line) { mvaddchstr(y, 0, line); }
};
