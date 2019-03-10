#pragma once

#include <SDL2/SDL.h>

template <size_t width, size_t height>
class BufferWindow {
 private:
  SDL_Window *window;
  SDL_Surface *screenSurface;

 public:
  uint32_t *pixels;
  uint64_t last = SDL_GetPerformanceCounter();

  BufferWindow() {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
      printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
    } else {
      window = SDL_CreateWindow("SDL Output", SDL_WINDOWPOS_UNDEFINED,
                                SDL_WINDOWPOS_UNDEFINED, width, height,
                                SDL_WINDOW_SHOWN);
      if (window == NULL) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
      } else {
        screenSurface = SDL_GetWindowSurface(window);
        pixels = (uint32_t *)screenSurface->pixels;
      }
    }
  }

  ~BufferWindow() {
    SDL_DestroyWindow(window);
    SDL_Quit();
  }

  auto checkExit() {
    SDL_Event e;
    int keysPointer = 0;
    while (SDL_PollEvent(&e) != 0) {
      if (e.type == SDL_QUIT) {
        return true;
      }
    }

    return false;
  }

  auto getDeltaTime() {
    auto now = SDL_GetPerformanceCounter();
    auto dt =
        (double)((now - last) * 1000 / (double)SDL_GetPerformanceFrequency());
    last = now;
    return dt;
  }

  void startDraw() { SDL_LockSurface(screenSurface); }

  void endDraw() {
    SDL_UnlockSurface(screenSurface);
    SDL_UpdateWindowSurface(window);
  }
};
