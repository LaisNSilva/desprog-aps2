#include <math.h>

#include "fourier.h"

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k] = 0;

        for (int j = 0; j < n; j++) {
            t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
        }
    }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    nft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    if (n == 1){
        t[0] = s[0];
        return;
    }

    double complex sp[n/2+1];
    double complex si[n/2+1];
    int p = 0;
    int i = 0;
    for (int e = 0; e <= n+1; e++){
        if (e%2 == 0){
            sp[p] = s[e];
            p++;
        }
        else {
            si[i]= s[e];
            i++;
        }
    }

    double complex tp[n/2+1];
    double complex ti[n/2+1];
    fft(sp, tp, n/2, sign);
    fft(si, ti, n/2, sign);

    for (int k= 0; k<n/2; k++){
        t[k] = tp[k] + ti[k]*cexp(sign * 2* PI * k * I / n);
        t[k+n/2] = tp[k] - ti[k]*cexp(sign * 2* PI * k * I / n);
    }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    fft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
  for (int l = 0; l<height; l++){
      double complex list[width];
      double complex t[width];
      for(int j=0; j<width; j++){
          list[j]=matrix[l][j];
 
    }
    fft_forward(list, t, width);
    for(int j=0; j<width; j++){
          matrix[l][j]=t[j];
 
    }

  }
  for (int coluna = 0; coluna<width; coluna++){
      double complex list_c[height];
      double complex t_c[height];
      for(int linha=0; linha<height; linha++){
          list_c[linha]=matrix[linha][coluna];
 
    }
    fft_forward(list_c, t_c, height);
    for(int c=0; c<height; c++){
          matrix[c][coluna]=t_c[c];
 
    }

  }


  

}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    for (int l = 0; l<height; l++){
      double complex list[width];
      double complex t[width];
      for(int j=0; j<width; j++){
          list[j]=matrix[l][j];
 
    }
    fft_inverse(list, t, width);
    for(int j=0; j<width; j++){
          matrix[l][j]=t[j];
 
    }

  }

  for (int coluna = 0; coluna<width; coluna++){
      double complex list_c[height];
      double complex t_c[height];
      for(int linha=0; linha<height; linha++){
          list_c[linha]=matrix[linha][coluna];
 
    }
    fft_inverse(list_c, t_c, height);
    for(int c=0; c<height; c++){
          matrix[c][coluna]=t_c[c];
 
    }

  }




}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;

            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x] = g * input[y][x];
        }
    }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
