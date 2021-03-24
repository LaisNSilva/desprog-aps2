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

    
    // acho que a bse vai ser pra quando n=1 pq ao dividir por 2 vai dar zero
    
    if (n == 1){
        t[0] = s[0]*cexp(sign * 2 * PI * 0 * I / n);
        return;
    }

    
    // dividir a entrada original em partes menores

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

    // resolver recursivamente o problema de cada parte
    double complex tp[n/2+1];
    double complex ti[n/2+1];
    fft(sp, tp, n/2, sign);
    tp[n/2-1]+= 0+ sp[n/2-1]*cexp(sign * 2 * PI * n/2-1 * I / n);
    fft(si, ti, n/2, sign);
    ti[n/2-1]+= 0+ si[n/2-1]*cexp(sign * 2 * PI * n/2-1 * I / n);
    // depois de chamar as recursivas não sei como aproveitar
    
    // combinar a solucao para resolver o problema para a entreda original

    for (int k= 0; k<n/2; k++){
        t[k] += tp[k] + ti[k]*cexp(sign * 2* PI * k * I / n);
        t[k+n/2] += tp[k] - ti[k]*cexp(sign * 2* PI * k * I / n);
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
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
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
