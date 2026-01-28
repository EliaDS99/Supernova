#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

// <fold> DICHIARAZIONE COST & STRUCT

#define C 299792458 // velocità luce nel vuoto [m/s]
#define M 5 // passi per l'errore

// Dichiarazione della struct per RK4
struct vector {
  double H; // ρ2 = ρ1 * H(ξ) * (γ + 1) / (γ - 1)
  double V; // v2 = 2 * v * V(ξ) / (γ + 1)
  double P; // p2 = 2 * ρ1 * v² * P(ξ) / (γ + 1)
  double dr; // passo d'integrazione
  double g; // grado politropica (p = cρ^gamma)
};

// </fold>
// <fold> DICHIARAZIONE FUNZIONI

// Dichiarazione fella funzione di addizione tra elementi di struct uguali
struct vector add(struct vector, struct vector);

// Dichiarazione della funzione f1(r, ρ, U) che restituisce dρ e dU
struct vector f(struct vector, double);

// Dichiarazione della funzione di integrazione RK4
struct vector RK4(struct vector, double);

// </fold>

// <fold> MAIN
int main(int argc, char** argv) {

  // <fold> VARIABILI & INPUT

  // Definizione delle variabili
  long int N, i, j;
  double r, rho1, v, I, Rs, t, rho2, v2, p2;
  struct vector y;

  // Controllo sul numero di argomenti inseriti da linea di comando
  if (argc != 3) {
    printf("\n### Inserire i 2 valori per: \n"
           "## Il grado γ della politropica \n"
           "## Selezione valori output: \n"
           "#    0: H, V e P \n"
           "#    1: ρ2, v2 e p2 \n"
           "#    2: E \n");
    return 1;
  }

  // </fold>
  // <fold> INIZIALIZZAZIONE VARIABILI

  y.g = atof(argv[1]);
  y.H = 1.;
  y.V = 1.;
  y.P = 1.;
  r = 1.;
  y.dr = r / 1000.;
  N = (long int) (r / y.dr);
  I = 0.;

  if (atoi(argv[2]) == 0) {
    rho1 = 1.;
    v = 1.;
    Rs = 1.;
    t = 1.;
  } else if (atoi(argv[2]) == 1) {
    printf("Inserisci i valori per la densità ρ1 dell'ambiente e la velocità v dell'onda d'urto:\n");
    printf("ρ1 = ");
    scanf("%lf", &rho1);
    printf("v = ");
    scanf("%lf", &v);
    Rs = 1.;
    t = 1.;
  } else if (atoi(argv[2]) == 2) {
    printf("Inserisci i valori per la densità ρ1 dell'ambiente, la velocità v dell'onda d'urto,\n"
           "il tempo t e il raggio Rs dell'onda d'urto al tempo t:\n");
    printf("ρ1 = ");
    scanf("%lf", &rho1);
    printf("v = ");
    scanf("%lf", &v);
    printf("Rs = ");
    scanf("%lf", &Rs);
    printf("t = ");
    scanf("%lf", &t);
  }

  rho2 = rho1*y.H*(y.g+1)/(y.g-1);
  v2 = v*y.V*2/(y.g+1);
  p2 = rho1*v*v*y.P*2/(y.g+1);

  // </fold>
  // <fold> STAMPA SU TERMINALE

  // Ciclo di integrazione e stampa su terminale (i<N && V'>0)
  for (i = 0; i < N && (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H) > 0.; i++) {
    I += (y.H*y.V*y.V+y.P)*r*r/(N*N*N*y.dr*y.dr*y.dr) * y.dr;
    printf("%.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n", r/(N*y.dr), y.H, y.V, y.P,
    rho1*y.H*(y.g+1)/(y.g-1) / rho2, v*y.V*2/(y.g+1) / v2, rho1*v*v*y.P*2/(y.g+1) / p2, I);
    // r, H, V, P, ρ2, v2, p2, I
    y = RK4(y, r/(N*y.dr));
    r -= y.dr;
  }

  printf("\n# ρ2(1) = %.2e\n"
         "# v2(1) = %.2e\n"
         "# p2(1) = %.2e\n", rho2, v2, p2);

  printf("\n# Energia dell'esplosione: %.2e J\n", 1.28*M_PI*rho1*pow(Rs,5)/(t*t*(y.g*y.g-1)) * I);
  printf("# Coefficiente moltiplicativo k: %.2e\n\n", powf(25.*(y.g*y.g-1) / (32*M_PI*I), 0.2));

  printf("\n# Velocità balistica dell'espansione libera: %.2e J\n", 2. * (1.28*M_PI*rho1*pow(Rs,5)/(t*t*(y.g*y.g-1)) * I) / (4./3.*M_PI * pow(Rs) * rho1)); // 2 * E / Mejected

  // </fold>
  // <fold> STAMPA SU FILE

  // Dichiarazione e definizione del file output.dat
  FILE* output;
  output = fopen("output.dat", "w");

  // Reset dei parametri
  y.H = 1.;
  y.V = 1.;
  y.P = 1.;
  r = N * y.dr;
  I = 0.;

  // Ciclo di integrazione e stampa su file
  for (i = 0; i < N && (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H) > 0.; i++) {
    I += (y.H*y.V*y.V+y.P)*r*r/(N*N*N*y.dr*y.dr*y.dr) * y.dr;
    fprintf(output, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", r/(N*y.dr), y.H, y.V, y.P,
    rho1*y.H*(y.g+1)/(y.g-1) / rho2, v*y.V*2/(y.g+1) / v2, rho1*v*v*y.P*2/(y.g+1) / p2, I);
    // r, H, V, P, ρ2, v2, p2, I
    y = RK4(y, r/(N*y.dr));
    r -= y.dr;
  }
  fclose(output);

  // </fold>
  // <fold> γ DIVERSI

  // Dichiarazione e definizione dei file output
  FILE* output1;
  FILE* output2;
  FILE* output3;
  FILE* output4;
  output1 = fopen("output1.dat", "w");
  output2 = fopen("output2.dat", "w");
  output3 = fopen("output3.dat", "w");
  output4 = fopen("output4.dat", "w");

  // Ciclo per i 4 file
  for (j = 1; j <= 4; j++) {

    // Reset dei parametri
    y.H = 1.;
    y.V = 1.;
    y.P = 1.;
    r = N * y.dr;
    I = 0.;

    // Cicli di integrazione, scelta di γ e stampa su j-esimo file
    if (j == 1) {
      y.g = 1.000000000000001;
      for (i = 0; i < N && (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H) > 0.; i++) {
        I += (y.H*y.V*y.V+y.P)*r*r/(N*N*N*y.dr*y.dr*y.dr) * y.dr;
        fprintf(output1, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", r/(N*y.dr), y.H, y.V, y.P, rho1*y.H*(y.g+1)/(y.g-1), v*y.V*2/(y.g+1), rho1*v*v*y.P*2/(y.g+1), I);
        // r, H, V, P, ρ2, v2, p2, I
        y = RK4(y, r/(N*y.dr));
        r -= y.dr;
      }
    } else if (j == 2) {
      y.g = 9./7.;
      for (i = 0; i < N && (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H) > 0.; i++) {
        I += (y.H*y.V*y.V+y.P)*r*r/(N*N*N*y.dr*y.dr*y.dr) * y.dr;
        fprintf(output2, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", r/(N*y.dr), y.H, y.V, y.P, rho1*y.H*(y.g+1)/(y.g-1), v*y.V*2/(y.g+1), rho1*v*v*y.P*2/(y.g+1), I);
        // r, H, V, P, ρ2, v2, p2, I
        y = RK4(y, r/(N*y.dr));
        r -= y.dr;
      }
    } else if (j == 3) {
      y.g = 7./5.;
      for (i = 0; i < N && (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H) > 0.; i++) {
        I += (y.H*y.V*y.V+y.P)*r*r/(N*N*N*y.dr*y.dr*y.dr) * y.dr;
        fprintf(output3, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", r/(N*y.dr), y.H, y.V, y.P, rho1*y.H*(y.g+1)/(y.g-1), v*y.V*2/(y.g+1), rho1*v*v*y.P*2/(y.g+1), I);
        // r, H, V, P, ρ2, v2, p2, I
        y = RK4(y, r/(N*y.dr));
        r -= y.dr;
      }
    } else if (j == 4) {
      y.g = 5./3.;
      for (i = 0; i < N && (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H) > 0.; i++) {
        I += (y.H*y.V*y.V+y.P)*r*r/(N*N*N*y.dr*y.dr*y.dr) * y.dr;
        fprintf(output4, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", r/(N*y.dr), y.H, y.V, y.P, rho1*y.H*(y.g+1)/(y.g-1), v*y.V*2/(y.g+1), rho1*v*v*y.P*2/(y.g+1), I);
        // r, H, V, P, ρ2, v2, p2, I
        y = RK4(y, r/(N*y.dr));
        r -= y.dr;
      }
    }
  }
  fclose(output1);
  fclose(output2);
  fclose(output3);
  fclose(output4);

  // </fold>
  // <fold> ERRORE

  // Dichiarazione e definizione del file error.dat
  FILE* error;
  error = fopen("error.dat", "w");

  // Reset indicie politropico γ
  y.g = atof(argv[1]);

  // Ciclo per gli M valori del passo di integrazione
  for (j = 1; j < M; j++) {

    // Reset dei parametri
    y.H = 1.;
    y.V = 1.;
    y.P = 1.;
    r = N * y.dr;
    I = 0.;
    y.dr = pow(10,-j); // cambio del passo di integrazione
    N = (long int) (r / y.dr); // conseguente cambio di N

    // Ciclo di integrazione e stampa su file
    for (i = 0; i < N && (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H) > 0.; i++) {
      I += (y.H*y.V*y.V+y.P)*r*r/(N*N*N*y.dr*y.dr*y.dr) * y.dr;
      y = RK4(y, r/(N*y.dr));
      r -= y.dr;
    }
    fprintf(error, "%.10e %.10e\n", y.dr, (4./3.*M_PI * pow(Rs,3) * rho1 * C*C - 1.28*M_PI*rho1*pow(Rs,5)/(t*t*(y.g*y.g-1)) * I) / (4./3.*M_PI * Rs*Rs*Rs * rho1 * C*C));
  }
  fclose(error);

  // </fold>
  // <fold> ISTRUZIONI GNUPLOT & END

  system ("gnuplot 'script.gp' -p");

  return 0;

  // </fold>
}
// </fold>

// <fold> FUNZIONI

// <fold> SOMMA

// Definizione della funzione di somma tra elementi di struct uguali
struct vector add(struct vector y, struct vector k) {
  struct vector sum;
  sum.H = y.H + k.H;
  sum.V = y.V + k.V;
  sum.P = y.P + k.P;
  sum.dr = y.dr;
  sum.g = y.g;
  return sum;
}

// </fold>
// <fold> VETTORE ACCELERAZIONE

// Definizione della funzione del vettore y derivato per dρ e dU
struct vector f(struct vector y, double r) {
  struct vector dy;


  dy.P = y.P * (3.*r*(y.g+1) - y.V * (6. + y.g - 8.*y.g*y.V/(r*(y.g+1.)))) / (2. * ((r*(y.g+1.)/2. - y.V) * (2.*y.V/(y.g+1.) - r) + y.g*y.P*(y.g-1.) / (y.H*(y.g+1.))));
  dy.H = y.H / y.g * (dy.P / y.P + 3. / (r - 2. * y.V / (y.g + 1.)));
  dy.V = dy.H / y.H * (0.5 * (y.g + 1.) * r - y.V) - 2. * y.V / r;


  /*
  dy.H = (y.g+1.) * (6. * y.P * (y.g-1.)/(2.*y.V - r*(y.g+1.)) - y.V*y.H * (3. + 4./r * (2.*y.V/(y.g+1.) - r))) / ((2.*y.V - r*(y.g+1.))*(2.*y.V - r*(y.g+1.)) - 2. * y.P/y.H * y.g*(y.g-1.));
  dy.P = y.P * (y.g * dy.H / y.H + 3. / (2. * y.V / (y.g + 1.) - r));
  dy.V = dy.H / y.H * (0.5 * (y.g + 1.) * r - y.V) - 2. * y.V / r;
  */

  /*
  dy.V = (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H);
  dy.H = y.H * (2. * y.V / r + dy.V) / (r * (y.g + 1.) / 2. - y.V);
  dy.P = y.P * (y.g * dy.H / y.H + 3. / (2. * y.V / (y.g + 1.) - r));
  */

  /*
  dy.H = (y.g+1.) * (6. * y.P * (y.g-1.)/(2.*y.V - r*(y.g+1.)) - y.V*y.H * (3. + 4./r * (2.*y.V/(y.g+1.) - r))) / ((2.*y.V - r*(y.g+1.))*(2.*y.V - r*(y.g+1.)) - 2. * y.P/y.H * y.g*(y.g-1.));
  dy.V = (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H);
  dy.P = y.P * (3.*r*(y.g+1) - y.V * (6. + y.g - 8.*y.g*y.V/(r*(y.g+1.)))) / (2. * ((r*(y.g+1.)/2. - y.V) * (2.*y.V/(y.g+1.) - r) + y.g*y.P*(y.g-1.) / (y.H*(y.g+1.))));
  */

  // H'(H,V,P)
  // H'(V') dy.H = y.H * (2. * y.V / r + dy.V) / (r * (y.g + 1.) / 2. - y.V);
  // H'(P') dy.H = y.H / y.g * (dy.P / y.P + 3. / (r - 2. * y.V / (y.g + 1.)));

  // V'(H,V,P) dy.V = (y.g+1.) * (y.P*(y.g-1.) * (4.*y.V*y.g/(r*(y.g+1.)) - 3.) / y.H + 3.*y.V * (y.V - 0.5*r*(y.g+1.))) / ((r*(y.g+1.)-2*y.V)*(r*(y.g+1.)-2.*y.V) - 2.*y.g*y.P*(y.g-1.)/y.H);
  // V'(H') dy.V = dy.H / y.H * (0.5 * (y.g + 1.) * r - y.V) - 2. * y.V / r;
  // V'(P') dy.V = (dy.P / y.H * (y.g - 1.) / (y.g + 1.) - 1.5 * y.V) / (r - 2. * y.V / (y.g + 1));

  // P'(H,V,P) dy.P = y.P * (3.*r*(y.g+1) - y.V * (6. + y.g - 8.*y.g*y.V/(r*(y.g+1.)))) / (2. * ((r*(y.g+1.)/2. - y.V) * (2.*y.V/(y.g+1.) - r) + y.g*y.P*(y.g-1.) / (y.H*(y.g+1.))));
  // P'(H') dy.P = y.P * (y.g * dy.H / y.H + 3. / (2. * y.V / (y.g + 1.) - r));
  // P'(V') dy.P = (y.g+1.)/(y.g-1.) * y.H * (dy.V * (r - 2.*y.V/(y.g+1.)) + 1.5*y.V);

  return dy;
}

// </fold>
// <fold> RK4

// Definizione della funzione di integrazione RK4
struct vector RK4(struct vector y, double r) {
  struct vector F, k1, k2, k3, k4, yn;

  F = f(y, r);
  k1.H = y.dr * F.H * 0.5;
  k1.V = y.dr * F.V * 0.5;
  k1.P = y.dr * F.P * 0.5;

  yn = add(y, k1);
  F = f(yn, r + y.dr*0.5);
  k2.H = y.dr * F.H * 0.5;
  k2.V = y.dr * F.V * 0.5;
  k2.P = y.dr * F.P * 0.5;

  yn = add(y, k2);
  F = f(yn, r + y.dr*0.5);
  k3.H = y.dr * F.H;
  k3.V = y.dr * F.V;
  k3.P = y.dr * F.P;

  yn = add(y, k3);
  F = f(yn, r + y.dr);
  k4.H = y.dr * F.H;
  k4.V = y.dr * F.V;
  k4.P = y.dr * F.P;

  y.H -= (2.*k1.H + 4.*k2.H + 2.*k3.H + k4.H) / 6.;
  y.V -= (2.*k1.V + 4.*k2.V + 2.*k3.V + k4.V) / 6.;
  y.P -= (2.*k1.P + 4.*k2.P + 2.*k3.P + k4.P) / 6.;

  return y;
}

// </fold>

// </fold>
