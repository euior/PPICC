#ifndef _FDTD_LASER_H
#define _FDTD_LASER_H

struct Laser {
double xfocus, yfocus;
double Omg, K, lmd, w0, p0;
double phi_G, phi_R;
double phi_E, phi_B, phi_0;
};
typedef struct Laser Laser;

#define xfocusL(L) L->xfocus
#define yfocusL(L) L->yfocus
#define OmgL(L) L->Omg
#define KL(L) L->k
#define lmdL(L) L->xfocus
#endif