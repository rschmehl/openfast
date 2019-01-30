#ifndef __KFC_H__
#define __KFC_H__

#ifdef __cplusplus
extern "C" {
#endif

// No major difference between DEBUG mode and Normal mode
// When DEBUG is enabled, a series of print statements will be enabled.
#define DEBUG 1 //Comment out to take out of debug mode

#include <stdio.h>
void controller_init(int *numFlaps, int *numPylons, double I_rot[], double *Requested_dT, int *errStat, char *errMsg); 
void controller_end(int *errStat, char *errMsg);
// void controller_step(double dcm_g2b_c[], double pqr_c[], double *acc_norm_c,
// 		double Xg_c[], double Vg_c[], double Vb_c[], double Ag_c[],
// 		double Ab_c[], double *rho_c, double apparent_wind_c[],
// 		double tether_force_c[], double wind_g_c[],
// 		double kFlapA_c[], double Motor_c[],
// 		int *errStat, char *errMsg);
void controller_step(double dcm_g2b_c[], double pqr_c[], double *acc_norm_c,
					 double Xg_c[], double Vg_c[], double Vb_c[], double Ag_c[],
					 double Ab_c[], double *rho_c, double apparent_wind_c[],
					 double tether_force_c[], double wind_g_c[],
					 double kFlapA_c[], double Motor_c[], double GenTorque[],
					 double RotorBladePitch[], double RotorAccel[], 
					 double RotorSpeed[], double AeroTorque[],
					 int *errStat, char *errMsg);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif

