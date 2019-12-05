#include <iostream>
#include <random>
#include <TMath.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPad.h>
#include <TPolyLine3D.h>
#include <TLatex.h>
#include <TRint.h>

using namespace std;
// Introduce a global random number engine
random_device r;
default_random_engine engine(r());

// Normalized 3d-vector
double nmvec(double xx,double yy,double zz,const char* axis){
	double sumsq = xx*xx + yy*yy + zz*zz;
	if (axis == "x"){
		return xx/TMath::Sqrt(sumsq);
	}else if (axis == "y"){
		return yy/TMath::Sqrt(sumsq);
	}else if (axis == "z"){
		return zz/TMath::Sqrt(sumsq);
	}else{
		return 0;
	}
}

// Solver for ax^2 + bx + c = 0 (0 <= x <= 1)
// Returns 0 if no solution exists
double solve_quad (double a, double b, double c){
	double positive = (-b + TMath::Sqrt(b*b - 4.0*a*c))/(2.0*a);
	double negative = (-b - TMath::Sqrt(b*b - 4.0*a*c))/(2.0*a);
	if ((negative < 0.0)&&(positive <= 1.0)){
		return positive;
	}else if ((negative >= 0.0)&&(positive > 1.0)){
		return negative;
	}else{
		return 0.0;
	}
}

// Start: Alpha rays emitted at...
// Returns trajectory as
// (x_M,y_M,z_M,a_x,a_y,a_z) depending on the
// emission surface

// Case 1: from the MCP-IN surface
double alpha_trajectory_1(const char* parameter,double* par){

	double R_MCP = par[0]; // Radius of MCP-IN
	double R_SSD = par[2];
	double x_0 = par[4];
	double z_0 = par[5]; // z_0

	// Based on SIMION simulation (online-3mmCentered-3kV-2810V.pdf)
	double centerX = -0.481228;
	double centerY = -0.433211;
	double stdevX = 2.88206;
	double stdevY = 3.24264;
//	// Based on SIMION simulation (online-2.5mmCentered-3kV-2810V.pdf)
//	double centerX = -0.511427;
//	double centerY = -0.320229;
//	double stdevX = 2.50386;
//	double stdevY = 2.74055;
//	// Based on SIMION simulation (online-2mmCentered-3kV-2810V.pdf)
//	double centerX = -0.453525;
//	double centerY = 0.0106807;
//	double stdevX = 2.13472;
//	double stdevY = 2.15487;
//	// Based on SIMION simulation (online-1.5mmCentered-3kV-2810V.pdf)
//	double centerX = -0.450175;
//	double centerY = 0.165801;
//	double stdevX = 1.79913;
//	double stdevY = 1.59912;
//	// Based on SIMION simulation (online-1mmCentered-3kV-2810V.pdf)
//	double centerX = -0.476828;
//	double centerY = 0.320813;
//	double stdevX = 1.55944;
//	double stdevY = 1.17039;
//	// Based on SIMION simulation (online-0.5mmCentered-3kV-2810V.pdf)
//	double centerX = -0.476919;
//	double centerY = 0.408884;
//	double stdevX = 1.41096;
//	double stdevY = 0.848104;
//	cout << "Fr distribution = Nx(" << centerX << ", " << stdevX << ") X Ny(" << centerY << ", " << stdevY << ")" << endl;


	normal_distribution<> randnorm(0.,1.);

	normal_distribution<> randy(centerX,stdevX);
	normal_distribution<> randz(centerY,stdevY);

	// Define Fr on MCP-IN surface
	double y_M, z_M;
	double M = 9999.;
	while (M > R_MCP){
		y_M = randy(engine);
		z_M = randz(engine);
		M = TMath::Sqrt(y_M*y_M + z_M*z_M);
	}

	// Define alpha emitted direction
	double a_x;
	while (a_x <= 0.0){
		a_x = randnorm(engine);
	}
	double a_y = randnorm(engine);
	double a_z = randnorm(engine);

	// Calculate approximate solid angle
	double theta = TMath::ATan(R_SSD/TMath::Sqrt(x_0*x_0 + z_0*z_0));
	double alpha = TMath::ATan(z_0/x_0);
	double sa = 2.*TMath::Pi()*(1. - TMath::Cos(theta))*TMath::Cos(alpha);

	if (parameter == "x_M"){
		return 0.0;
	}else if (parameter == "y_M"){
		return y_M;
	}else if (parameter == "z_M"){
		return z_M;
	}else if (parameter == "a_x"){
		return a_x;
	}else if (parameter == "a_y"){
		return a_y;
	}else if (parameter == "a_z"){
		return a_z;
	}else{
		return sa;
	}
}

// Case 2: from the MCP lid surface
double alpha_trajectory_2(const char *parameter,double preset_y,double *par){

	double R_MCP = par[0]; // Radius of MCP-IN
	double R_MCPlid = par[1];
	double R_SSD = par[2];
	double x_0 = par[4];
	double z_0 = par[5]; // z_0
	double T_lid = par[6];

	normal_distribution<double> randnorm(0.,1.);
	normal_distribution<double> randrad(0.0,R_MCP);

	if (parameter == "x_M"){
		return T_lid;
	}else if (parameter == "y_M"){
		return randrad(engine);
	}else if (parameter == "z_M"){
		double z_M;
		double M = 9999.;
		while ( (M > R_MCPlid) || (M < R_MCP) ){
			z_M = randrad(engine);
			M = TMath::Sqrt(preset_y*preset_y + z_M*z_M);
		}
		return z_M;
	}else if (parameter == "a_x"){
		double a_x = -1.0;
		while (a_x <= 0.0){
			a_x = randnorm(engine);
//			cout << "Trial a_x = " << a_x << ", ";
		}
//		cout << endl;
		return a_x;
	}else if (parameter == "a_y"){
		return randnorm(engine);
	}else if (parameter == "a_z"){
		return randnorm(engine);
	}else{
		// Calculate approximate solid angle
		double theta = TMath::ATan(R_SSD/TMath::Sqrt((x_0-T_lid)*(x_0-T_lid) + 2.*z_0*2.*z_0));
		double alpha = TMath::ATan(2.*z_0/(x_0-T_lid));
		double sa = 2.*TMath::Pi()*(1. - TMath::Cos(theta))*TMath::Cos(alpha);
		return sa;
	}
}

// Case 5: From the Am source
double alpha_trajectory_5 (const char* parameter, double *par){
	normal_distribution<double> randnorm(0.,1.);
	uniform_real_distribution<double> randx(par[4]-par[3],par[4]+par[3]);
	uniform_real_distribution<double> randy(-par[3],par[3]);

	// Define alpha on Am box lid
	double R = par[3] + 1.0;
	double x_M, y_M;
	while (R >= par[3]){
		x_M = randx(engine);
		y_M = randy(engine);
		R = TMath::Sqrt((x_M-par[4])*(x_M-par[4]) + y_M*y_M);
	}
	double z_M = par[5];
	double a_x = randnorm(engine);
	double a_y = randnorm(engine);
	double a_z;
	while (a_z >= 0.0){
		a_z = randnorm(engine);
	}
//			double a_x = 0.0;
//			double a_y = 0.0;

	// Calculate the simple Sr approximation
	double solid_angle = 2.*TMath::Pi() * (1. - TMath::Cos(TMath::ATan(par[2]/(2.*par[5]))));

	if (parameter == "x_M"){
		return x_M;
	}else if (parameter == "y_M"){
		return y_M;
	}else if (parameter == "z_M"){
		return z_M;
	}else if (parameter == "a_x"){
		return a_x;
	}else if (parameter == "a_y"){
		return a_y;
	}else if (parameter == "a_z"){
		return a_z;
	}else{
		return solid_angle*100./(4.*TMath::Pi());
	}
}

// Goal: SSD surface
// Returns t_SSD: where trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface z = -z_0
double t_SSD (double z_M, double a_z, double z_0){
	if (a_z == 0.0){
		return 0.0;
	}else{
		return -(z_0+z_M)/a_z;
	}
}

bool HitsSSD (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double x_0, double R_SSD){
	// vec{P} = (x_M,y_M,z_M) + t_SSD*(a_x,a_y,a_z)
	// is within (x-x_0)^2+y^2<R_SSD^2 ?
	double x = x_M + t_SSD*a_x;
	double y = y_M + t_SSD*a_y;
	return ((x-x_0)*(x-x_0) + y*y < R_SSD*R_SSD);
}

// Hit Component 1: MCP lid side
// Returns TRUE if trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface1 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double R_MCP, double T_lid){
	// vec{P} = (x_M,y_M,z_M) + t~*t_SSD*(a_x,a_y,a_z)
	// hits y^2+z^2 = R_MCP^2 at t~ = t_hitsurface1 ?
	// (y_M+t~*t_SSD*a_y)^2 + (z_M+t~*t_SSD*a_z)^2 = R_MCP^2
	// {(t_SSD*a_y)^2+(t_SSD*a_z)^2}*t~^2
	// + {2*y_M*t_SSD*a_y+2*z_M*t_SSD*a_z}*t~
	// + y_M^2 + z_M^2 - R_MCP^2

	double coef_a = ((a_y*a_y) + (a_z*a_z))*t_SSD*t_SSD;
	double coef_b = (y_M*a_y + z_M*a_z)*2.0*t_SSD;
	double coef_c = y_M*y_M + z_M*z_M - R_MCP*R_MCP;

	double t_hitsurface1 = solve_quad(coef_a,coef_b,coef_c);
	double x = x_M + t_hitsurface1*t_SSD*a_x;
	double y = y_M + t_hitsurface1*t_SSD*a_y;
	double z = z_M + t_hitsurface1*t_SSD*a_z;

	if ((t_hitsurface1 < 0.0)||(t_hitsurface1 > 1.0)){
		// The trajectory does not pass through here
		return false;
	}else{
		// is the intercept within 0 < x < T_lid ?
		return (x >= 0.0) && (x <= T_lid);
	}
}

// Hit Component 2: SSD holder back + inner side
// Returns TRUE if the trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface2 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double R_MCP, double SSDholder_front, double SSDholder_back, double W_SSDholder, double H_SSDholder){
	// part 1: Holder back

	// vec{P} = (x_M,y_M,z_M) + t~*t_SSD*(a_x,a_y,a_z)
	// hits x = SSDholder_back at t~ = t_hitsurface2_1 ?
	// x_M + t~*t_SSD*a_x = SSDholder_back
	// t~ = (SSDholder_back - x_M) / (t_SSD*a_x)

	bool holderback;
	if (t_SSD*a_x == 0.0){
		holderback = false;
	}else{
		double t_hitsurface2_1 = (SSDholder_back - x_M)/(t_SSD*a_x);
		if ((t_hitsurface2_1 <= 0.0) || (t_hitsurface2_1 >= 1.0)){
			holderback = false;
		}else{
			double x = x_M + t_hitsurface2_1*t_SSD*a_x;
			double y = y_M + t_hitsurface2_1*t_SSD*a_y;
			double z = z_M + t_hitsurface2_1*t_SSD*a_z;
			// -W_SSDholder/2 < y < W_SSDholder/2 ?
			bool holderback_y = (-W_SSDholder/2.0 <= y) && (y <= W_SSDholder/2.0);
			// -H_SSDholder/2 < z < H_SSDholder/2 ?
			bool holderback_z = (-H_SSDholder/2.0 <= z) && (z <= H_SSDholder/2.0);
			// y^2+z^2 > R_MCP^2 ?
			bool openingring = (y*y + z*z >= R_MCP*R_MCP);
			holderback = holderback_y && holderback_z && openingring;
		}
	}

	// part 2: Inner side

	// vec{P} = (x_M,y_M,z_M) + t~*t_SSD*(a_x,a_y,a_z)
	// hits y^2+z^2 = R_MCP^2 at t~ = t_hitsurface2_2 ?
	//(y_M+t~*t_SSD*a_y)^2 + (z_M+t~*t_SSD*a_z) = R_MCP^2
	// {(t_SSD*a_y)^2+(t_SSD*a_z)^2}*t~^2
	// + {2*y_M*t_SSD*a_y+2*z_M*t_SSD*a_z}*t~
	// + y_M^2 + z_M^2 - R_MCP^2

	double coef_a = ((a_y*a_y) + (a_z*a_z))*t_SSD*t_SSD;
	double coef_b = (y_M*a_y + z_M*a_z)*2.0*t_SSD;
	double coef_c = y_M*y_M + z_M*z_M - R_MCP*R_MCP;

	double t_hitsurface2_2 = solve_quad(coef_a,coef_b,coef_c);

	bool inner_side;
	if ((t_hitsurface2_2 <= 0.0) || (t_hitsurface2_2 >= 1.0)){
		// The trajectory does not pass through here
		inner_side = false;
	}else{
		// Is the intercept within SSDholder_back < x < SSDholder_front ?
		double x = x_M + t_hitsurface2_2*t_SSD*a_x;
		double y = y_M + t_hitsurface2_2*t_SSD*a_y;
		double z = z_M + t_hitsurface2_2*t_SSD*a_z;

		inner_side = (x >= SSDholder_back) && (x <= SSDholder_front);
	}

	return (holderback || inner_side);
}

// Hit Component 3: SSD box lid + side
// Returns TRUE if trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface3 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double R_SSD, double x_0, double z_0, double SSDholder_front, double W_SSDboxlid, double H_SSDboxlid){
	// part 1: Box lid

	// vec{P} = (x_M,y_M,z_M) + t~*t_SSD*(a_x,a_y,a_z)
	// hits z = -H_SSDboxlid at t~ = t_hitsurface3_1 ?
	// z_M + t~*t_SSD*a_z = -H_SSDboxlid
	// t~ = (-H_SSDboxlid - z_M) / (t_SSD*a_z)

	bool boxlid;
	if (t_SSD*a_z == 0.0){
		boxlid = false;
	}else{
		double t_hitsurface3_1 = (-H_SSDboxlid - z_M)/(t_SSD*a_z);
		if ((t_hitsurface3_1 <= 0.0) || (t_hitsurface3_1 >= 1.0)){
			boxlid = false;
		}else{
			double x = x_M + t_hitsurface3_1*t_SSD*a_x;
			double y = y_M + t_hitsurface3_1*t_SSD*a_y;
			double z = z_M + t_hitsurface3_1*t_SSD*a_z;
			// (x-x_0)^2+y^2 > R_SSD^2 ?
			bool openingring = ((x-x_0)*(x-x_0) + y*y >= R_SSD*R_SSD);
			// SSDholder_front < x < x_0 ?
			bool boxlid_x = (SSDholder_front <= x) && (x <= x_0);
			// -W_SSDboxlid/2 < y < W_SSDboxlid/2 ?
			bool boxlid_y = (-W_SSDboxlid/2.0 <= y) && (y <= W_SSDboxlid/2.0);
			bool boxlid_rect = boxlid_x && boxlid_y;
			// (x-x_0)^2+y^2 < (W_SSDboxlid/2)^2 ?
			bool boxlid_circ = (x-x_0)*(x-x_0)+y*y < (W_SSDboxlid/2.0)*(W_SSDboxlid/2.0);

			boxlid = (boxlid_rect || boxlid_circ) && openingring;
		}
	}

	// part 2: Inner side

	// vec{P} = (x_M,y_M,z_M) + t~*t_SSD*(a_x,a_y,a_z)
	// hits (x-x_0)^2+y^2 = R_SSD^2 at t~ = t_hitsurface3_2 ?
	// (x_M-x_0+t~*t_SSD*a_x)^2 + (y_M+t~*t_SSD*a_y) = R_SSD^2
	// {(t_SSD*a_x)^2+(t_SSD*a_y)^2}*t~^2
	// + {2*(x_M-x_0)*t_SSD*a_x+2*y_M*t_SSD*a_y}*t~
	// + (x_M-x_0)^2 + y_M^2 - R_SSD^2

	double coef_a = ((a_x*a_x) + (a_y*a_y))*t_SSD*t_SSD;
	double coef_b = ((x_M-x_0)*a_x + y_M*a_y)*2.0*t_SSD;
	double coef_c = (x_M-x_0)*(x_M-x_0) + y_M*y_M - R_SSD*R_SSD;

	double t_hitsurface3_2 = solve_quad(coef_a,coef_b,coef_c);
	bool inner_side;
	if ((t_hitsurface3_2 <= 0.0) || (t_hitsurface3_2 >= 1.0)){
		inner_side = false;
	}else{
		double x = x_M + t_hitsurface3_2*t_SSD*a_x;
		double y = y_M + t_hitsurface3_2*t_SSD*a_y;
		double z = z_M + t_hitsurface3_2*t_SSD*a_z;

		// -z_0 < z < -H_SSDboxlid ?
		inner_side = (-z_0 <= z) && (z <= -H_SSDboxlid);
	}

	return boxlid || inner_side;
}

// Hit Component 4: Am box lid + side
// Returns TRUE if trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface4 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double R_Am, double x_0, double z_0, double SSDholder_front, double W_SSDboxlid, double H_SSDboxlid){
	// part 1: Box lid

	// vec{P} = (x_M,y_M,z_M) + t~*t_SSD*(a_x,a_y,a_z)
	// hits z = H_SSDboxlid at t~ = t_hitsurface4_1 ?
	// z_M + t~*t_SSD*a_z = H_SSDboxlid
	// t~ = (H_SSDboxlid - z_M) / (t_SSD*a_z)

	bool boxlid;
	if (t_SSD*a_z == 0.0){
		boxlid = false;
	}else{
		double t_hitsurface4_1 = (H_SSDboxlid - z_M)/(t_SSD*a_z);
		if ((t_hitsurface4_1 <= 0.0) || (t_hitsurface4_1 >= 1.0)){
			boxlid = false;
		}else{
			double x = x_M + t_hitsurface4_1*t_SSD*a_x;
			double y = y_M + t_hitsurface4_1*t_SSD*a_y;
			double z = z_M + t_hitsurface4_1*t_SSD*a_z;
			// (x-x_0)^2+y^2 > R_SSD^2 ?
			bool openingring = ((x-x_0)*(x-x_0) + y*y >= R_Am*R_Am);
			// SSDholder_front < x < x_0 ?
			bool boxlid_x = (SSDholder_front <= x) && (x <= x_0);
			// -W_SSDboxlid/2 < y < W_SSDboxlid/2 ?
			bool boxlid_y = (-W_SSDboxlid/2.0 <= y) && (y <= W_SSDboxlid/2.0);
			bool boxlid_rect = boxlid_x && boxlid_y;
			// (x-x_0)^2+y^2 < (W_SSDboxlid/2)^2 ?
			bool boxlid_circ = (x-x_0)*(x-x_0)+y*y < (W_SSDboxlid/2.0)*(W_SSDboxlid/2.0);

			boxlid = (boxlid_rect || boxlid_circ) && openingring;
		}
	}

	// part 2: Inner side

	// vec{P} = (x_M,y_M,z_M) + t~*t_SSD*(a_x,a_y,a_z)
	// hits (x-x_0)^2+y^2 = R_Am^2 at t~ = t_hitsurface3_2 ?
	// (x_M-x_0+t~*t_SSD*a_x)^2 + (y_M+t~*t_SSD*a_y) = R_Am^2
	// {(t_SSD*a_x)^2+(t_SSD*a_y)^2}*t~^2
	// + {2*(x_M-x_0)*t_SSD*a_x+2*y_M*t_SSD*a_y}*t~
	// + (x_M-x_0)^2 + y_M^2 - R_Am^2

	double coef_a = ((a_x*a_x) + (a_y*a_y))*t_SSD*t_SSD;
	double coef_b = ((x_M-x_0)*a_x + y_M*a_y)*2.0*t_SSD;
	double coef_c = (x_M-x_0)*(x_M-x_0) + y_M*y_M - R_Am*R_Am;

	double t_hitsurface4_2 = solve_quad(coef_a,coef_b,coef_c);

	bool inner_side;
	if ((t_hitsurface4_2 <= 0.0) || (t_hitsurface4_2 >= 1.0)){
		inner_side = false;
	}else{
		double x = x_M + t_hitsurface4_2*t_SSD*a_x;
		double y = y_M + t_hitsurface4_2*t_SSD*a_y;
		double z = z_M + t_hitsurface4_2*t_SSD*a_z;

		// H_SSDboxlid < z < z_0 ?
		inner_side = (z > H_SSDboxlid) && (z < z_0);
	}

	return boxlid || inner_side;
}

// Does the trajectory reach the SSD?
bool reach_ssd (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double* par){
	// The geometry of the BPM
	double R_MCP = par[0]; // Radius of MCP-IN
	double R_MCPlid = par[1]; // Outer radius of MCP lid
	double R_SSD = par[2]; // Radius of hole on SSD box
	double R_Am = par[3]; // R_Am:Radius of hole on Am box
	double x_0 = par[4]; // SSDholder surface - R_SSD center
	double z_0 = par[5]; // MCP center - SSD box surface
	double T_lid = par[6]; // thickness of MCP lid
	double SSDholder_front = par[7]; // MCP-IN - SSDholder-upstream
	double SSDholder_back = par[8]; // MCP-IN - SSDholder-downstream
	double W_SSDholder = par[9]; // Width of SSD holder
	double H_SSDholder = par[10]; // Height of SSD holder
	double W_SSDboxlid = par[11]; // Width of SSD box
	double H_SSDboxlid = par[12]; // Height of SSD box

	// The particle position is defined as
	// vec{P} = (x_M+t~*t_SSD*a_x, y_M+t~*t_SSD*a_y, z_M+t~*t_SSD*a_z)
	// The particle hits the surface z = -z_0 when t~ = 1
	double t_s = t_SSD(z_M,a_z,z_0);

	if (t_s <= 0.0){
		// discard trajectories that do not hit the z = -z_0 surface
		return false;
	}else{
		// Hit Component 1
		bool hits_hitsurface1 = HitSurface1(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_MCP,T_lid);

		// Hit Component 2
		bool hits_hitsurface2 = HitSurface2(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_MCP,SSDholder_front,SSDholder_back,W_SSDholder,H_SSDholder);

		// Hit Component 3
		bool hits_hitsurface3 = HitSurface3(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_SSD,x_0,z_0,SSDholder_front,W_SSDboxlid,H_SSDboxlid);

		// Hit Component 4
		bool hits_hitsurface4 = HitSurface4(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_Am,x_0,z_0,SSDholder_front,W_SSDboxlid,H_SSDboxlid);

		// SSD Surface
		bool hits_SSD = HitsSSD(x_M,y_M,z_M,a_x,a_y,a_z,t_s,x_0,R_SSD);

		if (hits_hitsurface1){
			return false;
		}else if (hits_hitsurface2){
			return false;
		}else if (hits_hitsurface3){
			return false;
		}else if (hits_hitsurface4){
			return false;
		}else if (hits_SSD){
			return true;
		}else{
			return false;
		}
	}
}

// Draw objects on upper half of canvas
int draw_objects(double *par){

	double R_MCP = par[0]; // Radius of MCP-IN
	double R_MCPlid = par[1]; // Outer radius of MCP lid
	double R_SSD = par[2]; // Radius of hole on SSD box
	double R_Am = par[3]; // Radius of hole on Am box
	double x_0 = par[4]; // SSDholder surface - R_SSD center
	double z_0 = par[5]; // MCP center - SSD box surface
	double T_lid = par[6]; // thickness of MCP lid
	double SSDholder_front = par[7]; // MCP-IN - SSDholder-upstream
	double SSDholder_back = par[8]; // MCP-IN - SSDholder-downstream
	double W_SSDholder = par[9]; // Width of SSD holder
	double H_SSDholder = par[10]; // Height of SSD holder
	double W_SSDboxlid = par[11]; // Width of SSD box
	double H_SSDboxlid = par[12]; // Height of SSD box

	// Draw MCP
	int mcp_points = 1000;
	TPolyLine3D *mcp = new TPolyLine3D(mcp_points);
	for (int k = 0; k < mcp_points; ++k){
		mcp->SetPoint(k, 0.0, R_MCP*TMath::Cos(2.*double(k)*TMath::Pi()/double(mcp_points-1)), R_MCP*TMath::Sin(2.*double(k)*TMath::Pi()/double(mcp_points-1)));
	}
	mcp->SetLineWidth(3);
	mcp->SetLineColor(4);
	mcp->Draw();

	// Draw MCP lid side
	TPolyLine3D *mcp_lid_side = new TPolyLine3D(mcp_points);
	for (int k = 0; k < mcp_points; ++k){
		mcp_lid_side->SetPoint(k, T_lid, R_MCP*TMath::Cos(2.*double(k)*TMath::Pi()/double(mcp_points-1)), R_MCP*TMath::Sin(2.*double(k)*TMath::Pi()/double(mcp_points-1)));
	}
	mcp_lid_side->SetLineWidth(3);
	mcp_lid_side->SetLineColor(4);
	mcp_lid_side->Draw();
	TPolyLine3D *mcp_lid_outer = new TPolyLine3D(mcp_points);
	for (int k = 0; k < mcp_points; ++k){
		mcp_lid_outer->SetPoint(k, T_lid, R_MCPlid*TMath::Cos(2.*double(k)*TMath::Pi()/double(mcp_points-1)), R_MCPlid*TMath::Sin(2.*double(k)*TMath::Pi()/double(mcp_points-1)));
	}
	mcp_lid_outer->SetLineWidth(3);
	mcp_lid_outer->SetLineColor(4);
	mcp_lid_outer->Draw();

	// Draw SSD holder
	TPolyLine3D *holder_front = new TPolyLine3D(5);
	holder_front->SetPoint(0,SSDholder_front,W_SSDholder/2.0,H_SSDholder/2.0);
	holder_front->SetPoint(1,SSDholder_front,-W_SSDholder/2.0,H_SSDholder/2.0);
	holder_front->SetPoint(2,SSDholder_front,-W_SSDholder/2.0,-H_SSDholder/2.0);
	holder_front->SetPoint(3,SSDholder_front,W_SSDholder/2.0,-H_SSDholder/2.0);
	holder_front->SetPoint(4,SSDholder_front,W_SSDholder/2.0,H_SSDholder/2.0);
	holder_front->SetLineWidth(3);
	holder_front->SetLineColor(3);
	holder_front->Draw();
	TPolyLine3D *holder_back = new TPolyLine3D(5);
	holder_back->SetPoint(0,SSDholder_back,W_SSDholder/2.0,H_SSDholder/2.0);
	holder_back->SetPoint(1,SSDholder_back,-W_SSDholder/2.0,H_SSDholder/2.0);
	holder_back->SetPoint(2,SSDholder_back,-W_SSDholder/2.0,-H_SSDholder/2.0);
	holder_back->SetPoint(3,SSDholder_back,W_SSDholder/2.0,-H_SSDholder/2.0);
	holder_back->SetPoint(4,SSDholder_back,W_SSDholder/2.0,H_SSDholder/2.0);
	holder_back->SetLineWidth(3);
	holder_back->SetLineColor(3);
	holder_back->Draw();
	TPolyLine3D *holder_opening_front = new TPolyLine3D(mcp_points);
	for (int k = 0; k < mcp_points; ++k){
		holder_opening_front->SetPoint(k, SSDholder_front, R_MCP*TMath::Cos(2.*double(k)*TMath::Pi()/double(mcp_points-1)), R_MCP*TMath::Sin(2.*double(k)*TMath::Pi()/double(mcp_points-1)));
	}
	holder_opening_front->SetLineWidth(3);
	holder_opening_front->SetLineColor(3);
	holder_opening_front->Draw();
	TPolyLine3D *holder_opening_back = new TPolyLine3D(mcp_points);
	for (int k = 0; k < mcp_points; ++k){
		holder_opening_back->SetPoint(k, SSDholder_back, R_MCP*TMath::Cos(2.*double(k)*TMath::Pi()/double(mcp_points-1)), R_MCP*TMath::Sin(2.*double(k)*TMath::Pi()/double(mcp_points-1)));
	}
	holder_opening_back->SetLineWidth(3);
	holder_opening_back->SetLineColor(3);
	holder_opening_back->Draw();

	// Draw SSD Box contour
	int holder_points = 1000;
	TPolyLine3D *holder_c = new TPolyLine3D(holder_points);
	for (int k = 0; k < holder_points-3; ++k){
		holder_c->SetPoint(k,x_0+(W_SSDboxlid/2.0)*TMath::Sin(double(k)*TMath::Pi()/double(holder_points-4)),(W_SSDboxlid/2.0)*TMath::Cos(double(k)*TMath::Pi()/double(holder_points-4)),-H_SSDboxlid);
	}
	holder_c->SetPoint(holder_points-3,SSDholder_front,-W_SSDboxlid/2.0,-H_SSDboxlid);
	holder_c->SetPoint(holder_points-2,SSDholder_front,W_SSDboxlid/2.0,-H_SSDboxlid);
	holder_c->SetPoint(holder_points-1,x_0,W_SSDboxlid/2.0,-H_SSDboxlid);
	holder_c->SetLineWidth(3);
	holder_c->SetLineColor(3);
	holder_c->Draw();
	//Draw SSD Box Hole (upper surface)
	TPolyLine3D *holder_u = new TPolyLine3D(holder_points);
	for (int k = 0; k < holder_points-2; ++k){
		holder_u->SetPoint(k, x_0+R_SSD*TMath::Cos(2.*double(k)*TMath::Pi()/double(holder_points-1)), R_SSD*TMath::Sin(2.*double(k)*TMath::Pi()/double(holder_points-1)), -H_SSDboxlid);
	}
	holder_u->SetLineWidth(3);
	holder_u->SetLineColor(3);
	holder_u->Draw();
	//Draw SSD Box Hole (lower surface)
	TPolyLine3D *holder_l = new TPolyLine3D(holder_points);
	for (int k = 0; k < holder_points; ++k){
		holder_l->SetPoint(k, x_0+R_SSD*TMath::Cos(2.*double(k)*TMath::Pi()/double(holder_points-1)), R_SSD*TMath::Sin(2.*double(k)*TMath::Pi()/double(holder_points-1)), -z_0);
	}
	holder_l->SetLineWidth(3);
	holder_l->SetLineColor(3);
	holder_l->Draw();

	// Draw Am Box contour
	TPolyLine3D *Ambox_c = new TPolyLine3D(holder_points);
	for (int k = 0; k < holder_points-3; ++k){
		Ambox_c->SetPoint(k,x_0+(W_SSDboxlid/2.0)*TMath::Sin(double(k)*TMath::Pi()/double(holder_points-4)),(W_SSDboxlid/2.0)*TMath::Cos(double(k)*TMath::Pi()/double(holder_points-4)),H_SSDboxlid);
	}
	Ambox_c->SetPoint(holder_points-3,SSDholder_front,-W_SSDboxlid/2.0,H_SSDboxlid);
	Ambox_c->SetPoint(holder_points-2,SSDholder_front,W_SSDboxlid/2.0,H_SSDboxlid);
	Ambox_c->SetPoint(holder_points-1,x_0,W_SSDboxlid/2.0,H_SSDboxlid);
	Ambox_c->SetLineWidth(3);
	Ambox_c->SetLineColor(9);
	Ambox_c->Draw();
	//Draw Am Box Hole (upper surface)
	TPolyLine3D *Ambox_u = new TPolyLine3D(holder_points);
	for (int k = 0; k < holder_points-2; ++k){
		Ambox_u->SetPoint(k, x_0+R_Am*TMath::Cos(2.*double(k)*TMath::Pi()/double(holder_points-1)), R_Am*TMath::Sin(2.*double(k)*TMath::Pi()/double(holder_points-1)), H_SSDboxlid);
	}
	Ambox_u->SetLineWidth(3);
	Ambox_u->SetLineColor(9);
	Ambox_u->Draw();
	//Draw SSD Box Hole (lower surface)
	TPolyLine3D *Ambox_l = new TPolyLine3D(holder_points);
	for (int k = 0; k < holder_points; ++k){
		Ambox_l->SetPoint(k, x_0+R_Am*TMath::Cos(2.*double(k)*TMath::Pi()/double(holder_points-1)), R_Am*TMath::Sin(2.*double(k)*TMath::Pi()/double(holder_points-1)), z_0);
	}
	Ambox_l->SetLineWidth(3);
	Ambox_l->SetLineColor(9);
	Ambox_l->Draw();

	return 0;
}

int main(int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
	c1->Divide(1,2);

	int N_Fr = 100000; // 10^5 per sec.
//	int N_Fr = 1000; // for testing
	int N_Average = 180; // 3 min average
//	int N_Average = 2; // for testing
	cout << "Flying " << N_Fr << " alpha particles " << N_Average << " times." << endl;

	// BPM geometry parameters
	double geometry[13];
	geometry[0] = 14.0; // R_MCP:Radius of MCP-IN
	geometry[1] = 25.0; // R_MCPlid:Outer radius of MCP lid
	geometry[2] = 6.5; // R_SSD:Radius of hole on SSD box
	geometry[3] = 2.5; // R_Am:Radius of hole on Am box
	geometry[4] = 30.8; // x_0:SSDholder surface - SSD center
	geometry[5] = 24.0; // z_0:MCP center - SSD surface
	geometry[6] = 3.5; // T_lid:thickness of MCP lid
	geometry[7] = 16.8; // SSDholder_front:MCP-IN - SSDholder-upstream
	geometry[8] = 13.8; // SSDholder_back:MCP-IN - SSDholder-downstream
	geometry[9] = 50.0; // W_SSDholder:Width of SSD holder
	geometry[10] = 64.0; // H_SSDholder:Height of SSD holder
	geometry[11] = 28.0; // W_SSDboxlid:Width of SSD box
	geometry[12] = 23.0; // H_SSDboxlid:Height of SSD box

	// for CYRIC TOF BPM
//	double z_0 = 29.0;
//	double x_0 = 33.0;

//	cout << "z_0 = " << z_0 << " mm, x_0 = " << x_0 << " mm" << endl;

  c1->cd(1);

	draw_objects(geometry);
	double sa = 0.0; // calculated solid angle

	TGraph2D *g_traj = new TGraph2D();
	g_traj->SetName("traj");
	g_traj->SetTitle("Simulated Particle Trajectories; x (mm); y (mm); z (mm)");

	double detection = 0.0;
	double detect_sq = 0.0;

	int N_Detected;
	int tot_detected = 0;
	double x_mm = 0.0;
	double y_mm = 0.0;
	double z_mm = 0.0;
	double x_stdv = 0.0;
	double y_stdv = 0.0;
	double z_stdv = 0.0;
	double ax_mm = 0.0;
	double ay_mm = 0.0;
	double az_mm = 0.0;

	normal_distribution<double> randnorm(0.,1.);
	// for case 1
//	// Based on SIMION simulation (online-3mmCentered-3kV-2810V.pdf)
//	double centerX = -0.481228;
//	double centerY = -0.433211;
//	double stdevX = 2.88206;
//	double stdevY = 3.24264;
//	cout << "3 mm case" << endl;
//	// Based on SIMION simulation (online-2.5mmCentered-3kV-2810V.pdf)
//	double centerX = -0.511427;
//	double centerY = -0.320229;
//	double stdevX = 2.50386;
//	double stdevY = 2.74055;
//	cout << "2.5 mm case" << endl;
//	// Based on SIMION simulation (online-2mmCentered-3kV-2810V.pdf)
//	double centerX = -0.453525;
//	double centerY = 0.0106807;
//	double stdevX = 2.13472;
//	double stdevY = 2.15487;
//	cout << "2 mm case" << endl;
//	// Based on SIMION simulation (online-1.5mmCentered-3kV-2810V.pdf)
//	double centerX = -0.450175;
//	double centerY = 0.165801;
//	double stdevX = 1.79913;
//	double stdevY = 1.59912;
//	cout << "1.5 mm case" << endl;
//	// Based on SIMION simulation (online-1mmCentered-3kV-2810V.pdf)
//	double centerX = -0.476828;
//	double centerY = 0.320813;
//	double stdevX = 1.55944;
//	double stdevY = 1.17039;
//	cout << "1 mm case" << endl;
	// Based on SIMION simulation (online-0.5mmCentered-3kV-2810V.pdf)
	double centerX = -0.476919;
	double centerY = 0.408884;
	double stdevX = 1.41096;
	double stdevY = 0.848104;
	cout << "0.5 mm case" << endl;
//	cout << "Fr distribution = Nx(" << centerX << ", " << stdevX << ") X Ny(" << centerY << ", " << stdevY << ")" << endl;
	normal_distribution<double> randy(-centerX,stdevX);
	normal_distribution<double> randz(centerY,stdevY);
//	// for case 2
//	normal_distribution<double> randrad(0.0,geometry[0]);
//	// for case 3
//	normal_distribution<double> randrad(0.0,geometry[1]);
//	// for case 4
//	uniform_real_distribution<double> unirandx(geometry[7],geometry[4]+geometry[11]/2.);
//	uniform_real_distribution<double> unirandy(-geometry[11]/2.,geometry[11]/2.);
//	// for case 5
//	uniform_real_distribution<double> unirandx(geometry[4]-geometry[3],geometry[4]+geometry[3]);
//	uniform_real_distribution<double> unirandy(-geometry[3],geometry[3]);

	for (int j = 0; j < N_Average; ++j){
//		cout << "Start loop j = " << j+1 << endl;
		N_Detected = 0;
		double x_mean = 0.0;
		double x_sqmn = 0.0;
		double y_mean = 0.0;
		double y_sqmn = 0.0;
		double z_mean = 0.0;
		double z_sqmn = 0.0;
		double ax_mean = 0.0;
		double ay_mean = 0.0;
		double az_mean = 0.0;
		for (int i = 0; i < N_Fr; ++i){
//			cout << "Start loop i = " << i+1 << endl;


			// Particle generation for Case 1:
			double x_M = 0.0;
			double y_M, z_M;
			double M = 9999.;
			while (M > geometry[0]){
				y_M = randy(engine);
				z_M = randz(engine);
				M = TMath::Sqrt(y_M*y_M + z_M*z_M);
			}
			double a_x = -1.0;
			while (a_x <= 0.0){
				a_x = randnorm(engine);
				//cout << "Trial a_x = " << a_x << ", ";
			}
			//cout << endl;
			double a_y = randnorm(engine);
			double a_z = randnorm(engine);
//			// Calculate approximate solid angle
//			double theta = TMath::ATan(geometry[2]/TMath::Sqrt(geometry[4]*geometry[4] + geometry[5]*geometry[5]));
//			double alpha = TMath::ATan(geometry[5]/geometry[4]);
//			sa = 2.*TMath::Pi()*(1. - TMath::Cos(theta))*TMath::Cos(alpha);
			// End Case 1

//			// Particle generation for Case 2:
//			double x_M = geometry[6];
//			double y_M, z_M;
//			double M = 9999.;
//			while ( (M > geometry[1]) || (M < geometry[0]) ){
//				y_M = randrad(engine);
//				z_M = randrad(engine);
//				M = TMath::Sqrt(y_M*y_M + z_M*z_M);
//			}
//			double a_x = -1.0;
//			while (a_x <= 0.0){
//				a_x = randnorm(engine);
//				//cout << "Trial a_x = " << a_x << ", ";
//			}
//			//cout << endl;
//			double a_y = randnorm(engine);
//			double a_z = randnorm(engine);
//			// Calculate approximate solid angle
//			double theta = TMath::ATan(geometry[2]/TMath::Sqrt((geometry[4]-geometry[6])*(geometry[4]-geometry[6]) + 2.*geometry[5]*2.*geometry[5]));
//			double alpha = TMath::ATan(2.*geometry[5]/(geometry[4]-geometry[6]));
//			sa = 2.*TMath::Pi()*(1. - TMath::Cos(theta))*TMath::Cos(alpha);
			// End Case 2

//			// Particle generation for Case 3:
//			double x_M = geometry[7];
//			double y_M, z_M;
//			double M = 0.0;
//			while (M < geometry[0]){
//				y_M = randrad(engine);
//				z_M = randrad(engine);
//				if ((z_M < geometry[10]/2.0) && (z_M > -geometry[12])){
//					if ((y_M*y_M < geometry[9]*geometry[9]/4.0)){
//						M = TMath::Sqrt(y_M*y_M + z_M*z_M);
//					}else{
//						M = 0.0;
//					}
//				}else{
//					M = 0.0;
//				}
//			}
//			double a_x = -1.0;
//			while (a_x <= 0.0){
//				a_x = randnorm(engine);
//			}
//			double a_y = randnorm(engine);
//			double a_z = randnorm(engine);
//			// Calculate apprximate solid angle
//			double theta = TMath::ATan(geometry[2]/TMath::Sqrt((geometry[4]-geometry[7])*(geometry[4]-geometry[7]) + geometry[5]*geometry[5]));
//			double alpha = TMath::ATan(geometry[5]/(geometry[4]-geometry[7]));
//			sa = 2.*TMath::Pi()*(1. - TMath::Cos(theta))*TMath::Cos(alpha);
//			// End Case 3

//			// Particle generation for Case 4:
//			double x_M, y_M;
//			bool onlid = false;
//			while (!onlid){
//				x_M = unirandx(engine);
//				y_M = unirandy(engine);
//				bool on_hole = (TMath::Sqrt((x_M-geometry[4])*(x_M-geometry[4]) + y_M*y_M) < geometry[3]);
//				bool off_lid = false;
//				if (x_M > geometry[4]){
//					off_lid = (TMath::Sqrt((x_M-geometry[4])*(x_M-geometry[4]) + y_M*y_M) > geometry[11]/2.);
//				}
//				onlid = (!on_hole) && (!off_lid);
//			}
//			double z_M = geometry[12];
//			double a_x = randnorm(engine);
//			double a_y = randnorm(engine);
//			double a_z = 1.0;
//			while (a_z >= 0.0){
//				a_z = randnorm(engine);
//			}
//			// Calculate approximate solid angle
//			sa = 2.*TMath::Pi() * (1. - TMath::Cos(TMath::ATan(geometry[2]/(geometry[5]+geometry[12]))));
//			// End Case 4

//			// Particle generation for Case 5:
//			double x_M, y_M;
//			double M = 9999.;
//			while (M > geometry[3]){
//				x_M = unirandx(engine);
//				y_M = unirandy(engine);
//				M = TMath::Sqrt((x_M-geometry[4])*(x_M-geometry[4]) + y_M*y_M);
//			}
//			double z_M = geometry[5];
//			double a_x = randnorm(engine);
//			double a_y = randnorm(engine);
//			double a_z = 1.0;
//			while (a_z >= 0.0){
//				a_z = randnorm(engine);
//			}
//			// Calculate approximate solid angle
//			sa = 2.*TMath::Pi() * (1. - TMath::Cos(TMath::ATan(geometry[2]/(2.*geometry[5]))));
//			// End Case 5

//			double x_M = alpha_trajectory_2("x_M",0.0,geometry); // <-- gets stuck here
//			double x_M = 0.0;
//			cout << "x_M = " << x_M << ", ";
			x_mean += x_M;
			x_sqmn += x_M*x_M;
//			double y_M = alpha_trajectory_2("y_M",0.0,geometry);
//			y_M = randrad(engine);
//			z_M = randrad(engine);
//			cout << "y_M = " << y_M << ", ";
			y_mean += y_M;
			y_sqmn += y_M*y_M;
//			double z_M = alpha_trajectory_2("z_M",y_M,geometry);
//			cout << "z_M = " << z_M << endl;
			z_mean += z_M;
			z_sqmn += z_M*z_M;
//			double a_x = alpha_trajectory_2("a_x",0.0,geometry);
//			double a_y = alpha_trajectory_2("a_y",0.0,geometry);
//			double a_z = alpha_trajectory_2("a_z",0.0,geometry);
//			double a_x = 1.0;
//			double a_y = randnorm(engine);
//			double a_z = randnorm(engine);
			ax_mean += nmvec(a_x,a_y,a_z,"x");
			ay_mean += nmvec(a_x,a_y,a_z,"y");
			az_mean += nmvec(a_x,a_y,a_z,"z");
//			sa = alpha_trajectory_2("",0.0,geometry);

			bool alpha_detected = reach_ssd(x_M,y_M,z_M,a_x,a_y,a_z,geometry);
//			bool alpha_detected = false;
//			cout << "Hits surface 1?: " << HitSurface1(x_M,y_M,z_M,a_x,a_y,a_z,t_SSD(z_M,a_z,geometry[5]),geometry[0],geometry[6]) << endl;
//			cout << "Hits surface 2?: " << HitSurface2(x_M,y_M,z_M,a_x,a_y,a_z,t_SSD(z_M,a_z,geometry[5]),geometry[0],geometry[7],geometry[8],geometry[9],geometry[10]) << endl;
//			cout << "Hits surface 3?: " << HitSurface3(x_M,y_M,z_M,a_x,a_y,a_z,t_SSD(z_M,a_z,geometry[5]),geometry[2],geometry[4],geometry[5],geometry[7],geometry[11],geometry[12]) << endl;
//			cout << "Entered SSD ?: " << HitsSSD(x_M,y_M,z_M,a_x,a_y,a_z,t_SSD(z_M,a_z,geometry[5]),geometry[4],geometry[2]) << endl;
//			bool check_inregion = (TMath::Sqrt(y_M*y_M + z_M*z_M) > geometry[0]);
//			cout << "Particle generated within region : " << check_inregion << endl;
			// Detection
			if (alpha_detected){
				++N_Detected;
				++tot_detected;
				// Draw hit trajectories for all hits
        			if (j > -1){
					double t = t_SSD(z_M,a_z,geometry[5]);
					double x_t, y_t, z_t;
					if (t > 0){
						x_t = x_M + a_x*t;
						y_t = y_M + a_y*t;
						z_t = z_M + a_z*t;
					}else{
						x_t = x_M - a_x*t;
						y_t = y_M - a_y*t;
						z_t = z_M - a_z*t;
					}
//					double x_t = x_M + 50.*nmvec(a_x,a_y,a_z,"x");
//					double y_t = y_M + 50.*nmvec(a_x,a_y,a_z,"y");
//					double z_t = z_M + 50.*nmvec(a_x,a_y,a_z,"z");
					TPolyLine3D *trajectory = new TPolyLine3D(-1);
					trajectory->SetLineWidth(1);
					trajectory->SetLineColor(2);
					g_traj->SetPoint(2*(tot_detected-1),x_M,y_M,z_M);
					trajectory->SetPoint(0,x_M,y_M,z_M);
					g_traj->SetPoint(2*(tot_detected-1)+1,x_t,y_t,z_t);
					trajectory->SetPoint(1,x_t,y_t,z_t);
					trajectory->Draw();
//					cout << "Added trajectory " << tot_detected << endl;
				}
			}
		}
		detection += double(N_Detected)/double(N_Fr);
		detect_sq += double(N_Detected)*double(N_Detected)/(double(N_Fr)*double(N_Fr));
		x_mean /= double(N_Fr);
		x_sqmn /= double(N_Fr);
		x_mm += x_mean;
		x_stdv += TMath::Sqrt(x_sqmn - x_mean*x_mean);
		y_mean /= double(N_Fr);
		y_sqmn /= double(N_Fr);
		y_mm += y_mean;
		y_stdv += TMath::Sqrt(y_sqmn - y_mean*y_mean);
		z_mean /= double(N_Fr);
		z_sqmn /= double(N_Fr);
		z_mm += z_mean;
		z_stdv += TMath::Sqrt(z_sqmn - z_mean*z_mean);
		ax_mean /= double(N_Fr);
		ay_mean /= double(N_Fr);
		az_mean /= double(N_Fr);
		ax_mm += ax_mean;
		ay_mm += ay_mean;
		az_mm += az_mean;
		cout << j << " sets run" << endl;
	}

	detection /= double(N_Average);
	detect_sq /= double(N_Average);
	double dete_StDev = TMath::Sqrt(detect_sq - (detection*detection));
	double det_error = dete_StDev / TMath::Sqrt(N_Average);
	cout << 100.*detection << " +- " << 100.*det_error << "% entered the SSD holder." << endl;
	x_mm /= double(N_Average);
	x_stdv /= double(N_Average);
	y_mm /= double(N_Average);
	y_stdv /= double(N_Average);
	z_mm /= double(N_Average);
	z_stdv /= double(N_Average);
	ax_mm /= double(N_Average);
	ay_mm /= double(N_Average);
	az_mm /= double(N_Average);

	g_traj->Draw("SAME,P0,ah,fb,bb");
	g_traj->GetXaxis()->SetLimits(-10.,geometry[4]+geometry[2]+10.);
	g_traj->GetYaxis()->SetLimits(-geometry[0]-5.,geometry[0]+5.);
	g_traj->GetZaxis()->SetLimits(-5.,geometry[0]+5.);



	c1->cd(2);

	TLatex l;
	l.SetTextAlign(12);
	l.SetTextSize(0.05);
	l.DrawLatex(0.05,0.9,Form("SSD Holder position: (z_{0}, x_{0}) = (%g, %g) [mm]",geometry[5],geometry[4]));
	l.DrawLatex(0.05,0.8,Form("%d #alpha particles flown %d times",N_Fr,N_Average));
	l.DrawLatex(0.05,0.7,"#alpha initial position distribution:");
	l.DrawLatex(0.10,0.6,Form("N_{x}(%3.2f, %3.2f) #times N_{y}(%3.2f, %3.2f) #times N_{z}(%3.2f, %3.2f) [mm]",x_mm,x_stdv,y_mm,y_stdv,z_mm,z_stdv));
	l.DrawLatex(0.05,0.5,"#alpha average direction:");
	l.DrawLatex(0.10,0.4,Form("(%3.2f, %3.2f, %3.2f) (normalized)",ax_mm,ay_mm,az_mm));
	l.DrawLatex(0.05,0.3,Form("%3.3f #pm %3.3f %% of them reached the Si detector.",100.*detection,100.*det_error));
//	l.DrawLatex(0.05,0.2,Form("Based on calculated solid angle = %g %%",100.*sa/(2.*TMath::Pi())));


	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
