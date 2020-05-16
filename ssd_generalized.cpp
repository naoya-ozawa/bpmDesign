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



// Goal: SSD surface
// Returns t_SSD: where trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface z = 0
double t_SSD (double z_M, double a_z){
	if (a_z >= 0.0){
		return 0.0;
	}else{
		return -z_M/a_z;
	}
}

bool HitsSSD (double x_M, double y_M, double a_x, double a_y, double t_SSD, double R_SSD){
	// vec{P} = (x_M,y_M,z_M) + t_SSD*(a_x,a_y,a_z)
	// is within (x-x_0)^2+y^2<R_SSD^2 ?
	double x = x_M + t_SSD*a_x;
	double y = y_M + t_SSD*a_y;
	return (x*x + y*y < R_SSD*R_SSD);
}

// Hit Component 1: SSD case
// Returns TRUE if trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface1 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double R_SSD, double H_SSD){
	// Part 1
	//
	// Define t~_1 that
	// vec{P_1} = (x_M,y_M,z_M) + t~_1*t_SSD*(a_x,a_y,a_z)
	// is on surface z = H_SSD
	// ==> t~_1 = (H_SSD-z_M)/(t_SSD*a_z)
	// Check if x_1^2 + y_1^2 > R_SSD^2

	double t_tilde_1 = (H_SSD - z_M) / (t_SSD*a_z);
	double x_1 = x_M + t_tilde_1*t_SSD*a_x;
	double y_1 = y_M + t_tilde_1*t_SSD*a_y;

	bool H1P1 = (x_1*x_1+y_1*y_1-R_SSD*R_SSD > 0); // HIT if true


	// Part 2
	// 
	// Define t~'_1 that
	// vec{P'_1} = (x_M,y_M,z_M) + t~'_1*t_SSD*(a_x,a_y,a_z)
	// is on surface x'_1^2+y'_1^2 = R_SSD^2
	// Check if 0 < z'_1 < H_SSD

	double coef_a = ((a_x*a_x) + (a_y*a_y))*t_SSD*t_SSD;
	double coef_b = (x_M*a_x + y_M*a_y)*2.0*t_SSD;
	double coef_c = x_M*x_M + y_M*y_M - R_SSD*R_SSD;

	double t_tilde_prime_1 = solve_quad(coef_a,coef_b,coef_c);
	double z_1 = z_M + t_tilde_prime_1*t_SSD*a_z;

	bool H1P2 = (z_1 >= 0.0) && (z_1 <= H_SSD); // HIT if true

	if (H1P1){
		// The trajectory hits the surface of the SSD case
		return true;
	}else if (H1P2){
		// The trajectory hits the side of the SSD case
		return true;
	}else{
		// The trajectory does not hit this component
		return false;
	}
}

// Hit Component 2: SSD box lid
// Returns TRUE if the trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface2 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double R_Box, double H_SSD, double Z_Box){
	// Part 1
	//
	// Define t~_2 that
	// vec{P_2} = (x_M,y_M,z_M) + t~_2*t_SSD*(a_x,a_y,a_z)
	// is on surface z = Z_Box
	// ==> t~_1 = (Z_Box-z_M)/(t_SSD*a_z)
	// Check if x_2^2 + y_2^2 > R_Box^2

	double t_tilde_2 = (Z_Box - z_M) / (t_SSD*a_z);
	double x_2 = x_M + t_tilde_2*t_SSD*a_x;
	double y_2 = y_M + t_tilde_2*t_SSD*a_y;

	bool H2P1 = (x_2*x_2+y_2*y_2-R_Box*R_Box > 0); // HIT if true


	// Part 2
	// 
	// Define t~'_2 that
	// vec{P'_2} = (x_M,y_M,z_M) + t~'_2*t_SSD*(a_x,a_y,a_z)
	// is on surface x'_2^2+y'_2^2 = R_Box^2
	// Check if H_SSD < z'_2 < Z_Box

	double coef_a = ((a_x*a_x) + (a_y*a_y))*t_SSD*t_SSD;
	double coef_b = (x_M*a_x + y_M*a_y)*2.0*t_SSD;
	double coef_c = x_M*x_M + y_M*y_M - R_Box*R_Box;

	double t_tilde_prime_2 = solve_quad(coef_a,coef_b,coef_c);
	double z_2 = z_M + t_tilde_prime_2*t_SSD*a_z;

	bool H2P2 = (z_2 >= H_SSD) && (z_2 <= Z_Box); // HIT if true

	if (H2P1){
		// The trajectory hits the surface of the SSD case
		return true;
	}else if (H2P2){
		// The trajectory hits the side of the SSD case
		return true;
	}else{
		// The trajectory does not hit this component
		return false;
	}
}


// Does the trajectory reach the SSD?
bool reach_ssd (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double* par){
	// The geometry of the FC/SSD
	double R_FC = par[0]; // Radius of FC
	double z_0 = par[1]; // Position of FC bottom
	double z_FC = par[2]; // Height of FC inner side
	double x_Am = par[3]; // X position of Am
	double z_Am = par[4]; // Z position of Am
	double th_Am = par[5]; // tilt angle of Am
	double R_Am = par[6]; // Radius of Am source
	double R_SSD = par[7]; // Radius of SSD
	double H_SSD = par[8]; // Height of SSD case
	double R_Box = par[9]; // Radius of hole in SSD Box
	double Z_Box = par[10]; // Z position of SSD Box surface

	// The particle position is defined as
	// vec{P} = (x_M+t~*t_SSD*a_x, y_M+t~*t_SSD*a_y, z_M+t~*t_SSD*a_z)
	// The particle hits the surface z = 0 when t~ = 1
	double t_s = t_SSD(z_M,a_z);
	cout << "t_SSD = " << t_s << endl;
	if (t_s <= 0.0){
		// discard trajectories that do not hit the z = 0 surface
		return false;
	}else{
		// Hit Component 1
		bool hits_hitsurface1 = HitSurface1(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_SSD,H_SSD);
		cout << "Hit 1: " << hits_hitsurface1 << endl;
		// Hit Component 2
		bool hits_hitsurface2 = HitSurface2(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_Box,H_SSD,Z_Box);
		cout << "Hit 2: " << hits_hitsurface2 << endl;
		// SSD Surface
		bool hits_SSD = HitsSSD(x_M,y_M,a_x,a_y,t_s,R_SSD);
		cout << "Hit SSD: " << hits_SSD << endl;
		if (hits_hitsurface1){
			return false;
		}else if (hits_hitsurface2){
			return false;
		}else if (hits_SSD){
			return true;
		}else{
			cout << "Exception!" << endl;
			return false;
		}
	}
}

// Draw objects on upper half of canvas
int draw_objects(double *par){

	double R_FC = par[0]; // Radius of FC
	double z_0 = par[1]; // Position of FC bottom
	double z_FC = par[2]; // Height of FC inner side
	double x_Am = par[3]; // X position of Am
	double z_Am = par[4]; // Z position of Am
	double th_Am = par[5]; // tilt angle of Am
	double R_Am = par[6]; // Radius of Am source
	double R_SSD = par[7]; // Radius of SSD
	double H_SSD = par[8]; // Height of SSD case
	double R_Box = par[9]; // Radius of hole in SSD Box
	double Z_Box = par[10]; // Z position of SSD Box surface


	// Draw FC bottom
	int fc_points = 1000;
	TPolyLine3D *fcb = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		fcb->SetPoint(k, R_FC*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), R_FC*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)), z_0);
	}
	fcb->SetLineWidth(3);
	fcb->SetLineColor(4);
	fcb->Draw();

	// Draw FC top
	TPolyLine3D *fct = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		fct->SetPoint(k, R_FC*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), R_FC*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)), z_0-z_FC);
	}
	fct->SetLineWidth(3);
	fct->SetLineColor(4);
	fct->Draw();

	// Draw SSD 
	TPolyLine3D *ssd_surf = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		ssd_surf->SetPoint(k, R_SSD*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), R_SSD*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)), 0.0);
	}
	ssd_surf->SetLineWidth(3);
	ssd_surf->SetLineColor(3);
	ssd_surf->Draw();
	TPolyLine3D *ssd_top = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		ssd_top->SetPoint(k, R_SSD*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), R_SSD*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)), H_SSD);
	}
	ssd_top->SetLineWidth(3);
	ssd_top->SetLineColor(3);
	ssd_top->Draw();


	//Draw SSD Box Hole (upper surface)
	TPolyLine3D *holder_u = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		holder_u->SetPoint(k, R_Box*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), R_Box*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)), Z_Box);
	}
	holder_u->SetLineWidth(3);
	holder_u->SetLineColor(3);
	holder_u->Draw();
	//Draw SSD Box Hole (lower surface)
	TPolyLine3D *holder_l = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		holder_l->SetPoint(k, R_Box*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), R_Box*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)), H_SSD);
	}
	holder_l->SetLineWidth(3);
	holder_l->SetLineColor(3);
	holder_l->Draw();

	// Draw Am source
	TPolyLine3D *Amsrc = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		double rhox = R_Am / (TMath::Cos(th_Am) - TMath::Tan(th_Am)*TMath::Sin(th_Am));
		double rho0x = (x_Am*TMath::Cos(th_Am) - x_Am*TMath::Tan(th_Am)*TMath::Sin(th_Am)) / (TMath::Cos(th_Am) - TMath::Tan(th_Am)*TMath::Sin(th_Am));
		double rhoy = R_Am;
		double rhoz = rhox*TMath::Tan(th_Am);
		double rho0z = rho0x*TMath::Tan(th_Am) - x_Am*TMath::Tan(th_Am) + z_Am;
		double xx = rhox*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1))+rho0x;
		double yy = rhoy*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1));
		double zz = rhoz*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1))+rho0z;
		Amsrc->SetPoint(k, xx, yy, zz);
	}
	Amsrc->SetLineWidth(3);
	Amsrc->SetLineColor(9);
	Amsrc->Draw();

	return 0;
}

int main(int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
	c1->Divide(1,2);

//	int N_Fr = 100000; // 10^5 per sec.
	int N_Fr = 10; // for testing
//	int N_Average = 180; // 3 min average
	int N_Average = 1; // for testing
	cout << "Flying " << N_Fr << " alpha particles " << N_Average << " times." << endl;

	// BPM geometry parameters
	double geometry[11];
	geometry[0] = 13.0; // R_FC:Radius of FC
//	geometry[1] = 105.8; // z_0:Position of FC bottom
	geometry[1] = 18; // z_0:Position of FC bottom
	geometry[2] = 15.0; // z_FC:Height of FC inner side
	geometry[3] = 20.01; // x_Am:X position of Am
	geometry[4] = 7.48; // z_Am:Z position of Am
	geometry[5] = 15.0*2.0*TMath::Pi()/360.0; // th_Am:tilt angle of Am
	geometry[6] = 1.0; // R_Am:Radius of Am source
	geometry[7] = 13.8; // R_SSD:Radius of SSD
	geometry[8] = 2.0 ; // H_SSD:Height of SSD case
	geometry[9] = 20.0; // R_Box:Radius of hole in SSD Box
	geometry[10] = 3.1; // Z_Box:Z position of SSD Box surface

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
	double centerX = 0.0;
	double centerY = 0.0;
	double stdevX = 1.5;
	double stdevY = 1.5;
	normal_distribution<double> randx(centerX,stdevX);
	normal_distribution<double> randy(centerY,stdevY);
//	// for case 2
//	normal_distribution<double> randx(0.,2.*geometry[0]);
//	normal_distribution<double> randy(0.,2.*geometry[0]);
//	uniform_real_distribution<double> randz(geometry[1]-geometry[2],geometry[1]);
//	// for case 3
//	normal_distribution<double> randxi(0.,geometry[6]);
//	normal_distribution<double> randyeta(0.,geometry[6]);

	for (int j = 0; j < N_Average; ++j){
		cout << "Start loop j = " << j+1 << endl;
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
			cout << "Start loop i = " << i+1 << endl;

			// Particle generation for Case 1:
			double z_M = geometry[1];
			double x_M, y_M;
			double M = 9999.;
			while (M > geometry[0]){
				x_M = randx(engine);
				y_M = randy(engine);
				M = TMath::Sqrt(x_M*x_M + y_M*y_M);
			}
			double a_x = randnorm(engine);
			double a_y = randnorm(engine);
			double a_z = randnorm(engine);
//			// Calculate approximate solid angle
			double alpha = TMath::ATan(geometry[7]/geometry[1]);
			sa = 2.*TMath::Pi()*(1. - TMath::Cos(alpha));
			// End Case 1

//			// Particle generation for Case 2:
//			double x_M, y_M;
//			double M = 9999.;
//			while ( (M > geometry[0]) || (M < geometry[0]-0.1) ){
//				x_M = randx(engine);
//				y_M = randy(engine);
//				M = TMath::Sqrt(x_M*x_M + y_M*y_M);
//			}
//			double z_M = randz(engine);
//			double a_x = randnorm(engine);
//			double a_y = randnorm(engine);
//			double a_z = randnorm(engine);
//			// Calculate approximate solid angle
//			double alpha = TMath::ATan(geometry[7]/geometry[1]);
//			sa = 2.*TMath::Pi()*(1. - TMath::Cos(alpha));
			// End Case 2


//			// Particle generation for Case 3:
//			double x_M, y_M, z_M;
//			double M = 9999.;
//			while (M > geometry[6]){
//				double xi_M = randxi(engine);
//				double yeta_M = randyeta(engine);
//				double zeta_M = 0.0;
//				M = TMath::Sqrt(xi_M*xi_M + yeta_M*yeta_M);
//				x_M = xi_M*TMath::Cos(th_Am) + zeta_M*TMath::Sin(th_Am) + geometry[3];
//				y_M = yeta_M;
//				z_M = -xi_M*TMath::Sin(th_Am) + zeta_M*TMath::Cos(th_Am) + geometry[4];
//			}
//			double a_x = randnorm(engine);
//			double a_y = randnorm(engine);
//			double a_z = randnorm(engine);
//			// Calculate approximate solid angle
//			double alpha = 5.0 * 2.*TMath::Pi()/360.;
//			sa = 2.*TMath::Pi()*(1. - TMath::Cos(alpha))*TMath::Sin(21.0 * 2.*TMath::Pi()/360.);
//			// End Case 3

			x_mean += x_M;
			x_sqmn += x_M*x_M;
			y_mean += y_M;
			y_sqmn += y_M*y_M;
			z_mean += z_M;
			z_sqmn += z_M*z_M;
			ax_mean += nmvec(a_x,a_y,a_z,"x");
			ay_mean += nmvec(a_x,a_y,a_z,"y");
			az_mean += nmvec(a_x,a_y,a_z,"z");

			bool alpha_detected = reach_ssd(x_M,y_M,z_M,a_x,a_y,a_z,geometry);



			// Detection
			if (alpha_detected){
				++N_Detected;
				++tot_detected;
				// Draw hit trajectories for all hits
//        			if (j > -1){
//					double t = t_SSD(z_M,a_z);
//					double x_t, y_t, z_t;
//					if (t > 0){
//						x_t = x_M + a_x*t;
//						y_t = y_M + a_y*t;
//						z_t = z_M + a_z*t;
//					}else{
//						x_t = x_M - a_x*t;
//						y_t = y_M - a_y*t;
//						z_t = z_M - a_z*t;
//					}
//					TPolyLine3D *trajectory = new TPolyLine3D(-1);
//					trajectory->SetLineWidth(1);
//					trajectory->SetLineColor(2);
//					g_traj->SetPoint(2*(tot_detected-1),x_M,y_M,z_M);
//					trajectory->SetPoint(0,x_M,y_M,z_M);
//					g_traj->SetPoint(2*(tot_detected-1)+1,x_t,y_t,z_t);
//					trajectory->SetPoint(1,x_t,y_t,z_t);
//					trajectory->Draw();
//					cout << "Added trajectory " << tot_detected << endl;
//				}
			}
			if (true){
				double t = t_SSD(z_M,a_z);
				double x_t, y_t, z_t;
				x_t = x_M + a_x*t;
				y_t = y_M + a_y*t;
				z_t = z_M + a_z*t; // = 0
				double M = TMath::Sqrt(x_t*x_t + y_t*y_t);
				if (M > geometry[7]){
					double tR = solve_quad(a_x*a_x+a_y*a_y, 2.*x_M*a_x+2.*y_M*a_y, x_M*x_M+y_M*y_M-geometry[7]*geometry[7]);
					x_t = x_M + a_x*tR;
					y_t = y_M + a_y*tR;
					z_t = z_M + a_z*tR;
				}else{
					cout << "TRAJECTORY DRAWN" << endl;
				}
				TPolyLine3D *trajectory = new TPolyLine3D(-1);
				trajectory->SetLineWidth(1);
				trajectory->SetLineColor(2);
				g_traj->SetPoint(2*(tot_detected-1),x_M,y_M,z_M);
				trajectory->SetPoint(0,x_M,y_M,z_M);
				g_traj->SetPoint(2*(tot_detected-1)+1,x_t,y_t,z_t);
				trajectory->SetPoint(1,x_t,y_t,z_t);
				trajectory->Draw();
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
	g_traj->GetXaxis()->SetLimits(-30.,geometry[3]+10.);
	g_traj->GetYaxis()->SetLimits(-30.,30.);
	g_traj->GetZaxis()->SetLimits(-5.,geometry[1]+5.);


	c1->cd(2);

	TLatex l;
	l.SetTextAlign(12);
	l.SetTextSize(0.05);
//	l.DrawLatex(0.05,0.9,Form("SSD Holder position: (z_{0}, x_{0}) = (%g, %g) [mm]",geometry[5],geometry[4]));
	l.DrawLatex(0.05,0.8,Form("%d #alpha particles flown %d times",N_Fr,N_Average));
	l.DrawLatex(0.05,0.7,"#alpha initial position distribution:");
	l.DrawLatex(0.10,0.6,Form("N_{x}(%3.2f, %3.2f) #times N_{y}(%3.2f, %3.2f) #times N_{z}(%3.2f, %3.2f) [mm]",x_mm,x_stdv,y_mm,y_stdv,z_mm,z_stdv));
	l.DrawLatex(0.05,0.5,"#alpha average direction:");
	l.DrawLatex(0.10,0.4,Form("(%3.2f, %3.2f, %3.2f) (normalized)",ax_mm,ay_mm,az_mm));
	l.DrawLatex(0.05,0.3,Form("%3.3f #pm %3.3f %% of them reached the Si detector.",100.*detection,100.*det_error));
	l.DrawLatex(0.05,0.2,Form("Based on calculated solid angle = %g %%",100.*sa/(4.*TMath::Pi())));


	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
