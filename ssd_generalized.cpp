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
	// is within x^2+y^2<R_SSD^2 ?
	double x = x_M + t_SSD*a_x;
	double y = y_M + t_SSD*a_y;
	return (x*x + y*y < R_SSD*R_SSD);
}

// Hit Component 1: SSD case
// Returns TRUE if trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface1 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double R_SSD){
	// Part 1
	//
	// Define t~_1 that
	// vec{P_1} = (x_M,y_M,z_M) + t~_1*t_SSD*(a_x,a_y,a_z)
	// is on surface z = 0
	// ==> t~_1 = -z_M /(t_SSD*a_z)
	// Check if x_1^2 + y_1^2 > R_SSD^2

	double t_tilde_1 = -z_M / (t_SSD*a_z);
	double x_1 = x_M + t_tilde_1*t_SSD*a_x;
	double y_1 = y_M + t_tilde_1*t_SSD*a_y;

	bool H1P1 = (x_1*x_1+y_1*y_1-R_SSD*R_SSD > 0); // HIT if true


	// Part 2
	// 
	// Define t~'_1 that
	// vec{P'_1} = (x_M,y_M,z_M) + t~'_1*t_SSD*(a_x,a_y,a_z)
	// is on surface x'_1^2+y'_1^2 = R_SSD^2
	// Check if  < z'_1 < 2

	double coef_a = ((a_x*a_x) + (a_y*a_y))*t_SSD*t_SSD;
	double coef_b = (x_M*a_x + y_M*a_y)*2.0*t_SSD;
	double coef_c = x_M*x_M + y_M*y_M - R_SSD*R_SSD;

	double t_tilde_prime_1 = solve_quad(coef_a,coef_b,coef_c);
	double z_1 = z_M + t_tilde_prime_1*t_SSD*a_z;

	bool H1P2 = (z_1 >= 0.0) && (z_1 <= 2.0); // HIT if true

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
bool HitSurface2 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double R_Box, double Z_Box){
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
	// Check if 2 < z'_2 < Z_Box

	double coef_a = ((a_x*a_x) + (a_y*a_y))*t_SSD*t_SSD;
	double coef_b = (x_M*a_x + y_M*a_y)*2.0*t_SSD;
	double coef_c = x_M*x_M + y_M*y_M - R_Box*R_Box;

	double t_tilde_prime_2 = solve_quad(coef_a,coef_b,coef_c);
	double z_2 = z_M + t_tilde_prime_2*t_SSD*a_z;

	bool H2P2 = (z_2 >= 2.0) && (z_2 <= Z_Box); // HIT if true

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
	double R_C = par[0]; // Radius of catcher
	double x_C = par[1]; // x position of catcher center
	double z_C = par[2]; // z position of catcher center
	double R_Am = par[3]; // Radius of Am source
	double x_Am = par[4]; // X position of Am
	double z_Am = par[5]; // Z position of Am
	double R_SSD = par[6]; // Radius of SSD
	double R_Box = par[7]; // Radius of hole in SSD Box
	double Z_Box = par[8]; // Z position of SSD Box surface

	// The particle position is defined as
	// vec{P} = (x_M+t~*t_SSD*a_x, y_M+t~*t_SSD*a_y, z_M+t~*t_SSD*a_z)
	// The particle hits the surface z = 0 when t~ = 1
	double t_s = t_SSD(z_M,a_z);
//	cout << "t_SSD = " << t_s << endl;
	if (t_s <= 0.0){
		// discard trajectories that do not hit the z = 0 surface
		return false;
	}else{
		// Hit Component 1
		bool hits_hitsurface1 = HitSurface1(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_SSD);
//		cout << "Hit 1: " << hits_hitsurface1 << endl;
		// Hit Component 2
		bool hits_hitsurface2 = HitSurface2(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_Box,Z_Box);
//		cout << "Hit 2: " << hits_hitsurface2 << endl;
		// SSD Surface
		bool hits_SSD = HitsSSD(x_M,y_M,a_x,a_y,t_s,R_SSD);
//		cout << "Hit SSD: " << hits_SSD << endl;
		if (hits_hitsurface1){
			return false;
		}else if (hits_hitsurface2){
			return false;
		}else if (hits_SSD){
			return true;
		}else{
//			cout << "Exception! (t_SSD = " << t_s << ")" << endl;
			return false;
		}
	}
}

// Draw objects on upper half of canvas
int draw_objects(double *par){

	double R_C = par[0]; // Radius of catcher
	double x_C = par[1]; // x position of catcher center
	double z_C = par[2]; // z position of catcher center
	double R_Am = par[3]; // Radius of Am source
	double x_Am = par[4]; // X position of Am
	double z_Am = par[5]; // Z position of Am
	double R_SSD = par[6]; // Radius of SSD
	double R_Box = par[7]; // Radius of hole in SSD Box
	double Z_Box = par[8]; // Z position of SSD Box surface

	// Draw catcher
	int fc_points = 1000;
	TPolyLine3D *cat = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		cat->SetPoint(k, x_C, R_C*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), z_C+R_C*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)));
	}
	cat->SetLineWidth(3);
	cat->SetLineColor(4);
	cat->Draw();

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
		ssd_top->SetPoint(k, R_SSD*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), R_SSD*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)), 2.0);
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

	// Draw Am source
	TPolyLine3D *Amsrc = new TPolyLine3D(fc_points);
	for (int k = 0; k < fc_points; ++k){
		Amsrc->SetPoint(k, x_Am+R_Am*TMath::Cos(2.*double(k)*TMath::Pi()/double(fc_points-1)), R_Am*TMath::Sin(2.*double(k)*TMath::Pi()/double(fc_points-1)), z_Am);
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

	int N_Fr = 100000; // 10^5 per sec.
//	int N_Fr = 10; // for testing
	int N_Average = 180; // 3 min average
//	int N_Average = 2; // for testing
	cout << "Flying " << N_Fr << " alpha particles " << N_Average << " times." << endl;

	// BPM geometry parameters
	double geometry[9];
	geometry[0] = 10./2.; // R_C: Radius of catcher
	geometry[1] = -18.08; // x_C: x position of catcher center
	geometry[2] = 11.6; // z_C: z position of catcher center
	geometry[3] = 3.0; // R_Am: Radius of Am source
	geometry[4] = 1.0; // x_Am: X position of Am
	geometry[5] = 11.5+11.6; // z_Am: Z position of Am
	geometry[6] = 13.8/2.; // R_SSD: Radius of SSD
	geometry[7] = 20./2.; // R_Box: Radius of hole in SSD Box
	geometry[8] = 3.1; // Z_Box: Z position of SSD Box surface

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
	// // for case 1 (from catcher)
	// double centerX = geometry[2];
	// double centerY = 0.0;
	// double stdevX = 3.0;
	// double stdevY = 3.0;
	// normal_distribution<double> randx(centerX,stdevX);
	// normal_distribution<double> randy(centerY,stdevY);
	// for case 4 (from Am)
	normal_distribution<double> randx(geometry[4],geometry[3]/4.);
	normal_distribution<double> randy(0.,geometry[3]/4.);

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

			// // Particle generation for Case 1 (catcher):
			// double x_M, y_M, z_M;
			// double M = 9999.;
			// bool offtarget = true;
			// while (offtarget){
			// 	x_M = geometry[1];
			// 	y_M = randy(engine);
			// 	z_M = randx(engine);
			// 	M = TMath::Sqrt((z_M-geometry[2])*(z_M-geometry[2]) + y_M*y_M);
			// 	offtarget = (M >= geometry[0]);
			// }
			// double a_x = randnorm(engine);
			// double a_y = randnorm(engine);
			// double a_z = randnorm(engine);
			// // Calculate approximate solid angle
			// double alpha = TMath::ATan(geometry[6]/TMath::Sqrt(geometry[1]*geometry[1]+geometry[2]*geometry[2]));
			// double beta = TMath::ATan(-geometry[2]/geometry[1]);
			// sa = 2.*TMath::Pi()*(1. - TMath::Cos(alpha))*TMath::Sin(beta);
			// // End Case 1

			// Particle generation for Case 4:
			double x_M, y_M;
			double z_M = geometry[5];
			double M = 9999.;
			while (M > geometry[3]){
				x_M = randx(engine);
				y_M = randy(engine);
				M = TMath::Sqrt((x_M-geometry[4])*(x_M-geometry[4]) + y_M*y_M);
			}
			double a_x = randnorm(engine);
			double a_y = randnorm(engine);
			double a_z = randnorm(engine);
			// Calculate approximate solid angle
			double alpha = TMath::ATan(geometry[6]/(TMath::Sqrt(geometry[5]*geometry[5] + geometry[4]*geometry[4])));
			sa = 2.*TMath::Pi()*(1. - TMath::Cos(alpha));
			// End Case 3

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
				double t = t_SSD(z_M,a_z);
				double x_t, y_t, z_t;
				x_t = x_M + a_x*t;
				y_t = y_M + a_y*t;
				z_t = z_M + a_z*t; // = 0
				// double M = TMath::Sqrt(x_t*x_t + y_t*y_t);
// 				if (M > geometry[6]){
// 					double tR = solve_quad(a_x*a_x+a_y*a_y, 2.*(x_M-geometry[11])*a_x+2.*y_M*a_y, (x_M-geometry[11])*(x_M-geometry[11])+y_M*y_M-geometry[7]*geometry[7]);
// 					x_t = x_M + a_x*tR;
// 					y_t = y_M + a_y*tR;
// 					z_t = z_M + a_z*tR;
// 				}else{
// //					cout << "TRAJECTORY DRAWN" << endl;
// 				}
				TPolyLine3D *trajectory = new TPolyLine3D(-1);
				trajectory->SetLineWidth(1);
				trajectory->SetLineColor(2);
				g_traj->SetPoint(2*(tot_detected-1),x_M,y_M,z_M);
				trajectory->SetPoint(0,x_M,y_M,z_M);
				g_traj->SetPoint(2*(tot_detected-1)+1,x_t,y_t,z_t);
				trajectory->SetPoint(1,x_t,y_t,z_t);
				trajectory->Draw();			
			}
			// Progress
			cout << "\r" << 100.*double(j*N_Fr+(i+1))/double(N_Fr*N_Average) << "% completed...";
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
	}

	detection /= double(N_Average);
	detect_sq /= double(N_Average);
	double dete_StDev = TMath::Sqrt(detect_sq - (detection*detection));
	double det_error = dete_StDev / TMath::Sqrt(N_Average);
	cout << endl;
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
