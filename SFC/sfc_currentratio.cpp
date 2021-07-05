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

// Draw objects on upper half of canvas
int draw_objects(double *par){

    double R_ou = par[0]; // Outer radius of FP up/down
    double R_in = par[1]; // Inner radius of FP up/down
    double C_fp = par[2]; // Cut of FP up/down
    double R_ce = par[3]; // Radius of center FP
    double z_ce = par[4]; // z position of center FP
	double theta = 0.0;

	// Draw FP up
	TPolyLine3D *fpup = new TPolyLine3D(205);
	fpup->SetPoint(0,-R_ou,C_fp,0.0);
	fpup->SetPoint(1,-R_in,C_fp,0.0);
	theta = TMath::ATan(C_fp/R_in);
	for (int k = 0; k < 100; ++k){
		double theta_k = (TMath::Pi()-theta) - (double(k)/100)*((TMath::Pi()-theta) - (theta));
		double x_k = R_in*TMath::Cos(theta_k);
		double y_k = R_in*TMath::Sin(theta_k);
		fpup->SetPoint(k+2, x_k, y_k, 0.0);
	}
	fpup->SetPoint(102,R_in,C_fp,0.0);
	fpup->SetPoint(103,R_ou,C_fp,0.0);
	theta = TMath::ATan(C_fp/R_ou);
	for (int k = 0; k < 100; ++k){
		double theta_k = theta + (double(k)/100)*((TMath::Pi()-theta) - (theta));
		double x_k = R_ou*TMath::Cos(theta_k);
		double y_k = R_ou*TMath::Sin(theta_k);
		fpup->SetPoint(k+104, x_k, y_k, 0.0);
	}
	fpup->SetPoint(204,-R_ou,C_fp,0.0);
	fpup->SetLineWidth(3);
	fpup->SetLineColor(4);
	fpup->Draw();

	// Draw FP down
	TPolyLine3D *fpdo = new TPolyLine3D(205);
	fpdo->SetPoint(0,-R_ou,-C_fp,0.0);
	fpdo->SetPoint(1,-R_in,-C_fp,0.0);
	theta = TMath::ATan(C_fp/R_in);
	for (int k = 0; k < 100; ++k){
		double theta_k = (TMath::Pi()-theta) - (double(k)/100)*((TMath::Pi()-theta) - (theta));
		double x_k = R_in*TMath::Cos(theta_k);
		double y_k = -R_in*TMath::Sin(theta_k);
		fpdo->SetPoint(k+2, x_k, y_k, 0.0);
	}
	fpdo->SetPoint(102,R_in,-C_fp,0.0);
	fpdo->SetPoint(103,R_ou,-C_fp,0.0);
	theta = TMath::ATan(C_fp/R_ou);
	for (int k = 0; k < 100; ++k){
		double theta_k = theta + (double(k)/100)*((TMath::Pi()-theta) - (theta));
		double x_k = R_ou*TMath::Cos(theta_k);
		double y_k = -R_ou*TMath::Sin(theta_k);
		fpdo->SetPoint(k+104, x_k, y_k, 0.0);
	}
	fpdo->SetPoint(204,-R_ou,-C_fp,0.0);
	fpdo->SetLineWidth(3);
	fpdo->SetLineColor(4);
	fpdo->Draw();

	// Draw FP center
	TPolyLine3D *fpce = new TPolyLine3D(100);
	for (int k=0; k<100; ++k){
		fpce->SetPoint(k,R_ce*TMath::Cos(2.*TMath::Pi()*(double(k)/99)),R_ce*TMath::Sin(2.*TMath::Pi()*(double(k)/99)),z_ce);
	}
	fpce->SetLineWidth(3);
	fpce->SetLineColor(4);
	fpce->Draw();

	return 0;
}

// Hit plate?
bool hit_plate(double x_P, double y_P, double z_P, double a_x_P, double a_y_P, double a_z_P, double *par, const char* plate){
    double R_ou = par[0]; // Outer radius of FP up/down
    double R_in = par[1]; // Inner radius of FP up/down
    double C_fp = par[2]; // Cut of FP up/down
    double R_ce = par[3]; // Radius of center FP
    double z_ce = par[4]; // z position of center FP
	if (plate == "up"){
		bool od = (TMath::Sqrt(x_P*x_P + y_P*y_P) < R_ou);
		bool id = (TMath::Sqrt(x_P*x_P + y_P*y_P) > R_in);
		bool doughnut = (od && id);
		bool cut = (y_P > C_fp);
		return (doughnut && cut);
	}else if (plate == "center"){
		bool id = (TMath::Sqrt(x_P*x_P + y_P*y_P) < R_in);
		bool cut = (y_P < C_fp) && (y_P > -C_fp);
		bool thru = id || cut;
		// r_P + t a_P = r_C
		// --> t = (z_C-z_P)/a_z
		double t = (z_ce - z_P) / a_z_P;
		double x_C = x_P + t * a_x_P;
		double y_C = y_P + t * a_y_P;
		return (TMath::Sqrt(x_C*x_C + y_C*y_C) < R_ce) && thru;
	}else if (plate == "down"){
		bool od = (TMath::Sqrt(x_P*x_P + y_P*y_P) < R_ou);
		bool id = (TMath::Sqrt(x_P*x_P + y_P*y_P) > R_in);
		bool doughnut = (od && id);
		bool cut = (y_P < -C_fp);
		return (doughnut && cut);
	}else{
		return false;
	}
}

int main(int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
	c1->Divide(1,2);

	int N_Fr = 5000;
	int N_Average = 100;

	// FC geometry parameters
	double geometry[5];
    geometry[0] = 30.0; // R_ou: Outer radius of FP up/down
    geometry[1] = 7.5; // R_in: Inner radius of FP up/down
    geometry[2] = 1.5; // C_fp: Cut of FP up/down
    geometry[3] = 5.0; // R_ce: Radius of center FP
    geometry[4] = -56.08; // z_ce: z position of center FP

	// Generate particles with given beam profile and focal distance
	double centerX = 0.0;
	double stdevX = 5.0;
	double centerY = 0.0;
	double stdevY = 5.0;
	double f_beam = -1000.;
	normal_distribution<double> randx(centerX,stdevX);
	normal_distribution<double> randy(centerY,stdevY);

	c1->cd(1);

	TGraph2D *g_traj = new TGraph2D();
	g_traj->SetName("traj");
	g_traj->SetTitle("Simulated Particle Trajectories; x (mm); y (mm); z (mm)");
	g_traj->SetMarkerStyle(20);
	int g_traj_pt = 0;

	draw_objects(geometry);

	double N_up_me = 0;
	double N_ce_me = 0;
	double N_do_me = 0;
	double N_up_sq = 0;
	double N_ce_sq = 0;
	double N_do_sq = 0;
	for (int j = 0; j < N_Average; ++j){
//		cout << "Start loop j = " << j+1 << endl;
		int N_up = 0;
		int N_ce = 0;
		int N_do = 0;
		for (int i = 0; i < N_Fr; ++i){
//			cout << "Start loop i = " << i+1 << endl;
			// Create normalized trajectories
			double x_M = randx(engine);
			double y_M = randy(engine);
			double z_M = 0.0;
			double a_x = centerX - x_M;
			double a_y = centerY - y_M;
			double a_z = f_beam - z_M;
			a_x = nmvec(a_x,a_y,a_z,"x");
			a_y = nmvec(a_x,a_y,a_z,"y");
			a_z = nmvec(a_x,a_y,a_z,"z");
			double t = (geometry[4] - z_M) / a_z; // r_P + t a_P = r_C --> t = (z_C-z_P)/a_z
			double x_t = x_M+t*a_x;
			double y_t = y_M+t*a_y;
			double z_t = z_M+t*a_z;
			if (i==0){
				TPolyLine3D *trajectory = new TPolyLine3D(-1);
				trajectory->SetLineWidth(1);
				trajectory->SetLineColor(2);
				trajectory->SetPoint(0,x_M,y_M,z_M);
				trajectory->SetPoint(1,x_t,y_t,z_t);
				trajectory->Draw();
			}
			if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"up")){
				// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Up" << endl;
				++N_up;
				if (i==0){
					g_traj->SetPoint(g_traj_pt,x_M,y_M,z_M);
					++g_traj_pt;
				}
			}
			if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"center")){
				// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Center" << endl;
				++N_ce;
				if (i==0){
					g_traj->SetPoint(g_traj_pt,x_t,y_t,z_t);
					++g_traj_pt;
				}				
			}
			if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"down")){
				// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Down" << endl;
				++N_do;
				if (i==0){
					g_traj->SetPoint(g_traj_pt,x_M,y_M,z_M);
					++g_traj_pt;
				}
			}
			// Progress
			cout << "\r" << 100.*double(j*N_Fr+(i+1))/double(N_Fr*N_Average) << "% completed...";
		}
		N_up_me += double(N_up)/double(N_Average);
		N_up_sq += double(N_up*N_up)/double(N_Average);
		N_ce_me += double(N_ce)/double(N_Average);
		N_ce_sq += double(N_ce*N_ce)/double(N_Average);
		N_do_me += double(N_do)/double(N_Average);
		N_do_sq += double(N_do*N_do)/double(N_Average);
	}
	cout << endl;
	double N_up_sd = TMath::Sqrt(N_up_me*N_up_me - N_up_sq);
	double N_ce_sd = TMath::Sqrt(N_ce_me*N_ce_me - N_ce_sq);
	double N_do_sd = TMath::Sqrt(N_do_me*N_do_me - N_do_sq);
	if (TMath::Max(N_up_sd/N_up_me,TMath::Max(N_ce_sd/N_ce_me,N_do_sd/N_do_me)) > 0.0001) cout << "^^^^^ WARNING: LARGE ERROR ^^^^^" << endl;

	g_traj->Draw("SAME,PCOL,ah,fb,bb");
	g_traj->GetXaxis()->SetLimits(-geometry[0]-5.,geometry[0]+5.);
	g_traj->GetYaxis()->SetLimits(-geometry[0]-5.,geometry[0]+5.);
	g_traj->GetZaxis()->SetLimits(geometry[4]-5.,5.);



	c1->cd(2);

	double rat_up = 100. * N_up_me / (N_up_me + N_ce_me + N_do_me);
	double rat_ce = 100. * N_ce_me / (N_up_me + N_ce_me + N_do_me);
	double rat_do = 100. * N_do_me / (N_up_me + N_ce_me + N_do_me);

	TLatex l;
	l.SetTextAlign(12);
	l.SetTextSize(0.05);
	l.DrawLatex(0.05,0.9,Form("Beam center : (#mu[x], #mu[y]) = (%g, %g) [mm]",centerX,centerY));
	l.DrawLatex(0.05,0.8,Form("Beam spread : (#sigma[x], #sigma[y]) = (%g, %g) [mm]",stdevX,stdevY));
	l.DrawLatex(0.05,0.7,Form("Beam focal distance : f = %g [mm]",f_beam));
	l.DrawLatex(0.05,0.6,Form("%d particles flown %d times",N_Fr,N_Average));
	l.DrawLatex(0.05,0.5,"Detected current ratio :");
	l.DrawLatex(0.10,0.4,Form("FP_{up} : FP_{center} : FP_{down} = (%3.2f, %3.2f, %3.2f) [%%]",rat_up,rat_ce,rat_do));

	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
