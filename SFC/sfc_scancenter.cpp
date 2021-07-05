// Assuming a centralized symmetric gaussian beam,
// scan the beam spread (sigma) and focal point f to plot the dependence of FP_ce on them

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

	int N_Fr = 1000;
	int N_Average = 100;

	// FC geometry parameters
	double geometry[5];
    geometry[0] = 30.0; // R_ou: Outer radius of FP up/down
    geometry[1] = 7.5; // R_in: Inner radius of FP up/down
    geometry[2] = 1.5; // C_fp: Cut of FP up/down
    geometry[3] = 5.0; // R_ce: Radius of center FP
    geometry[4] = -56.08; // z_ce: z position of center FP

	// Generate particles with given beam profile and focal distance
	double centerXY = 0.0;
	double stdevXY = 0.0;
	double f_beam = 0.0;


	c1->cd(1);

	TGraph2D *g_traj = new TGraph2D();
	g_traj->SetName("traj");
	g_traj->SetTitle("FP center hit ratio #frac{FP_{center}}{FP_{up}+FP_{center}+FP_{down}};Beam spread #sigma at FP up/down surface (mm);Beam focal point f wrt FP up/down surface (mm);FP center hit ratio (%)");
	g_traj->SetMarkerStyle(20);
	int g_traj_pt = 0;

	for (int m = 0; m < 20; ++m){
		stdevXY = 2.0 + double(m)*1.0;
		normal_distribution<double> randxy(centerXY,stdevXY);
		for (int n = 0; n < 20; ++n){
			f_beam = double(TMath::Power(n-10,3))*0.1;
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
					double x_M = randxy(engine);
					double y_M = randxy(engine);
					double z_M = 0.0;
					double a_x = centerXY - x_M;
					double a_y = centerXY - y_M;
					double a_z = f_beam - z_M;
					a_x = nmvec(a_x,a_y,a_z,"x");
					a_y = nmvec(a_x,a_y,a_z,"y");
					a_z = nmvec(a_x,a_y,a_z,"z");
					double t = (geometry[4] - z_M) / a_z; // r_P + t a_P = r_C --> t = (z_C-z_P)/a_z
					double x_t = x_M+t*a_x;
					double y_t = y_M+t*a_y;
					double z_t = z_M+t*a_z;
					if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"up")){
						// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Up" << endl;
						++N_up;
					}
					if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"center")){
						// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Center" << endl;
						++N_ce;
					}
					if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"down")){
						// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Down" << endl;
						++N_do;
					}
				}
				N_up_me += double(N_up)/double(N_Average);
				N_up_sq += double(N_up*N_up)/double(N_Average);
				N_ce_me += double(N_ce)/double(N_Average);
				N_ce_sq += double(N_ce*N_ce)/double(N_Average);
				N_do_me += double(N_do)/double(N_Average);
				N_do_sq += double(N_do*N_do)/double(N_Average);
			}
			double N_up_sd = TMath::Sqrt(N_up_me*N_up_me - N_up_sq);
			double N_ce_sd = TMath::Sqrt(N_ce_me*N_ce_me - N_ce_sq);
			double N_do_sd = TMath::Sqrt(N_do_me*N_do_me - N_do_sq);
			double N_count_err = TMath::Max(N_up_sd/N_up_me,TMath::Max(N_ce_sd/N_ce_me,N_do_sd/N_do_me));
			if (N_count_err > 0.0001) cout << "^^^^^ WARNING: LARGE ERROR " << 100.*N_count_err << "% ^^^^^" << endl;
			double rat_ce = 100. * N_ce_me / (N_up_me + N_ce_me + N_do_me);
			g_traj->SetPoint(g_traj_pt,stdevXY,f_beam,rat_ce);
			++g_traj_pt;

			// Progress
			cout << "\rCalculated " << g_traj_pt << " points...";
		}
	}
	cout << endl;

	g_traj->Draw("SURF4");
	g_traj->GetXaxis()->SetTitleOffset(1.8);


	c1->cd(2);

	TGraphErrors *g_spread = new TGraphErrors();
	g_spread->SetName("spread");
	g_spread->SetTitle("FP center hit ratio #frac{FP_{center}}{FP_{up}+FP_{center}+FP_{down}} for laminar beam (f=-1000);Beam spread #sigma at FP up/down surface (mm);FP center hit ratio (%)");
	int g_spread_pt = 0;

	f_beam = -1000.;

	for (int m = 0; m < 40; ++m){
		stdevXY = 2.0 + double(m)*0.5;
		normal_distribution<double> randxy(centerXY,stdevXY);
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
				double x_M = randxy(engine);
				double y_M = randxy(engine);
				double z_M = 0.0;
				double a_x = centerXY - x_M;
				double a_y = centerXY - y_M;
				double a_z = f_beam - z_M;
				a_x = nmvec(a_x,a_y,a_z,"x");
				a_y = nmvec(a_x,a_y,a_z,"y");
				a_z = nmvec(a_x,a_y,a_z,"z");
				double t = (geometry[4] - z_M) / a_z; // r_P + t a_P = r_C --> t = (z_C-z_P)/a_z
				double x_t = x_M+t*a_x;
				double y_t = y_M+t*a_y;
				double z_t = z_M+t*a_z;
				if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"up")){
					// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Up" << endl;
					++N_up;
				}
				if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"center")){
					// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Center" << endl;
					++N_ce;
				}
				if (hit_plate(x_M,y_M,z_M,a_x,a_y,a_z,geometry,"down")){
					// cout << "Particle #" << j*N_Fr+(i+1) << " Hit Down" << endl;
					++N_do;
				}
			}
			N_up_me += double(N_up)/double(N_Average);
			N_up_sq += double(N_up*N_up)/double(N_Average);
			N_ce_me += double(N_ce)/double(N_Average);
			N_ce_sq += double(N_ce*N_ce)/double(N_Average);
			N_do_me += double(N_do)/double(N_Average);
			N_do_sq += double(N_do*N_do)/double(N_Average);
		}
		double N_up_sd = TMath::Sqrt(N_up_me*N_up_me - N_up_sq);
		double N_ce_sd = TMath::Sqrt(N_ce_me*N_ce_me - N_ce_sq);
		double N_do_sd = TMath::Sqrt(N_do_me*N_do_me - N_do_sq);
		double N_tot_me = N_up_me + N_ce_me + N_do_me;
		double N_tot_sd = TMath::Sqrt(N_up_sd*N_up_sd + N_ce_sd*N_ce_sd + N_do_sd*N_do_sd);
		double rat_ce = 100. * N_ce_me / N_tot_me;
		double rat_ce_err = rat_ce * TMath::Sqrt(TMath::Power(N_ce_sd/N_ce_me,2.) + TMath::Power(N_tot_sd/N_tot_me,2.));
		g_spread->SetPoint(g_spread_pt,stdevXY,rat_ce);
		g_spread->SetPointError(g_spread_pt,0.,rat_ce_err);
		++g_spread_pt;
	}

	g_spread->Draw("ALP");

	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
