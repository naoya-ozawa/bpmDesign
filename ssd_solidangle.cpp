#include <iostream>
#include <random>
#include <TMath.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TRint.h>

using namespace std;


// Simple solid angle for reference
double solidangle(double z_0, double x_0, double R_SSD){

	double adj = TMath::Sqrt(z_0*z_0 + x_0*x_0);
	double opp = R_SSD/TMath::Sqrt(2.0);
	double hyp = TMath::Sqrt(adj*adj + opp*opp);

	double cosine = adj/hyp;
	double Sr = 2.0 * TMath::Pi() *  (1.0 - cosine);

	return Sr;
}


// Detection efficiency of the SSD
double efficiency(int N_Fr, int N_Average, double z_l, double R_SSD, double M_H, double M_W, double centerX, double stdevX, double centerY, double stdevY, double z_0, double x_0){

	double detection = 0.0;

	double y_M_ll = -M_W/2.0;
	double y_M_ul = M_W/2.0;
	double z_M_ll = z_0 - (M_H/2.0);
	double z_M_ul = z_0 + (M_H/2.0);

	random_device rnd;
	default_random_engine engine(rnd());

	uniform_real_distribution<> randhalf(0.,1.);
	uniform_real_distribution<> randall(-1.,1.);

	normal_distribution<> randy(centerX,stdevX);
	normal_distribution<> randz(centerY+z_0,stdevY);

	for (int j = 0; j < N_Average; j++){

		int N_Detected = 0;

		for (int i = 0; i < N_Fr; ++i){

			// Define Fr on MESH
			double y_M = randy(engine);
			while ((y_M < y_M_ll)||(y_M_ul < y_M)){
				y_M = randy(engine);
			}
			double z_M = randz(engine);
			while ((z_M < z_M_ll)||(z_M_ul < z_M)){
				z_M = randz(engine);
			}

			// Define alpha emitted direction
			double a_x = randhalf(engine);
			double a_y = randall(engine);
			double a_z = randall(engine);
			while (a_z == 0.0){
				a_z = randall(engine);
			}

			// Check if alpha hits upper lid
			double t_u = (z_l - z_M)/a_z;
			double x_u = t_u * a_x;
			double y_u = y_M + t_u * a_y;

			double lid_upper = (x_u-x_0)*(x_u-x_0) + y_u*y_u - R_SSD*R_SSD;


			// Check if alpha hits lower lid
			double t_l = -z_M/a_z;
			double x_l = t_l * a_x;
			double y_l = y_M + t_l * a_y;

			double lid_lower = (x_l-x_0)*(x_l-x_0) + y_l*y_l - R_SSD*R_SSD;


			// Detection
			if ((lid_upper<0.0)&&(lid_lower<0.0)){
				++N_Detected;
			}
		}
		detection += 100.*double(N_Detected)/double(N_Fr);
	}	
	return detection/double(N_Average);
}





int main(int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
	c1->Divide(2,1);

	int N_Fr = 10000;
	int N_Average = 1000;

	cout << "Testing with " << N_Fr << " Fr ions on the mesh." << endl;
	cout << "Taking average of " << N_Average << " trials." << endl;

	double z_l = 1.0;
	double R_SSD = 6.5;
	double M_H = 70.0;
	double M_W = 70.0;

	// Based on SIMION simulation	
	double centerX = 0.0;
	double centerY = 0.0;
	double stdevX = 5.0;
	double stdevY = 1.5;
	cout << "Fr distribution = Nx(" << centerX << ", " << stdevX << ") X Ny(" << centerY << ", " << stdevY << ")" << endl;

	// Geometrical variable of the BPM
	double z_0_mean = 26.0;
	double x_0_mean = 26.0;

	// For (z_0,x_0) optimization
	double halfwidth = 5.0;
	double delta = 5.0;
	int step1d = round(2*halfwidth/delta + 1);
	TGraph2D *zx = new TGraph2D(step1d*step1d);
	zx->SetName("zx");
	zx->SetTitle("Simulated Geometrical Detection Efficiency; z_{0} (mm); x_{0} (mm); #varepsilon_{geometry} (%)");

	TGraph2D *sa = new TGraph2D(step1d*step1d);
	sa->SetName("sa");
	sa->SetTitle("Approx. Acceptance based on Solid Angle; z_{0} (mm); x_{0} (mm); #Omega/(2#pi) (%)");

	double z_0 = z_0_mean - halfwidth;
	double x_0 = x_0_mean - halfwidth;

	for (int k = 0; k < step1d; ++k){
		for (int l = 0; l < step1d; ++l){
			double eff = efficiency(N_Fr, N_Average, z_l, R_SSD, M_H, M_W, centerX, stdevX, centerY, stdevY, z_0, x_0);
			double geoeff = 100.*solidangle(z_0,x_0,R_SSD)/(2.0*TMath::Pi());
//			printf("(z_0, x_0) = (%f, %f) mm : %f %% detection\r",z_0,x_0,eff);
			cout << "(z_0, x_0) = (" << z_0 << ", " << x_0 << ") mm: " << eff << "% detection / Approx. (Sol. Ang./2pi) = " << geoeff << "%" << endl;

			zx->SetPoint(k*step1d+l,z_0,x_0,eff);
			sa->SetPoint(k*step1d+l,z_0,x_0,geoeff);
			x_0 += delta;
		}
		x_0 = x_0_mean - halfwidth;
		z_0 += delta;
	}


	// Draw

	c1->cd(1);
	// The simulated eff_geometry based on the 2D-Gaussian emittance on the mesh
	zx->Draw("SURF1");


	c1->cd(2);
	sa->Draw("SURF1");

	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
