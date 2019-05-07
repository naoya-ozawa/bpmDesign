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

double t_sph (double ym, double zm, double ax, double ay, double az, double r){
	double a = ax*ax + ay*ay + az*az;
	double b = ay*ym + az*zm;
	double c = ym*ym + zm*zm - r*r;
	return (-b + TMath::Sqrt(b*b - a*c)) / a;
}

int main(int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
	c1->Divide(2,1);

	int N_Fr = 100000;

	double R = 500.0; // Unit test sphere

	double M_H = 70.0;
	double M_W = 70.0;

	// Based on SIMION simulation	
	double centerX = -1.06;
	double centerY = 0.643;
	double stdevX = 3.90;
	double stdevY = 2.37;

	double y_M_ll = -M_W/2.0;
	double y_M_ul = M_W/2.0;
	double z_M_ll = -M_H/2.0;
	double z_M_ul = M_H/2.0;


	random_device rnd;
	default_random_engine engine(rnd());



        c1->cd(1);

	TGraph2D *traj_u = new TGraph2D();
	traj_u->SetName("traj_u");
	traj_u->SetTitle("Uniformly Distributed Direction; x (mm); y (mm); z (mm)");

	uniform_real_distribution<> randhalf(0.,1.);
	uniform_real_distribution<> randall(-1.,1.);

	normal_distribution<> randy(centerX,stdevX);
	normal_distribution<> randz(centerY,stdevY);

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
		double t = t_sph(y_M,z_M,a_x,a_y,a_z,R);
//		traj_u->SetPoint(2*i,0.0,y_M,z_M);
		traj_u->SetPoint(2*i+1,t*a_x,y_M+t*a_y,z_M+t*a_z);
	}

	traj_u->Draw("P,ah,fb,bb");
	traj_u->GetXaxis()->SetLimits(-R-5.,R+5.);
	traj_u->GetYaxis()->SetLimits(-R-5.,R+5.);
	traj_u->GetZaxis()->SetLimits(-R-5.,R+5.);




        c1->cd(2);

	TGraph2D *traj_n = new TGraph2D();
	traj_n->SetName("traj_n");
	traj_n->SetTitle("Normally Distributed Direction; x (mm); y (mm); z (mm)");

	normal_distribution<> normall(0.,R);

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
		double a_x = normall(engine);
		while (a_x <= 0.0){
			a_x = normall(engine);
		}
		double a_y = normall(engine);
		double a_z = normall(engine);
		while (a_z == 0.0){
			a_z = randall(engine);
		}
		double t = t_sph(y_M,z_M,a_x,a_y,a_z,R);
//		traj_n->SetPoint(2*i,0.0,y_M,z_M);
		traj_n->SetPoint(2*i+1,t*a_x,y_M+t*a_y,z_M+t*a_z);
	}

	traj_n->Draw("P,ah,fb,bb");
	traj_n->GetXaxis()->SetLimits(-R-5.,R+5.);
	traj_n->GetYaxis()->SetLimits(-R-5.,R+5.);
	traj_n->GetZaxis()->SetLimits(-R-5.,R+5.);



	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
