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
#include <TRint.h>

using namespace std;

int main(int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
//	c1->Divide(2,1);

	int N_Fr = 1000;
	cout << "Flying " << N_Fr << " alpha particles from the mesh." << endl;

	double z_l = 1.0;
	double R_SSD = 6.5;
	double M_H = 70.0;
	double M_W = 70.0;
	double R_MCP = 21.0;

	// Based on SIMION simulation	
	double centerX = 0.0;
	double centerY = 0.0;
	double stdevX = 5.0;
	double stdevY = 1.5;
	cout << "Fr distribution = Nx(" << centerX << ", " << stdevX << ") X Ny(" << centerY << ", " << stdevY << ")" << endl;

	// Geometrical variable of the BPM
	double z_0 = 26.0;
	double x_0 = 26.0;
	cout << "z_0 = " << z_0 << " mm, x_0 = " << x_0 << " mm" << endl;


	// Draw MESH
	TPolyLine3D *mesh = new TPolyLine3D(5);
	mesh->SetPoint(0, 0., -M_W/2.0, z_0+M_H/2.0);
	mesh->SetPoint(1, 0., M_W/2.0, z_0+M_H/2.0);
	mesh->SetPoint(2, 0., M_W/2.0, z_0-M_H/2.0);
	mesh->SetPoint(3, 0., -M_W/2.0, z_0-M_H/2.0);
	mesh->SetPoint(4, 0., -M_W/2.0, z_0+M_H/2.0);
	mesh->SetLineWidth(3);
	mesh->SetLineColor(4);
	mesh->Draw();

	// Draw MCP
	int mcp_points = 1000;
	TPolyLine3D *mcp = new TPolyLine3D(mcp_points);
	for (int k = 0; k < mcp_points; ++k){
		mcp->SetPoint(k, 0.0, R_MCP*TMath::Cos(2.*double(k)*TMath::Pi()/double(mcp_points-1)), z_0+R_MCP*TMath::Sin(2.*double(k)*TMath::Pi()/double(mcp_points-1)));
	}
	mcp->SetLineWidth(3);
	mcp->SetLineColor(4);
	mcp->Draw("SAME");

	//Draw SSD Holder Hole (upper surface)
	int holder_points = 1000;
	TPolyLine3D *holder_u = new TPolyLine3D(holder_points);
	for (int k = 0; k < holder_points; ++k){
		holder_u->SetPoint(k, x_0+R_SSD*TMath::Cos(2.*double(k)*TMath::Pi()/double(holder_points-1)), R_SSD*TMath::Sin(2.*double(k)*TMath::Pi()/double(holder_points-1)), z_l);
	}
	holder_u->SetLineWidth(3);
	holder_u->SetLineColor(3);
	holder_u->Draw("SAME");

	//Draw SSD Holder Hole (lower surface)
	TPolyLine3D *holder_l = new TPolyLine3D(holder_points);
	for (int k = 0; k < holder_points; ++k){
		holder_l->SetPoint(k, x_0+R_SSD*TMath::Cos(2.*double(k)*TMath::Pi()/double(holder_points-1)), R_SSD*TMath::Sin(2.*double(k)*TMath::Pi()/double(holder_points-1)), 0.);
	}
	holder_l->SetLineWidth(3);
	holder_l->SetLineColor(3);
	holder_l->Draw("SAME");


//	TPad *p1 = new TPad("p1","p1",0.0,0.0,1.0,1.0,0,0,0);
//	p1->Draw();
//	p1->cd();
	
//	TView *view = TView::CreateView(1);
//	view->SetRange(-5.,-M_W/2.0 - 5.,-5.,x_0+R_SSD+5.,M_W/2.0 + 5.,M_H + 5.);

	TGraph2D *traj = new TGraph2D();
	traj->SetName("traj");
	traj->SetTitle("Simulated Particle Trajectories; x (mm); y (mm); z (mm)");
//	traj->Draw("P0");

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


		// Draw all trajectories
		double t = TMath::Sqrt(z_0*z_0 + x_0*x_0); // mm
		double x_t = a_x*t;
		double y_t = a_y*t + y_M;
		double z_t = a_z*t + z_M;

		TPolyLine3D *trajectory = new TPolyLine3D(-1);
		trajectory->SetPoint(0,0.0,y_M,z_M);
		trajectory->SetPoint(1,x_t,y_t,z_t);
		trajectory->SetLineWidth(1);
		trajectory->SetLineStyle(2);
		trajectory->Draw("SAME");
	
		// Detection
		if ((lid_upper<0.0)&&(lid_lower<0.0)){
			traj->SetPoint(2*N_Detected,0.0,y_M,z_M);
			traj->SetPoint(2*N_Detected+1,x_l,y_l,0.0);
			trajectory->SetLineColor(2);
			++N_Detected;
		}
	}
	detection = 100.*double(N_Detected)/double(N_Fr);

	cout << detection << "% entered the SSD holder." << endl;

	traj->Draw("SAME P0");
	traj->GetXaxis()->SetLimits(-5.,x_0+R_SSD+5.);
	traj->GetYaxis()->SetLimits(-M_W/2.0-5.,M_W/2.0+5.);
	traj->GetZaxis()->SetLimits(-5.,M_H+5.);



	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
