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

// Simple solid angle for reference
double solidangle(double z_0, double x_0, double R_SSD){

	double adj = TMath::Sqrt(z_0*z_0 + x_0*x_0);
	double opp = R_SSD/TMath::Sqrt(2.0);
	double hyp = TMath::Sqrt(adj*adj + opp*opp);

	double cosine = adj/hyp;
	double Sr = 2.0 * TMath::Pi() *  (1.0 - cosine);

	return Sr;
}

int main(int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
	c1->Divide(1,2);

	int N_Fr = 10000000; // 10^7 per sec.
        int N_Average = 60; // 1 min average
	cout << "Flying " << N_Fr << " alpha particles from the mesh." << endl;

	double z_l = 1.0; // Thickness of the SSD lid
	double R_SSD = 6.5; // Radius of the hole on the SSD lid
	double R_MCP = 14.0; //Radius of the MCP-IN

	// Based on SIMION simulation (20190814_01)
	double centerX = -0.42909;
	double centerY = 1.29343;
	double stdevX = 2.44142;
	double stdevY = 1.63478;
	cout << "Fr distribution = Nx(" << centerX << ", " << stdevX << ") X Ny(" << centerY << ", " << stdevY << ")" << endl;

	// Geometrical variable of the BPM
	double z_0 = 24.0;
	double x_0 = 30.4;
	// for CYRIC TOF BPM
//	double z_0 = 29.0;
//	double x_0 = 33.0;

	cout << "z_0 = " << z_0 << " mm, x_0 = " << x_0 << " mm" << endl;

        c1->cd(1);

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
	double detect_sq = 0.0;

	random_device rnd;
	default_random_engine engine(rnd());

	normal_distribution<> randnorm(0.,1.);

	normal_distribution<> randy(centerX,stdevX);
	normal_distribution<> randz(centerY+z_0,stdevY);

        int N_Detected = 0;
	for (int j = 0; j < N_Average; ++j){
		for (int i = 0; i < N_Fr; ++i){

			// Define Fr on MCP-IN surface
			double y_M, z_M;
			double M = R_MCP + 1.0;
			while (M > R_MCP){
				y_M = randy(engine);
				z_M = randz(engine);
				double y_mm = y_M;
				double z_mm = z_M - z_0;
				M = TMath::Sqrt(y_mm*y_mm + z_mm*z_mm);
			}

			// Define alpha emitted direction
			double a_x = randnorm(engine);
			while (a_x <= 0.0){
				a_x = randnorm(engine);
			}
			double a_y = randnorm(engine);
			double a_z = randnorm(engine);
			while (a_z == 0.0){
				a_z = randnorm(engine);
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



			// Draw hit trajectories for the 1st sample
                        if (j == 0){
				double t = TMath::Sqrt(z_0*z_0 + x_0*x_0); // mm
				double x_t = a_x*t;
				double y_t = a_y*t + y_M;
				double z_t = a_z*t + z_M;

				TPolyLine3D *trajectory = new TPolyLine3D(-1);
				trajectory->SetPoint(0,0.0,y_M,z_M);
				trajectory->SetPoint(1,x_t,y_t,z_t);
				trajectory->SetLineWidth(1);
//				trajectory->SetLineStyle(2);
//				trajectory->Draw("SAME");
				if ((lid_upper<0.0)&&(lid_lower<0.0)){
					trajectory->SetLineColor(2);
					traj->SetPoint(2*N_Detected,0.0,y_M,z_M);
					traj->SetPoint(2*N_Detected+1,x_l,y_l,0.0);
					trajectory->Draw("SAME");
				}
			}
	

			// Detection
			if ((lid_upper<0.0)&&(lid_lower<0.0)){
				++N_Detected;
			}
		}
		detection += double(N_Detected)/double(N_Fr);
		detect_sq += double(N_Detected)*double(N_Detected)/(double(N_Fr)*double(N_Fr));
		N_Detected = 0;
	}

        detection /= double(N_Average);
	detect_sq /= double(N_Average);

	double dete_StDev = TMath::Sqrt(detect_sq - (detection*detection));
	cout << 100.*detection << " +- " << 100.*dete_StDev << "% entered the SSD holder." << endl;

	traj->Draw("SAME,P0,ah,fb,bb");
	traj->GetXaxis()->SetLimits(-5.,x_0+R_SSD+5.);
	traj->GetYaxis()->SetLimits(-R_MCP-5.,R_MCP+5.);
	traj->GetZaxis()->SetLimits(-5.,R_MCP+5.);



	c1->cd(2);

	double geoeff = 100.*solidangle(z_0,x_0,R_SSD)/(2.0*TMath::Pi());

	TLatex l;
	l.SetTextAlign(12);
	l.SetTextSize(0.05);
	l.DrawLatex(0.15,0.9,Form("SSD Holder position: (z_{0}, x_{0}) = (%g, %g) [mm]",z_0,x_0));
	l.DrawLatex(0.15,0.8,Form("%d #alpha particles flown %d times in upstream direction",N_Fr,N_Average));
	l.DrawLatex(0.15,0.7,"#alpha initial position distribution on MCP:");
	l.DrawLatex(0.25,0.6,Form("N_{x}(%g, %g) #times N_{y}(%g, %g)",centerX,stdevX,centerY,stdevY));
	l.DrawLatex(0.15,0.5,Form("%g #pm %g %% of them reached the Si detector.",100.*detection,100.*dete_StDev));
	l.DrawLatex(0.15,0.4,Form("For reference: (Solid angle)/(2#pi) = %g %%",geoeff));



	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
