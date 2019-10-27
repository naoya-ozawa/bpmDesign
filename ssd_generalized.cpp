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
// Returns t_SSD: where trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface z = -z_0
double t_SSD (double z_M, double a_z, double z_0){
	if (a_z == 0.0){
		return 0.0;
	}else{
		return -(z_0+z_M)/a_z;
	}
}

bool HitsSSD (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double z_0, double x_0){
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
	// (y_M+t~*t_SSD*a_y)^2 + (z_M+t~*t_SSD*a_z) = R_MCP^2
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
		return (0.0 <= x) && (x <= T_lid);
	}
}

// Hit Component 2: SSD holder back + inner side
// Returns TRUE if the trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface2 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double SSDholder_back, double W_SSDholder, double H_SSDholder, double R_MCP, double SSDholder_front){
	
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

	double x = x_M + t_hitsurface2_2*t_SSD*a_x;
	double y = y_M + t_hitsurface2_2*t_SSD*a_y;
	double z = z_M + t_hitsurface2_2*t_SSD*a_z;

	bool inner_side = (SSDholder_back <= x) && (x <= SSDholder_front);

	return holderback || inner_side;
}


// Hit Component 3: SSD box lid + side
// Returns TRUE if trajectory (x_M,y_M,z_M,a_x,a_y,a_z) penetrates the surface
bool HitSurface3 (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double t_SSD, double SSDholder_front, double x_0, double W_SSDboxlid, double R_MCP, double R_SSD, double H_SSDboxlid){

	// part 1: Box lid
	
	// vec{P} = (x_M,y_M,z_M) + t~*t_SSD*(a_x,a_y,a_z)
	// hits z = H_SSDboxlid at t~ = t_hitsurface3_1 ?
	// z_M + t~*t_SSD*a_z = H_SSDboxlid
	// t~ = (H_SSDboxlid - z_M) / (t_SSD*a_z)

	bool boxlid;
	if (t_SSD*a_z == 0.0){
		boxlid = false;
	}else{
		double t_hitsurface3_1 = (H_SSDboxlid - z_M)/(t_SSD*a_z);
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

	double x = x_M + t_hitsurface3_2*t_SSD*a_x;
	double y = y_M + t_hitsurface3_2*t_SSD*a_y;
	double z = z_M + t_hitsurface3_2*t_SSD*a_z;

	// -z_0 < z < H_SSDboxlid ?
	bool inner_side = (-z_0 <= z) && (z <= H_SSDboxlid);

	return boxlid || inner_side;
}


// Does the trajectory reach the SSD?
bool reach_ssd (double x_M, double y_M, double z_M, double a_x, double a_y, double a_z, double* par){
	// The geometry of the BPM
	double z_0 = par[0];
	double R_MCP = par[1];
	double SSDholder_back = par[2];
	double W_SSDholder = par[3];
	double H_SSDholder = par[4];
	double SSDholder_front = par[5];
	double W_SSDboxlid = par[6];
	double R_SSD = par[7];
	double H_SSDboxlid = par[8];

	// The particle position is defined as
	// vec{P} = (x_M+t~*t_SSD*a_x, y_M+t~*t_SSD*a_y, z_M+t~*t_SSD*a_z)
	// The particle hits the surface z = -z_0 when t~ = 1
	double t_s = t_SSD(z_M,a_z,z_0);

	// Hit Component 1
	bool hits_hitsurface1 = HitSurface1(x_M,y_M,z_M,a_x,a_y,a_z,t_s,R_MCP,T_lid);

	// Hit Component 2
	bool hits_hitsurface2 = HitSurface2(x_M,y_M,z_M,a_x,a_y,a_z,t_s,SSDholder_back,W_SSDholder,H_SSDholder,R_MCP,SSDholder_front);

	// Hit Component 3
	bool hits_hitsurface3 = HitSurface3(x_M,y_M,z_M,a_x,a_y,a_z,t_s,SSDholder_front,x_0,W_SSDboxlid,R_MCP,R_SSD,H_SSDboxlid);

	// SSD Surface
	bool hits_SSD = HitsSSD(x_M,y_M,z_M,a_x,a_y,a_z,t_s,z_0,x_0);

	if (hits_hitsurface1){
		return false;
	}else if (hits_hitsurface2){
		return false;
	}else if (hits_hitsurface3){
		return false;
	}else if (hits_SSD){
		return true;
	}else{
		return false;
	}
}



int main(int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
	c1->Divide(1,2);

//	int N_Fr = 10000000; // 10^7 per sec.
	int N_Fr = 100; // for testing
        int N_Average = 60; // 1 min average
	cout << "Flying " << N_Fr << " alpha particles from the mesh." << endl;

	double z_l = 1.0; // Thickness of the SSD lid
	double R_SSD = 6.5; // Radius of the hole on the SSD lid
	double R_MCP = 14.0; // Radius of the MCP-IN
	double T_lid = 3.5; // Thickness of the MCP lid from the MCP surface

	// BPM geometry parameters
	double geometry[13];
	geometry[0] = 14.0; // R_MCP:Radius of MCP-IN
	geometry[1] = 25.0; // R_MCPlid:Outer radius of MCP lid
	geometry[2] = 6.5; // R_SSD:Radius of hole on SSD box
	geometry[3] = 2.5; // R_Am:Radius of hole on Am box
	geometry[4] = 30.8; // x_0:SSDholder surface - R_SSD center
	geometry[5] = 24.0; // z_0:MCP center - SSD box surface
	geometry[6] = 3.5; // T_lid:thickness of SSD box surface
	geometry[7] = 16.8; // SSDholder_front:MCP-IN - SSDholder-upstream
	geometry[8] = 13.8; // SSDholder_back:MCP-IN - SSDholder-downstream
	geometry[9] = 50.0; // W_SSDholder:Width of SSD holder
	geometry[10] = 64.0; // H_SSDholder:Height of SSD holder
	geometry[11] = 28.0; // W_SSDboxlid:Width of SSD box
	geometry[12] = 23.0; // H_SSDboxlid:Height of SSD box



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
	mcp->Draw("SAME F");

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


	TPad *p1 = new TPad("p1","p1",0.0,0.0,1.0,1.0,0,0,0);
	p1->Draw();
	p1->cd();
	
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
