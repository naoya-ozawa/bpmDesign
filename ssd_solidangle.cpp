#include <iostream>
#include <random>

using namespace std;


// Detection efficiency of the SSD
double efficiency(int N_Fr, int N_Average, double z_l, double R_SSD, double M_H, double M_W, double z_0, double x_0){

	double detection = 0.0;

	double y_M_ll = -M_W/2.0;
	double y_M_ul = M_W/2.0;
	double z_M_ll = z_0 - (M_H/2.0);
	double z_M_ul = z_0 + (M_H/2.0);

	random_device rnd;
	default_random_engine engine(rnd());

	uniform_real_distribution<> randhalf(0.,1.);
	uniform_real_distribution<> randall(-1.,1.);

	uniform_real_distribution<> randy(y_M_ll,y_M_ul);
	uniform_real_distribution<> randz(z_M_ll,z_M_ul);

	for (int j = 0; j < N_Average; j++){

		int N_Detected = 0;

		for (int i = 0; i < N_Fr; ++i){

			// Define Fr on MESH
			double y_M = randy(engine);
			double z_M = randz(engine);

			// Define alpha emitted direction
			double a_x = 0.0;
			while (a_x == 0.0){
				a_x = randhalf(engine);
			}
			double a_y = randall(engine);
			double a_z = randall(engine);
		

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





int main(){


	int N_Fr = 10000;
	int N_Average = 1000;

	cout << "Testing with " << N_Fr << " Fr ions on the mesh." << endl;
	cout << "Taking average of " << N_Average << " trials." << endl;

	double z_l = 1.0;
	double R_SSD = 6.5;
	double M_H = 70.0;
	double M_W = 70.0;
	
	double z_0 = 26.0;
	double x_0 = 26.0;

//	cout << "Geometry z_0 = " << z_0 << " mm, x_0 = " << x_0 << " mm." << endl;

	cout << "For z_0 = " << z_0 << " mm, x_0 = " << x_0 << " mm: ";

	double eff = efficiency(N_Fr, N_Average, z_l, R_SSD, M_H, M_W, z_0, x_0);

	cout << eff << "% detection." << endl;


	return 0;
}
