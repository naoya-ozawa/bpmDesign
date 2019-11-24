void sigma_dep(){

	Double_t target_stdev[6] = {0.5,1.0,1.5,2.0,2.5,3.0}; // mm
	Double_t target_stdev_err[6] = {0.,0.,0.,0.,0.,0.};

	Double_t detect_eff[6] = {0.52,0.51,0.49,0.47,0.43,0.42}; // %
	Double_t detect_eff_err[6] = {0.02,0.02,0.02,0.02,0.02,0.02};

	TGraphErrors *ge = new TGraphErrors(6,target_stdev,detect_eff,target_stdev_err,detect_eff_err);
	ge->SetTitle("SSD Detection Efficiency Depending on the Beam StDev at the Target Surface;StDev at Target (mm);Simulated Detection Efficiency (%)");
	ge->Draw("ALP");
}
