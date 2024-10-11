#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


using namespace std;



void analyzeHEPMCEvents_MC(){

	int numHEPMCFiles = 1;	

	int elecE = 18;
	int hadronE = 275;
	
	TString inputFileString = "lager-vmp-00mrad.jpsi-18on275.4pi.disp-jpsi-00-electron.run00001-lumi1.hepmc";
	
	TString dateString = "10_11_2024";
	
	TString fileType = "DEMP_JPsi";

	string line;

	string a1, a2, a3, a4, a5;

	string firstStr;
    //for particle listing
	int particleIdx[300];
	int motherIdx[300];
	int pdgCode[300];
	double px[300], py[300], pz[300], E[300], mass[300]; int status[300];

	string lineIdentifier = "";
	double vtx_x, vtx_y, vtx_z;
	

	int evLineCounter = 0;
	int particleCounter = 0;
	int numParticlesInEvent = -999;
	int eventCounter = 0;
	int globalLineCounter = 0;
	int lundIdx = 1;

	bool foundInternalPhoton = false;
    bool foundScatteredElectron = false;
		
	int internalPhotonIdx = -999;
	int scatteredElectronIdx = -999;
	int beamElectronIdx = -999;
	int beamProtonIdx = -999;

	TH1D * hist_qSquared_gammaStar    = new TH1D("qSquared_gammaStar", "qSquared_gammaStar", 1000, 1.0, 10.0);
	TH1D * hist_qSquared_ScatElectron = new TH1D("qSquared_ScatElectron", "qSquared_ScatElectron", 1000, 1.0, 100.0);
	TH2D * hist_eta_vs_qSquared_reco  = new TH2D("eta_vs_qSquared_reco", "eta_vs_qSquared_reco", 100, 1.0, 10.0, 500, -4.0, 1.0);
	TH1D * hist_eta_scat_electron     = new TH1D("eta_scat_electron", "eta_scat_electron", 100, -4.0, 4.0);
	TH1D * hist_pt_scat_electron      = new TH1D("pt_scat_electron", "pt_scat_electron", 100, 0.0, 4.0);
	TH1D * hist_px_scat_electron      = new TH1D("px_scat_electron", "px_scat_electron", 100, -4.0, 4.0);
	TH1D * hist_py_scat_electron      = new TH1D("py_scat_electron", "py_scat_electron", 100, -4.0, 4.0);
	TH1D * hist_pz_scat_electron      = new TH1D("pz_scat_electron", "pz_scat_electron", 100, 0.0, 20.0);

	TH1D * hist_t_distribution        = new TH1D("proton_pt_squared", "proton_pt_squared", 100, 0.0, 1.7);
	TH1D * hist_pt_scat_proton        = new TH1D("pt_scat_proton", "pt_scat_proton", 100, 0.0, 4.0);
    TH1D * hist_px_scat_proton        = new TH1D("px_scat_proton", "px_scat_proton", 100, -4.0, 4.0);
    TH1D * hist_py_scat_proton        = new TH1D("py_scat_proton", "py_scat_proton", 100, -4.0, 4.0);
    TH1D * hist_pz_scat_proton        = new TH1D("pz_scat_proton", "pz_scat_proton", 100, hadronE - hadronE*0.2, hadronE + hadronE*0.05);		
	
	TH1D * hist_j_psi_invMass = new TH1D("j_psi_invMass", "j_psi_invMass", 100, 0.0, 6.0);


	for(int ihepmcFile = 1; ihepmcFile < numHEPMCFiles+1; ihepmcFile++){

		ifstream inputTextFile(inputFileString.Data());
		
		cout << "processing HEPMC file: " << Form("%s", inputFileString.Data()) << endl;
				
		particleCounter = 0;
        numParticlesInEvent = -999;
		lundIdx = 1;
			
		TLorentzVector scatProton(0.,0.,0.,0.);
		TLorentzVector scatElectron(0.,0.,0.,0.);
		TLorentzVector jPsiDaughter_one(0.,0.,0.,0.);
		TLorentzVector jPsiDaughter_two(0.,0.,0.,0.);
		TLorentzVector reco_JPsi(0.,0.,0.,0.);

		TLorentzVector beamElectron(0.,0.,0.,0.);
		TLorentzVector beamProton(0.,0.,0.,0.);		
		TLorentzVector gammaStar(0.,0.,0.,0.);	

		while (getline(inputTextFile, line) ){
		
			istringstream ss(line);
			ss >> lineIdentifier;
			
			if(lineIdentifier == "E"){ 

				cout << "Processing event " << eventCounter << endl; 

				ss >> a1 >> a2 >> numParticlesInEvent >> a4;
                ss >> vtx_x >> vtx_y >> vtx_z;

                vtx_x = 0.0;
                vtx_y = 0.0;
                vtx_z = 0.0;
				
				scatProton.SetPxPyPzE(0.,0.,0.,0.);
				scatElectron.SetPxPyPzE(0.,0.,0.,0.);
				jPsiDaughter_one.SetPxPyPzE(0.,0.,0.,0.);
				jPsiDaughter_two.SetPxPyPzE(0.,0.,0.,0.);
				reco_JPsi.SetPxPyPzE(0.,0.,0.,0.);
				beamElectron.SetPxPyPzE(0.,0.,0.,0.);
				beamProton.SetPxPyPzE(0.,0.,0.,0.);	
				gammaStar.SetPxPyPzE(0.,0.,0.,0.);	
			}
				
			if(lineIdentifier == "P"){ 

				ss >> particleIdx[particleCounter] >> motherIdx[particleCounter] >> pdgCode[particleCounter]
				>> px[particleCounter] >> py[particleCounter] >> pz[particleCounter] >> E[particleCounter] >> mass[particleCounter] >> status[particleCounter];

				//not really being used right now, just for later use so you can see how it works.
				if(status[particleCounter] == 4 && pdgCode[particleCounter] == 11){beamElectronIdx = particleCounter;}	
				if(status[particleCounter] == 4 && pdgCode[particleCounter] == 2212){beamProtonIdx = particleCounter;}	
				if(status[particleCounter] == 1 && pdgCode[particleCounter] == 11){scatteredElectronIdx = particleCounter;}

				particleCounter++;
			}
			
			if(particleCounter == numParticlesInEvent){
				
				for(int idx = 0; idx < particleCounter; idx++){				
					if(status[idx] == 1 && pdgCode[idx] == 2212){

                        scatProton.SetPxPyPzE(px[idx], py[idx], pz[idx], E[idx]);
					}
					if(status[idx] == 1 && pdgCode[idx] == 11 && particleIdx[idx] == 3){ // Scattered electron -- we are cheating, for now!

                        scatElectron.SetPxPyPzE(px[idx], py[idx], pz[idx], E[idx]);

                    }
					if(status[idx] == 1 && TMath::Abs(pdgCode[idx]) == 11 && particleIdx[idx] == 8){ // jPsi -> ee daughter number 1 -- we are cheating, for now!

                        jPsiDaughter_one.SetPxPyPzE(px[idx], py[idx], pz[idx], E[idx]);

                    }
					if(status[idx] == 1 && TMath::Abs(pdgCode[idx]) == 11 && particleIdx[idx] == 9){ // jPsi -> ee daughter number 1 -- we are cheating, for now!

                        jPsiDaughter_two.SetPxPyPzE(px[idx], py[idx], pz[idx], E[idx]);

                    }
					if(status[idx] == 4 && pdgCode[idx] == 2212){

                        beamProton.SetPxPyPzE(px[idx], py[idx], pz[idx], E[idx]);

                    }
					if(status[idx] == 4 && pdgCode[idx] == 11){

                        beamElectron.SetPxPyPzE(px[idx], py[idx], pz[idx], E[idx]);

                    }
					if(status[idx] == 13 && pdgCode[idx] == 22){

                        gammaStar.SetPxPyPzE(px[idx], py[idx], pz[idx], E[idx]);

                    }
				}

				
				double qSquared_gammaStar = -1*gammaStar.Mag2();
				hist_qSquared_gammaStar->Fill(qSquared_gammaStar);
				
				/*
					Exercise 1: calculate Q2 using the scatter electron information, and plot it.
				
					Exercise 2: calculate the momentum transfer, t, and plot it
				
					Exercise 3: reconstruct the JPsi vector, and plot the invariant mass.
				
					Exercise 4: What looks strange about the J/Psi invariant mass distribution? Why?
				*/
				
				internalPhotonIdx = -999;
				scatteredElectronIdx = -999;
				beamElectronIdx = -999;
				beamProtonIdx = -999;
			
				eventCounter++; particleCounter = 0; numParticlesInEvent = -999;
				
			}// end loop over event particles stored in memory

		}//while loop to read event

		inputTextFile.close();
	}

	cout << "total events = " << eventCounter << endl;


	TCanvas * canvas1 = new TCanvas("can1", "can1", 1600, 800);
	canvas1->Divide(2,1);

	canvas1->cd(1)->SetLogy();
	canvas1->cd(1)->SetLogx();

	hist_qSquared_gammaStar->SetLineColor(kBlack);
	hist_qSquared_ScatElectron->SetLineColor(kRed);

	hist_qSquared_gammaStar->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    hist_qSquared_ScatElectron->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");

	hist_qSquared_gammaStar->SetStats(0);
	hist_qSquared_ScatElectron->SetStats(0);
	hist_qSquared_gammaStar->Draw();
	//qSquared_ScatElectron->Draw("SAME");

	canvas1->cd(2);
    canvas1->cd(2)->SetLogz();
	canvas1->cd(2)->SetLeftMargin(0.15);
	canvas1->cd(2)->SetRightMargin(0.15);

	hist_j_psi_invMass->GetXaxis()->SetTitle("J/#Psi #rightarrow e^{+}e^{-} invariant mass [GeV/c^{2}]");
	
	hist_j_psi_invMass->Draw();

	
	TFile * outputFile = new TFile(Form("Yug_analysis_%s_MC_%dx%d_GeV_%s.root",fileType.Data(), elecE, hadronE, dateString.Data()), "RECREATE");

	hist_qSquared_gammaStar->Write();
	hist_qSquared_ScatElectron->Write();
	hist_eta_vs_qSquared_reco->Write();
	hist_t_distribution->Write();


	hist_eta_scat_electron->Write(); 
    hist_pt_scat_electron->Write();
    hist_px_scat_electron->Write();
    hist_py_scat_electron->Write();
    hist_pz_scat_electron->Write(); 

   	hist_pt_scat_proton->Write(); 
    hist_px_scat_proton->Write(); 
    hist_py_scat_proton->Write(); 
    hist_pz_scat_proton->Write(); 
	
	hist_j_psi_invMass->Write();

	return;

}

