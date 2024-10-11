#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


using namespace std;

/*
E 0 4 9
U GEV CM
P 1 0 11 0.0000000000000000e+00 0.0000000000000000e+00 -1.7999999992746666e+01 1.8000000000000000e+01 5.1099904266560590e-04 4
P 2 1 22 1.0579540327082353e-01 1.3003692159396107e+00 -2.7210330144669115e-01 2.2416038112714887e-01 -1.3137522960237802e+00 13
P 3 1 11 -1.0579540327082353e-01 -1.3003692159396107e+00 -1.7727896691299975e+01 1.7775839618872851e+01 5.4163501256537891e-04 1
P 4 0 2212 0.0000000000000000e+00 0.0000000000000000e+00 2.7499839935073504e+02 2.7500000000000000e+02 9.3827210000933003e-01 4
P 5 4 2212 0.0000000000000000e+00 0.0000000000000000e+00 2.7499839934916560e+02 2.7499999999843055e+02 9.3827210000933003e-01 2
V -3 0 [2,5]
P 6 -3 443 -3.5551862194740669e-01 1.4958487072197826e+00 1.1811206391896121e+01 1.2306919129834569e+01 3.0970478186295050e+00 2
P 7 -3 2212 4.6131402512181968e-01 -1.9547949246518961e-01 2.6291508940720183e+02 2.6291724100064704e+02 9.3827210000157535e-01 1
P 8 6 -11 1.1797602583060831e+00 1.4803344190948167e+00 7.3075640979819489e+00 7.5487560277617902e+00 5.1099902180814746e-04 1
P 9 6 11 -1.5352788802540114e+00 1.5514288127159673e-02 4.5036422939314944e+00 4.7581631020908288e+00 5.1099900790317480e-04 1
*/

void analyzeHEPMCEvents_MC(){

	int numHEPMCFiles = 1;	

	int elecE = 18;
	int hadronE = 275;

	TString inputFileString = "lager-vmp-00mrad.jpsi-18on275.4pi.disp-jpsi-00-electron.run00001-lumi1.hepmc";
	
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

	TH1D * qSquared_gammaStar = new TH1D("qSquared_gammaStar", "qSquared_gammaStar", 1000, 1.0, 100.0);
	TH1D * qSquared_ScatElectron = new TH1D("qSquared_ScatElectron", "qSquared_ScatElectron", 1000, 1.0, 100.0);
	TH2D * eta_vs_qSquared_reco = new TH2D("eta_vs_qSquared_reco", "eta_vs_qSquared_reco", 100, 1.0, 10.0, 500, -4.0, 1.0);
	TH1D * eta_scat_electron = new TH1D("eta_scat_electron", "eta_scat_electron", 100, -4.0, 4.0);
	TH1D * pt_scat_electron = new TH1D("pt_scat_electron", "pt_scat_electron", 100, 0.0, 4.0);
	TH1D * px_scat_electron = new TH1D("px_scat_electron", "px_scat_electron", 100, -4.0, 4.0);
	TH1D * py_scat_electron = new TH1D("py_scat_electron", "py_scat_electron", 100, -4.0, 4.0);
	TH1D * pz_scat_electron = new TH1D("pz_scat_electron", "pz_scat_electron", 100, 0.0, 20.0);

	TH1D * t_distribution = new TH1D("proton_pt_squared", "proton_pt_squared", 100, 0.0, 1.7);
	TH1D * pt_scat_proton = new TH1D("pt_scat_proton", "pt_scat_proton", 100, 0.0, 4.0);
    TH1D * px_scat_proton = new TH1D("px_scat_proton", "px_scat_proton", 100, -4.0, 4.0);
    TH1D * py_scat_proton = new TH1D("py_scat_proton", "py_scat_proton", 100, -4.0, 4.0);
    TH1D * pz_scat_proton = new TH1D("pz_scat_proton", "pz_scat_proton", 100, hadronE - hadronE*0.2, hadronE + hadronE*0.05);		
	
	TH1D * j_psi_invMass = new TH1D("j_psi_invMass", "j_psi_invMass", 100, 0.0, 6.0);


	for(int ihepmcFile = 1; ihepmcFile < numHEPMCFiles+1; ihepmcFile++){

			ifstream inputTextFile(inputFileString.Data());
			//ifstream inputTextFile("ab_output.hepmc.hepmc");
			//ifstream inputTextFile(Form("../../../../sznajder/eic_analysis_alex/job_24_05_01_0000/events_1.txt.HepMC3")); //, elecE, hadronE, ihepmcFile));
	
			//ifstream inputTextFile(Form("../../../../sznajder/eic_analysis/FINAL/%d_%d_plus/events_%d.txt.HepMC3",elecE, hadronE, ihepmcFile));
			
			if(eventCounter > 0 && eventCounter%1000 == 0){ cout << "1000 events processed" << endl; }

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

			while (getline(inputTextFile, line) ){
		
				istringstream ss(line);
				ss >> lineIdentifier;
			
				if(lineIdentifier == "E"){ 

					//E 0 9 15 @ -5.0868376913709906e-02 -2.5795220595253325e-03 7.7078711664698085e+00 -1.0678171802415379e+01				

					cout << "Event " << eventCounter << " start..." << endl;

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
				}
				

				if(lineIdentifier == "P"){ 

					ss >> particleIdx[particleCounter] >> motherIdx[particleCounter] >> pdgCode[particleCounter]
				   	>> px[particleCounter] >> py[particleCounter] >> pz[particleCounter] >> E[particleCounter] >> mass[particleCounter] >> status[particleCounter];


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
					}

					
					reco_JPsi = jPsiDaughter_one + jPsiDaughter_two;
					
					
				}

	
				if(particleCounter == numParticlesInEvent){
						
					//first, do something with the jPsi vector
					
					j_psi_invMass->Fill(reco_JPsi.M());
						
						
					for(int idx = 0; idx < particleCounter; idx++){
	            
						if(status[idx] == 1 && pdgCode[idx] == 2212){

							TLorentzVector scatProton(px[idx], py[idx], pz[idx], E[idx]);
							
							/*
							
							Exercise for the student!!
							
							*/
							
							double t = scatProton.Pt()*scatProton.Pt();
							t_distribution->Fill(t);

						}
						
						if(status[idx] == 13 && idx == internalPhotonIdx && E[idx] != -1){

							TLorentzVector gammaStar(px[idx], py[idx], pz[idx], E[idx]);
							
							double qSquared = -1*gammaStar.Mag2();
							qSquared_gammaStar->Fill(qSquared);

                        }	

						if(status[idx] == 1 && idx == scatteredElectronIdx && E[idx] != -1){
						
							TLorentzVector scatElectron(px[idx], py[idx], pz[idx], E[idx]);

							/*
							
							
							Exercise for the student!!
							
							
							*/

							double electronEnergy = elecE;

                            double theta = scatElectron.Theta();
                            theta = theta - TMath::Pi();
                            double qSquared = 2*electronEnergy*E[idx]*(1 - TMath::Cos(theta));

							qSquared_ScatElectron->Fill(qSquared);

							eta_vs_qSquared_reco->Fill(qSquared, scatElectron.Eta());

							eta_scat_electron->Fill(scatElectron.Eta());
    						pt_scat_electron->Fill(scatElectron.Pt());
    						px_scat_electron->Fill(scatElectron.Px());
    						py_scat_electron->Fill(scatElectron.Py());
    						pz_scat_electron->Fill(scatElectron.Pz());


						}
						
					}
			
				
					internalPhotonIdx = -999;
					scatteredElectronIdx = -999;
					beamElectronIdx = -999;
					beamProtonIdx = -999;
				
					eventCounter++; particleCounter = 0; numParticlesInEvent = -999;
					//if(eventCounter == 50000){break;}
				}


			} //while loop to read event

		inputTextFile.close();

	}//loop for number of input HEMPC files

	cout << "total events = " << eventCounter << endl;


	TCanvas * canvas1 = new TCanvas("can1", "can1", 1600, 800);
	canvas1->Divide(2,1);

	canvas1->cd(1)->SetLogy();
	canvas1->cd(1)->SetLogx();

	qSquared_gammaStar->SetLineColor(kBlack);
	qSquared_ScatElectron->SetLineColor(kRed);

	qSquared_gammaStar->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    qSquared_ScatElectron->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");

	qSquared_gammaStar->SetStats(0);
	qSquared_ScatElectron->SetStats(0);
	qSquared_gammaStar->Draw();
	//qSquared_ScatElectron->Draw("SAME");

	canvas1->cd(2);
    canvas1->cd(2)->SetLogz();
	canvas1->cd(2)->SetLeftMargin(0.15);
	canvas1->cd(2)->SetRightMargin(0.15);

    //eta_vs_qSquared_reco->SetStats(0);
	
	//eta_vs_qSquared_reco->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    //eta_vs_qSquared_reco->GetYaxis()->SetTitle("Scat. Electron Pseudorapidity, #eta");
    //eta_vs_qSquared_reco->Draw("COLZ");
	
	j_psi_invMass->Draw();

	
	TFile * outputFile = new TFile(Form("DVCS_analysis_%s_MC_%dx%d_GeV_5_29_2024.root",fileType.Data(), elecE, hadronE), "RECREATE");

	qSquared_gammaStar->Write();
	qSquared_ScatElectron->Write();
	eta_vs_qSquared_reco->Write();
	t_distribution->Write();


	eta_scat_electron->Write(); 
    pt_scat_electron->Write();
    px_scat_electron->Write();
    py_scat_electron->Write();
    pz_scat_electron->Write(); 

    pt_scat_proton->Write(); 
    px_scat_proton->Write(); 
    py_scat_proton->Write(); 
    pz_scat_proton->Write(); 
	
	j_psi_invMass->Write();

	return;

}

