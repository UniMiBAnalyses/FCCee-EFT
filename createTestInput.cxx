//
// create a test input for likelihood.cxx
//
//   root -l -b -q createTestInput.cxx
//
// 3 histograms: SM, C=+1, C=-1
// built from a falling SM spectrum with linear and quadratic EFT terms growing with the bin index:
//
//    N(C) = SM + C * L + C^2 * Q
//


void createTestInput(TString name_file = "input.root") {

  TFile* f_output = new TFile (name_file, "RECREATE");

  int n_bins = 10;
  TH1F* h_SM  = new TH1F ("histo_SM",  "SM",   n_bins, 0, 10);
  TH1F* h_cP1 = new TH1F ("histo_cP1", "C=+1", n_bins, 0, 10);
  TH1F* h_cM1 = new TH1F ("histo_cM1", "C=-1", n_bins, 0, 10);

  for (int iBin = 1; iBin <= n_bins; iBin++) {
    double sm   = 1000. * exp (-0.4 * iBin);
    double lin  = 0.02  * iBin * sm;
    double quad = 0.001 * iBin * iBin * sm;

    h_SM ->SetBinContent(iBin, sm);
    h_cP1->SetBinContent(iBin, sm + lin + quad);
    h_cM1->SetBinContent(iBin, sm - lin + quad);
  }

  std::cout << " SM    = " << h_SM->Integral()  << std::endl;
  std::cout << " C=+1  = " << h_cP1->Integral() << std::endl;
  std::cout << " C=-1  = " << h_cM1->Integral() << std::endl;

  f_output->Write();
  f_output->Close();

  std::cout << " written: " << name_file << std::endl;

}
