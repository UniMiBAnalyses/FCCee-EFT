//
// draw the 3 input histograms: SM, C=+1, C=-1
//
//   root -l drawInput.cxx
//   root -l 'drawInput.cxx("myfile.root", "histo_SM", "histo_cP1", "histo_cM1")'
//


void drawInput(TString name_file = "input.root",
               TString name_SM   = "histo_SM",
               TString name_cP1  = "histo_cP1",
               TString name_cM1  = "histo_cM1"
              ) {

  gStyle->SetOptStat(0);

  TFile* f_input = new TFile (name_file, "READ");
  if (!f_input || f_input->IsZombie()) {
    std::cout << " file not found: " << name_file << std::endl;
    return;
  }

  TH1* h_SM  = (TH1*) f_input->Get(name_SM);
  TH1* h_cP1 = (TH1*) f_input->Get(name_cP1);
  TH1* h_cM1 = (TH1*) f_input->Get(name_cM1);

  if (!h_SM || !h_cP1 || !h_cM1) {
    std::cout << " histogram not found: " << name_SM << " " << name_cP1 << " " << name_cM1 << std::endl;
    return;
  }

  //---- table

  std::cout << " bin      SM     C=+1     C=-1" << std::endl;
  for (int iBin = 1; iBin <= h_SM->GetNbinsX(); iBin++) {
    printf (" %3d %8.2f %8.2f %8.2f\n", iBin, h_SM->GetBinContent(iBin), h_cP1->GetBinContent(iBin), h_cM1->GetBinContent(iBin));
  }

  //---- plot

  TCanvas* cc = new TCanvas ("cc", "", 800, 600);
  cc->SetLogy();

  h_SM ->SetLineColor(kBlack);
  h_cP1->SetLineColor(kRed);
  h_cM1->SetLineColor(kBlue);
  h_SM ->SetLineWidth(2);
  h_cP1->SetLineWidth(2);
  h_cM1->SetLineWidth(2);
  h_cM1->SetLineStyle(2);

  h_cP1->SetTitle("");
  h_cP1->GetXaxis()->SetTitle("bin");
  h_cP1->GetYaxis()->SetTitle("events");
  h_cP1->Draw("histo");
  h_cM1->Draw("histo same");
  h_SM ->Draw("histo same");

  TLegend* leg = new TLegend (0.65, 0.70, 0.88, 0.88);
  leg->AddEntry(h_SM,  "SM",   "l");
  leg->AddEntry(h_cP1, "C=+1", "l");
  leg->AddEntry(h_cM1, "C=-1", "l");
  leg->Draw();

  cc->SaveAs("input.png");

}
