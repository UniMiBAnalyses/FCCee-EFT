//
// likelihood scan for one EFT Wilson coefficient C
//
//   root -l likelihood.cxx
//   root -l 'likelihood.cxx("input.root", "histo_SM", "histo_cP1", "histo_cM1", -1, 1, 201)'
//
// Model, bin by bin (as in AnalyticAnomalousCoupling):
//
//    N(C) = SM + C * L + C^2 * Q
//
// with L and Q obtained from the histograms at C = +1 and C = -1:
//
//    N(+1) = SM + L + Q
//    N(-1) = SM - L + Q
//    -->  L = ( N(+1) - N(-1) ) / 2
//    -->  Q = ( N(+1) + N(-1) ) / 2 - SM
//
// Data = SM histogram (Asimov)
//
// Likelihood = product of Poisson, normalised to the saturated model:
//
//    -2 lnL(C) = 2 * sum_i [ N_i(C) - n_i + n_i * ln( n_i / N_i(C) ) ]
//


void likelihood(TString name_file = "input.root",
                TString name_SM   = "histo_SM",
                TString name_cP1  = "histo_cP1",
                TString name_cM1  = "histo_cM1",
                float c_min = -1.,
                float c_max =  1.,
                int n_points = 201
               ) {

  gStyle->SetOptStat(0);

  //---- input

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

  int n_bins = h_SM->GetNbinsX();
  if (h_cP1->GetNbinsX() != n_bins || h_cM1->GetNbinsX() != n_bins) {
    std::cout << " different number of bins! " << std::endl;
    return;
  }


  //---- algebra: linear and quadratic components

  TH1* h_lin = (TH1*) h_cP1->Clone("h_lin");
  h_lin->Add(h_cM1, -1.);
  h_lin->Scale(0.5);

  TH1* h_quad = (TH1*) h_cP1->Clone("h_quad");
  h_quad->Add(h_cM1);
  h_quad->Scale(0.5);
  h_quad->Add(h_SM, -1.);

  std::cout << " SM    = " << h_SM->Integral()   << std::endl;
  std::cout << " C=+1  = " << h_cP1->Integral()  << std::endl;
  std::cout << " C=-1  = " << h_cM1->Integral()  << std::endl;
  std::cout << " lin   = " << h_lin->Integral()  << std::endl;
  std::cout << " quad  = " << h_quad->Integral() << std::endl;

  for (int iBin = 1; iBin <= n_bins; iBin++) {
    if (h_quad->GetBinContent(iBin) < 0) {
      std::cout << " WARNING: negative quadratic term in bin " << iBin << " = " << h_quad->GetBinContent(iBin) << std::endl;
    }
  }


  //---- -2 lnL
  //
  // Where is the Poisson? And the factorial?
  //
  //    Pois(n|mu) = mu^n * exp(-mu) / n!
  //
  // We use the ratio to the saturated model, i.e. the same Poisson evaluated at mu = n:
  //
  //    Pois(n|mu) / Pois(n|n) = [ mu^n exp(-mu) / n! ] / [ n^n exp(-n) / n! ]
  //                           = (mu/n)^n * exp(-(mu-n))
  //
  //    --> the n! cancels (same n = data above and below)
  //
  //    -2 ln [ Pois(n|mu) / Pois(n|n) ] = 2 * [ mu - n + n * ln(n/mu) ]
  //
  // and the -2lnL of the product of Poisson is the sum over bins of this.
  //
  // The denominator, Pois(n|n), depends only on data, not on C:
  // it is a constant shift of -2lnL and it cancels in -2DlnL = -2lnL(C) - min(-2lnL).
  // --> scan, best fit and intervals are identical to the explicit Poisson version (commented below).
  //
  // Why this form:
  //   - no n! --> fine also with non-integer n (Asimov data = SM yields, e.g. 670.32)
  //   - sum of logs, no product of small numbers --> no numerical underflow with many bins / large n
  //   - -2lnL = 0 when mu = n
  //

  auto minus2LogL = [&] (double c) {
    double value = 0.;
    for (int iBin = 1; iBin <= n_bins; iBin++) {
      double n  = h_SM->GetBinContent(iBin);
      double mu = h_SM->GetBinContent(iBin) + c * h_lin->GetBinContent(iBin) + c * c * h_quad->GetBinContent(iBin);
      if (mu <= 0) {
        if (n > 0) return 1.e9;   // expectation <= 0 with observed events: not allowed
        continue;
      }
      value += 2. * (mu - n);
      if (n > 0) value += 2. * n * log (n / mu);

      //
      // explicit Poisson version, equivalent for -2DlnL:
      //   TMath::Poisson(n, mu) = mu^n exp(-mu) / Gamma(n+1), i.e. n! --> Gamma(n+1), works for non-integer n
      //   NB: -2lnL is shifted by a constant (min not at 0 any more), and for large n the Poisson
      //       can become very small: the log of the product is safe only as sum of logs, as done here
      //
//       value += -2. * log (TMath::Poisson(n, mu));

    }
    return value;
  };


  //---- scan

  std::vector<double> v_c;
  std::vector<double> v_m2lnL;

  for (int iPoint = 0; iPoint < n_points; iPoint++) {
    double c = c_min + (c_max - c_min) * iPoint / (n_points - 1);
    v_c.push_back(c);
    v_m2lnL.push_back(minus2LogL (c));
  }


  //---- minimum of -2lnL in [c_min, c_max]

  TF1* f_m2lnL = new TF1 ("f_m2lnL", [&] (double* x, double* p) { return minus2LogL (x[0]); }, c_min, c_max, 0);

  double c_best    = f_m2lnL->GetMinimumX (c_min, c_max);
  double min_m2lnL = minus2LogL (c_best);

  std::cout << " min -2lnL = " << min_m2lnL << " at C = " << c_best << std::endl;


  //---- shift: -2DlnL = -2lnL - min(-2lnL), so the minimum of the curve is at 0

  TGraph* gr_scan = new TGraph ();
  for (int iPoint = 0; iPoint < n_points; iPoint++) {
    gr_scan->SetPoint(iPoint, v_c[iPoint], v_m2lnL[iPoint] - min_m2lnL);
  }
  gr_scan->SetPoint(n_points, c_best, 0.);   // add the minimum itself
  gr_scan->Sort();

  std::cout << " best fit C = " << c_best << std::endl;


  //---- intervals (Wilks, 1 dof): 68% CL -> 1.00, 95% CL -> 3.84

  std::vector<double> thresholds = {1.00, 3.84};
  std::vector<TString> names_CL = {"68% CL", "95% CL"};

  for (int iThreshold = 0; iThreshold < thresholds.size(); iThreshold++) {
    std::cout << " " << names_CL[iThreshold] << " (-2DlnL = " << thresholds[iThreshold] << ") crossings: ";
    int n_crossings = 0;
    for (int iPoint = 0; iPoint < n_points-1; iPoint++) {
      double y1 = v_m2lnL[iPoint]   - min_m2lnL - thresholds[iThreshold];
      double y2 = v_m2lnL[iPoint+1] - min_m2lnL - thresholds[iThreshold];
      if (y1 * y2 < 0) {
        double c_cross = v_c[iPoint] + (v_c[iPoint+1] - v_c[iPoint]) * (-y1) / (y2 - y1);
        std::cout << "  " << c_cross;
        n_crossings++;
      }
    }
    if (n_crossings == 0) std::cout << "  none in [" << c_min << ", " << c_max << "] --> enlarge the range";
    std::cout << std::endl;
  }


  //---- plot: scan

  TCanvas* cc_scan = new TCanvas ("cc_scan", "scan", 800, 600);
  gr_scan->SetLineColor(kBlue);
  gr_scan->SetLineWidth(2);
  gr_scan->SetTitle("");
  gr_scan->GetXaxis()->SetTitle("C");
  gr_scan->GetYaxis()->SetTitle("-2 #Delta lnL");
  gr_scan->Draw("AL");

  for (auto threshold : thresholds) {
    TLine* line = new TLine (c_min, threshold, c_max, threshold);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw();
  }


  //---- plot: inputs

  TCanvas* cc_inputs = new TCanvas ("cc_inputs", "inputs", 800, 600);
  h_SM->SetLineColor(kBlack);
  h_cP1->SetLineColor(kRed);
  h_cM1->SetLineColor(kBlue);
  h_SM->SetLineWidth(2);
  h_cP1->SetLineWidth(2);
  h_cM1->SetLineWidth(2);

  h_cP1->Draw("histo");
  h_cM1->Draw("histo same");
  h_SM->Draw("histo same");

  TLegend* leg = new TLegend (0.70, 0.70, 0.88, 0.88);
  leg->AddEntry(h_SM,  "SM",   "l");
  leg->AddEntry(h_cP1, "C=+1", "l");
  leg->AddEntry(h_cM1, "C=-1", "l");
  leg->Draw();

  std::cout << " I'm done " << std::endl;

}
