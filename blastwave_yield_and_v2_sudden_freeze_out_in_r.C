//
// Implementation of the blast wave v2 formula for an elliptic freeze-out hyper-surface
// Klaus Reygers, October 2019
//

class blastwave_yield_and_v2 {

    // function for v2 numerator and denominator
    // fos = freeze-out surface
    // fos1: freeze-out time tf in the lab system: tf = sqrt(tau^2 + z^2) (no radial dependence)
    static double v2_fos1_numerator(const double *x, const double *p);
    static double v2_fos1_denominator(const double *x, const double *p);

    // wrapper functions for v2 numerator and denominator
    ROOT::Math::WrappedParamFunction<> w_v2_fos1_num;
    ROOT::Math::WrappedParamFunction<> w_v2_fos1_den;

    ROOT::Math::AdaptiveIntegratorMultiDim ig;

  public:
    blastwave_yield_and_v2()
        : w_v2_fos1_num(&blastwave_yield_and_v2::v2_fos1_numerator, 2, 7),
          w_v2_fos1_den(&blastwave_yield_and_v2::v2_fos1_denominator, 2, 7) {}

    void calc_blastwave_yield_and_v2_fos1(const double &pt, const double &m, const double &T, const double &rho0,
                                          const double &rho2, const double &RxOverRy, double &inv_yield, double &v2);

    ClassDef(blastwave_yield_and_v2, 1)
};

// numerator of the blastwave v2 formula
double blastwave_yield_and_v2::v2_fos1_numerator(const double *x, const double *p) {

    // integration variables
    double rHat = x[0];
    double PhiHat = x[1];

    // parameters
    double pt = p[0];
    double m = p[1];
    double T = p[2];
    double rho0 = p[3];
    double rho2 = p[4];
    double RxOverRy = p[5];
    double A = p[6]; // arbitrary factor, can be adjusted to improve numerical stability

    double PhiB = TMath::ATan(RxOverRy * TMath::Tan(PhiHat)); // boost angle
    double mt = TMath::Sqrt(m * m + pt * pt);
    double rho = rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)); // transverse rapidity
    double xip = pt * TMath::SinH(rho) / T;
    double xim = mt * TMath::CosH(rho) / T;

    return A * rHat * TMath::BesselI(2, xip) * TMath::BesselK(1, xim) * TMath::Cos(2 * PhiB);
}

// denominator of the blastwave v2 formula
double blastwave_yield_and_v2::v2_fos1_denominator(const double *x, const double *p) {

    // integration variables
    double rHat = x[0];
    double PhiHat = x[1];

    // parameters
    double pt = p[0];
    double m = p[1];
    double T = p[2];
    double rho0 = p[3];
    double rho2 = p[4];
    double RxOverRy = p[5];
    double A = p[6]; // arbitrary factor, can be adjusted to improve numerical stability

    double PhiB = TMath::ATan(RxOverRy * TMath::Tan(PhiHat)); // boost angle
    double mt = TMath::Sqrt(m * m + pt * pt);
    double rho = rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)); // transverse rapidity
    double xip = pt * TMath::SinH(rho) / T;
    double xim = mt * TMath::CosH(rho) / T;

    return A * rHat * TMath::BesselI(0, xip) * TMath::BesselK(1, xim);
}

void blastwave_yield_and_v2::calc_blastwave_yield_and_v2_fos1(const double &pt, const double &m, const double &T,
                                                              const double &rho0, const double &rho2,
                                                              const double &RxOverRy, double &inv_yield, double &v2) {

    // blast wave parameters:
    // the last number (par[6]) is an arbitrary normalization which we will adjust
    // in order to have a good numerical stability for different masses and pt
    double pars[7] = {pt, m, T, rho0, rho2, RxOverRy, 1.};

    // determine scale factor which ensures good numerical stability
    const double xsf[2] = {1., 0.};
    double sf = v2_fos1_denominator(xsf, pars);
    pars[6] = 1. / sf;

    // set parameters
    w_v2_fos1_num.SetParameters(pars);
    w_v2_fos1_den.SetParameters(pars);

    // define integrator
    // ROOT::Math::AdaptiveIntegratorMultiDim ig;
    ig.SetRelTolerance(1e-6);

    // integration range
    double xmin[2] = {0., 0.};
    double xmax[2] = {1., 2. * TMath::Pi()};

    // integrate
    ig.SetFunction(w_v2_fos1_num);
    double v2_num = ig.Integral(xmin, xmax);
    // if (ig_num.Status() != 0) cout << ig_num.Status() << endl;

    ig.SetFunction(w_v2_fos1_den);
    double v2_den = ig.Integral(xmin, xmax);
    // if (ig_den.Status() != 0) cout << ig_den.Status() << endl;

    if (v2_den != 0) {
        v2 = v2_num / v2_den;
    } else {
        cout << "WARNING: v2 denominator zero!!!" << endl;
    }

    inv_yield = sf * TMath::Sqrt(m * m + pt * pt) * v2_den;
}

// main function:
// show usage of blast wave v2 function
void blastwave_yield_and_v2_sudden_freeze_out_in_r() {

    blastwave_yield_and_v2 bw;

    //
    // Plot v2 vs pt
    //
    const double T_BW = 0.10;
    const double Rho0_BW = 1.1;
    const double Rho2_BW = 0.065;
    const double RxOverRy_BW = 0.8;

    // const double m_BW = 9.460; // Upsilon
    // const double m_BW = 3.096; // J/psi
    // const double m_BW = 0.14; // pion
    const double m_BW = 0.938; // proton

    const double pt_min = 0;
    const double pt_max = 16.;
    const double dpt = 0.1;

    const int i_max = (pt_max - pt_min) / dpt + 1;
    TGraph giy(i_max);
    TGraph gv2(i_max);

    TStopwatch t;
    t.Start();

    for (int i = 0; i < i_max; ++i) {
        double pt_BW = pt_min + i * dpt + 0.2 * dpt;
        double inv_yield_BW = 0;
        double v2_BW = 0;
        bw.calc_blastwave_yield_and_v2_fos1(pt_BW, m_BW, T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
        giy.SetPoint(i, pt_BW, inv_yield_BW);
        gv2.SetPoint(i, pt_BW, v2_BW);
    }

    t.Stop();
    t.Print();

    // plot invariant Yield and v2 vs pt
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1");
    TH2F *fr = new TH2F("frame", "blast wave #it{v}_{2}", 1, pt_min, pt_max + 0.5, 1, 0., 0.8);
    fr->SetXTitle("p_{T} (GeV/#it{c})");
    fr->SetYTitle("#it{v}_{2}");
    fr->Draw();
    gv2.SetMarkerStyle(kFullCircle);
    gv2.SetMarkerSize(0.5);
    gv2.DrawClone("l");

    TCanvas *c2 = new TCanvas("c2");
    c2->SetLogy();
    double xmin, ymin, xmax, ymax;
    giy.ComputeRange(xmin, ymin, xmax, ymax);
    TH2F *fr2 = new TH2F("frame", "blast wave inv. yield", 1, pt_min, pt_max + 0.5, 1, ymin, ymax);
    fr2->SetXTitle("p_{T} (GeV/#it{c})");
    fr2->SetYTitle("inv. yield");
    fr2->Draw();
    giy.SetMarkerStyle(kFullCircle);
    giy.SetMarkerSize(0.5);
    giy.DrawClone("l");
}