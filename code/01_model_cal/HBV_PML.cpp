#include <Rcpp.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <vector>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
DataFrame HBV_PML_Rcpp(DataFrame pars, 
                       DataFrame input, 
                       int ifPML = 1,
                       double area = 1.0,
                       int ifRoute = 1) {
  
  double parAlpha = as<NumericVector>(pars["Alpha"])[0], parThelta = as<NumericVector>(pars["Thelta"])[0];
  double par_m = as<NumericVector>(pars["m"])[0], parAm_25 = as<NumericVector>(pars["Am_25"])[0];
  double parD0 = as<NumericVector>(pars["D0"])[0], parkQ = as<NumericVector>(pars["kQ"])[0];
  double parkA = as<NumericVector>(pars["kA"])[0], parS_sls = as<NumericVector>(pars["S_sls"])[0];
  double parfER0 = as<NumericVector>(pars["fER0"])[0], parVPDmin = as<NumericVector>(pars["VPDmin"])[0];
  double parVPDmax = as<NumericVector>(pars["VPDmax"])[0], parhc = as<NumericVector>(pars["hc"])[0];
  double parTs = as<NumericVector>(pars["Ts"])[0], parCFMAX = as<NumericVector>(pars["CFMAX"])[0];
  double parCFR = as<NumericVector>(pars["CFR"])[0], parCWH = as<NumericVector>(pars["CWH"])[0]; 
  double parBETA = as<NumericVector>(pars["BETA"])[0], parLP = as<NumericVector>(pars["LP"])[0];
  double parFC = std::max(1e-10, (double)as<NumericVector>(pars["FC"])[0]);
  double parPERC = as<NumericVector>(pars["PERC"])[0], parK0 = as<NumericVector>(pars["K0"])[0];
  double parK1 = as<NumericVector>(pars["K1"])[0], parK2 = as<NumericVector>(pars["K2"])[0];
  double parUZL = as<NumericVector>(pars["UZL"])[0], parMAXBAS = as<NumericVector>(pars["MAXBAS"])[0];
  
  NumericVector P = input["P_MSWEP"], Temp = input["Temp"];
  NumericVector LAI, Pres, Wind, SWd, CO2, Rn, VPD, Eeq, PET_FAO;
  
  if (ifPML == 1) {
    LAI = input["LAI_GLASS"]; Pres = input["Pres"]; Wind = input["Wind"];
    SWd = input["SWd"]; CO2 = input["CO2"]; Rn = input["Rn"];
    VPD = input["VPD"]; Eeq = input["Eeq"];
  } else { PET_FAO = input["PET_FAO"]; }
  
  int n = P.size();
  NumericVector v(n, 0.0), vl(n, 0.0), SM(n, 0.0), UZ(n, 0.0), LZ(n, 0.0);
  NumericVector P_eff(n, 0.0), EA(n, 0.0), Q0(n, 0.0), Q1(n, 0.0), Qsim_mm(n, 0.0), Qsim_m3s(n, 0.0), QMOD(n, NA_REAL);
  NumericVector Ei_out(n, 0.0), Es_out(n, 0.0), Et_out(n, 0.0), Est_out(n, 0.0), EPT_out(n, NA_REAL);
  
  std::vector<double> P_minus_Ei_hist(n, 0.0), Es_eq_hist(n, 0.0);
  double prev_v = 0.0, prev_vl = 0.0, prev_SM = 0.0, prev_UZ = 0.0, prev_LZ = 0.0;
  
  for (int i = 0; i < n; ++i) {
    double temp = Temp[i], prec = P[i];
    double rain = (temp >= parTs) ? prec : 0.0, snow = (temp < parTs) ? prec : 0.0;
    double Ta = (temp > parTs) ? (temp - parTs) : 0.0, Tn = (temp < parTs) ? (parTs - temp) : 0.0;
    
    double melt = std::min(parCFMAX * Ta, prev_v);
    double rfz = std::min(parCFR * parCFMAX * Tn, prev_vl);
    
    v[i] = prev_v - melt + snow + rfz;
    vl[i] = prev_vl + melt + rain - rfz;
    
    if (vl[i] > parCWH * v[i]) {
      P_eff[i] = vl[i] - (parCWH * v[i]);
      vl[i] = parCWH * v[i];
    } else { P_eff[i] = 0.0; }
    
    double ept = 0.0, Ei = 0.0, Es = 0.0, Et = 0.0, Est = 0.0, Prcp = P_eff[i];
    
    if (ifPML == 1) {
      double tmean = Temp[i], Pa = Pres[i] / 1000.0, PAR_mol = (0.45 * SWd[i]) * 4.57;
      double f_VPD_gc = 1.0 / (1.0 + VPD[i] / parD0);
      double f_VPD = (VPD[i] > parVPDmax) ? 0.0 : ((VPD[i] >= parVPDmin) ? (parVPDmax - VPD[i]) / (parVPDmax - parVPDmin) : 1.0);
      double fT2 = std::min(std::exp(0.031 * (tmean - 25.0)) / (1.0 + std::exp(0.115 * (tmean - 41.0))), 1.0);
      
      double Am = parAm_25 * fT2;
      double P1 = Am * parAlpha * parThelta * PAR_mol, P2 = Am * parAlpha * PAR_mol;
      double P3 = Am * parThelta * CO2[i], P4 = parAlpha * PAR_mol * parThelta * CO2[i];
      
      double log_term = (P2 + P3 + P4) / (P2 + P3 * std::exp(parkQ * LAI[i]) + P4);
      double Ags = CO2[i] * P1 / (P2 * parkQ + P4 * parkQ) * (parkQ * LAI[i] + std::log(log_term));
      double Gc = std::max((par_m / CO2[i]) * (Ags * f_VPD) * 1.6 * f_VPD_gc * 1e-2 / (0.446 * (273.0 / (273.0 + tmean)) * (Pa / 101.3)), 1e-6);
      
      double Tou = std::exp(-parkA * LAI[i]), d = 0.667 * parhc, zom = 0.125 * parhc, zoh = 0.135 * zom;
      double Cp = 4.2 * 0.242, lamada = 2500.0 - 2.4 * tmean, gama = Cp * Pa / (0.622 * lamada);
      double slop = 4098.0 * 0.6108 * std::exp((17.27 * tmean) / (tmean + 237.3)) / std::pow(tmean + 237.3, 2.0);
      double epsilon = slop / gama;
      
      double uz = std::log(67.8 * 15.0 - 5.42) / 4.87 * Wind[i];
      double Ga = uz * 0.16 / (std::log((15.0 - d) / zom) * std::log((15.0 - d) / zoh));
      
      double LEcr = epsilon * Rn[i] * (1.0 - Tou) / (epsilon + 1.0 + Ga / Gc);
      double LEca = ((3.846 * 1000.0 * Pa / (tmean + 273.15)) * Cp * Ga * VPD[i] / gama) / (epsilon + 1.0 + Ga / Gc);
      double Ecr = LEcr / lamada * 86.4, Eca = LEca / lamada * 86.4;
      
      double fveg = 1.0 - std::exp(-LAI[i] / 5.0);
      double Pwet = -std::log(1.0 - parfER0) / parfER0 * (parS_sls * LAI[i]) / fveg;
      Ei = (Prcp < Pwet) ? (fveg * Prcp) : (fveg * Pwet + (parfER0 * fveg) * (Prcp - Pwet));
      if (std::isnan(Ei) || std::isinf(Ei)) Ei = 0.0;
      
      P_minus_Ei_hist[i] = Prcp - Ei;
      Es_eq_hist[i] = Eeq[i] * Tou;
      
      double sum_PEi = 0.0, sum_Es_eq = 0.0;
      // MATCH MATLAB EXACTLY: [32, 0] means 33 days total window.
      for (int k = 0; k <= 32; ++k) {
        if (i - k >= 0) {
          sum_PEi += P_minus_Ei_hist[i - k];
          sum_Es_eq += Es_eq_hist[i - k];
        }
      }
      Es = std::min(1.0, std::max(0.0, (sum_Es_eq > 0.0) ? (sum_PEi / sum_Es_eq) : 0.0)) * Es_eq_hist[i];
      
      Est = Es + Eca + Ecr;
      Et = Eca + Ecr;
      
      double Pei = (Prcp >= Ei) ? (Prcp - Ei) : 0.0;
      if (Prcp < Ei) Ei = Prcp;
      
      ept = (Est <= Pei + prev_UZ + prev_LZ) ? (Est + Ei) : (Pei + prev_UZ + prev_LZ + Ei);
      
      Ei_out[i] = Ei; Es_out[i] = Es; Et_out[i] = Et; Est_out[i] = Est; EPT_out[i] = ept;
    } else { ept = PET_FAO[i]; }
    
    double R = Prcp * std::pow(std::max(0.0, prev_SM) / parFC, parBETA);
    double SM_dummy = std::max(0.0, std::min(prev_SM + Prcp - R, parFC));
    R = R + std::max(0.0, prev_SM + Prcp - R - parFC) + std::min(0.0, prev_SM + Prcp - R);
    
    EA[i] = ept * std::min(SM_dummy / (parFC * parLP), 1.0);
    SM[i] = std::max(0.0, std::min(SM_dummy - EA[i], parFC));
    EA[i] = EA[i] + std::max(0.0, SM_dummy - EA[i] - parFC) + std::min(0.0, SM_dummy - EA[i]);
    
    Q0[i] = std::max(0.0, std::min(parK1 * prev_UZ + parK0 * std::max(prev_UZ - parUZL, 0.0), prev_UZ));
    double RL = std::max(0.0, std::min(prev_UZ - Q0[i], parPERC));
    
    UZ[i] = prev_UZ + R - Q0[i] - RL;
    Q1[i] = std::max(0.0, std::min(parK2 * prev_LZ, prev_LZ));
    LZ[i] = prev_LZ + RL - Q1[i];
    
    Qsim_mm[i] = Q0[i] + Q1[i];
    Qsim_m3s[i] = Qsim_mm[i] * area / 86.4;
    
    prev_v = v[i]; prev_vl = vl[i]; prev_SM = SM[i]; prev_UZ = UZ[i]; prev_LZ = LZ[i];
  }
  
  if (ifRoute == 1) {
    int maxb = std::max(1, static_cast<int>(std::round(parMAXBAS)));
    std::vector<double> c_list(maxb, 0.0);
    double sum_c = 0.0;
    for (int j = 1; j <= maxb; ++j) {
      double p1 = 0.0, p2 = (maxb + 1) / 2.0, p3 = maxb + 1.0;
      if (j > p1 && j <= p2) c_list[j-1] = (j - p1) / (p2 - p1);
      else if (j > p2 && j <= p3) c_list[j-1] = (p3 - j) / (p3 - p2);
      sum_c += c_list[j-1];
    }
    

    if (sum_c > 0) {
      for (int j = 0; j < maxb; ++j) {
        c_list[j] /= sum_c;
      }
    }
    
    for (int t = 0; t < n; ++t) {
      if (t >= maxb - 1) {
        double temp_q = 0.0;
        for (int k = 0; k < maxb; ++k) temp_q += c_list[k] * Qsim_m3s[t - maxb + 1 + k];
        QMOD[t] = temp_q;
      }
    }
  } else {
    for (int t = 0; t < n; ++t) QMOD[t] = Qsim_m3s[t];
  }
  
  return DataFrame::create(
    Named("QMOD") = QMOD, Named("QT_m3") = Qsim_m3s, Named("QT_mm") = Qsim_mm,
          Named("ET_sim") = EA, Named("P_after_s") = P_eff,
          Named("Ei") = Ei_out, Named("Es") = Es_out, Named("Et") = Et_out,
                Named("Est") = Est_out, Named("EPT") = EPT_out,
                Named("SM") = SM, Named("UZ") = UZ, Named("LZ") = LZ,
                      Named("v") = v, Named("vl") = vl
  );
}