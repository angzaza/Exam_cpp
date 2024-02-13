#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "VBFHCC_Ana.h"
#include <cmath>

//CSV_READER
int csvFile_reader::open(const std::string &file_name)
{
  infile.open(file_name.c_str());
  if (!infile)
  {
    return -1;
  }
  return 0;
}

std::vector<std::string> csvFile_reader::read_row(){
    std::vector<std::string> row;
    std::string column;
    std::string line_break = "\n";
    
    char c;
    while(infile.get(c)){
        if(c==','){
            row.push_back(column);
            column.clear();
        }
        else{
            column += c;
            if(column.find(line_break)!= std::string::npos) {
                row.push_back(column);
                column.clear();
                break;
            }
        }
    }
    return row;
}


double Jet::mass(){
    double P = pt_ * cosh(eta_); //momentum
    double M2 = pow(en_,2.0) - pow(P,2.0); //4-vector magnitude
    return M2 < 0.0 ? -sqrt(-M2) : sqrt(M2);
}

Jet& Jet::operator+=(const Jet& rjet) {
   //sum 3-momentum components
   double px_tot = px() + rjet.px();
   double py_tot = py() + rjet.py();
   double pz_tot = pz() + rjet.pz();
   //convert rapresentation from (px,py,pz) to (pt, eta, phi)
   double pt_tot = sqrt( pow(px_tot,2.0) + pow(py_tot,2.0) );
   double theta = atan2(pt_tot,pz_tot);
   double eta_tot = -log(tan(theta/2));
   double phi_tot = atan2(py_tot,px_tot);

   this->pt_  = pt_tot;
   this->eta_ = eta_tot;
   this->phi_ = phi_tot;
   this->ctag_  = 0.;
   this->en_  = rjet.en()+en();
   return *this;
}


Jet operator+(const Jet& j1, const Jet& j2) {
    Jet temp{j1}; //using copy constructor
    temp+=j2;
    return temp;
}

bool compareByCvsAll(const Jet& a, const Jet& b){
//  return a.ctag() > b.ctag();
  if (std::isnan(a.ctag())) {
        // If a is NaN, b should come before a
        return false;
    } else if (std::isnan(b.ctag())) {
        // If b is NaN, a should come before b
        return true;
    } else {
        // Neither a nor b is NaN, compare their ctag values
        return a.ctag() > b.ctag();
    }


}


bool EventCandidate::isGoodEvent() {
    //std::cout<<"jet1 eta: "<<candidate[0].eta()<<"   CvsAll: "<<candidate[0].ctag()<<std::endl;
    //std::cout<<"jet2 eta: "<<candidate[1].eta()<<"   CvsAll: "<<candidate[1].ctag()<<std::endl;
    if(fabs(candidate[0].eta())<2.5 && fabs(candidate[1].eta())<2.5 && candidate[0].ctag()>0.5 && candidate[1].ctag()>0.2){
      double mVBF = (candidate[2]+candidate[3]).mass();
      double dEtaVBF = fabs(candidate[2].eta()-candidate[3].eta());
      //std::cout<<"mVBF: "<<mVBF<<"   dEtaVBF: "<<dEtaVBF<<std::endl;
      if(mVBF > 500. && dEtaVBF > 3.8){
        return true;
      }
      else return false;
    }
    else return false;
}

void EventCandidate::print(){
    std::cout << "Higgs candidate mass = " << (candidate[0]+candidate[1]).mass() << std::endl;
}

EventCandidate EventCandidate::analyseRow(std::vector<std::string> row) {
    double jetpt[4]={0.};
    double jeteta[4]={0.};
    double jetphi[4]={0.};
    double jeten[4]={0.};
    double jetctag[4]={0.}; 
    std::array<Jet, 4> incandidate;
    //File-specific!! //take muon info from row
    for(int i = 0; i<4; i++){
        //string to double conversion done using stod from <string> which requires std=c++11
        jetpt[i]  = std::stod(row.at(i*5+1));
        jeteta[i] = std::stod(row.at(i*5+2));
        jetphi[i] = std::stod(row.at(i*5+3));
        jeten[i]  = std::stod(row.at(i*5+4));
        jetctag[i]  = std::stod(row.at(i*5+5));
        //std::cout<<"c tag: "<<jetctag[i]<<std::endl;

        Jet jet_i(jetpt[i], jeteta[i], jetphi[i], jeten[i], jetctag[i]);
        incandidate[i]=jet_i;
    }
    //std::cout<<"before sorting c-tag: "<<incandidate[0].ctag()<<std::endl;
    //std::sort(incandidate.begin(),incandidate.end(),compareByCvsAll );
    //std::cout<<"highest c-tag: "<<incandidate[0].ctag()<<"   eta: "<<incandidate[0].eta()<<std::endl;
    //File-specific!! //take event info from row
    Event inevt(std::stod(row.at(0)));

    //Build EventCandidate to be analysed
    EventCandidate output(inevt , incandidate);
    return output;
}


int main() { 

    csvFile_reader input_jet;
    std::vector<Jet> candidate_list;

    const std::string input_file_name("VBFHCC.csv");
    if (input_jet.open(input_file_name) < 0)
    {
        std::cout << "Cannot open file " << input_file_name.c_str() << std::endl;
        return 1;
    } else std::cout << "Reading file " << input_file_name.c_str() << std::endl;

    std::vector<std::string> row;

    int line_index = 0;
    int header_index = 0; //header position
    while(true){
        if(line_index == header_index) {
            row = input_jet.read_row();
            line_index++; //skip header
        } else {
            //read current line
            //std::cout<<"reading line"<<std::endl;
            row = input_jet.read_row();
            if(!(row.size() == 0)) {
                line_index++;
            }
            else break;

            //expected 3 columns for event info, plus 3 * 4 columns for 3mu info
            //std::cout<<"row size: "<<row.size()<<std::endl;
            if(row.size()<21) break;
             
            EventCandidate current_candidate = current_candidate.analyseRow(row);
            //std::cout<<"current cand good event: "<<current_candidate.isGoodEvent()<<std::endl;
            if(current_candidate.isGoodEvent()) current_candidate.print(); 
                
        }
    }
    return 0;
}

