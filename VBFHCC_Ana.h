#include <fstream>
#include <string>
#include <vector>
#include <array>

class csvFile_reader
{
public:
  csvFile_reader(){};

  int open(const std::string &file_name);
  void read_file(int header);
  std::vector<std::string> read_row();

private:
  std::ifstream infile;
};


class Jet
{
public:
    Jet(){};

    Jet(double i_pt, double i_eta, double i_phi, double i_en, double i_ctag)
        : pt_(i_pt), eta_(i_eta), phi_(i_phi), en_(i_en), ctag_(i_ctag)
    {}

    Jet(const Jet & j){
        pt_  = j.pt_;
        eta_ = j.eta_;
        phi_ = j.phi_;
        en_  = j.en_;
        ctag_ = j.ctag_;
    }

    double pt()  const { return pt_; }
    double eta() const { return eta_; }
    double phi() const { return phi_; }
    double en()  const { return en_; }
    double ctag()  const { return ctag_; }
    double mass();
    double px() const { return (pt_ * cos(phi_ )); }
    double py() const { return (pt_ * sin(phi_ )); }
    double pz() const { return (pt_ * sinh(eta_)); }
    void print();


    Jet& operator+=(const Jet& rjet);

private:
    double pt_ {0};
    double eta_{0};
    double phi_{0};
    double en_ {0};
    double ctag_ {0};

};


class Event
{
public:
    Event(){};

    Event(int inevt)
        : evt_(inevt)
    {}

    int evt()  const { return evt_; }
    // Overload == operator to compare events
    //bool operator==(const Event);

private:
    int evt_ {0};
};


class EventCandidate
{
public:
    EventCandidate(){}

    EventCandidate(Event inevent, std::array<Jet, 4>& incandidate)
        : event(inevent), candidate(incandidate)
    {}

    bool isGoodEvent();
    EventCandidate analyseRow(std::vector<std::string> row);
    void print();

private:
    Event event; //call Event()
    std::array<Jet, 4> candidate; //call Jet()
    double chi2{0};
};

