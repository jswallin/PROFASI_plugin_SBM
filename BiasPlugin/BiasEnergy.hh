#ifndef BiasEnergy_HH
#define BiasEnergy_HH
#include <Energy/Energy.hh>

namespace prf
{
  class BiasEnergy : public Energy	// inherit the Energy class
  {
  public:
    BiasEnergy();	// constructor
    ~BiasEnergy() = default;	//destructor
    void init(); //initialize variables
    void set_pars(std::string s) override;
    double evaluate() override; // perform actual calculation
    
  private:
    std::string filename;
    size_t MAXCONT=1000;
    size_t atm1[1000],atm2[1000],ncont;
    size_t atm1_2[1000],atm2_2[1000],ncont2;
    size_t ntd1,ntd2,ctd1,ctd2;
    int doublet[1000];
    double atm_dist[1000],lamSC1[1000];
    double atm_dist2[1000],lamSC2[1000];
    double lambda_SC,lambda_SC_loc,ksi_SC;
    double kappa1,kappa2,kappa_d;
  };
}

#endif



