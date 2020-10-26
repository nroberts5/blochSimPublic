#include <iostream>
#include <string>
#include "blochSim.h"
#include "params.h"

#include <armadillo>
#include <ctime>

#define EPS 1e-6
#define minNEX 100

using namespace std;
using namespace arma;

void standardSGRE(blochSim& p, paramsStruct& params, acqParamsStruct& acqParams, simParamsStruct& simParams, int save0 = 0)
{
    int k = 0; double theta = 0;
    vec  ss_sig(acqParams.necho, arma::fill::zeros);
    vec nss_sig(save0 + acqParams.necho, arma::fill::zeros);
    for (int nex = 0; nex < simParams.blochNEX; nex++)
    {
        // rf spoiling
        theta = fmod((theta + k * acqParams.seed), 360.0); k++;

        // TIP
        p.tip(params.beta*acqParams.flip*datum::pi/180.0, datum::pi/180.0*theta);

        double t = 0;
        for (int kt = 0; kt < acqParams.necho; kt++)
        {
            // Free Precess
            p.freeprecess(acqParams.tes(kt)-t);
            // Measure
            p.save(save0+kt);
            t = acqParams.tes(kt);
        }

        //Decay
        p.freeprecess(acqParams.tr-t);

        //Spoil
        if (!simParams.forceSpoil)
        {
            p.spoil();
        } else{
            p.forcespoil();
        }

        nss_sig = arma::abs(p.steadyStateSignal());
        ss_sig = arma::abs(nss_sig(arma::span(save0,save0+acqParams.necho-1)) - ss_sig);
        if (arma::norm(ss_sig)<EPS && nex >= minNEX)
        {
           cout<<"Steady State Achieved after "<<nex<<" NEX"<<endl;
            break;
        } else
        {
            ss_sig = nss_sig(arma::span(save0,save0+acqParams.necho-1));
        }
    }
}

void initFWSpins(blochSim& p, paramsStruct& params, acqParamsStruct& acqParams, simParamsStruct& simParams)
{
    p.addSpin(params.rhoW, params.T1W,1.0/params.R2s, params.psi, simParams.BlochSimReps);
    for (unsigned int k = 0; k < acqParams.relAmps.n_elem; k++)
    {
        double m0 = params.rhoF * acqParams.relAmps(k);
        if (m0>0.0)
        {
            p.addSpin(m0, params.T1F, 1.0/params.R2s,
                          params.psi + acqParams.fp(k),
                          floor(simParams.BlochSimReps/acqParams.relAmps.n_elem));
        }
    }
}

void simpleDemo()
{
    paramsStruct params{};
    params.init();
    acqParamsStruct acqParams{};
    acqParams.init();
    simParamsStruct simParams{};
    simParams.BlochSimReps = 10000;
    simParams.forceSpoil = false;
    simParams.init();

    blochSim mySim{};
    initFWSpins(mySim, params, acqParams, simParams);
    standardSGRE(mySim, params, acqParams, simParams);
    cx_rowvec ss_sig = mySim.steadyStateSignal().st();
    arma::abs(ss_sig).print("Steady State Magnitude:");
    arma::atan(imag(ss_sig)/real(ss_sig)).print("Steady State Phase:");
}

int main(int argc, char const *argv[])
{
    time_t tstart, tend; 
    tstart = time(0);
    simpleDemo();

    tend = time(0); 
    cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    return 0;
}
