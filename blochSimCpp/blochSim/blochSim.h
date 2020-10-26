#ifndef BLOCHSIM_BLOCHSIM_H
#define BLOCHSIM_BLOCHSIM_H

#pragma once

#include <armadillo>
#include <complex>
#include <map>


using namespace arma;
using namespace std;

#define PI datum::pi

struct blochSimData
{
    vec time;
    cx_vec transSignal;
    vec longSignal;
};

class blochSim
{
public:
    blochSim();
    ~blochSim(){};
    void addSpin(double M0, double T1, double T2, double df, int reps);
    void tip(double alpha, double phase);
    void spoil(double low_ext=-2.0*PI, double up_ext=2.0*PI);
    void forcespoil();
    void freeprecess(double timestep);
    void setReceiverPhase(double phase);
    complex<double> transverseSignal();
    double longitudinalSignal();
    void save(unsigned int key);
    cx_vec steadyStateSignal(int lastNpts = 16);
    blochSimData operator [](int i);
    mat operator ()(int i);
    void print();
    mat getM();
private:
    mat xrot(double theta);
    mat yrot(double theta);
    mat zrot(double theta);
    mat throt(double phi, double theta);
    void AB(mat& A, vec& B, double T, int packet);

    vec M0, T1, T2, df;
    uvec repsPerPacket, packetData0, repsNormalizer;
    mat M;

    int spin_count, tip_count, PrecessionDirection;
    double time, last_tip_time;
    double exciterPhase, receiverPhase;

    std::map<unsigned int,blochSimData> memory;
    uvec mem_keys;
};





#endif //BLOCHSIMCPP_BLOCHSIM_H