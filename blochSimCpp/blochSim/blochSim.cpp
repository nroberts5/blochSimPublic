#include "blochSim.h"

#include <iostream>
#include <armadillo>
#include <cmath>

using namespace arma;
using namespace std;

/*************************************************** PUBLIC METHODS ***************************************************/
blochSim::blochSim(){
    this->PrecessionDirection = 1;
    this->spin_count = this->tip_count = 0;
    this->time = this->last_tip_time = 0.0;
    this->exciterPhase = this->receiverPhase = 0.0;
}

void blochSim::addSpin(double M0, double T1, double T2, double df, int reps){
    this->M0.resize(this->M0.n_elem+1); this->M0(this->M0.n_elem-1) = M0;
    this->T1.resize(this->T1.n_elem+1); this->T1(this->T1.n_elem-1) = T1;
    this->T2.resize(this->T2.n_elem+1); this->T2(this->T2.n_elem-1) = T2;
    this->df.resize(this->df.n_elem+1); this->df(this->df.n_elem-1) = df;
    this->repsPerPacket.resize(this->repsPerPacket.n_elem+1); this->repsPerPacket(this->repsPerPacket.n_elem-1) = reps;
    this->packetData0.resize(this->packetData0.n_elem+1); this->packetData0(this->packetData0.n_elem-1) = this->M.n_cols;


    this->M.resize(3,this->M.n_cols+reps); vec m(3); m << 0 << 0 << M0;
    this->M(arma::span::all,arma::span(this->M.n_cols-reps,this->M.n_cols-1)).each_col() = m;

    this->repsNormalizer.resize(this->repsNormalizer.n_elem + reps);
    this->repsNormalizer(arma::span(this->repsNormalizer.n_elem-reps,this->repsNormalizer.n_elem-1)).fill(reps);

    this->spin_count++;
}

void blochSim::tip(double alpha, double phase){
    this->exciterPhase = phase;
    this->receiverPhase = phase;
    this->M = this->throt(alpha,this->exciterPhase)*this->M;
    this->tip_count++;
    this->last_tip_time = this->time;
}

 void blochSim::spoil(double low_ext, double up_ext){
     // Each species recieves the same spoiling
     for (unsigned int packet = 0; packet < this->repsPerPacket.n_elem; packet++)
     {
        int packet_start = this->packetData0(packet);
        int packet_stop = packet_start + repsPerPacket(packet) - 1;
        vec phases = arma::linspace(low_ext, up_ext, packet_stop - packet_start + 1);
        for (unsigned int k = 0; k < phases.n_elem; k++)
        {
            this->M(arma::span::all,packet_start+k) = this->zrot(phases(k))*this->M(arma::span::all,packet_start+k);
        }
     }
 }

void blochSim::forcespoil(){
    this->M(arma::span(0,1),arma::span::all).fill(0.0);
}

void blochSim::freeprecess(double timestep){
    for (unsigned int packet = 0; packet < this->repsPerPacket.n_elem; packet++)
    {
        mat A(3,3,arma::fill::zeros); vec B(3,arma::fill::zeros);
        int packet_start = this->packetData0(packet);
        int packet_stop = packet_start + repsPerPacket(packet) - 1;
        this->AB(A, B, timestep, packet);
        this->M(arma::span::all,arma::span(packet_start,packet_stop)).each_col([&](vec& m){ m = A*m + B;});
    }
    this->time += timestep;
};

void blochSim::setReceiverPhase(double phase){
    this->receiverPhase = phase;
}

complex<double> blochSim::transverseSignal(){
    arma::span all = arma::span::all;
    cx_rowvec sig = cx_rowvec(M(0,all),this->M(1,all));
    sig *= std::exp(-(std::complex<double>)1i*this->receiverPhase);
    sig = sig / this->repsNormalizer.st();
    return arma::as_scalar(arma::sum(sig));
}

double blochSim::longitudinalSignal(){
    return arma::as_scalar(sum(this->M(2,arma::span::all) / this->repsNormalizer.st()));
}

void blochSim::save(unsigned int key){
    bool newKey = !(bool)(this->memory.count(key));
    int numEls = this->memory[key].time.n_elem;
    this->memory[key].time.resize(numEls+1);
    this->memory[key].transSignal.resize(numEls+1);
    this->memory[key].longSignal.resize(numEls+1);

    this->memory[key].time(numEls) = this->time;
    this->memory[key].transSignal(numEls) = this->transverseSignal();
    this->memory[key].longSignal(numEls) = this->longitudinalSignal();

    if (newKey)
    {
        this->mem_keys.resize(this->mem_keys.n_elem+1);
        this->mem_keys(this->mem_keys.n_elem-1) = key;
    }
}

cx_vec blochSim::steadyStateSignal(int lastNpts){
    cx_vec ss_sig(mem_keys.n_elem,arma::fill::zeros);
    for (unsigned int k = 0; k < this->mem_keys.n_elem; ++k) {
        unsigned int key = this->mem_keys(k);
        int index1 = (int)this->memory[key].transSignal.n_elem-lastNpts > 0 ? this->memory[key].transSignal.n_elem-lastNpts : 0;
        ss_sig(k) = arma::median(this->memory[key].transSignal(arma::span(index1,this->memory[key].transSignal.n_elem-1)));
    }
    return ss_sig;
}

blochSimData blochSim::operator[](int i){
    return this->memory[i];
}

mat blochSim::operator()(int i){
    arma::span all = arma::span::all;
    mat data(4,this->memory[i].time.n_elem);
    data(0,all) = real(this->memory[i].transSignal).st();
    data(1,all) = imag(this->memory[i].transSignal).st();
    data(2,all) = this->memory[i].longSignal.st();
    data(3,all) = this->memory[i].time.st();
    return data;
}

void blochSim::print(){
    vec meanvec = arma::mean(this->M,1);
    cout<<"transverseSignal:"<<this->transverseSignal()<<endl;
    cout<<"longitudinalSignal:"<<this->longitudinalSignal()<<endl;
    meanvec.print("M: ");
};

mat blochSim::getM(){
    return this->M;
}

/*************************************************** PRIVATE METHODS ***************************************************/
mat blochSim::xrot(double theta){
    mat Rx(3, 3, arma::fill::zeros);
    Rx << 1.0 << 0.0 << 0.0 << endr
       << 0.0 << cos(theta) << -sin(theta) << endr 
       << 0.0 << sin(theta) << cos(theta) << endr;
    return Rx;
}

mat blochSim::yrot(double theta){
    mat Ry(3, 3, arma::fill::zeros);
    Ry << cos(theta) << 0 << sin(theta) << endr 
       << 0 << 1 << 0 << endr 
       << -sin(theta) << 0 << cos(theta) << endr;
    return Ry;
}

mat blochSim::zrot(double theta){
    mat Rz(3, 3, arma::fill::zeros);
    Rz << cos(theta) << -sin(theta) << 0 << endr 
       << sin(theta) << cos(theta) << 0 << endr 
       << 0 << 0 << 1 << endr;
    return Rz;
}

mat blochSim::throt(double phi, double theta){
    mat Rz = this->zrot(-theta);
    mat Ry = this->yrot(phi);
    mat Rth = Rz.i()*Ry*Rz;
    return Rth;
}

void blochSim::AB(mat& A, vec& B, double T, int packet){
    double phi = this->PrecessionDirection*2*PI*this->df(packet)*T;
    double E1 = exp(-T/this->T1(packet));
    double E2 = exp(-T/this->T2(packet));

    mat A1(3, 3, arma::fill::zeros);
    A1 << E2 << 0 << 0 << endr 
       << 0 << E2 << 0 << endr 
       << 0 << 0 << E1 << endr;

    A = A1 * this->zrot(phi);
    B << 0 << 0 << this->M0(packet)*(1-E1);
}