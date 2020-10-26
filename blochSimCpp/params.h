#ifndef BLOCHSIMCPP_MYPARAMS_H
#define BLOCHSIMCPP_MYPARAMS_H

#pragma once

#include <iostream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

struct speciesStruct{
    std::string name;
    vec frequency;
    vec relAmps;

    speciesStruct(std::string name, const vec& frequency, const vec& relAmps)
    {
        this->name = name;
        this->frequency = frequency;
        this->relAmps = relAmps;
    }
};

struct paramsStruct{
    paramsStruct() : 
    beta(1.0), 
    M0(100.0),
    T1W(576e-3),
    T1F(280e-3),
    R2s(40),
    psi(25.0),
    PDFF(0.0),
    rhoW(0.0),
    rhoF(0.0)
    {}
    double beta;
    double M0;
    double T1W;
    double T1F;
    double R2s;
    double psi;
    double PDFF;
    double rhoW;
    double rhoF;

    void init() {
        this->rhoW = (1-this->PDFF) * this->M0;
        this->rhoF = this->PDFF * this->M0;
    }
};

struct acqParamsStruct{
    acqParamsStruct() :
    isPhantom(false),
    field_strength(1.5),
    flip(5.0),
    TE1(1.2e-3),
    DTE(2.0e-3),
    necho(6),
    tr(15.2e-3),
    seed(117.0),
    precessionDirection(1.0)
    {}
    bool isPhantom;
    double field_strength;
    double flip;
    double TE1;
    double DTE;
    int necho;
    double tr;
    double seed;
    int precessionDirection;
    vec tes;
    std::vector<speciesStruct> species;
    vec relAmps;
    vec ppm;
    vec fp;

    void init() {
        this->tes = linspace(this->TE1,this->TE1+(this->necho-1)*this->DTE,this->necho);

        vec water_frequency = {0};
        vec water_amps = {1};
        speciesStruct water("water", water_frequency, water_amps);
        this->species.emplace_back(water);

        vec fat_frequency = {-3.80, -3.40, -2.60, -1.94, -0.39, 0.60};
        vec fat_amps = {0.087, 0.693, 0.128, 0.004, 0.039, 0.048};
        speciesStruct fat("fat", fat_frequency, fat_amps);
        this->species.emplace_back(fat);

        vec fatPhantom_frequency = fat_frequency - 0.1096;
        vec fatPhantom_amps = fat_amps;
        speciesStruct fatPhantom("fatPhantom", fatPhantom_frequency, fatPhantom_amps);
        this->species.emplace_back(fatPhantom);

        if (isPhantom)
        {
            this->ppm = species[2].frequency;
            this->relAmps = species[2].relAmps;
        } else{
            this->ppm = species[1].frequency;
            this->relAmps = species[1].relAmps;
        }
        this->fp = 42.577 * this->ppm * this->field_strength;
    }

};

struct simParamsStruct{
    simParamsStruct() :
    BlochSimReps(1000),
    blochNEX(1000),
    forceSpoil(false)
    {}

    int BlochSimReps;
    int blochNEX;
    bool forceSpoil;

    void init() {}
};

#endif //BLOCHSIMCPP_MYPARAMS_H
