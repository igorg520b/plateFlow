#ifndef PARAMSWRAPPER_H
#define PARAMSWRAPPER_H

#include <QObject>
#include "parameters_sim.h"

class ParamsWrapper : public QObject
{
    Q_OBJECT

public:
    icy::Params p;

    Q_PROPERTY(bool f_EnableFracture READ getEnableFracture WRITE setEnableFracture NOTIFY propertyChanged)
    // EnableCollisions
    bool getEnableEnableFracture() {return p.EnableFracture;}
    void setEnableFracture(bool val) {p.EnableFracture=val; Q_EMIT propertyChanged();}

Q_SIGNALS:
    void propertyChanged();
};


struct icy::SimParams
{

public:
    SimParams() {p.Reset();}


    // EnableCZs
    bool getEnableCZs() {return p.EnableCZs;}
    void setEnableCZs(bool val) {p.EnableCZs=val; Q_EMIT propertyChanged();}

    // LinearCollisions
    bool getLinearCollisions() {return p.LinearCollisions;}
    void setLinearCollisions(bool val) {p.LinearCollisions=val; Q_EMIT propertyChanged();}

    // LinearCollisions
    bool getLinearIndenter() {return p.LinearIndenter;}
    void setLinearIndenter(bool val) {p.LinearIndenter=val; Q_EMIT propertyChanged();}

    // MaxSteps
    int getMaxSteps() {return p.MaxSteps;}
    void setMaxSteps(int val) {p.MaxSteps=val; Q_EMIT propertyChanged();}

    // Initial Time Step
    double getInitialTimeStep() {return p.InitialTimeStep;}
    void setInitialTimeStep(double val) {p.InitialTimeStep=val; Q_EMIT propertyChanged();}

    // Convergence Epsilon
    double getConvergenceEpsilon() {return p.ConvergenceEpsilon;}
    void setConvergenceEpsilon(double val) {p.ConvergenceEpsilon=val; Q_EMIT propertyChanged();}

    // ConvergenceCutoff
    double getConvergenceCutoff() {return p.ConvergenceCutoff;}
    void setConvergenceCutoff(double val) {p.ConvergenceCutoff=val; Q_EMIT propertyChanged();}


    // Video Time Step
    double getVideoTimeStep() {return p.VideoTimeStep;}
    void setVideoTimeStep(double val) {p.VideoTimeStep=val; Q_EMIT propertyChanged();}

    // Initial Time Step
    double getGravity() {return p.Gravity;}
    void setGravity(double val) {p.Gravity=val; Q_EMIT propertyChanged();}

    // Initial Time Step
    double getDensity() {return p.Density;}
    void setDensity(double val) {p.Density=val; Q_EMIT propertyChanged();}

    double getYoungsModulus() {return p.YoungsModulus;}
    void setYoungsModulus(double ym) { p.YoungsModulus=ym; p.RecomputeLamdaAndMu(); Q_EMIT propertyChanged(); }

    double getPoissonsRatio() {return p.PoissonsRatio;}
    void setPoissonsRatio(double nu) { p.PoissonsRatio=nu; p.RecomputeLamdaAndMu(); Q_EMIT propertyChanged(); }

    // Kappa
    double getKappa() {return p.Kappa;}
    void setKappa(double val) {p.Kappa=val; Q_EMIT propertyChanged();}

    // Dhat
    double getDhat() {return p.dhat;}
    void setDhat(double val) {p.dhat=val; Q_EMIT propertyChanged();}

    // Indenter Radius
    double getIndenterRadius() {return p.R_indenter;}
    void setIndenterRadius(double val) {p.R_indenter=val; Q_EMIT propertyChanged();}

    // Indentation rate
    double getIndentationRate() {return p.indentation_rate;}
    void setIndentationRate(double val) {p.indentation_rate=val; Q_EMIT propertyChanged();}

    // Indenter Kappa
    double getIndKappa() {return p.ind_Kappa;}
    void setIndKappa(double val) {p.ind_Kappa=val; Q_EMIT propertyChanged();}

    // Indenter Dhat
    double getIndDhat() {return p.ind_dhat;}
    void setIndDhat(double val) {p.ind_dhat=val; Q_EMIT propertyChanged();}


    void set_cz_alpha(double value) { p.cz_alpha=value; p.RecomputeCZParams(); Q_EMIT propertyChanged(); }
    double get_cz_alpha() {return p.cz_alpha;}

    void set_cz_beta(double value) { p.cz_beta=value; p.RecomputeCZParams(); Q_EMIT propertyChanged(); }
    double get_cz_beta() {return p.cz_beta;}

    void set_cz_lambda_n(double value) {p.cz_lambda_n = value; p.RecomputeCZParams(); Q_EMIT propertyChanged(); }
    double get_cz_lambda_n() {return p.cz_lambda_n;}

    void set_cz_lambda_t(double value) {p.cz_lambda_t = value; p.RecomputeCZParams(); Q_EMIT propertyChanged(); }
    double get_cz_lambda_t() {return p.cz_lambda_t;}

    void set_cz_phi_n(double value) {p.cz_phi_n = value; p.RecomputeCZParams(); Q_EMIT propertyChanged(); }
    double get_cz_phi_n() {return p.cz_phi_n;}

    void set_cz_phi_t(double value) {p.cz_phi_t = value; p.RecomputeCZParams(); Q_EMIT propertyChanged(); }
    double get_cz_phi_t() {return p.cz_phi_t;}

    void set_cz_sigma_max(double value) {p.cz_sigma_max = value; p.RecomputeCZParams(); Q_EMIT propertyChanged(); }
    double get_cz_sigma_max() {return p.cz_sigma_max;}

    void set_cz_tau_max(double value) {p.cz_tau_max = value; p.RecomputeCZParams(); Q_EMIT propertyChanged(); }
    double get_cz_tau_max() {return p.cz_tau_max;}
    double get_cz_del_n() {return p.cz_del_n;}
    double get_cz_del_t() {return p.cz_del_t;}

    void Reset() {p.Reset(); Q_EMIT propertyChanged();}


};

#endif // PARAMSWRAPPER_H
