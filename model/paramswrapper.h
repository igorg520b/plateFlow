#ifndef PARAMSWRAPPER_H
#define PARAMSWRAPPER_H

#include <QObject>
#include "parameters_sim.h"

class ParamsWrapper : public QObject
{
    Q_OBJECT

public:
    icy::Params *p;

    Q_PROPERTY(bool f_EnableFracture READ getEnableFracture WRITE setEnableFracture NOTIFY propertyChanged)
//    Q_PROPERTY(double i_ READ get WRITE set NOTIFY propertyChanged)
//    Q_PROPERTY(int i_ READ get WRITE set NOTIFY propertyChanged)

    Q_PROPERTY(double i_InitialTimeStep READ getInitialTimeStep WRITE setInitialTimeStep NOTIFY propertyChanged)
    Q_PROPERTY(int i_MaxSteps READ getMaxSteps WRITE setMaxSteps NOTIFY propertyChanged)
    Q_PROPERTY(int i_MinIter READ getMinIter WRITE setMinIter NOTIFY propertyChanged)
    Q_PROPERTY(int i_MaxIter READ getMaxIter WRITE setMaxIter NOTIFY propertyChanged)
    Q_PROPERTY(double i_ConvergenceEpsilon READ getConvergenceEpsilon WRITE setConvergenceEpsilon NOTIFY propertyChanged)
    Q_PROPERTY(double i_ConvergenceCutoff READ getConvergenceCutoff WRITE setConvergenceCutoff NOTIFY propertyChanged)


    Q_PROPERTY(double p_Gravity READ getGravity WRITE setGravity NOTIFY propertyChanged)
    Q_PROPERTY(double p_WaterDensity READ getWaterDensity WRITE setWaterDensity NOTIFY propertyChanged)
    Q_PROPERTY(double p_IceDensity READ getIceDensity WRITE setIceDensity NOTIFY propertyChanged)
    Q_PROPERTY(double p_PoissonsRatio READ getPoissonsRatio WRITE setPoissonsRatio NOTIFY propertyChanged)
    Q_PROPERTY(double p_YoungsModulus READ getYoungsModulus WRITE setYoungsModulus NOTIFY propertyChanged)
    Q_PROPERTY(double p_Thickness READ getThickness WRITE setThickness NOTIFY propertyChanged)


    Q_PROPERTY(double f_FractureTractionThreshold READ getFractureTractionThreshold WRITE setFractureTractionThreshold NOTIFY propertyChanged)
    Q_PROPERTY(double f_FractureWeakeningCoeff READ getFractureWeakeningCoeff WRITE setFractureWeakeningCoeff NOTIFY propertyChanged)
    Q_PROPERTY(double f_CutoffCoefficient READ getCutoffCoefficient WRITE setCutoffCoefficient NOTIFY propertyChanged)
    Q_PROPERTY(double f_FractureAngleThreshold READ getFractureAngleThreshold WRITE setFractureAngleThreshold NOTIFY propertyChanged)
    Q_PROPERTY(double f_FractureAreaThreshold READ getFractureAreaThreshold WRITE setFractureAreaThreshold NOTIFY propertyChanged)
    Q_PROPERTY(double f_SubsteppingTimestepFactor READ getSubsteppingTimestepFactor WRITE setSubsteppingTimestepFactor NOTIFY propertyChanged)
    Q_PROPERTY(int f_SubstepIterations READ getSubstepIterations WRITE setSubstepIterations NOTIFY propertyChanged)
    Q_PROPERTY(int f_SubstepRadius READ getSubstepRadius WRITE setSubstepRadius NOTIFY propertyChanged)
    Q_PROPERTY(int f_FractureMaxSubsteps READ getFractureMaxSubsteps WRITE setFractureMaxSubsteps NOTIFY propertyChanged)

    Q_PROPERTY(double m_VideoTimeStep READ getVideoTimeStep WRITE setVideoTimeStep NOTIFY propertyChanged)

    // EnableCollisions
    bool getEnableFracture() {return p->EnableFracture;}
    void setEnableFracture(bool val) {p->EnableFracture=val; Q_EMIT propertyChanged();}


    // InitialTimeStep
    double getInitialTimeStep() {return p->InitialTimeStep;}
    void setInitialTimeStep(double val) {p->InitialTimeStep=val; Q_EMIT propertyChanged();}


    // MaxSteps
    int getMaxSteps() {return p->MaxSteps;}
    void setMaxSteps(int val) {p->MaxSteps=val; Q_EMIT propertyChanged();}

    // MinIter
    int getMinIter() {return p->MinIter;}
    void setMinIter(int val) {p->MinIter=val; Q_EMIT propertyChanged();}

    // MaxIter
    int getMaxIter() {return p->MaxIter;}
    void setMaxIter(int val) {p->MaxIter=val; Q_EMIT propertyChanged();}


    // ConvergenceEpsilon
    double getConvergenceEpsilon() {return p->ConvergenceEpsilon;}
    void setConvergenceEpsilon(double val) {p->ConvergenceEpsilon=val; Q_EMIT propertyChanged();}

    // ConvergenceCutoff
    double getConvergenceCutoff() {return p->ConvergenceCutoff;}
    void setConvergenceCutoff(double val) {p->ConvergenceCutoff=val; Q_EMIT propertyChanged();}



    // Gravity
    double getGravity() {return p->Gravity;}
    void setGravity(double val) {p->Gravity=val; Q_EMIT propertyChanged();}

    // WaterDensity
    double getWaterDensity() {return p->WaterDensity;}
    void setWaterDensity(double val) {p->WaterDensity=val; Q_EMIT propertyChanged();}

    // IceDensity
    double getIceDensity() {return p->IceDensity;}
    void setIceDensity(double val) {p->IceDensity=val; Q_EMIT propertyChanged();}

    // PoissonsRatio
    double getPoissonsRatio() {return p->PoissonsRatio;}
    void setPoissonsRatio(double val) {p->PoissonsRatio=val; Q_EMIT propertyChanged();}

    // YoungsModulus
    double getYoungsModulus() {return p->YoungsModulus;}
    void setYoungsModulus(double val) {p->YoungsModulus=val; Q_EMIT propertyChanged();}

    // Thickness
    double getThickness() {return p->Thickness;}
    void setThickness(double val) {p->Thickness=val; Q_EMIT propertyChanged();}



    // fracture

    // FractureTractionThreshold
    double getFractureTractionThreshold() {return p->FractureTractionThreshold;}
    void setFractureTractionThreshold(double val) {p->FractureTractionThreshold=val; Q_EMIT propertyChanged();}

    // FractureWeakeningCoeff
    double getFractureWeakeningCoeff() {return p->FractureWeakeningCoeff;}
    void setFractureWeakeningCoeff(double val) {p->FractureWeakeningCoeff=val; Q_EMIT propertyChanged();}

    // CutoffCoefficient
    double getCutoffCoefficient() {return p->CutoffCoefficient;}
    void setCutoffCoefficient(double val) {p->CutoffCoefficient=val; Q_EMIT propertyChanged();}

    // FractureAngleThreshold
    double getFractureAngleThreshold() {return p->FractureAngleThreshold;}
    void setFractureAngleThreshold(double val) {p->FractureAngleThreshold=val; Q_EMIT propertyChanged();}

    // FractureAreaThreshold
    double getFractureAreaThreshold() {return p->FractureAreaThreshold;}
    void setFractureAreaThreshold(double val) {p->FractureAreaThreshold=val; Q_EMIT propertyChanged();}

    // SubsteppingTimestepFactor
    double getSubsteppingTimestepFactor() {return p->SubsteppingTimestepFactor;}
    void setSubsteppingTimestepFactor(double val) {p->SubsteppingTimestepFactor=val; Q_EMIT propertyChanged();}


    // SubstepIterations
    int getSubstepIterations() {return p->SubstepIterations;}
    void setSubstepIterations(int val) {p->SubstepIterations=val; Q_EMIT propertyChanged();}

    // SubstepRadius
    int getSubstepRadius() {return p->SubstepRadius;}
    void setSubstepRadius(int val) {p->SubstepRadius=val; Q_EMIT propertyChanged();}

    // FractureMaxSubsteps
    int getFractureMaxSubsteps() {return p->FractureMaxSubsteps;}
    void setFractureMaxSubsteps(int val) {p->FractureMaxSubsteps=val; Q_EMIT propertyChanged();}


    // misc
    // VideoTimeStep
    double getVideoTimeStep() {return p->VideoTimeStep;}
    void setVideoTimeStep(double val) {p->VideoTimeStep=val; Q_EMIT propertyChanged();}


    //
    //double get() {return p->;}
    //void set(double val) {p->=val; Q_EMIT propertyChanged();}

    //    int get() {return p->;}
    //    void set(int val) {p->=val; Q_EMIT propertyChanged();}



Q_SIGNALS:
    void propertyChanged();
};


#endif // PARAMSWRAPPER_H
