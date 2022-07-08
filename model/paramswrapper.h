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


#endif // PARAMSWRAPPER_H
