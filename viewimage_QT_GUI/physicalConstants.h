#ifndef PhysicalConstants_H
#define PhysicalConstants_H
#include "logger.h"
class PhysicalConstants
{
public:
    PhysicalConstants();
    static PhysicalConstants *getInstance();
    double  getGravitationant(){return gravitationant;}
    double  getAstronomicalUnit(){return astronomicalUnit;}
    double  getParsec(){return parsec;}
    double  getProtonMass(){return protonMass;}
    double  getElectronMass(){return electronMass;}
    double  getbolzmenant(){return bolzmenant;}
    double  getPlanckconstant(){return planckConstant;}
    double  getUnitCharge(){return unitCharge;}
    double  getLightSpeed(){return lightSpeed;}
    double  getThomsonCrossSection(){return thomsonCrossSection;}
    double  getStefanBolzmanant(){return stefanBolzmanant;}
    double  getRadiationDensityant(){return radiationDensityant;}
    double  getMeanMolecWeight(){return meanMolecWeight;}
    double  getElectronVolt(){return electronVolt;}
    double  getKiloElectronVolt(){return kiloElectronVolt;}
    double  getMicron(){return micron;}
    double  getKilometer(){return kilometer;}
    double  getAngstroem(){return angstroem;}
    double  getSolarLuminosity(){return solarLuminosity;}
    double  getSolarRadius(){return solarRadius;}
    double  getSolarMass(){return solarMass;}
    double  getSolarTemperature(){return solarTemperature;}
    double  getEarthMass(){return earthMass;}
    double  getEarthEquatorialRadius(){return earthEquatorialRadius;}
    double  getMoonMass(){return moonMass;}
    double  getMmoonRadius(){return moonRadius;}
    double  getDistanceEarthMoon(){return distanceEarthMoon;}
    double  getJupiterMass(){return jupiterMass;}
    double  getJupiterRadius(){return jupiterRadius;}
    double  getDistanceSunJupiter(){return distanceSunJupiter;}
    double  getYear(){return year;}
    double  getHour(){return hour;}
    double  getDay(){return day;}
    double  getPi(){return pi;}
private:
    static PhysicalConstants *instance;
    //Natural ants in CGS
    double  gravitationant; // Gravitational ant
    double  astronomicalUnit;      // Astronomical Unit       [cm]
    double  parsec ;                    // Parsec                  [cm]
    double   protonMass;                         // Mass of proton          [g]
    double   electronMass;                    // Mass of electron        [g]
    double  bolzmenant;                         // Bolzmann's ant     [erg/K]
    double   planckConstant;                        // Planck's constant       [erg.s]
    double   unitCharge;                       // Unit charge
    double  lightSpeed;               // Light speed             [cm/s]
    double  thomsonCrossSection;                         // Thompson cross-section  [cm^2]
    double  stefanBolzmanant;                        // Stefan-Boltzmann   [erg/cm^2/K^4/s]
    double  radiationDensityant;                       // 4 ss / cc               [erg/cm^3/K^4]
    //
    //     Gas ants
    //
    double  meanMolecWeight;                           // Mean molec weight H2+He+Metals

    //
    //     Alternative units
    //
    double  electronVolt;                          // Electronvolt            [erg]
    double  kiloElectronVolt;                         // Kilo electronvolt       [erg]
    double  micron;                               // Micron                  [cm]
    double  kilometer;                             // Kilometer               [cm]
    double  angstroem;                              // angstroem               [cm]
    //
    //     Astronomy ants
    //
    double  solarLuminosity ;                         // Solar luminosity        [erg/s]
    double     solarRadius;                           // Solar radius            [cm]
    //sm  = 1.99d33                          // Solar mass              [g]
    double      solarMass;                        // Solar mass              [g]
    double  solarTemperature;                            // Solar temperature       [K]
    double     earthMass;                          // Mass of Earth           [g]
    double    earthEquatorialRadius ;                          // Equatorial radius Earth [cm]
    double   moonMass;                           // Mass of moon            [g]
    double     moonRadius;                           // Radius of moon          [cm]
    double     distanceEarthMoon;                          // Distance earth-moon (center-to-center)  [cm]
    double     jupiterMass ;                           // Mass of Jupiter         [g]
    double  jupiterRadius ;                           // Equat. radius Jupiter   [rm]
    double   distanceSunJupiter ;                        // Distance Jupiter-Sun    [cm]
    //
    //     Time units
    //
    //y= 3.1536d7                            // Year                    [s]
    double    year;                          // Year                    [s]
    double     hour;                             // Hour                    [s]
    double    day ;                              // Day                     [s]
    //
    //     Math ants
    //
    double   pi  ;
};

#endif // PhysicalConstants_H
