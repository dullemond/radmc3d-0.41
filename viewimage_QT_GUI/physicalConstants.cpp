#include "physicalConstants.h"
PhysicalConstants *PhysicalConstants ::instance=0;

/**
 *  @brief This method gives a uniq instance of PhysicalConstants class.
 *
 *
 */
PhysicalConstants *PhysicalConstants::getInstance(){
    if(!instance){
        Logger *logger=Logger ::getInstance();
        instance=new PhysicalConstants();
        logger->writeToLogFile("physical constants are initialized");
    }
    return instance;

}
PhysicalConstants::PhysicalConstants(){
    this->gravitationant=6.673e-8; // Gravitational ant
    this->astronomicalUnit = 1.49598e13;      // Astronomical Unit       [cm]
    this->parsec = 3.08572e18;                    // Parsec                  [cm]
    this->protonMass = 1.6726e-24;                         // Mass of proton          [g]
    this->electronMass  = 9.1095e-28;                    // Mass of electron        [g]
    this->bolzmenant = 1.3807e-16;                         // Bolzmann's ant     [erg/K]
    this->planckConstant  = 6.6262e-27;                        // Planck's ant       [erg.s]
    this->unitCharge  = 4.8032e-10;                       // Unit charge
    this->lightSpeed  = 2.9979245800000e10 ;               // Light speed             [cm/s]
    this->thomsonCrossSection= 6.6524e-25;                         // Thompson cross-section  [cm^2]
    this->stefanBolzmanant  = 5.6703e-5;                        // Stefan-Boltzmann   [erg/cm^2/K^4/s]
    this->radiationDensityant  = 7.5657e-15;                       // 4 ss / cc               [erg/cm^3/K^4]
    //
    //     Gas ants
    //
    this->meanMolecWeight= 2.3000;                           // Mean molec weight H2+He+Metals

    //
    //     Alternative units
    //
    this->electronVolt= 1.6022e-12;                          // Electronvolt            [erg]
    this->kiloElectronVolt = 1.6022e-9;                         // Kilo electronvolt       [erg]
    this->micron= 1.e-4;                               // Micron                  [cm]
    this->kilometer   = 1.e5;                             // Kilometer               [cm]
    this->angstroem= 1.e-8;                              // angstroem               [cm]
    //
    //     Astronomy ants
    //
    this->solarLuminosity  = 3.8525e33;                         // Solar luminosity        [erg/s]
    this->solarRadius  = 6.96e10;                           // Solar radius            [cm]
    //sm  = 1.99d33                          // Solar mass              [g]
    this->solarMass  = 1.98892e33;                        // Solar mass              [g]
    this->solarTemperature  = 5.78e3;                            // Solar temperature       [K]
    this->earthMass = 5.9736e27;                          // Mass of Earth           [g]
    this->earthEquatorialRadius = 6.375e08;                          // Equatorial radius Earth [cm]
    this->moonMass = 7.347e25;                           // Mass of moon            [g]
    this->moonRadius = 1.738e08;                           // Radius of moon          [cm]
    this->distanceEarthMoon = 3.844e10;                          // Distance earth-moon (center-to-center)  [cm]
    this->jupiterMass = 1.899e30;                           // Mass of Jupiter         [g]
    this->jupiterRadius = 7.1492e9;                           // Equat. radius Jupiter   [rm]
    this->distanceSunJupiter = 7.78412e13;                        // Distance Jupiter-Sun    [cm]
    //
    //     Time units
    //
    //y= 3.1536d7                            // Year                    [s]
    this->year= 31557600;                          // Year                    [s]
    this->hour= 3.6000e3;                             // Hour                    [s]
    this->day = 8.64e4;                              // Day                     [s]
    //
    //     Math ants
    //
    this->pi  = 3.1415926535897932385;
}
