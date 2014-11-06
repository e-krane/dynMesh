// Defining constants and variables
const scalar Pi = 3.14159265359;
const vector CofR = dict_.subDict("coordinateSystem").lookup("origin");
const scalar rhoInf = readScalar( dict_.subDict("forcesDict").lookup("rhoInf") );
const wordList patches = dict_.subDict("forcesDict").lookup("patches");
const scalar Iz = readScalar( dict_.subDict("variableRpm").lookup("Iz") );
const scalar cutoffRpm = readScalar( dict_.subDict("variableRpm").lookup("cutoffRpm") );

scalar Pnom = readScalar( dict_.subDict("variableRpm").lookup("P") );
scalar rpmNom = readScalar( dict_.subDict("variableRpm").lookup("rpmNom") );
scalar P;
scalar alpha;
scalar Mflow;
scalar Mpow;
scalar transitionRpm(15);
scalar K;

// Calculate forces on patches and moment about axis
dictionary forcesDict;
forcesDict.add("patches", patches);
forcesDict.add("rhoName", "rhoInf");
forcesDict.add("rhoInf", rhoInf);
forcesDict.add("CofR", CofR);

forces F("forces", *this, forcesDict);
forces::forcesMoments fm = F.calcForcesMoment();
scalar Mpress = fm.second().first().z();
scalar Mtau = fm.second().second().z();
Mflow = fm.second().first().z()+fm.second().second().z();

// Apply moment if the rotational speed is high enough (because rpm=0 => Mpow=Inf)
if (abs(rpmOld) > abs(cutoffRpm))
{
    P = Pnom;
}
// Linear increase in the transition region
else if ( abs(rpmOld) > (cutoffRpm-transitionRpm) && abs(rpmOld) < cutoffRpm )
{
    K = Pnom / (cutoffRpm-transitionRpm);
    P = rpmOld * K;
}
// P=0 for low rpm:s
else 
{
    P = 0;
}

Mpow = Pnom/(rpmNom*Pi/30);

// Calculate angular acceleration and rpm (with accleration limiter)
alpha = ( Mflow-Mpow )/Iz;
if (abs(alpha) > 5)
{
    if (pos(alpha))
    {
        alpha = 5;
    }
    else if (neg(alpha))
    {
        alpha = -5;
    }
}

rpm = rpmOld + ( (30/Pi)*alpha )*time().deltaT().value();
