package ini.cx3d.spatialOrganization.factory;

import ini.cx3d.spatialOrganization.RationalJava;
import ini.cx3d.spatialOrganization.interfaces.Rational;
import ini.cx3d.swig.spatialOrganization.spatialOrganization;

import java.math.BigInteger;

/**
 * Factory that generates Rational object and enables quick switching between Java and native CPP implementation
 */
public class RationalFactory {

    private static boolean NATIVE = spatialOrganization.useNativeRational;

    public Rational create(long numerator, long denominator) {
        if (NATIVE) {
            return ini.cx3d.swig.spatialOrganization.Rational.create(numerator, denominator);
        } else {
            return new RationalJava(numerator, denominator);
        }
    }

    public Rational create(BigInteger numerator, BigInteger denominator) {
        if(NATIVE) {
            return ini.cx3d.swig.spatialOrganization.Rational.create(numerator, denominator);
        }else {
            return new RationalJava(numerator, denominator);
        }
    }

    public Rational create(double value) {
        if(NATIVE) {
            return ini.cx3d.swig.spatialOrganization.Rational.create(value);
        }else {
            return new RationalJava(value);
        }
    }
}
