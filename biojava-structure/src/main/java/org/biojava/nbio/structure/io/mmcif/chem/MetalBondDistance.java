package org.biojava.nbio.structure.io.mmcif.chem;

/** A bean that contains cutoffs for correctly detecting metal bonds.
 * Definitions are in file bond_distance_limits.cif.gz
 *
 * Created by andreas on 6/9/16.
 */
public class MetalBondDistance {

    private String atomType1;
    private String atomType2;
    private float lowerLimit;
    private float upperLimit;

    public String getAtomType1() {
        return atomType1;
    }

    public void setAtomType1(String atomType1) {
        this.atomType1 = atomType1;
    }

    public String getAtomType2() {
        return atomType2;
    }

    public void setAtomType2(String atomType2) {
        this.atomType2 = atomType2;
    }

    public float getLowerLimit() {
        return lowerLimit;
    }

    public void setLowerLimit(float lowerLimit) {
        this.lowerLimit = lowerLimit;
    }

    public float getUpperLimit() {
        return upperLimit;
    }

    public void setUpperLimit(float upperLimit) {
        this.upperLimit = upperLimit;
    }

    @Override
    public String toString() {
        return "MetalBindDistance{" +
                "atomType1='" + atomType1 + '\'' +
                ", atomType2='" + atomType2 + '\'' +
                ", lowerLimit=" + lowerLimit +
                ", upperLimit=" + upperLimit +
                '}';
    }
}
