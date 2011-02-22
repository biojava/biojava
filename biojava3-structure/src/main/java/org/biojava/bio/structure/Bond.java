/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.bio.structure;

/**
 * Simple bond  - it's just an edge/vertex for two {@link Atom} nodes.
 * @author Jules Jacobsen <jacobsen@ebi.ac.uk>
 */
public class Bond {

    private double length;
    private BondType type;
    private Atom atomA;
    private Atom atomB;

    public Bond() {
    }

    public Bond(double length, BondType type, Atom atomA, Atom atomB) {
        this.length = length;
        this.type = type;
        this.atomA = atomA;
        this.atomB = atomB;
    }

    @Override
    public String toString() {
        return "Bond{" + "length=" + length + " type=" + type + " atomA=" + atomA + " atomB=" + atomB + '}';
    }
}
