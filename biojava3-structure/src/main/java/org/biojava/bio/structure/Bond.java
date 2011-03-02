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
    private Group groupA;
    private Atom atomA;
    private Group groupB;
    private Atom atomB;

    public Bond() {
    }

    public Atom getAtomA() {
        return atomA;
    }

    public Atom getAtomB() {
        return atomB;
    }

    public Group getGroupA() {
        return groupA;
    }

    public Group getGroupB() {
        return groupB;
    }

    public double getLength() {
        return length;
    }

    public BondType getType() {
        return type;
    }

    
    public Bond(double length, BondType type, Group resA, Atom atomA, Group resB, Atom atomB) {
        this.length = length;
        this.type = type;
        this.groupA = resA;
        this.atomA = atomA;
        this.groupB = resB;
        this.atomB = atomB;
    }

    @Override
    public String toString() {
        return "Bond{" + "length=" + length + " type=" + type + " resA=" + groupA.getPDBName() + " " + groupA.getResidueNumber() + " atomA=" + atomA + " resB=" + groupB.getPDBName() + " " + groupB.getResidueNumber() + " atomB=" + atomB + '}';
    }


}
