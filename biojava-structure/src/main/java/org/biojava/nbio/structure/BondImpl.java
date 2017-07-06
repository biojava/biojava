/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure;

import java.util.ArrayList;
import java.util.List;

/**
 * A simple bond -- it stores information about two atoms as well as information
 * about its bond order.
 *
 * @author Jules Jacobsen <jacobsen@ebi.ac.uk>
 * @author Ulysse Carion
 */
public class BondImpl implements Bond {

	private static final long serialVersionUID = 8836120946858134380L;
	private Atom atomA;
	private Atom atomB;
	private int bondOrder;

	/**
	 * Constructs a new bond from a pair of atoms and the bond order of the bond
	 * between them.
	 * <p>
	 * Note that by forming a bond between atoms 'A' and 'B' with this
	 * constructor, atoms 'A' and 'B' will be updated to have this bond in their
	 * list of bonds. If you do not want this automatic updating, instead use
	 * {@link #Bond(Atom, Atom, int, boolean)} with the
	 * <code>addSelfToAtoms</code> flag set to <code>false</code>.
	 *
	 * @param atomA one of the atoms in this bond
	 * @param atomB the other atom in this bond
	 * @param bondOrder the bond order of this bond
	 */
	public BondImpl(Atom atomA, Atom atomB, int bondOrder) {
		this(atomA, atomB, bondOrder, true);
	}

	/**
	 * Constructs a new bond from a pair of atoms and the bond order of the bond
	 * between them.
	 *
	 * @param atomA one of the atoms in this bond
	 * @param atomB the other atom in this bond
	 * @param bondOrder the bond order of this bond
	 * @param addSelfToAtoms if set to true, this bond, once created, will
	 * automatically add itself to atomA and atomB's bond lists. (If this
	 * argument is set to false, the list returned from {@link Atom#getBonds()}
	 * will not contain this bond.)
	 */
	public BondImpl(Atom atomA, Atom atomB, int bondOrder, boolean addSelfToAtoms) {
		this.atomA = atomA;
		this.atomB = atomB;
		this.bondOrder = bondOrder;

		if (addSelfToAtoms) {
			addSelfToAtoms();
		}
	}

	/**
	 * Adds this Bond to its atoms bond lists. If this method is not called,
	 * then the list returned from calling {@link Atom#getBonds()} will not
	 * include this bond.
	 * <p>
	 * If you created your Bond with the constructor
	 * {@link #Bond(Atom, Atom, int)}, this method has already been called for
	 * you and should not be called again.
	 */
	// TODO first check if those bonds haven't been made already
	private void addSelfToAtoms() {
		List<Bond> bonds = atomA.getBonds();
		if (bonds==null) {
			bonds = new ArrayList<Bond>(AtomImpl.BONDS_INITIAL_CAPACITY);
			atomA.setBonds(bonds);
		}

		boolean exists = false;
		for (Bond bond : bonds) {
			// TODO is it equals() what we want here, or is it == ? - JD 2016-03-02
			if (bond.getOther(atomA).equals(atomB)) {
				exists = true;
				break;
			}
		}
		if (!exists) {
			atomA.addBond(this);
			atomB.addBond(this);
		}
	}

	/**
	 * Gets atom 'A' of this bond. There is no meaning to which atom is 'A' and
	 * which is 'B'; the atoms are labeled 'A' or 'B' based on the order in
	 * which they are passed to this class's constructor.
	 *
	 * @see #getAtomB()
	 * @return one of the two atoms in this bond
	 */
	@Override
	public Atom getAtomA() {
		return atomA;
	}

	/**
	 * Gets atom 'B' of this bond. There is no meaning to which atom is 'A' and
	 * which is 'B'; the atoms are labeled 'A' or 'B' based on the order in
	 * which they are passed to this class's constructor.
	 *
	 * @see #getAtomA()
	 * @return one of the two atoms in this bond
	 */
	@Override
	public Atom getAtomB() {
		return atomB;
	}

	/**
	 * A utility method to get the other atom in a bond, given one of its atoms.
	 * If the atom passed is one of the atoms in this bond, then this method is
	 * essentially equivalent to saying
	 * <code>atom == bond.getAtomA() ? bond.getAtomB() : bond.getAtomA()</code>.
	 * <p>
	 * <i>Note:</i> Comparison of atoms in this method is done with
	 * <code>==</code>, not <code>equals</code>.
	 *
	 * @param exclude the atom of the bond to not return
	 * @throws IllegalArgumentException if the passed atom is not in this bond
	 * @return the atom in this bond that was not passed as an argument
	 */
	@Override
	public Atom getOther(Atom exclude) {
		if (exclude != atomA && exclude != atomB) {
			throw new IllegalArgumentException(
					"Atom to exclude is not in bond.");
		}

		if (exclude == atomA) {
			return atomB;
		} else {
			return atomA;
		}
	}

	/**
	 * Gets the bond order of this bond. A return value of '1' corresponds to a
	 * single bond, '2' to a double bond, etc.
	 *
	 * @return this bond's bond order
	 */
	@Override
	public int getBondOrder() {
		return bondOrder;
	}

	/**
	 * Gets the distance between the two atoms of this bond.
	 * <p>
	 * This distance is calculated by {@link Calc#getDistance(Atom, Atom)}, but
	 * this method will suppress the empty threat of a
	 * {@link StructureException} that method makes.
	 *
	 * @return the distance between the two atoms of this bond.
	 */
	@Override
	public double getLength() {

		return Calc.getDistance(atomA, atomB);

	}

	@Override
	public String toString() {
		return "Bond [atomA=" + atomA + ", atomB=" + atomB + ", bondOrder="
				+ bondOrder + "]";
	}



}
