/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.bio.structure;

/**
 * A simple bond -- it stores information about two atoms as well as information
 * about the type of bond it is and its bond order.
 * 
 * @see BondType
 * @author Jules Jacobsen <jacobsen@ebi.ac.uk>
 * @author Ulysse Carion
 */
public class Bond {
	private BondType type;
	private Atom atomA;
	private Atom atomB;
	private int bondOrder;

	/**
	 * Constructs a new bond from a pair of atoms, the type of bond this is, and
	 * the bond order of the bond between them.
	 * <p>
	 * Note that when using this constructor, creating a Bond with some atoms A
	 * and B is not equivalent to forming a bond between A and B; to do this,
	 * one must either call {@link #addSelfToAtoms()} or use the constructor
	 * {@link #Bond(Atom, Atom, BondType, int, boolean)}, with the last argument
	 * set to 'true'.
	 * 
	 * @param atomA
	 *            one of the atoms in this bond
	 * @param atomB
	 *            the other atom in this bond
	 * @param type
	 *            this bond's type
	 * @param bondOrder
	 *            the bond order of this bond
	 */
	public Bond(Atom atomA, Atom atomB, BondType type, int bondOrder) {
		this(atomA, atomB, type, bondOrder, false);
	}

	/**
	 * Constructs a new bond from a pair of atoms, the type of bond this is, and
	 * the bond order of the bond between them.
	 * 
	 * @param atomA
	 *            one of the atoms in this bond
	 * @param atomB
	 *            the other atom in this bond
	 * @param type
	 *            this bond's type
	 * @param bondOrder
	 *            the bond order of this bond
	 * @param addSelfToAtoms
	 *            if set to true, this bond, once created, will automatically
	 *            add itself to atomA and atomB's bond lists.
	 */
	public Bond(Atom atomA, Atom atomB, BondType type, int bondOrder,
			boolean addSelfToAtoms) {
		this.type = type;
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
	 */
	// TODO first check if those bonds haven't been made already
	public void addSelfToAtoms() {
		atomA.getBonds().add(this);
		atomB.getBonds().add(this);
	}

	/**
	 * Gets atom 'A' of this bond. There is no meaning to which atom is 'A' and
	 * which is 'B'; the atoms are labeled 'A' or 'B' based on the order in
	 * which they are passed to this class's constructor.
	 * 
	 * @see #getAtomB()
	 * @return one of the two atoms in this bond
	 */
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
	 * @param exclude
	 *            the atom of the bond to not return
	 * @throws IllegalArgumentException
	 *             if the passed atom is not in this bond
	 * @return the atom in this bond that was not passed as an argument
	 */
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
	 * Gets the BondType of this Bond.
	 * 
	 * @return this bond's BondType
	 */
	public BondType getType() {
		return type;
	}

	/**
	 * Gets the bond order of this bond. A return value of '1' corresponds to a
	 * single bond, '2' to a double bond, etc.
	 * 
	 * @return this bond's bond order
	 */
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
	public double getLength() {
		try {
			return Calc.getDistance(atomA, atomB);
		} catch (StructureException e) {
			return -1; // this will never happen.
		}
	}

	@Override
	public String toString() {
		return "Bond [type=" + type + ", atomA=" + atomA + ", atomB=" + atomB
				+ ", bondOrder=" + bondOrder + "]";
	}
}
