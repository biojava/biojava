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

import java.io.Serializable;

/**
 * A simple bond -- it stores information about two atoms as well as information
 * about its bond order.
 *
 * @author Jules Jacobsen <jacobsen@ebi.ac.uk>
 * @author Ulysse Carion
 */
public interface Bond extends Serializable {

	/**
	 * Gets atom 'A' of this bond. There is no meaning to which atom is 'A' and
	 * which is 'B'; the atoms are labeled 'A' or 'B' based on the order in
	 * which they are passed to this class's constructor.
	 *
	 * @see #getAtomB()
	 * @return one of the two atoms in this bond
	 */
	public Atom getAtomA();

	/**
	 * Gets atom 'B' of this bond. There is no meaning to which atom is 'A' and
	 * which is 'B'; the atoms are labeled 'A' or 'B' based on the order in
	 * which they are passed to this class's constructor.
	 *
	 * @see #getAtomA()
	 * @return one of the two atoms in this bond
	 */
	public Atom getAtomB();

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
	public Atom getOther(Atom exclude);

	/**
	 * Gets the bond order of this bond. A return value of '1' corresponds to a
	 * single bond, '2' to a double bond, etc.
	 *
	 * @return this bond's bond order
	 */
	public int getBondOrder();

	/**
	 * Gets the distance between the two atoms of this bond.
	 * <p>
	 * This distance is calculated by {@link Calc#getDistance(Atom, Atom)}, but
	 * this method will suppress the empty threat of a
	 * {@link StructureException} that method makes.
	 *
	 * @return the distance between the two atoms of this bond.
	 */
	public double getLength();
}
