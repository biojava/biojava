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
package org.biojava.nbio.structure.cluster;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureTools;

/**
 * A Subunit consists of a set of residues from a Structure, which may
 * correspond to an entire Chain, a Domain, or any subset or combination of
 * residues from them.
 * <p>
 * There are two fundamental requirements. First, the residues have to be
 * sequential and connected in the original Structure (with the exception of
 * missing residues or loops in between). All the residues are of the same type.
 * These requirements are not checked when constructing the Object, but they are
 * assumed to always hold.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class Subunit {

	private Atom[] reprAtoms;
	private ProteinSequence sequence = null;

	/**
	 * A Subunit is uniquely defined by the coordinates of the representative
	 * Atoms of its residues, in sequential order.
	 * 
	 * @param repAtoms
	 *            representative Atoms
	 */
	public Subunit(Atom[] reprAtoms) {

		this.reprAtoms = reprAtoms;
	}

	/**
	 * Get all the representative Atoms of the Subunit. These Atoms are used for
	 * clustering and displaying the Subunit.
	 * 
	 * @return representative Atom[]
	 */
	public Atom[] getRepresentativeAtoms() {
		return reprAtoms;
	}

	/**
	 * The size of a Subunit is the number of residues it contains.
	 * 
	 * @return the size of the Subunit
	 */
	public int size() {
		return reprAtoms.length;
	}

	/**
	 * Get the protein sequence of the Subunit as String.
	 * 
	 * @return protein sequence String
	 */
	public String getProteinSequenceString() {

		if (sequence != null)
			return sequence.toString();

		StringBuilder builder = new StringBuilder();
		for (Atom a : reprAtoms)
			// This method preferred over getChemComp.getOneLetterCode because
			// it returns always X for Unknown residues
			builder.append(StructureTools.get1LetterCode(a.getGroup()
					.getPDBName()));

		return builder.toString();
	}

	/**
	 * Get the protein sequence of the Subunit.
	 * 
	 * @return sequence ProteinSequence
	 * @throws CompoundNotFoundException
	 */
	public ProteinSequence getProteinSequence()
			throws CompoundNotFoundException {

		if (sequence == null)
			sequence = new ProteinSequence(getProteinSequenceString());

		return sequence;
	}

	@Override
	public String toString() {
		return "Subunit [Size=" + size() + ", Sequence=" + sequence + "]";
	}

}
