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
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.StructureTools;

/**
 * A Subunit consists of a set of residues from a Structure, which may
 * correspond to an entire Chain, a Domain, or any subset or combination of
 * residues from them. All the residues of a Subunit are of the same type.
 * <p>
 * The Subunit object can contain additional fields for identification and
 * annotation.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class Subunit {

	// Optional fields for Subunit annotation
	private String name;
	private Structure structure;
	private StructureIdentifier identifier;

	// Required fields for Subunit definition
	private Atom[] reprAtoms;
	private ProteinSequence sequence = null;

	/**
	 * A Subunit is solely defined by the coordinates of the representative
	 * Atoms of its residues. It can be identified with a StructureIdentifier
	 * and/or a name and stores a reference to the Structure from which the
	 * Atoms were obtained.
	 * 
	 * @param reprAtoms
	 *            representative Atoms. It cannot be null or empty
	 * @param name
	 *            String field that identifies the Subunit. It can be null
	 * @param identifier
	 *            StructureIdentifier. It can be null
	 * @param structure
	 *            parent Structure object. It can be null
	 */
	public Subunit(Atom[] reprAtoms, String name,
			StructureIdentifier identifier, Structure structure) {

		if (reprAtoms == null)
			throw new IllegalArgumentException(
					"Representative Atom Array of the Subunit is null");
		if (reprAtoms.length==0)
			throw new IllegalArgumentException(
					"Representative Atom Array of the Subunit has 0 length");

		this.reprAtoms = reprAtoms;
		this.name = name;
		this.identifier = identifier;
		this.structure = structure;
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

	/**
	 * The Name of a Subunit is a free-text field, user defined.
	 * 
	 * @return the Subunit name
	 */
	public String getName() {
		return name;
	}

	/**
	 * The Name of a Subunit is a free-text field, user defined.
	 * 
	 * @param name
	 *            of the Subunit
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * The parent Structure from which the Subunit atoms were obtained.
	 * 
	 * @return Structure object
	 */
	public Structure getStructure() {
		return structure;
	}

	/**
	 * The standard identifier of the Subunit.
	 * 
	 * @return StructureIdentifier object
	 */
	public StructureIdentifier getIdentifier() {
		return identifier;
	}

	@Override
	public String toString() {
		return "Subunit [Name: " + name + ", Identifier: " + identifier
				+ ", Size:" + size() + ", Sequence:" + sequence + "]";
	}

}
