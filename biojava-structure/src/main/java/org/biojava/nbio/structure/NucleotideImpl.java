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
 * Created on 05.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.nbio.structure;

/**
 * A nucleotide group is almost the same as a Hetatm group.
 * @see HetatomImpl
 * @see AminoAcidImpl
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */
public class NucleotideImpl extends HetatomImpl {

	private static final long serialVersionUID = -7467726932980288712L;
	/** this is a "nucleotide", a special occurance of a Hetatom. */
	public static final GroupType type = GroupType.NUCLEOTIDE;

	/*
	 * inherits most from Hetero and has just a few extensions.
	 */
	public NucleotideImpl() {
		super();

	}

	@Override
	public GroupType getType(){ return type;}


	@Override
	public String toString(){

		String str = "PDB: "+ pdb_name + " " + residueNumber +  " "+ pdb_flag;
		if (pdb_flag) {
			str = str + "atoms: "+atoms.size();
		}
		return str ;

	}

	/**
	 * Returns the O3' atom if present, otherwise null
	 * @return O3' atom or null
	 */
	public Atom getO3Prime() {

		return getAtom("O3'");

	}

	/**
	 * Returns the O5' atom if present, otherwise null
	 * @return O5' atom or null
	 */
	public Atom getO5Prime() {

		return getAtom("O5'");

	}

	/**
	 * Returns the P atom if present, otherwise null
	 * @return P atom or null
	 */
	public Atom getP() {

		return getAtom("P");

	}

	// note we need to implement a clone here, despite there's one in super class already,
	// that's due to issue https://github.com/biojava/biojava/issues/631 - JD 2017-01-21
	@Override
	public Object clone() {

		NucleotideImpl n = new NucleotideImpl();
		n.setPDBFlag(has3D());
		n.setResidueNumber(getResidueNumber());

		n.setPDBName(getPDBName());

		//clone atoms and bonds.
		cloneAtomsAndBonds(n);
		
		// copying the alt loc groups if present, otherwise they stay null
		if (getAltLocs()!=null && !getAltLocs().isEmpty()) {
			for (Group altLocGroup:this.getAltLocs()) {
				Group nAltLocGroup = (Group)altLocGroup.clone();
				n.addAltLoc(nAltLocGroup);
			}
		}
		
		if (chemComp!=null)
			n.setChemComp(chemComp);


		return n;
	}
}
