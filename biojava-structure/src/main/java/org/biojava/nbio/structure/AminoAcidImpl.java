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
 * AminoAcid inherits most from Hetatom.  Adds a few AminoAcid
 * specific methods.
 * @author Andreas Prlic
 * @author Jules Jacobsen
 * @since 1.4
 * @version %I% %G%
 *
 */
public class AminoAcidImpl extends HetatomImpl implements AminoAcid {

	private static final long serialVersionUID = -6018854413829044230L;

	/** this is an Amino acid. type is "amino". */
	public static final GroupType type = GroupType.AMINOACID;

	/** IUPAC amino acid residue names
	 */
	private Character amino_char ;

	private String recordType; // allows to distinguish between AAs that have been created from SEQRES records and ATOM records

	/**
	 * inherits most from Hetero and has just a few extensions.
	 */
	public AminoAcidImpl() {
		super();

		amino_char = null;
		recordType = ATOMRECORD;
	}

	@Override
	public GroupType getType(){ return type;}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Atom getN()    {return getAtom("N");  }

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Atom getCA()   {
		// note CA can also be Calcium, but that can't happen in a standard aminoacid, so this should be safe
		return getAtom("CA");
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Atom getC()    {return getAtom("C");  }

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Atom getO()    {return getAtom("O");  }

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Atom getCB()   {return getAtom("CB"); }


	/**
	 * {@inheritDoc}
	 */
	@Override
	public  Character getAminoType() {
		return amino_char;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setAminoType(Character aa){
		amino_char  = aa ;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setRecordType(String recordName) {
		recordType = recordName;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getRecordType() {
		return recordType;
	}

	/** string representation. */
	@Override
	public String toString(){

		String str = "AminoAcid "+ recordType + ":"+ pdb_name + " " + amino_char +
				" " + residueNumber +  " "+ pdb_flag + " " + recordType  ;
		if (pdb_flag) {
			str = str + " atoms: "+atoms.size();
		}
		if (!getAltLocs().isEmpty())
			str += " has altLocs :" + getAltLocs().size();

		return str ;

	}
	/** set three character name of AminoAcid.
	 *
	 * @param s  a String specifying the PDBName value
	 * @see #getPDBName()
	 */
	@Override
	public void setPDBName(String s) {

		pdb_name =s ;

	}


	/** returns and identical copy of this Group object .
	 * @return  and identical copy of this Group object
	 */
	@Override
	public Object clone() {

		AminoAcidImpl n = new AminoAcidImpl();
		n.setPDBFlag(has3D());
		n.setResidueNumber(getResidueNumber());

		n.setPDBName(getPDBName());

		n.setAminoType(getAminoType());
		n.setRecordType(recordType);

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
