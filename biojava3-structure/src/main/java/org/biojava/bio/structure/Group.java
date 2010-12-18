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
package org.biojava.bio.structure;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.io.PDBParseException;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;


/**
 *
 * This is the data structure for a single Group of atoms.  A protein
 * sequence ({@link Chain} in PDB file) is represented as a list of Groups.
 * There are 3 types of Groups:
 *
 * <ul>
 * <li>{@link AminoAcid}</li>
 * <li>{@link HetatomImpl Hetatom}</li>
 * <li>{@link NucleotideImpl Nucleotide}</li>
 * </ul>
 *
 *
 * @see HetatomImpl
 * @see AminoAcidImpl
 * @see NucleotideImpl
 * @author Andreas Prlic
 * @author Horvath Tamas
 * @since 1.4
 * @version %I% %G%
 */
public interface Group {


	/* returns and identical copy of this Group .
     public Object clone() ;
	 */

	/**
	 * Return the PDBcode (residue number + insertion code ) of this group.
	 *
	 * The residue number is treated as a String for the following reasons:
	 *
	 * Every amino acid in a PDB file is
	 * identified uniquely by 3 things: The chain ID, the residue number and
	 * the insertion code. To make sure one does not forget about the
	 * insertion code, in BioJava it is appended to the residue number.
	 *
	 * To add to this, residue numbers can be negative, non-consecutive and
	 * also non-sequential. As such it is often easiest, to treat them as
	 * public identifiers and within your own code work with  the internal
	 * atom or group positions...
	 *
	 * @see #setPDBCode
	 * @return a String representing the PDBCode value
	 * @deprecated replaced by {@link #getResidueNumber()}
	 */
	@Deprecated
	public String getPDBCode();

	/**
	 * Specifies the PDBCode (residue number + insertion code) value.
	 *
	 * @param pdbcode  a String specifying the PDBCode value
	 * @see #getPDBCode
	 * @deprecated replaced by {@link #setResidueNumber(ResidueNumber)}
	 */
	@Deprecated
	public void setPDBCode(String pdbcode);


	/** getnumber of atoms.
	 *  @return number of atoms of this Group
	 */
	public int size();

	/**
	 *  returns true or false, depending if this group has 3D coordinates or not.
	 *
	 * @return true if Group has 3D coordinates
	 */
	public boolean has3D ();

	/** flag if group has 3D data .
	 *
	 * @param flag  true to set flag that this Group has 3D coordinates
	 */
	public void setPDBFlag(boolean flag);

	/**
	 * get Type of group, e.g. amino, hetatom, nucleotide.
	 *
	 *
	 * @return a String representing the type value
	 */
	public String getType();

	/** add an atom to this group.
	 *
	 * @param atom  an Atom object
	 */
	public void addAtom(Atom atom);

	/** Get list of atoms.
	 *
	 * @return an List object representing the atoms
	 * @see #setAtoms(List)
	 */
	public List<Atom> getAtoms() ;


	/** Set the atoms of this group.
	 * @see org.biojava.bio.structure.Atom
	 * @param atoms a list of atoms
	 */
	public void setAtoms(List<Atom> atoms);

	/** Remove all atoms from this group.
	 *
	 */
	public void clearAtoms();

	/** Get an atom.  Throws StructureException if atom not found.
	 *
	 * @param name  a String
	 * @return an Atom object
	 * @throws StructureException if atom not found.
	 */
	public Atom getAtom(String name) throws StructureException;


	/**  Get an atom by the full PDB name e.g. " N  " for N. Throws StructureException if atom not found.
	 * @param pdbName  a String
	 * @return an Atom object
	 * @throws StructureException ...
	 */
	public  Atom getAtomByPDBname(String pdbName) throws StructureException;

	/** Get at atom by position.
	 *
	 * @param position  an int
	 * @return an Atom object
	 * @throws StructureException if not atom at this position
	 */
	public Atom getAtom(int position) throws StructureException;

	/** Teturns flag whether a particular atom is existing within this group .
	 *
	 * @param name  a String ...
	 * @return true if Atom with name is existing within this group
	 */
	public boolean hasAtom(String name);

	/** Get the PDB 3 character name for this group. (e.g. ALA)
	 *
	 * @return a String representing the PDBName value
	 * @see #setPDBName
	 */
	public String getPDBName();

	/** Set the PDB 3 letter name for this group. (e.g. ALA)
	 *
	 * @param s  a String specifying the PDBName value
	 * @throws PDBParseException ...
	 * @see #getPDBName
	 */
	public void setPDBName(String s) throws PDBParseException;


	/** calculate if a groups has all atoms required for an amino acid.
     this allows to include chemically modified amino acids that
     are labeled hetatoms into some computations ... the usual way
     to identify if a group is an amino acid is getType() !

     <p>
     amino atoms are : N, CA, C, O, CB
     GLY does not have CB (unless we would calculate some artificially
     </p>

     Example: 1DW9 chain A first group is a Selenomethionine, provided as HETATM, but here returns true.
     <pre>
     HETATM    1  N   MSE A   1      11.720  20.973   1.584  0.00  0.00           N
     HETATM    2  CA  MSE A   1      10.381  20.548   1.139  0.00  0.00           C
     HETATM    3  C   MSE A   1       9.637  20.037   2.398  0.00  0.00           C
     HETATM    4  O   MSE A   1      10.198  19.156   2.985  0.00  0.00           O
     HETATM    5  CB  MSE A   1      10.407  19.441   0.088  0.00  0.00           C
     </pre>
	 *
	 * @return true if all Atoms required for an AminoAcid are available (N, CA, C, O, CB)
     @see #getType
	 */
	public boolean hasAminoAtoms() ;



	/** properties of this amino acid. currerntly available properties.
	 * are:
	 * phi
	 * psi
	 *
	 *
	 * @param properties  a Map object specifying the properties value
	 * @see #getProperties

	 */

	public void setProperties(Map<String,Object> properties) ;

	/** return properties.
	 * @see #setProperties
	 *
	 * @return a HashMap object representing the properties value
	 */
	public Map<String,Object> getProperties() ;

	/** set a single property .
	 *
	 * @param key    a String
	 * @param value  an Object
	 * @see #getProperty

	 */
	public void setProperty(String key, Object value) ;

	/** get a single property .
	 *
	 * @param key  a String
	 * @return an Object
	 * @see #setProperty
	 */
	public Object getProperty(String key) ;

	/** get an Atom Iterator.
	 *
	 * @return an Iterator object
	 */
	public Iterator<Atom> iterator() ;


	/** returns and identical copy of this Group object .
	 * @return  and identical copy of this Group object
	 */
	public Object clone();


	/** Set the back-reference (to its parent Chain).
	 *
	 * @param parent the parent Chain
	 * @see #setChain(Chain)
	 * @see #getChain()
	 */
	@Deprecated
	public void setParent(Chain parent) ;

	/** Returns the parent Chain of the Group.
	 *
	 * @return Chain the Chain object that contains the Group
	 * @see #setChain(Chain)
	 * @deprecated replaced by {@link #getChain()}
	 */
	@Deprecated
	public Chain getParent() ;

	/**
	 * Sets the back-reference to its parent Chain.
	 * @param chain the parent Chain
	 * @see #getChain()
	 * @since 3.0
	 */
	public void setChain(Chain chain);

	/**
	 * Returns the parent Chain of the Group.
	 *
	 * @return Chain the Chain object that contains the Group
	 * @see #setChain(Chain)
	 * @since 3.0
	 */
	public Chain getChain() ;

	/**
	 * returns a dynamically created ResidueNumber for the group - this
	 * contains the chainId, resNum and insCode of the group.
	 * @see ResidueNumber
	 * @return ResidueNumber for the group.
	 * @since 3.0
	 */
	public ResidueNumber getResidueNumber();

	
	/** sets the ResidueNumber for this Group
	 * 
	 * @param residueNumber the PDB residueNumber
	 */
	public void setResidueNumber(ResidueNumber residueNumber);

	/** Utility method to temporarily set a chainID for a group, if a parent chain object does not exist yet.
	 * Not recommended for general use other than parsing.
	 * 
	 * @param chainId
	 * @param residueNumber
	 * @param iCode
	 */
	public void setResidueNumber(String chainId, Integer residueNumber, Character iCode);

	/**
	 * Utility method for returning the chainId of the Group or null if no
	 * Chain has been set. This replaces the need to use the expression
	 * group.getChain().getId() 
	 * @since 3.0
	 * @return  the ID of the chain
	 */
	public String getChainId();

	/** Set the Chemical Component that closer describes this group.
	 * 
	 * @param cc the chemical component
	 */
	public void setChemComp(ChemComp cc);

	/** Get the chemical component that closer describes this group. If the information does not exist yet, fetches the information from PDB web site.
	 *  
	 * @return the Chemical Component definition for this Group.
	 */
	public ChemComp getChemComp();

	
	/** Test if this group has alternate locations.
	 * 
	 * @return boolean flag if there are alternate locations.
	 */
	public boolean hasAltLoc();
	
	
	/** Get the list of alternate locations.
	 * 
	 * @return List of other groups that are on alternate locations
	 */
	public List<Group> getAltLocs();
	
	/** Add a group that is an alternate location for this group.
	 * 
	 */
	public void addAltLoc(Group g);

}
