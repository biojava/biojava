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

import org.biojava.nbio.structure.io.mmcif.model.ChemComp;

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

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
 *
 */
public interface Group extends Serializable {

	/** Group property key for secondary structure annotation */
	public static final String SEC_STRUC = "secstruc";

	/**
	 * Get number of atoms.
	 * @return number of atoms of this Group
	 */
	public int size();

	/**
	 * Return true or false, depending if this group has 3D coordinates or not.
	 *
	 * @return true if Group has 3D coordinates
	 */
	public boolean has3D ();

	/**
	 * Flag if group has 3D data .
	 *
	 * @param flag  true to set flag that this Group has 3D coordinates
	 */
	public void setPDBFlag(boolean flag);

	/**
	 * Get Type of group, one of {@link GroupType#AMINOACID}, {@link GroupType#HETATM}
	 * or {@link GroupType#NUCLEOTIDE}
	 *
	 * @return a String representing the type value
	 */
	public GroupType getType();

	/**
	 * Add an atom to this group.
	 *
	 * @param atom  an Atom object
	 */
	public void addAtom(Atom atom);

	/**
	 * Get list of atoms.
	 *
	 * @return a List object representing the atoms
	 * @see #setAtoms(List)
	 */
	public List<Atom> getAtoms() ;


	/**
	 * Set the atoms of this group.
	 * @see {@link Atom}
	 * @param atoms a list of atoms
	 */
	public void setAtoms(List<Atom> atoms);

	/** Remove all atoms from this group.
	 *
	 */
	public void clearAtoms();

	/**
	 * Get an atom given its PDB name.
	 * Beware that some PDB atom names are ambiguous (e.g. CA, which means C-alpha or Calcium),
	 * ambiguities should not occur within the same group though. To solve these ambiguities
	 * one would need to check the atom returned for the required element with {@link Atom#getElement()}
	 *
	 * @param name  a trimmed String representing the atom's PDB name, e.g. "CA"
	 * @return an Atom object or null if no such atom exists within this group
	 */
	public Atom getAtom(String name) ;


	/**
	 * Get at atom by position.
	 *
	 * @param position  an int
	 * @return an Atom object or null if no Atom exists for given position
	 */
	public Atom getAtom(int position) ;

	/**
	 * Tell whether a particular atom exists within this group.
	 * Beware that some PDB atom names are ambiguous (e.g. CA, which means C-alpha or Calcium),
	 * ambiguities should not occur within the same group though.
	 *
	 * @param name  a trimmed String representing the atom's PDB name, e.g. "CA"
	 * @return true if Atom with name exists within this group
	 */
	public boolean hasAtom(String name);

	/**
	 * Get the PDB 3-letter name for this group. (e.g. ALA)
	 *
	 * @return a String representing the PDBName value
	 * @see #setPDBName
	 */
	public String getPDBName();

	/**
	 * Set the PDB 3-letter name for this group. (e.g. ALA)
	 *
	 * @param s  a String specifying the PDBName value
	 * @see #getPDBName
	 */
	public void setPDBName(String s) ;


	/**
	 * Calculate if this group has all atoms required for an amino acid backbone.
	 * This allows to include chemically modified amino acids that
	 * are labeled hetatoms into some computations, the usual way
	 * to identify if a group is an amino acid is {@link #getType()}
	 * <p>
	 * amino atoms are : N, CA, C, O
	 * </p>
	 *
	 * Example: 1DW9 chain A first group is a Selenomethionine, provided as HETATM, but here returns true.
	 * <pre>
	 * HETATM    1  N   MSE A   1      11.720  20.973   1.584  0.00  0.00           N
	 * HETATM    2  CA  MSE A   1      10.381  20.548   1.139  0.00  0.00           C
	 * HETATM    3  C   MSE A   1       9.637  20.037   2.398  0.00  0.00           C
	 * HETATM    4  O   MSE A   1      10.198  19.156   2.985  0.00  0.00           O
	 * HETATM    5  CB  MSE A   1      10.407  19.441   0.088  0.00  0.00           C
	 * </pre>
	 *
	 * @return true if all Atoms required for an AminoAcid are available (N, CA, C, O)
	 * @see #getType
	 */
	public boolean hasAminoAtoms() ;


	/**
	 * Check if this group is a polymeric group, from the definition in Chemical Component Dictionary
	 *
	 * @return true if a polymeric group
	 */
	public boolean isPolymeric();


	/**
	 * Check if this group is an aminoacid group, from the definition in Chemical Component Dictionary
	 *
	 * @return true if an amino acid
     */
	public boolean isAminoAcid();


	/**
	 * Check if this group is a nucleotide group, from the definition in Chemical Component Dictionary
	 *
	 * @return true if a nucleotide
     */
	public boolean isNucleotide();



	/**
	 * Properties of this amino acid. Currently available properties are:
	 * phi
	 * psi
	 * secstruc
	 *
	 * @param properties  a Map object specifying the properties value
	 * @see #getProperties
	 */
	public void setProperties(Map<String,Object> properties) ;

	/**
	 * Return properties.
	 * @see #setProperties
	 *
	 * @return a HashMap object representing the properties value
	 */
	public Map<String,Object> getProperties() ;

	/**
	 * Set a single property .
	 *
	 * @param key    a String
	 * @param value  an Object
	 * @see #getProperty
	 */
	public void setProperty(String key, Object value) ;

	/**
	 * Get a single property .
	 *
	 * @param key  a String
	 * @return an Object
	 * @see #setProperty
	 */
	public Object getProperty(String key) ;

	/**
	 * Get an Atom Iterator.
	 *
	 * @return an Iterator object
	 */
	public Iterator<Atom> iterator() ;


	/**
	 * Returns and identical copy of this Group object .
	 * @return  and identical copy of this Group object
	 */
	public Object clone();

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
	 * Returns a dynamically created ResidueNumber for the group - this
	 * contains the chainId, resNum and insCode of the group.
	 * @see ResidueNumber
	 * @return ResidueNumber for the group.
	 * @since 3.0
	 */
	public ResidueNumber getResidueNumber();


	/**
	 * Sets the ResidueNumber for this Group
	 *
	 * @param residueNumber the PDB residueNumber
	 */
	public void setResidueNumber(ResidueNumber residueNumber);

	/**
	 * Utility method to temporarily set a chainID for a group, if a parent chain object does not exist yet.
	 * Not recommended for general use other than parsing.
	 *
	 * @param chainId
	 * @param residueNumber
	 * @param iCode
	 */
	public void setResidueNumber(String chainId, Integer residueNumber, Character iCode);

	/**
	 * Utility method for returning the chainId of the Group or null if no
	 * Chain has been set. This is equivalent to calling getChain().getId()
	 * 
	 * Prior to version 5.0 this method returned the chain name.
	 * @since 3.0
	 * @return  the ID of the chain
	 */
	public String getChainId();

	/**
	 * Set the Chemical Component that closer describes this group.
	 *
	 * @param cc the chemical component
	 */
	public void setChemComp(ChemComp cc);

	/**
	 * Get the chemical component that closer describes this group. If the information does not exist yet, fetches the information from PDB web site.
	 *
	 * @return the Chemical Component definition for this Group.
	 */
	public ChemComp getChemComp();


	/**
	 * Check if this group has alternate location groups.
	 *
	 * @return boolean flag if there are alternate locations.
	 * @see #getAltLocs()
	 */
	public boolean hasAltLoc();


	/**
	 * Get the list of other alternate location groups.
	 * <p>
	 * The main group (this group) will contain the first altloc (be it the default '.' or 'A' or a mix of '.' and 'A').
	 * <p>
	 * This method will return the altloc groups that are not the main group, e.g.:
	 *
	 * <li> if '.' (default), 'A' and 'B' altlocs are present in file, the main group will contain
	 * the default '.' and this method will return 2 altloc groups
	 * </li>
	 *
	 * <li> if 'A' and 'B' are present in file without a default '.' group, then the main group will contain the 'A'
	 * location whilst this method will return only 1 altloc group with the 'B' location
	 * </li>
	 *
	 * <p>
	 * Note that atoms with the default altloc (.) are included in all groups. Atoms with other altlocs (typically A, B, etc)
	 * will be sorted into groups by altloc.
	 * <p>
	 * Thus it can happen that an altloc group duplicate the contents of the main group.
	 *
	 * @return List of other groups that are on alternate locations
	 */
	public List<Group> getAltLocs();

	/**
	 * Add a group that is an alternate location for this group.
	 *
	 * @param g the altloc group to add
	 */
	public void addAltLoc(Group g);

	/**
	 * Determines if this group is water.
	 *
	 * @see GroupType#WATERNAMES
	 * @return true if it's water, false otherwise.
	 */
	public boolean isWater();

	/**
	 * Gets the alternate location group to this group that has the alt-loc character code passed.
	 *
	 * @param altLoc the alternate location code of the group desired
	 * @return the alternate location group if found, or null otherwise
	 */
	public Group getAltLocGroup(Character altLoc);


	/**
	 * Attempts to reduce the memory imprint of this group by trimming
	 * all internal Collection objects to the required size.
	 *
	 */
	public void trimToSize();

	/**
	 * Function to get the Group as an MDL molblock
	 * @return the string of the MDL molblock
	 */
	public String toSDF();
	
	/**
	 * Tells whether the group is annotated as HETATM in the file.
	 * To be used only at parsing time to be able to infer that a
	 * polymeric group is in a ligand chain or not.
	 * @return
	 */
	public boolean isHetAtomInFile();
	
	/**
	 * Sets the field isHetAtomInFile which is intented only for 
	 * helping in infering if a polymeric group is in a ligand chain
	 * or in a polymeric chain.
	 * @param isHetAtomInFile
	 */
	public void setHetAtomInFile(boolean isHetAtomInFile);
}
