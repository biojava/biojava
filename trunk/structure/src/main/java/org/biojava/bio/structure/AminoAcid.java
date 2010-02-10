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

import java.util.Map;


/**
 * <p>
 * A {@link Group} that represents an AminoAcid.
 * </p>
 * 
 * <p>
 * In PDB files information on AminoAcids can be observed in the SEQRES and in the ATOM records.
 * Since frequently coordinates for some of the amino acids are missing, during parsing of the PDB
 * files the SEQRES and the ATOM records are aligned and a whenever possible the AminoAcid objects 
 * are combined. Access to the SEQRES and ATOM sequence is possible through the {@link Chain} object.
 * It is possible to distinguish between SEQRES and ATOM derived AminoAcids by {@link #getRecordType()}.
 * </p>
 * 
 *<p>
 *  AminoAcid inherits most from {@link HetatomImpl }.  Adds a few AminoAcid
 *  specific methods.
 *  </p>
 *  
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 * 
 */
public interface AminoAcid extends Group {
	
	/** Field to distribguish AminoAcids that have been created from SEQRES records and ATOM records.
	 * 
	 */
	public static final String ATOMRECORD = "ATOM";
	
	/** Field to distribguish AminoAcids that have been created from SEQRES records and ATOM records.
	 * 
	 */
	public static final String SEQRESRECORD = "SEQRES";
	
    /**
     * Specifies the secondary structure as a Map.
     *
     * @param secstr  a Map object specifying the sec struc 
     * @see #getSecStruc
     */
    public void setSecStruc(Map<String,String> secstr) ;
    
    /** get secondary structure data .
     *
     * @return a Map object representing the sec struc value
     *
     * @see #setSecStruc
     */
    public Map<String,String> getSecStruc() ;

    /** get N atom.
     *
     * @return an Atom object
     * @throws StructureException ...
    */
    public Atom getN()  throws StructureException ;

    /** get CA atom.
     * @return an Atom object
     * @throws StructureException ...
     */
    public Atom getCA() throws StructureException ;

    /** get C atom.
     * @return an Atom object
     * @throws StructureException ...
     */
    public Atom getC()  throws StructureException ;

    /** get O atom.
     * @return an Atom object
     * @throws StructureException ...
     */
    public Atom getO()  throws StructureException ;

    /** get CB atom.
     * @return an Atom object
     * @throws StructureException ...
     */
    public Atom getCB() throws StructureException ;

    


    /** returns the name of the AA, in single letter code.
     *
     * @return a Character object representing the amino type value
     * @see #setAminoType
     */
    public  Character getAminoType() ;

    /** set the name of the AA, in single letter code .
     *
     * @param aa  a Character object specifying the amino type value
     * @see #getAminoType
     */
    public void setAminoType(Character aa) ;

	/** Allows to distinguish between amino acids that are provided
		as ATOM records and a SEQRES records.
	 * @param recordName either ATOMRECORD or SEQRESRECORD
	 * @see #getRecordType()
	*/
	public void setRecordType(String recordName);

	/** Allows to distinguish between amino acids that are provided
		as ATOM records and a SEQRES records.
	 * @return the origin of this amino acid (ATOM or SEQRES records)
	 * @see #setRecordType(String)
	*/
    public String getRecordType();

    /** string representation. */
    public String toString() ;

}
