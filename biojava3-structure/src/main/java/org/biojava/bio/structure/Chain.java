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
 * Created on 25.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure;

import java.util.List;

import org.biojava3.core.sequence.template.Sequence;

/**
 * <p>
 * Defines the interface for a Chain. A Chain corresponds to a Chain in a PDB file.
 * A chain consists out of a list of {@link Group} objects. A Group can either be 
 * an {@link AminoAcid}, {@link HetatomImpl Hetatom} or {@link NucleotideImpl Nucleotide}.
 * </p>
 * 
 * <p>
 * The BioJava API provides access to both the ATOM and SEQRES records in a PDB file.
 * During parsing of a PDB file it aligns the ATOM and SEQRES groups and joins them.
 * The SEQRES sequence can be accessed via  {@link #getSeqResGroups()} and the
 * ATOM groups via {@link #getAtomGroups()}. Groups that have been observed
 * (i.e. they are in the ATOM records) can be detected by {@link Group}.has3D()
 *  </p>
 *
 * @author Andreas Prlic
 * @version %I% %G%
 * @since 1.4
 */
public interface Chain {
	
    /** returns an identical copy of this Chain. 
     * @return  an identical copy of this Chain 
     */
    public Object clone();

    /** add a group to the list of ATOM record group of this chain.
     * To add SEQRES records a more complex alignment between ATOM and SEQRES residues
     * is required, please see SeqRes2AtomAligner for more details on that.
     * @param group  a Group object    
     */
    public void addGroup(Group group);
    
    /** Get the ID used by Hibernate.
     * 
     * @return the ID used by Hibernate
     * @see #setId(Long)
     */
    public Long getId() ;

    /** Set the ID used by Hibernate.
     * 
     * @param id assigned by Hibernate
     * @see #getId()
     */ 
    public void setId(Long id) ;

	
    /** return the amino acid at position X.
     * @param position  an int
     * @return a Group object
     * @deprecated use getAtomGroup or getSeqResGroup instead
     * @see #getAtomGroup(int)
     * @see #getSeqResGroup(int)
     */
    public Group getGroup (int position);
	
    /** 
     * Return the Group at given position, 
     * from within Groups with observed density in the chain, i.e.
     * those with coordinates in ATOM and HETATMS (including waters) records.
     * @param position  an int
     * @return a Group object
     * @see #getAtomLength()
     * @see #getAtomGroups()
     * @see #getSeqResGroup(int)
     */
    public Group getAtomGroup (int position);
    
    /** 
     * Return the Group at given position, 
     * from within groups in the SEQRES records of the chain, i.e.
     * the aminoacids/nucleotides in the construct. 
     * @param position  an int
     * @return a Group object
     * @see #getSeqResLength()
     * @see #getSeqResGroups()
     * @see #getAtomGroup(int)
     */
    public Group getSeqResGroup (int position);
    
    
    /** return a List of all groups of a special type (e.g. amino,
     * hetatm, nucleotide).
     * @param type  a String
     * @return a List object
     * @deprecated use getAtomGroups or getSeqResGroups instead
     */
    public List<Group> getGroups (String type);

    /** return all groups of this chain.
     * @return a List of all Group objects of this chain
     * @deprecated use getAtomGroups or getSeqResGroups instead
     */
    public List<Group> getGroups ();
    
    /** 
     * Return all Groups with observed density in the chain, i.e.
     * those with coordinates in ATOM and HETATMS (including waters) records.
     * 
     * @return a List object representing the Groups of this Chain.
     * @see #setAtomGroups(List)
     * @see #getAtomLength()
     * @see #getSeqResGroups()
     */
    public List<Group> getAtomGroups();
    
    /** 
     * Set all Groups with observed density in the chain, i.e. 
     * those with coordinates in ATOM and HETATMs (including waters) records.
     * @param groups a List object representing the Groups of this Chain.
     * @see #getAtomGroups()
     */
    public void setAtomGroups(List<Group> groups);
    
    /** 
     * Return a List of all (observed) Groups of a special type, one of: {@link GroupType#AMINOACID},
     * {@link GroupType#HETATM} or {@link GroupType#NUCLEOTIDE}.
     * Note that if a standard aminoacid appears as a HETATM (because it is part of a ligand) then
     * it is still considered as {@link GroupType#AMINOACID} and not as {@link GroupType#HETATM}.
     * @param type  a String
     * @return a List object
     * @see #setAtomGroups(List)
     * @see #getSeqResGroups(String)
     */
    public List<Group> getAtomGroups (String type);

    /** get a group by its PDB residue numbering. if the PDB residue number is not know,
     * throws a StructureException.
     * 
     * @param pdbresnum the PDB residue number of the group
     * @return the matching group
     * @throws StructureException
     * @deprecated replaced by {@link #getGroupByPDB(ResidueNumber)}
     */
    @Deprecated
    public Group getGroupByPDB(String pdbresnum) throws StructureException;
    

    /** 
     * Get a group by its PDB residue numbering. If the PDB residue number is not known,
     * throws a StructureException.
     * 
     * @param resNum the PDB residue number of the group
     * @return the matching group
     * @throws StructureException
     */
    public Group getGroupByPDB(ResidueNumber resNum) throws StructureException;
    
    /** Get all groups that are located between two PDB residue numbers.
     * 
     * @param pdbresnumStart PDB residue number of start
     * @param pdbresnumEnd PDB residue number of end
     * @return Groups in between. or throws a StructureException if either start or end can not be found,
     * @throws StructureException
     * @deprecated replaced by {@link #getGroupsByPDB(ResidueNumber, ResidueNumber)}
     */
    @Deprecated
    public Group[] getGroupsByPDB(String pdbresnumStart, String pdbresnumEnd) throws StructureException;

    /** Get all groups that are located between two PDB residue numbers.
     * 
     * @param pdbresnumStart PDB residue number of start
     * @param pdbresnumEnd PDB residue number of end
     * @return Groups in between. or throws a StructureException if either start or end can not be found,
     * @throws StructureException
     */
    public Group[] getGroupsByPDB(ResidueNumber pdbresnumStart, ResidueNumber pdbresnumEnd) throws StructureException;

    
    
    /** Get all groups that are located between two PDB residue numbers. In contrast to getGroupsByPDB
     * this method call ignores if the exact outer groups are not found. This is useful e.g. when requesting the range
     * of groups as specified by the DBREF records - these frequently are rather inaccurate.
     * 
     * 
     * @param pdbresnumStart PDB residue number of start
     * @param pdbresnumEnd PDB residue number of end
     * @param ignoreMissing ignore missing groups in this range.
     * @return Groups in between. or throws a StructureException if either start or end can not be found,
     * @throws StructureException
     * @deprecated replaced by #{@link #getGroupsByPDB(ResidueNumber, ResidueNumber, boolean)}
     */
    @Deprecated
    public Group[] getGroupsByPDB(String pdbresnumStart, String pdbresnumEnd,boolean ignoreMissing) throws StructureException;

    
    /** Get all groups that are located between two PDB residue numbers. In contrast to getGroupsByPDB
     * this method call ignores if the exact outer groups are not found. This is useful e.g. when requesting the range
     * of groups as specified by the DBREF records - these frequently are rather inaccurate.
     * 
     * 
     * @param pdbresnumStart PDB residue number of start
     * @param pdbresnumEnd PDB residue number of end
     * @param ignoreMissing ignore missing groups in this range.
     * @return Groups in between. or throws a StructureException if either start or end can not be found,
     * @throws StructureException
     * 
     */
    public Group[] getGroupsByPDB(ResidueNumber pdbresnumStart, ResidueNumber pdbresnumEnd,boolean ignoreMissing) throws StructureException;
    
    
    /** get total length of chain, including HETATMs..
     * @return an int representing the length of the whole chain including HETATMs
     * @deprecated please use getAtomLength or getLengthSeqRes instead
     * @see #getAtomLength()
     * @see #getSeqResLength()
     */
    public int getLength();
    
    
    /** 
     * Return the number of Groups with observed density in the chain, i.e. 
     * those with coordinates in ATOM and HETATMs (including waters) records
     * 
     * @return the length
     * @see #getAtomGroup(int)
     * @see #getAtomGroups()
     * @see #getSeqResLength()) 
     */
    public int getAtomLength();
    
    /** 
     * Return the number of groups in the SEQRES records of the chain, i.e.
     * the number of aminoacids/nucleotides in the construct 
     * 
     * @return the length
     * @see #getSeqResGroup(int)
     * @see #getSeqResGroups()
     * @see #getAtomLength()
     */
    public int getSeqResLength();
    
    /** returns the length of the AminoAcids in the ATOM records of this chain.
     * note: not all amino acids need to have 3D coords, in fact in could be that none
     * has!    
     * @return an int representing the length of the AminoAcids in the ATOM records of the chain.
     * @deprecated use getAtomGroups("amino").size() instead.
     */

    public int getLengthAminos();

  
    /** 
     * Set the Compound
     * @param compound the Compound 
     * @see #getCompound()
    */
    public void setCompound(Compound compound);

    /** 
     * Return the Compound for this chain.
     * 
     * @return the Compound object 
     * @see #setCompound(Compound)
     */
    public Compound getCompound();
    
    /** 
     * Set the name of this chain (Chain id in PDB file ).
     * @param name  a String specifying the name value
     * @see #getChainID()
     */
    public void setChainID(String name);	

    
    
    /** 
     * Get the name of this chain (Chain id in PDB file ).
     * @return a String representing the name value
     * @see #setChainID(String) 
     */
    public String getChainID();
    
    
    /** If available, returns the internal chain ID that is used in mmCif files, otherwise null
     * 
     * @return String or null
     * @since 3.0.5
     */
    public String getInternalChainID();
    
    /** 
     * Set the internal chain ID that is used in mmCif files
     * 
     * @param internalChainID
     * @since 3.0.5
     */
    public void setInternalChainID(String internalChainID);
    
    /** string representation.  */
    public String toString();
	
    
    /** Convert the SEQRES groups of a Chain to a Biojava Sequence object.
     * 
     * @return the SEQRES groups of the Chain as a Sequence object.
     * @throws IllegalSymbolException 
     */
    public Sequence<?> getBJSequence()  ;
       
    /** 
     * Return the sequence of amino acids as it has been provided in the ATOM records.
     * 
     * @return amino acid sequence as string
     * @see #getSeqResSequence()
     */
    public String getAtomSequence();
    
    /**
     * Return the PDB SEQRES sequence as a one-letter sequence string.
	 * Non-standard residues are represented by an "X".
	 * @return one-letter PDB SEQRES sequence as string 
     * @see #getAtomSequence() 
     */
    public String getSeqResSequence();
    
    /** Set the Swissprot id of this chain.
     * @param sp_id  a String specifying the swissprot id value
     * @see #getSwissprotId()
     */
    public void setSwissprotId(String sp_id);

    /** Get the Swissprot id of this chain.
     * @return a String representing the swissprot id value
     * @see #setSwissprotId(String sp_id)
     */
    public String getSwissprotId() ;
    
    
    /** 
     * Return a List of all SEQRES groups of a special type, one of: {@link GroupType#AMINOACID},
     * {@link GroupType#HETATM} or {@link GroupType#NUCLEOTIDE}.
     * @param type  a String
     * @return an List object
     * @see #setSeqResGroups(List)
     * @see #getAtomGroups(String)
     */
    public List<Group> getSeqResGroups (String type);

    /** 
     * Return a list of all groups in SEQRES records of the chain, i.e.
     * the aminoacids/nucleotides in the construct. 
     * @return a List of all Group objects of this chain
     * @see #setSeqResGroups(List)
     * @see #getSeqResLength()
     * @see #getAtomGroups()
     */
    public List<Group> getSeqResGroups ();
 
    /** 
     * Set the list of SeqResGroups for this chain.
     * 
     * @param seqResGroups a List of Group objects that from the SEQRES groups of this chain.
     * @see #getSeqResGroups()
     */
    public void setSeqResGroups(List<Group> seqResGroups);

    /** Set the back-reference to its parent Structure.
     * @param parent the parent Structure object for this Chain
     * @see #getParent()
     *  
     */
    public void setParent(Structure parent) ; 
    
    /** Returns the parent Structure of this chain.
     * 
     * @return the parent Structure object
     * @see #setParent(Structure)
     */
    
    public Structure getParent() ;
    
    /** Get all groups that are not polymer groups and that are not solvent groups.
     *  Will automatically fetch Chemical Component files from the PDB web site, even if
     *  PDBFileReader.setLoadChemCompInfo(flag) has not been set to true. Otherwise the Ligands could not
     *  correctly be identified.
     * @return list of Groups that are ligands
     */
    public List<Group> getAtomLigands();
    
    
   
    
}
