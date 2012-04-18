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
 * Created on 26.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure;

import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.io.FileConvert;
import org.biojava.bio.structure.io.PDBFileReader;


/**
 *
 * Interface for a structure object. Provides access to the data of a PDB file.
 *
 * A structure object allows to access the PDB header information as well
 * as to the data from the ATOM records. The header information is
 * currently available through the following objects:
 * <ul>
 * <li>{@link PDBHeader}</li>
 * <li>{@link DBRef}</li>
 * <li>{@link Compound}</li>
 * </ul>
 *
 * The structure object provides access to the data from the ATOM records through
 * a hierarchy of sub-object:
 * <pre>
 * Structure
 *         |
 *         {@link Chain}
 *             |
 *             {@link Group}
 *                 |
 *                 {@link Atom}
 * </pre>
 *
 * For more documentation on how to work with the Structure API please
 * see <a href="http://biojava.org/wiki/BioJava:CookBook#Protein_Structure" target="_top">
 * http://biojava.org/wiki/BioJava:CookBook#Protein_Structure</a>
 *
 *
 *
 *
<hr/>
</hr>
 * <p>
 * Q: How can I get a Structure object from a PDB file?
 * </p>
 * <p>
 * A:
 </p>
 * <pre>
public {@link Structure} loadStructure(String pathToPDBFile){
		{@link PDBFileReader} pdbreader = new {@link PDBFileReader}();

		{@link Structure} structure = null;
		try{
			structure = pdbreader.getStructure(pathToPDBFile);
			System.out.println(structure);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return structure;
	}
 </pre>

<hr>
</hr>
<p>
Q: How can I calculate Phi and Psi angles of AminoAcids?
</p>
<p>
A:
</p>
<pre>
public void calcPhiPsi({@link Structure} structure){


		// get the first chain from the structure

		{@link Chain} chain  = structure.getChain(0);

		// A protein chain consists of a number of groups. These can be either
		// {@link AminoAcid}, {@link HetatomImpl Hetatom} or {@link NucleotideImpl Nucleotide} groups.
		//
		// Note: BioJava provides access to both the ATOM and SEQRES data in a PDB file.
		// since we are interested in doing calculations here, we only request the groups
		// from the ATOM records

		//  get the Groups of the chain that are AminoAcids.
		List<Group> groups = chain.getAtomGroups("amino");

		{@link AminoAcid} a;
		{@link AminoAcid} b;
		{@link AminoAcid} c ;

		for ( int i=0; i < groups.size(); i++){

			// since we requested only groups of type "amino" they will always be amino acids
			// Nucleotide and Hetatom groups will not be present in the groups list.

			b = ({@link AminoAcid})groups.get(i);

			double phi =360.0;
			double psi =360.0;

			if ( i > 0) {
				a = ({@link AminoAcid})groups.get(i-1) ;
				try {

					// the Calc class provides utility methods for various calculations on
					// structures, groups and atoms

					phi = {@link Calc}.getPhi(a,b);
				} catch ({@link StructureException} e){
					e.printStackTrace();
					phi = 360.0 ;
				}
			}
			if ( i < groups.size()-1) {
				c = ({@link AminoAcid})groups.get(i+1) ;
				try {
					psi = {@link Calc}.getPsi(b,c);
				}catch ({@link StructureException} e){
					e.printStackTrace();
					psi = 360.0 ;
				}
			}

			System.out.print(b.getPDBCode() + " " + b.getPDBName() + ":"  );

			System.out.println(String.format("\tphi: %+7.2f psi: %+7.2f", phi, psi));

		}
</pre>
<hr>
</hr>

 *
 *
 *
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */
public interface Structure extends Cloneable{


	/** returns an identical copy of this Structure object
	 *
	 * @return identical copy of this Structure object
	 */
	public Structure clone();

    /**
     * String representation of object.
     */
    public String toString();

    /**
     *
     * set PDB code of structure .
     *
     * @param pdb_id  a String specifying the PDBCode
     * @see #getPDBCode
     *
     */
    public void setPDBCode (String pdb_id) ;

    /**
     *
     * get PDB code of structure.
     *
     * @return a String representing the PDBCode value
     * @see #setPDBCode
     */
    public String  getPDBCode () ;

    /** set biological name of Structure .
     *
     * @param name  a String specifying the biological name of the Structure
     * @see #getName
     */
    public void setName(String name);

    /** get biological name of Structure.
     *
     * @return a String representing the biological name of the Structure
     * @see #setName
     */
    public String getName();

    /** set the Header data .
     *
     * @param h  a Map object specifying the header
     * @see #getHeader
     * @deprecated
     */
    public void setHeader(Map<String,Object> h) ;

    /** get Header data .
     *
     * @return a Map object representing the header value
     * @see #setHeader
     * @deprecated use getPDBHeader instead
     * @see #getPDBHeader
     */
    public Map<String,Object> getHeader() ;

    /**
       sets/gets an List of  Maps which corresponds to the CONECT lines in the PDB file:

       <pre>
       COLUMNS         DATA TYPE        FIELD           DEFINITION
       ---------------------------------------------------------------------------------
        1 -  6         Record name      "CONECT"
        7 - 11         Integer          serial          Atom serial number
       12 - 16         Integer          serial          Serial number of bonded atom
       17 - 21         Integer          serial          Serial number of bonded atom
       22 - 26         Integer          serial          Serial number of bonded atom
       27 - 31         Integer          serial          Serial number of bonded atom
       32 - 36         Integer          serial          Serial number of hydrogen bonded
       atom
       37 - 41         Integer          serial          Serial number of hydrogen bonded
       atom
       42 - 46         Integer          serial          Serial number of salt bridged
       atom
       47 - 51         Integer          serial          Serial number of hydrogen bonded
       atom
       52 - 56         Integer          serial          Serial number of hydrogen bonded
       atom
       57 - 61         Integer          serial          Serial number of salt bridged
       atom
       </pre>

       the HashMap for a single CONECT line contains the following fields:

       <li> atomserial (mandatory) : Atom serial number</li>
       <li> bond1 .. bond4 (optional): Serial number of bonded atom</li>
       <li> hydrogen1 .. hydrogen4 (optional):Serial number of hydrogen bonded atom</li>
       <li> salt1 .. salt2 (optional): Serial number of salt bridged atom</li>

       *
       * @param connections  a List object specifying the connections
       * @see #getConnections
    */
    public void setConnections(List<Map<String,Integer>> connections);

    /**
     * Returns the connections value.
     * @return a List object representing the connections value
     * @see #setConnections
     */
    public List<Map<String,Integer>> getConnections();

    /** return number of Chains in this Structure.
     * @return an int representing the number of Chains in this Structure
     */
    public int size() ;

    /** return number of chains of model.
    *
    * @param modelnr  an int specifying the number of the Model that should be used
    * @return an int representing the number of Chains in this Model
    */
    public int size(int modelnr);

    /** return number of models .
     * in this implementation also XRAY structures have "1 model", since
     * model is the container for the chains.
     * to test if a Structure is an NMR structure use @see isNMR ,
     * since this is based on the info in the PDB file header.
     *
     * @return an int representing the number of models in this Structure
     */
    public int nrModels() ;

    /** test if this structure is an nmr structure.
     *
     * @return true if this Structure has been resolved by NMR
     */
    public boolean isNmr() ;

    /** set NMR flag.
     *
     * @param nmr  true to declare that this Structure has been solved by NMR.
     */
    public void setNmr(boolean nmr);


    /** add a new model.
     *
     * @param model  a List object containing the Chains of the new Model
     */
    public void addModel(List<Chain> model);


    /** a convenience function if one wants to edit and replace the
     * models in a structure. allows to set (replace) the model at position
     * with the new List of Chains.
     * @param position starting at 0
     * @param model
     */
    public void setModel(int position, List<Chain> model);

    /** retrieve all Chains belonging to a model .
     * @see #getChains(int modelnr)
     *
     * @param modelnr  an int
     * @return a List object containing the Chains of Model nr. modelnr

     */
    public List<Chain> getModel(int modelnr);

    /** retrieve all chains - if it is a NMR structure will return the chains of the first model.
     * This is the same as getChains(0);
     * @see #getModel(int modelnr)
     * @see #getChains(int modelnr)
     *
     * @return a List object containing the Chains of Model nr. modelnr
     */
    public List<Chain> getChains();


    /** set the chains of a structure, if this is a NMR structure,
     * this will only set model 0.
     *
     * @see #setChains(int, List)
     *
     * @param chains the list of chains for this structure.
     */
    public void setChains(List<Chain> chains);

    /** retrieve all chains of a model.
     * @see #getModel
     *
     * @param modelnr  an int
     * @return a List object containing the Chains of Model nr. modelnr
     */
    public List<Chain> getChains(int modelnr);

    /** set the chains for a model
     * @param chains
     * @param modelnr
     */
    public void setChains( int modelnr, List<Chain> chains);

    /** add a new chain.
     *
     * @param chain  a Chain object
     */
    public void addChain(Chain chain);

    /** add a new chain, if several models are available.
     *
     * @param chain    a Chain object
     * @param modelnr  an int specifying to which model the Chain should be added
     */
    public void addChain(Chain chain, int modelnr);

    /** retrieve a chain by it's position within the Structure .
     *
     * @param pos  an int for the position in the List of Chains.
     * @return a Chain object
    */
    public Chain getChain(int pos);

    /** retrieve a chain by it's position within the Structure and model number.
     *
     * @param pos      an int
     * @param modelnr  an int
     * @return a Chain object
    */
    public Chain getChain( int modelnr, int pos);



    /** request a particular chain from a structure.
     * by default considers only the first model.
     * @param chainId the ID of a chain that should be returned
     * @return Chain the requested chain
     * @throws StructureException
     */
    public Chain findChain(String chainId)
    throws StructureException;


    /** check if a chain with the id chainId is contained in this structure.
     *
     * @param chainId the name of the chain
     * @return true if a chain with the id (name) chainId is found
     */
    public boolean hasChain(String chainId);

    /** request a particular chain from a particular model
     * @param modelnr the number of the model to use
     * @param chainId the ID of a chain that should be returned
     * @return Chain the requested chain
     * @throws StructureException
     */
    public Chain findChain(String chainId, int modelnr)
    throws StructureException;

    /** request a particular group from a structure.
    * by default considers only the first model in the structure.
    * @param chainId the ID of the chain to use
    * @param pdbResnum the PDB residue number of the requested group
    * @return Group the requested Group
    * @throws StructureException
    *
    */
    public  Group findGroup(String chainId, String pdbResnum)
    throws StructureException;

    /** request a particular group from a structure.
     * considers only model nr X. count starts with 0.
     * @param chainId the ID of the chain to use
     * @param pdbResnum the PDB residue number of the requested group
     * @param modelnr the number of the model to use
     * @return Group the requested Group
     * @throws StructureException
     */
     public  Group findGroup(String chainId, String pdbResnum, int modelnr)
     throws StructureException;


     /** request a chain by it's PDB code
      * by default takes only the first model
      *
      * @param chainId the chain identifier
      * @return the Chain that matches the chainID
      * @throws StructureException
      */
     public Chain getChainByPDB(String chainId)
         throws StructureException;

     /** request a chain by it's PDB code
      * by default takes only the first model
      *
      * @param chainId the chain identifier
      * @param modelnr request a particular model;
      * @return the Chain that matches the chainID in the model
      * @throws StructureException
      */
     public Chain getChainByPDB(String chainId, int modelnr)
         throws StructureException;


    /** create a String that contains the contents of a PDB file .
     *
     * @return a String that looks like a PDB file
     * @see FileConvert
     */
    public String toPDB();

    /** set the compounts
     *
     * @param molList
     */
    public void setCompounds(List<Compound>molList);

    /** get all the Compounds that are defined in the PDB Header
     *
     * @return a list of compound
     */
    public List<Compound> getCompounds();

    /** set the list of database references for this structure
     * @param dbrefs list of DBRef objects
     *
     *
     */
    public void setDBRefs(List<DBRef> dbrefs);

    /** get the list of database references
     *
     * @return list of DBRef objects
     */
    public List<DBRef> getDBRefs();

    /** request a particular compound by its id
     *
     * @param molId
     * @return a compound
     */
    public Compound getCompoundById(String molId);


    /** return the header information for this PDB file
     *
     * @return the PDBHeader object
     */
    public PDBHeader getPDBHeader();

    /**
     * return whether or not the entry has an associated journal article
     * or publication. The JRNL section is not mandatory and thus may not be
     * present.
     * @return flag if a JournalArticle has been found.
     */
    public boolean hasJournalArticle();

    /**
     * get the associated publication as defined by the JRNL records in a PDB
     * file.
     * @return a JournalArticle
     */
    public JournalArticle getJournalArticle();

    /**
     * set the associated publication as defined by the JRNL records in a PDB
     * file.
     * @param journalArticle
     */
    public void setJournalArticle(JournalArticle journalArticle);

    /** get the list of SSBonds as they have been defined in the PDB files
     *
     * @return a list of SSBonds
     */
    public List<SSBond> getSSBonds();

    /** set the list of SSBonds for this structure
     *
     * @param ssbonds
     */
    public void setSSBonds(List<SSBond> ssbonds);

    /** add a single SSBond to this structure
     *
     * @param ssbond
     */
    public void addSSBond(SSBond ssbond);

    /** the the header information for this PDB file
     *
     * @param header the PDBHeader object
     */
    public void setPDBHeader(PDBHeader header);

    /** get the ID used by Hibernate
     *
     * @return the ID used by Hibernate
     */
    public Long getId() ;

    /** set the ID used by Hibernate
     *
     * @param id
     */
    public void setId(Long id) ;

    /**
     * @param sites the sites to set in the structure
     */
    public void setSites(List<Site> sites);

    /**
     * @return the sites contained in this structure
     */
    public List<Site> getSites();

    public List<Group> getHetGroups();
    
    /**
     * Sets a flag to indicate if this structure is a biological assembly
     * @param biologicalAssembly true if biological assembly, otherwise false
     * @since 3.2
     */
    public void setBiologicalAssembly(boolean biologicalAssembly);

    /**
     * Gets flag that indicates if this structure is a biological assembly
     * @return  true if biological assembly, otherwise false
     * @since 3.2
     */
    public boolean isBiologicalAssembly();

    /**
     * Sets crystallographic information for this structure
     * @param PDBCrystallographicInfo crystallographic information
     * @since 3.2
     */
    
    public void setCrystallographicInfo(PDBCrystallographicInfo crystallographicInfo);
    /**
     * Gets crystallographic information for this structure
     * @return PDBCrystallographicInfo crystallographic information
     * @since 3.2
     */
    public PDBCrystallographicInfo getCrystallographicInfo();
    
}
