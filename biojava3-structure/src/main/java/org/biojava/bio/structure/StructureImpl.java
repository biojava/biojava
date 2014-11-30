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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.structure.io.EntityFinder;
import org.biojava.bio.structure.io.FileConvert;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implementation of a PDB Structure. This class
 * provides the data contained in a PDB file.
 * to get structure objects from different sources
 * see io package.
 *
 * @author Andreas Prlic
 * @author Jules Jacobsen
 * @since 1.4
 * @version %I% %G%
 */
public class StructureImpl implements Structure, Serializable {

	private static final long serialVersionUID = -8344837138032851347L;
	
	private static final Logger logger = LoggerFactory.getLogger(StructureImpl.class);

	private String pdb_id ;
	/* models is an ArrayList of ArrayLists */
	private List<List<Chain>> models;
	//List<Chain> seqResList;
	private List<Map <String,Integer>> connections ;
	private List<Compound> compounds;
	private List<DBRef> dbrefs;
	private List<SSBond> ssbonds;
	private List<Site> sites;
	private List<Group> hetAtoms;
	private String name ;
	
	private PDBHeader pdbHeader;

	private Long id;
	private boolean biologicalAssembly;

	/**
	 * A map of chain identifiers to entities
	 */
	private TreeMap<String,Entity> chainIds2entities;
	
	/**
	 *  Constructs a StructureImpl object.
	 */
	public StructureImpl() {
		super();

		models         = new ArrayList<List<Chain>>();
		name           = "";
		connections    = new ArrayList<Map<String,Integer>>();
		compounds      = new ArrayList<Compound>();
        dbrefs         = new ArrayList<DBRef>();
        pdbHeader      = new PDBHeader();
        ssbonds        = new ArrayList<SSBond>();
        sites          = new ArrayList<Site>();
        hetAtoms          = new ArrayList<Group>();
	}

	/** get the ID used by Hibernate
     *
     * @return the ID used by Hibernate
     */
    public Long getId() {
        return id;
    }

    /** set the ID used by Hibernate
     *
     * @param id
     */
    public void setId(Long id) {
        this.id = id;
    }


	/** construct a Structure object that only contains a single group
	 *
	 * @param g
	 */
	public StructureImpl(Group g){
		this();

		Chain c = new ChainImpl();
		c.addGroup(g);

		addChain(c);
	}

	/** construct a Structure object that contains a particular chain
	 *
	 * @param c
	 */
	public StructureImpl(Chain c){
		this();
		addChain(c);
	}

	/** returns an identical copy of this structure .
	 * @return an identical Structure object
	 */
	public Structure clone() {

		Structure n = new StructureImpl();
		// go through whole substructure and clone ...

		// copy structure data

		n.setPDBCode(getPDBCode());
		n.setName(getName());
		//TODO: do deep copying of data!
		n.setPDBHeader(pdbHeader);
		n.setDBRefs(this.getDBRefs());
		n.setConnections(getConnections());
		n.setSites(getSites());
		n.setCrystallographicInfo(getCrystallographicInfo());

		// go through each chain and clone chain
		for (int i=0;i<nrModels();i++){
			List<Chain> cloned_model = new ArrayList<Chain>();

			for (int j=0;j<size(i);j++){

				Chain current_chain = (Chain) getChain(i,j);
				Chain cloned_chain  = (Chain) current_chain.clone();

				cloned_model.add(cloned_chain);
			}
			n.addModel(cloned_model);

		}

		for (SSBond ssbond: ssbonds){
			n.addSSBond(ssbond.clone());
		}

		//n.setSeqRes(this.getSeqRes());
		return n ;
	}


	public Group findGroup(String chainId, String pdbResnum, int modelnr)
	throws StructureException {


		// if structure is xray there will be only one "model".
		if ( modelnr > models.size())
			throw new StructureException(" no model nr " + modelnr +
					" in this structure. (contains "+models.size()+")");


		Chain c = findChain(chainId,modelnr);

		List<Group> groups = c.getAtomGroups();

		// now iterate over all groups in this chain.
		// in order to find the amino acid that has this pdbRenum.

		Iterator<Group> giter = groups.iterator();
		while (giter.hasNext()){
			Group g =  giter.next();
			String rnum = g.getResidueNumber().toString();
			//System.out.println(g + " >" + rnum + "< >" + pdbResnum + "<");
			// we only mutate amino acids
			// and ignore hetatoms and nucleotides in this case
			if (rnum.equals(pdbResnum))
				return g;
		}

		throw new StructureException("could not find group " + pdbResnum +
				" in chain " + chainId);
	}


	public Group findGroup(String chainName, String pdbResnum) throws StructureException
	{
		return findGroup(chainName, pdbResnum, 0);

	}




	public Chain findChain(String chainId, int modelnr) throws StructureException {

		List<Chain> chains = getChains(modelnr);

		// iterate over all chains.
		Iterator<Chain> iter = chains.iterator();
		while (iter.hasNext()){
			Chain c = iter.next();

			if (c.getChainID().equals(chainId)) {
				return c;
			}
		}
		throw new StructureException("could not find chain \"" + chainId + "\" for PDB id " + pdb_id);
	}


	public Chain findChain(String chainId) throws StructureException {

		return findChain(chainId,0);
	}


	/**
	 *
	 * set PDB code of structure .
	 * @see #getPDBCode
	 *
	 */
	public void setPDBCode (String pdb_id_) {
		pdb_id = pdb_id_ ;
	}
	/**
	 *
	 * get PDB code of structure .
	 *
	 * @return a String representing the PDBCode value
	 * @see #setPDBCode
	 */
	public String  getPDBCode () {
		return pdb_id ;
	}



	/** set biological name of Structure.
	 *
	 * @see #getName
	 *
	 */
	public void   setName(String nam) { name = nam; }

	/** get biological name of Structure.
	 *
	 * @return a String representing the name
	 * @see #setName
	 */
	public String getName()           { return name;  }


	
	public void      setConnections(List<Map<String,Integer>> conns) { connections = conns ; }
	
	/**
	 * Returns the connections value.
	 *
	 * @return a List object representing the connections value
	 * @see Structure interface
	 * @see #setConnections
	 */
	public List<Map<String,Integer>> getConnections()                { return connections ;}

	public void addChain(Chain chain) {
		int modelnr = 0 ;
		addChain(chain,modelnr);
	}
	
	public void addChain(Chain chain, int modelnr) {
		// if model has not been initialized, init it!
		chain.setParent(this);
		if ( models.size() == 0  ) {
			List<Chain> model = new ArrayList<Chain>() ;
			model.add(chain);
			models.add(model);

		} else {
			List<Chain> model = models.get(modelnr);
			model.add(chain);
		}



	}


	
	public Chain getChain(int number) {

		int modelnr = 0 ;

		return getChain(modelnr,number);
	}

	
	public Chain getChain(int modelnr,int number) {

		List<Chain> model  =  models.get(modelnr);

		Chain chain =   model.get (number );

		return chain ;
	}


	
	public void addModel(List<Chain> model){
		for (Chain c: model){
    		c.setParent(this);
    	}
		models.add(model);
	}


	public void setChains(List<Chain> chains){

		setModel(0,chains);
	}



    public void setModel(int position, List<Chain> model){
    	if (model == null)
    		throw new IllegalArgumentException("trying to set model to null!");

    	for (Chain c: model)
    		c.setParent(this);

    	//System.out.println("model size:" + models.size());

    	if (models.size() ==0){
    		models.add(model);
    	} else {
    		models.set(position, model);
    	}
    }

	/** string representation.
	 *
	 */
	public String toString(){
		String newline = System.getProperty("line.separator");
		StringBuffer str = new StringBuffer();
		str.append("structure ");
		str.append(name);
		str.append(" ");
		str.append(pdb_id);
		str.append(" ");

		if ( nrModels()>1 ){
			str.append( " models: ");
			str.append(nrModels());
			str.append(newline) ;
		}

        str.append(pdbHeader.toString());
        str.append(newline) ;

        for (int i=0;i<nrModels();i++){
			if ( nrModels()>1 ) {
				str.append(" model[");
				str.append(i);
				str.append("]:");
				str.append(newline);
			}
			str.append(" chains:");
			str.append(newline);

			for (int j=0;j<size(i);j++){

				Chain cha = (Chain)getChain(i,j);
				List<Group> agr = cha.getAtomGroups("amino");
				List<Group> hgr = cha.getAtomGroups("hetatm");
				List<Group> ngr = cha.getAtomGroups("nucleotide");

				str.append("chain " + j + ": >"+cha.getChainID()+"< ");
				if ( cha.getHeader() != null){
					Compound comp = cha.getHeader();
					String molName = comp.getMolName();
					if ( molName != null){
						str.append(molName);
					}
				}


				str.append(newline);
                str.append(" length SEQRES: ").append(cha.getSeqResLength());
                str.append(" length ATOM: ").append(cha.getAtomLength());
                str.append(" aminos: ").append(agr.size());
                str.append(" hetatms: ").append(hgr.size());
				str.append(" nucleotides: "+ngr.size() + newline);
			}

		}
        str.append("DBRefs: "+ dbrefs.size()+ newline);
        for (DBRef dbref: dbrefs){
            str.append(dbref.toPDB()).append(newline);
        }
        str.append("Molecules: ").append(newline);
		Iterator<Compound> iter = compounds.iterator();
		while (iter.hasNext()){
			Compound mol = iter.next();
            str.append(mol).append(newline);
		}


		return str.toString() ;
	}

	/** return number of chains , if NMR return number of chains of first model .
	 *
	 */
	public int size() {
		int modelnr = 0 ;

		if ( models.size() > 0) {
			return models.get(modelnr).size();
		}
		else {
			return 0 ;
		}

	}

	/** return number of chains  of model.
	 *
	 */
	public int size(int modelnr) { return getChains(modelnr).size();   }

	// some NMR stuff :

	/** return number of models. */
	public int nrModels() {
		return models.size() ;
	}

	/**
	 * Whether this Structure is a crystallographic structure or not.
	 * It will first check the experimental technique and if not present it will try
	 * to guess from the presence of a space group and sensible cell parameters  
	 * 
	 * @return true if crystallographic, false otherwise
	 */
	public boolean isCrystallographic() {
		if (pdbHeader.getExperimentalTechniques()!=null) {
			return ExperimentalTechnique.isCrystallographic(pdbHeader.getExperimentalTechniques());
		} else {
			// no experimental technique known, we try to guess...
			if (pdbHeader.getCrystallographicInfo().getSpaceGroup()!=null) {
				return pdbHeader.getCrystallographicInfo().getCrystalCell().isCellReasonable();
			}
		}
		return false;
	}
	
	/**
	 * Whether this Structure is a NMR structure or not.
	 * It will first check the experimental technique and if not present it will try
	 * to guess from the presence of more than 1 model and from b-factors being 0 in first chain of first model
	 * @return true if NMR, false otherwise
	 */
	public boolean isNmr() {
		
		// old implementation was:
		//return nmrflag;
		
		if (pdbHeader.getExperimentalTechniques()!=null) {
			return ExperimentalTechnique.isNmr(pdbHeader.getExperimentalTechniques());
		} else {
			// no experimental technique known, we try to guess...
			if (nrModels()>1) {
				if (pdbHeader.getCrystallographicInfo().getSpaceGroup()!=null) {
					// multi-model and cell unreasonable: must be NMR
					if (!pdbHeader.getCrystallographicInfo().getCrystalCell().isCellReasonable())
						return true;
				} else { 
					// multi-model and missing space group: must be NMR
					return true; 
				}
			}
		}
		return false;
	}
	
	@Deprecated
	public void setNmr(boolean nmr) {	
		// old implementation was:
		// this.nmrflag = nmr;
	}


	/** retrieve all chains of a model.
	 *
	 * @param modelnr  an int
	 * @return a List object
	 */
	public List<Chain> getChains(int modelnr){
		return getModel(modelnr);
	}

	public List<Chain> getChains(){
		return getModel(0);
	}

    public void setChains(int modelnr, List<Chain> chains){
    	for (Chain c: chains){
    		c.setParent(this);
    	}
        models.remove(modelnr);
        models.add(modelnr, chains);

    }

	/** retrieve all Chains belonging to a model .
	 *
	 * @param modelnr  an int
	 * @return a List object
	 */
	public List<Chain> getModel(int modelnr) {

		List<Chain> model = models.get(modelnr);
		return model;
	}




	public Chain getChainByPDB(String chainId, int modelnr)
	throws StructureException{

		List<Chain> chains = getChains(modelnr);
		Iterator<Chain> iter = chains.iterator();
		while ( iter.hasNext()){
			Chain c = iter.next();
			if ( c.getChainID().equals(chainId))
				return c;
		}
		throw new StructureException("did not find chain with chainId \"" + chainId + "\"" + " for PDB id " + pdb_id);

	}


	public Chain getChainByPDB(String chainId)
	throws StructureException{
		return getChainByPDB(chainId,0);
	}


	/** create a String that contains the contents of a PDB file.
	 *
	 * @return a String that represents the structure as a PDB file.
	 */
	public String toPDB() {
		FileConvert f = new FileConvert(this) ;

		String str = f.toPDB();


		return str ;

	}


	public boolean hasChain(String chainId) {
		int modelnr = 0;

		List<Chain> chains = getChains(modelnr);
		Iterator<Chain> iter = chains.iterator();
		while ( iter.hasNext()){
			Chain c = iter.next();
			// we check here with equals because we might want to distinguish between upper and lower case chains!
			if ( c.getChainID().equals(chainId))
				return true;
		}
		return false;
	}

	public void setCompounds(List<Compound>molList){
		this.compounds = molList;
	}

	public List<Compound> getCompounds() {
		return compounds;
	}

	public Compound getCompoundById(String molId) {
		for (Compound mol : this.compounds){
			if (mol.getMolId().equals(molId)){
				return mol;
			}
		}
		return null;
	}


    public List<DBRef> getDBRefs() {
       return dbrefs;
    }


    public void setDBRefs(List<DBRef> dbrefs) {
    	if ( dbrefs == null)
    		throw new IllegalArgumentException("trying to set dbrefs to null!");

    	for( DBRef ref : dbrefs){
    		ref.setParent(this);
    	}
        this.dbrefs = dbrefs;
    }


	public PDBHeader getPDBHeader() {
		return pdbHeader;
	}

	public void setPDBHeader(PDBHeader pdbHeader){
		this.pdbHeader = pdbHeader;
	}

    /** get the list of SSBonds as they have been defined in the PDB files
     *
     * @return a list of SSBonds
     */
    public List<SSBond> getSSBonds(){
    	return ssbonds;

    }
    /** set the list of SSBonds for this structure
     *
     * @param ssbonds
     */
    public void setSSBonds(List<SSBond> ssbonds){
    	this.ssbonds = ssbonds;
    }

    /** add a single SSBond to this structure
     *
     * @param ssbond the SSBond.
     */
    public void addSSBond(SSBond ssbond){
    	ssbonds.add(ssbond);
    	ssbond.setSerNum(ssbonds.size());
    }

    /**
     * Return whether or not the entry has an associated journal article
     * or publication. The JRNL section is not mandatory and thus may not be
     * present.
     * @return flag if a JournalArticle could be found.
     */
    public boolean hasJournalArticle() {
    	return this.pdbHeader.hasJournalArticle();
    }

    /**
     * get the associated publication as defined by the JRNL records in a PDB
     * file.
     * @return a JournalArticle
     */
    public JournalArticle getJournalArticle() {
        return this.pdbHeader.getJournalArticle();
    }

    /**
     * set the associated publication as defined by the JRNL records in a PDB
     * file.
     * @param journalArticle the article
     */
    public void setJournalArticle(JournalArticle journalArticle) {
        this.pdbHeader.setJournalArticle(journalArticle);
    }

    /**
     * @return the sites contained in this structure
     */
  
    public List<Site> getSites() {
            return sites;
    }

    /**
     * @param sites the sites to set in the structure
     */
 
    public void setSites(List<Site> sites) {
            this.sites = sites;
    }

    /** Caution: we should probably remove this to avoid confusion. Currently this is always an empty list!
     *
     * @return a list of Groups listed in the HET records - this will not
     * include any waters.
     */
    
    public List<Group> getHetGroups() {
        return hetAtoms;
    }
    
    /**
     * Sets a flag to indicate if this structure is a biological assembly
     * @param biologicalAssembly true if biological assembly, otherwise false
     * @since 3.2
     */
    public void setBiologicalAssembly(boolean biologicalAssembly) {
    	this.biologicalAssembly = biologicalAssembly;
    }

    /**
     * Gets flag that indicates if this structure is a biological assembly
     * @return the sites contained in this structure
     * @since 3.2
     */
    public boolean isBiologicalAssembly() {
    	return biologicalAssembly;
    }
    
    /**
     * Sets crystallographic information for this structure
     * @param PDBCrystallographicInfo crystallographic information
     * @since 3.2
     */
    
    public void setCrystallographicInfo(PDBCrystallographicInfo crystallographicInfo) {
    	this.pdbHeader.setCrystallographicInfo(crystallographicInfo);
    }
    
    /**
     * Gets crystallographic information for this structure
     * @return PDBCrystallographicInfo crystallographic information
     * @since 3.2
     */
    public PDBCrystallographicInfo getCrystallographicInfo() {
    	return pdbHeader.getCrystallographicInfo();
    }

	@Override
	public String getIdentifier() {
		return pdb_id;
	}

	@Override
	public String getPdbId() {
		return pdb_id;
	}

	@Override
	public List<ResidueRange> getResidueRanges() {
		List<ResidueRange> range = new ArrayList<ResidueRange>();
		for (Chain chain : getChains()) {
			range.add(ResidueRange.parse(pdb_id + "." + chain.getChainID()));
		}
		return range;
	}

	@Override
	public List<String> getRanges() {
		return ResidueRange.toStrings(getResidueRanges());
	}

	@Override
	public List<Entity> getEntities() {

		if (chainIds2entities !=null) {
			return findUniqueEntities();
		}
		
		// if null, it hasn't been initialised yet, let's init it:		
		chainIds2entities = new TreeMap<String,Entity>();

		// finding out whether we have SEQRES: if at least 1 chain has seqres groups, we will consider true
		boolean hasSeqRes = false;
		for (Chain chain: getChains()) {
			if (chain.getSeqResLength()>0) { 
				hasSeqRes = true;
				break;
			}
		}
		
		if (hasSeqRes) { // we have SEQRES

			logger.debug("Getting entities from SEQRES sequences");
			// map of sequences to list of chain identifiers
			Map<String, List<String>> uniqSequences = new HashMap<String, List<String>>();
			// finding the entities (groups of identical chains)
			for (Chain chain:getChains()) {

				String seq = chain.getSeqResSequence();
					
				if (uniqSequences.containsKey(seq)) {
					uniqSequences.get(seq).add(chain.getChainID());
				} else {
					List<String> list = new ArrayList<String>();
					list.add(chain.getChainID());
					uniqSequences.put(seq, list);
				}		

			}

			for (List<String> chainIds:uniqSequences.values()) {
				// sorting ids in alphabetic order
				Collections.sort(chainIds);
				List<Chain> chains = new ArrayList<Chain>();
				for (String chainId:chainIds) {
					// chains will be sorted in ids' alphabetic order
					try {
						chains.add(this.getChainByPDB(chainId));
					} catch (StructureException e) {
						// this basically can't happen, if it does it is some kind of bug
						logger.error("Unexpected exception!",e);
					}
				}
				// the representative will be the one with first chain id in alphabetic order 
				Entity entity = new Entity(chains.get(0), chains);
				for (Chain member:entity.getMembers()) {
					chainIds2entities.put(member.getChainID(), entity);
				}
			}


		} else {
			logger.debug("Getting entities from aligning ATOM sequences. If you have SEQRES in your file make sure you are using the right FileParsingParams");

			EntityFinder ef = new EntityFinder(this);
			
			chainIds2entities = ef.findEntities();
		}

		return findUniqueEntities();
	}
	
	@Override
	public Entity getEntity(String chainId) {
		if (chainIds2entities == null) 
			getEntities();
		
		return chainIds2entities.get(chainId);
	}
	
	/**
	 * Utility method to obtain a list of unique entities from the chainIds2entities map
	 * @return
	 */
	private List<Entity> findUniqueEntities() {
		
		List<Entity> list = new ArrayList<Entity>();
		
		for (Entity cluster:chainIds2entities.values()) {
			boolean present = false;
			for (Entity cl:list) {
				if (cl==cluster) {
					present = true;
					break;
				}
			}
			if (!present) list.add(cluster);
		}
		return list;
	} 
	
}
