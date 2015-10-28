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
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.CompoundFinder;
import org.biojava.nbio.structure.io.FileConvert;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

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
		hetAtoms       = new ArrayList<Group>();
	}

	/** get the ID used by Hibernate
	 *
	 * @return the ID used by Hibernate
	 */
	@Override
	public Long getId() {
		return id;
	}

	/** set the ID used by Hibernate
	 *
	 * @param id
	 */
	@Override
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
	@Override
	public Structure clone() {

		Structure n = new StructureImpl();
		// go through whole substructure and clone ...

		// copy structure data

		n.setPDBCode(getPDBCode());
		n.setName(getName());
		//TODO the header data is not being deep-copied, that's a minor issue since it is just some static metadata, but we should recheck this if needed - JD 2014-12-11
		n.setPDBHeader(pdbHeader);
		n.setDBRefs(this.getDBRefs());
		n.setConnections(getConnections());
		n.setSites(getSites());
				

		// go through each chain and clone chain
		for (int i=0;i<nrModels();i++){
			List<Chain> cloned_model = new ArrayList<Chain>();

			for (int j=0;j<size(i);j++){

				Chain cloned_chain  = (Chain) getChain(i,j).clone();

				// setting the parent: can only be done from the parent
				cloned_chain.setStructure(n);

				cloned_model.add(cloned_chain);
				
			}
			n.addModel(cloned_model);

		}

		// deep-copying of Compounds is tricky: there's cross references also in the Chains
		// beware: if we copy the compounds we would also need to reset the references to compounds in the individual chains
		List<Compound> newCompoundList = new ArrayList<Compound>();
		for (Compound compound:this.compounds) {
			Compound newCompound = new Compound(compound); // this sets everything but the chains
			for (String chainId:compound.getChainIds()) {
				
					for (int modelNr=0;modelNr<n.nrModels();modelNr++) {
						try {
							Chain newChain = n.getChainByPDB(chainId,modelNr);
							newChain.setCompound(newCompound);
							newCompound.addChain(newChain);
						} catch (StructureException e) {
							// this actually happens for structure 1msh, which has no chain B for model 29 (clearly a deposition error)
							logger.warn("Could not find chain id "+chainId+" of model "+modelNr+" while cloning compound "+compound.getMolId()+". Something is wrong!");
						}
					}
			}
			newCompoundList.add(newCompound);
		}
		n.setCompounds(newCompoundList); 


		for (SSBond ssbond: ssbonds){
			n.addSSBond(ssbond.clone());
		}

		return n ;
	}


	@Override
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

		for (Group g : groups) {
			String rnum = g.getResidueNumber().toString();
			//System.out.println(g + " >" + rnum + "< >" + pdbResnum + "<");
			// we only mutate amino acids
			// and ignore hetatoms and nucleotides in this case
			if (rnum.equals(pdbResnum)) {
				return g;
			}
		}

		throw new StructureException("could not find group " + pdbResnum +
				" in chain " + chainId);
	}


	@Override
	public Group findGroup(String chainName, String pdbResnum) throws StructureException
	{
		return findGroup(chainName, pdbResnum, 0);

	}




	@Override
	public Chain findChain(String chainId, int modelnr) throws StructureException {

		List<Chain> chains = getChains(modelnr);

		// iterate over all chains.
		for (Chain c : chains) {
			if (c.getChainID().equals(chainId)) {
				return c;
			}
		}
		throw new StructureException("Could not find chain \"" + chainId + "\" for PDB id " + pdb_id);
	}


	@Override
	public Chain findChain(String chainId) throws StructureException {

		return findChain(chainId,0);
	}


	@Override
	public void setPDBCode (String pdb_id_) {
		pdb_id = pdb_id_ ;
	}

	@Override
	public String  getPDBCode () {
		return pdb_id ;
	}



	@Override
	public void   setName(String nam) { name = nam; }

	@Override
	public String getName()           { return name;  }



	@Override
	public void      setConnections(List<Map<String,Integer>> conns) { connections = conns ; }

	/**
	 * Return the connections value.
	 *
	 * @return a List object representing the connections value
	 * @see Structure interface
	 * @see #setConnections
	 */
	@Override
	public List<Map<String,Integer>> getConnections()                { return connections ;}

	@Override
	public void addChain(Chain chain) {
		int modelnr = 0 ;
		addChain(chain,modelnr);
	}

	@Override
	public void addChain(Chain chain, int modelnr) {
		// if model has not been initialized, init it!
		chain.setStructure(this);
		if (models.isEmpty()) {
			List<Chain> model = new ArrayList<Chain>() ;
			model.add(chain);
			models.add(model);

		} else {
			List<Chain> model = models.get(modelnr);
			model.add(chain);
		}



	}



	@Override
	public Chain getChain(int number) {

		int modelnr = 0 ;

		return getChain(modelnr,number);
	}


	@Override
	public Chain getChain(int modelnr,int number) {

		List<Chain> model  =  models.get(modelnr);

		return model.get (number );
	}



	@Override
	public void addModel(List<Chain> model){
		for (Chain c: model){
			c.setStructure(this);
		}
		models.add(model);
	}


	@Override
	public void setChains(List<Chain> chains){

		setModel(0,chains);
	}



	@Override
	public void setModel(int position, List<Chain> model){
		if (model == null)
			throw new IllegalArgumentException("trying to set model to null!");

		for (Chain c: model)
			c.setStructure(this);

		//System.out.println("model size:" + models.size());

		if (models.isEmpty()){
			models.add(model);
		} else {
			models.set(position, model);
		}
	}

	/** string representation.
	 *
	 */
	@Override
	public String toString(){
		String newline = System.getProperty("line.separator");
		StringBuilder str = new StringBuilder();
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

		str.append(pdbHeader);
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

				Chain cha = getChain(i,j);
				List<Group> agr = cha.getAtomGroups(GroupType.AMINOACID);
				List<Group> hgr = cha.getAtomGroups(GroupType.HETATM);
				List<Group> ngr = cha.getAtomGroups(GroupType.NUCLEOTIDE);

				str.append("chain ").append(j).append(": >").append(cha.getChainID()).append("< ");
				if ( cha.getCompound() != null){
					Compound comp = cha.getCompound();
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
				str.append(" nucleotides: ").append(ngr.size()).append(newline);
			}

		}
		str.append("DBRefs: ").append(dbrefs.size()).append(newline);
		for (DBRef dbref: dbrefs){
			str.append(dbref.toPDB()).append(newline);
		}
		str.append("Molecules: ").append(newline);
		for (Compound mol : compounds) {
			str.append(mol).append(newline);
		}


		return str.toString() ;
	}

	/** return number of chains , if NMR return number of chains of first model .
	 *
	 */
	@Override
	public int size() {
		int modelnr = 0 ;

		if (!models.isEmpty()) {
			return models.get(modelnr).size();
		}
		else {
			return 0 ;
		}

	}

	/** return number of chains  of model.
	 *
	 */
	@Override
	public int size(int modelnr) { return getChains(modelnr).size();   }

	// some NMR stuff :

	/** return number of models. */
	@Override
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
	@Override
	public boolean isCrystallographic() {
		if (pdbHeader.getExperimentalTechniques()!=null) {
			return ExperimentalTechnique.isCrystallographic(pdbHeader.getExperimentalTechniques());
		} else {
			// no experimental technique known, we try to guess...
			if (pdbHeader.getCrystallographicInfo().getSpaceGroup()!=null) {
				if (pdbHeader.getCrystallographicInfo().getCrystalCell()==null) {
					return false; // space group defined but no crystal cell: incomplete info, return false
				} else {
					return pdbHeader.getCrystallographicInfo().getCrystalCell().isCellReasonable();
				}
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
	@Override
	public boolean isNmr() {

		// old implementation was:
		//return nmrflag;

		if (pdbHeader.getExperimentalTechniques()!=null) {
			return ExperimentalTechnique.isNmr(pdbHeader.getExperimentalTechniques());
		} else {
			// no experimental technique known, we try to guess...
			if (nrModels()>1) {
				if (pdbHeader.getCrystallographicInfo().getSpaceGroup()!=null) {
					// multimodel, sg defined, but missing cell: must be NMR
					if (pdbHeader.getCrystallographicInfo().getCrystalCell()==null) 
						return true;					
					// multi-model, sg defined and cell unreasonable: must be NMR
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

	@Override
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
	@Override
	public List<Chain> getChains(int modelnr){
		return getModel(modelnr);
	}

	@Override
	public List<Chain> getChains(){
		return getModel(0);
	}

	@Override
	public void setChains(int modelnr, List<Chain> chains){
		for (Chain c: chains){
			c.setStructure(this);
		}
		models.remove(modelnr);
		models.add(modelnr, chains);

	}

	/** retrieve all Chains belonging to a model .
	 *
	 * @param modelnr  an int
	 * @return a List object
	 */
	@Override
	public List<Chain> getModel(int modelnr) {

		return models.get(modelnr);
	}




	@Override
	public Chain getChainByPDB(String chainId, int modelnr)
			throws StructureException{

		List<Chain> chains = getChains(modelnr);
		for (Chain c : chains) {
			if (c.getChainID().equals(chainId)) {
				return c;
			}
		}
		throw new StructureException("did not find chain with chainId \"" + chainId + "\"" + " for PDB id " + pdb_id);

	}


	@Override
	public Chain getChainByPDB(String chainId)
			throws StructureException{
		return getChainByPDB(chainId,0);
	}


	@Override
	public String toPDB() {
		FileConvert f = new FileConvert(this) ;
		return f.toPDB();
	}

	@Override
	public String toMMCIF() {
		FileConvert f = new FileConvert(this);
		return f.toMMCIF();
	}

	@Override
	public boolean hasChain(String chainId) {
		int modelnr = 0;

		List<Chain> chains = getChains(modelnr);
		for (Chain c : chains) {
			// we check here with equals because we might want to distinguish between upper and lower case chains!
			if (c.getChainID().equals(chainId)) {
				return true;
			}
		}
		return false;
	}

	@Override
	public void setCompounds(List<Compound> molList){
		this.compounds = molList;
	}

	@Override
	public void addCompound(Compound compound) {
		this.compounds.add(compound);
	}

	@Override
	public List<Compound> getCompounds() {
		// compounds are parsed from the PDB/mmCIF file normally
		// but if the file is incomplete, it won't have the Compounds information and we try 
		// to guess it from the existing seqres/atom sequences
		if (compounds==null || compounds.isEmpty()) {
			CompoundFinder cf = new CompoundFinder(this);
			this.compounds = cf.findCompounds();

			// now we need to set references in chains:
			for (Compound compound:compounds) {
				for (Chain c:compound.getChains()) {
					c.setCompound(compound);
				}
			}
		}
		return compounds;
	}

	@Override
	public Compound getCompoundById(int molId) {
		for (Compound mol : this.compounds){
			if (mol.getMolId()==molId){
				return mol;
			}
		}
		return null;
	}


	@Override
	public List<DBRef> getDBRefs() {
		return dbrefs;
	}


	@Override
	public void setDBRefs(List<DBRef> dbrefs) {
		if ( dbrefs == null)
			throw new IllegalArgumentException("trying to set dbrefs to null!");

		for( DBRef ref : dbrefs){
			ref.setParent(this);
		}
		this.dbrefs = dbrefs;
	}


	@Override
	public PDBHeader getPDBHeader() {
		return pdbHeader;
	}

	@Override
	public void setPDBHeader(PDBHeader pdbHeader){
		this.pdbHeader = pdbHeader;
	}

	/** get the list of SSBonds as they have been defined in the PDB files
	 *
	 * @return a list of SSBonds
	 */
	@Override
	public List<SSBond> getSSBonds(){
		return ssbonds;

	}
	/** set the list of SSBonds for this structure
	 *
	 * @param ssbonds
	 */
	@Override
	public void setSSBonds(List<SSBond> ssbonds){
		this.ssbonds = ssbonds;
	}

	/** add a single SSBond to this structure
	 *
	 * @param ssbond the SSBond.
	 */
	@Override
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
	@Override
	public boolean hasJournalArticle() {
		return this.pdbHeader.hasJournalArticle();
	}

	/**
	 * get the associated publication as defined by the JRNL records in a PDB
	 * file.
	 * @return a JournalArticle
	 */
	@Override
	public JournalArticle getJournalArticle() {
		return this.pdbHeader.getJournalArticle();
	}

	/**
	 * set the associated publication as defined by the JRNL records in a PDB
	 * file.
	 * @param journalArticle the article
	 */
	@Override
	public void setJournalArticle(JournalArticle journalArticle) {
		this.pdbHeader.setJournalArticle(journalArticle);
	}

	/**
	 * @return the sites contained in this structure
	 */

	@Override
	public List<Site> getSites() {
		return sites;
	}

	/**
	 * @param sites the sites to set in the structure
	 */

	@Override
	public void setSites(List<Site> sites) {
		this.sites = sites;
	}

	/** Caution: we should probably remove this to avoid confusion. Currently this is always an empty list!
	 *
	 * @return a list of Groups listed in the HET records - this will not
	 * include any waters.
	 */

	@Override
	public List<Group> getHetGroups() {
		return hetAtoms;
	}

	/**
	 * Sets a flag to indicate if this structure is a biological assembly
	 * @param biologicalAssembly true if biological assembly, otherwise false
	 * @since 3.2
	 */
	@Override
	public void setBiologicalAssembly(boolean biologicalAssembly) {
		this.biologicalAssembly = biologicalAssembly;
	}

	/**
	 * Gets flag that indicates if this structure is a biological assembly
	 * @return the sites contained in this structure
	 * @since 3.2
	 */
	@Override
	public boolean isBiologicalAssembly() {
		return biologicalAssembly;
	}

	/**
	 * Sets crystallographic information for this structure
	 * @param crystallographicInfo crystallographic information
	 * @since 3.2
	 */

	@Override
	public void setCrystallographicInfo(PDBCrystallographicInfo crystallographicInfo) {
		this.pdbHeader.setCrystallographicInfo(crystallographicInfo);
	}

	/**
	 * Gets crystallographic information for this structure
	 * @return PDBCrystallographicInfo crystallographic information
	 * @since 3.2
	 */
	@Override
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
	public void resetModels() {
		models = new ArrayList<List<Chain>>();
	}

}
