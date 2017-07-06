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

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

import org.biojava.nbio.structure.io.FileConvert;
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
public class StructureImpl implements Structure {

	private static final long serialVersionUID = -8344837138032851348L;

	private static final Logger logger = LoggerFactory.getLogger(StructureImpl.class);

	private String pdb_id ;

	private List<Model> models;

	private List<EntityInfo> entityInfos;
	private List<DBRef> dbrefs;
	private List<Bond> ssbonds;
	private List<Site> sites;
	private String name ;
	private StructureIdentifier structureIdentifier;

	private PDBHeader pdbHeader;

	private Long id;
	private boolean biologicalAssembly;

	/**
	 *  Constructs a StructureImpl object.
	 */
	public StructureImpl() {
		super();

		models         = new ArrayList<>();
		name           = "";
		entityInfos      = new ArrayList<>();
		dbrefs         = new ArrayList<>();
		pdbHeader      = new PDBHeader();
		ssbonds        = new ArrayList<>();
		sites          = new ArrayList<>();
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
	 * @param id the hibernate ID
	 */
	@Override
	public void setId(Long id) {
		this.id = id;
	}



	/** Construct a Structure object that only contains a single group
	 *
	 * @param g group object
	 */
	public StructureImpl(Group g){
		this();

		Chain c = new ChainImpl();
		c.addGroup(g);

		addChain(c);
	}

	/** construct a Structure object that contains a particular chain
	 *
	 * @param c chain
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
		// Note: structures are also cloned in SubstructureIdentifier.reduce().
		// Changes might need to be made there as well

		Structure n = new StructureImpl();
		// go through whole substructure and clone ...

		// copy structure data

		n.setPDBCode(getPDBCode());
		n.setName(getName());
		//TODO the header data is not being deep-copied, that's a minor issue since it is just some static metadata, but we should recheck this if needed - JD 2014-12-11
		n.setPDBHeader(pdbHeader);
		n.setDBRefs(this.getDBRefs());
		n.setSites(getSites());


		// go through each chain and clone chain
		for (int i=0;i<nrModels();i++){
			List<Chain> cloned_model = new ArrayList<Chain>();

			for (int j=0;j<size(i);j++){

				Chain cloned_chain  = (Chain) getChainByIndex(i,j).clone();

				// setting the parent: can only be done from the parent
				cloned_chain.setStructure(n);

				cloned_model.add(cloned_chain);

			}
			n.addModel(cloned_model);

		}

		// deep-copying of entityInfofos is tricky: there's cross references also in the Chains
		// beware: if we copy the entityInfos we would also need to reset the references to entityInfos in the individual chains
		List<EntityInfo> newEntityInfoList = new ArrayList<EntityInfo>();
		for (EntityInfo entityInfo : this.entityInfos) {
			EntityInfo newEntityInfo = new EntityInfo(entityInfo); // this sets everything but the chains
			for (String asymId:entityInfo.getChainIds()) {

				for (int modelNr=0;modelNr<n.nrModels();modelNr++) {
					Chain newChain = n.getChain(asymId,modelNr);
					if (newChain==null) {
						// this actually happens for structure 1msh, which has no chain B for model 29 (clearly a deposition error)
						logger.warn("Could not find chain asymId "+asymId+" of model "+modelNr+" while cloning entityInfo "+entityInfo.getMolId()+". Something is wrong!");
						continue;
					}
					newChain.setEntityInfo(newEntityInfo);
					newEntityInfo.addChain(newChain);
				}
			}
			newEntityInfoList.add(newEntityInfo);
		}
		n.setEntityInfos(newEntityInfoList);
		// TODO ssbonds are complicated to clone: there are deep references inside Atom objects, how would we do it? - JD 2016-03-03

		return n ;
	}


	/** {@inheritDoc} */
	@Override
	public Group findGroup(String chainName, String pdbResnum, int modelnr)
			throws StructureException {


		// if structure is xray there will be only one "model".
		if ( modelnr > models.size())
			throw new StructureException(" no model nr " + modelnr +
					" in this structure. (contains "+models.size()+")");


		// first we need to gather all groups with the author id chainName: polymers, non-polymers and waters
		Chain polyChain = getPolyChainByPDB(chainName, modelnr);
		if(polyChain != null) {
			List<Group> groups = new ArrayList<>();

			groups.addAll(polyChain.getAtomGroups());


			// there can be more than one non-poly chain for a given author id
			for (Chain chain: getNonPolyChainsByPDB(chainName, modelnr)) {
				groups.addAll(chain.getAtomGroups());
			}

			Chain water = getWaterChainByPDB(chainName, modelnr);

			if (water!=null)
				groups.addAll(water.getAtomGroups());



			// now iterate over all groups 
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
		}
		throw new StructureException("could not find group " + pdbResnum +
				" in chain " + chainName);
	}


	/** {@inheritDoc} */
	@Override
	public Group findGroup(String chainName, String pdbResnum) throws StructureException
	{
		return findGroup(chainName, pdbResnum, 0);

	}




	/** {@inheritDoc} */
	@Override
	public Chain findChain(String chainName, int modelnr) throws StructureException {

		return getChainByPDB(chainName, modelnr);
		
	}


	/** {@inheritDoc} */
	@Override
	public Chain findChain(String chainId) throws StructureException {

		return findChain(chainId,0);
	}


	/** {@inheritDoc} */
	@Override
	public void setPDBCode (String pdb_id_) {
		pdb_id = pdb_id_ ;
	}

	/** {@inheritDoc} */
	@Override
	public String  getPDBCode () {
		return pdb_id ;
	}



	/** {@inheritDoc} */
	@Override
	public void   setName(String nam) { name = nam; }

	/** {@inheritDoc} */
	@Override
	public String getName()           { return name;  }



	/**
	 * @return The StructureIdentifier used to create this structure
	 */
	@Override
	public StructureIdentifier getStructureIdentifier() {
		return structureIdentifier;
	}

	/**
	 * @param structureIdentifier the structureIdentifier corresponding to this structure
	 */
	@Override
	public void setStructureIdentifier(StructureIdentifier structureIdentifier) {
		this.structureIdentifier = structureIdentifier;
	}

	/** {@inheritDoc} */
	@Override
	public void addChain(Chain chain) {
		int modelnr = 0 ;
		addChain(chain,modelnr);
	}

	/** {@inheritDoc} */
	@Override
	public void addChain(Chain chain, int modelnr) {
		// if model has not been initialized, init it!
		chain.setStructure(this);
		if (models.isEmpty()) {
			Model model = new Model();
			List<Chain> modelChains = new ArrayList<>() ;
			modelChains.add(chain);
			model.setChains(modelChains);
			models.add(model);

		} else {
			Model model = models.get(modelnr);
			model.addChain(chain);
		}



	}



	/** {@inheritDoc} */
	@Override
	public Chain getChainByIndex(int number) {

		int modelnr = 0 ;

		return getChainByIndex(modelnr,number);
	}


	/** {@inheritDoc} */
	@Override
	public Chain getChainByIndex(int modelnr,int number) {

		Model model = models.get(modelnr);

		return model.getChains().get(number);
	}



	/** {@inheritDoc} */
	@Override
	public void addModel(List<Chain> modelChains){
		for (Chain c: modelChains){
			c.setStructure(this);
		}
		Model model = new Model();
		model.setChains(modelChains);
		models.add(model);
	}


	/** {@inheritDoc} */
	@Override
	public void setChains(List<Chain> chains){

		setModel(0,chains);
	}



	/** {@inheritDoc} */
	@Override
	public void setModel(int position, List<Chain> modelChains){
		if (modelChains == null)
			throw new IllegalArgumentException("trying to set model to null!");

		for (Chain c: modelChains)
			c.setStructure(this);

		//System.out.println("model size:" + models.size());


		Model model = new Model();
		model.setChains(modelChains);

		if (models.isEmpty()){
			models.add(model);
		} else {
			models.set(position, model);
		}
	}

	/** String representation.
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

				Chain cha = getChainByIndex(i,j);
				List<Group> agr = cha.getAtomGroups(GroupType.AMINOACID);
				List<Group> hgr = cha.getAtomGroups(GroupType.HETATM);
				List<Group> ngr = cha.getAtomGroups(GroupType.NUCLEOTIDE);




				str.append("chain ")
						.append(j).append(": asymId:")
						.append(cha.getId())
						.append(" authId:")
						.append(cha.getName()).append(" ");


				if ( cha.getEntityInfo() != null){
					EntityInfo comp = cha.getEntityInfo();
					String molName = comp.getDescription();
					if ( molName != null){
						str.append(molName);
					}
					String type =  comp.getType().toString();
					str.append(" (")
							.append(type)
							.append(")");
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
		for (EntityInfo mol : entityInfos) {
			str.append(mol).append(newline);
		}


		return str.toString() ;
	}

	@Override
	public int size() {
		int modelnr = 0 ;

		if (!models.isEmpty()) {
			return models.get(modelnr).getPolyChains().size();
		}
		else {
			return 0 ;
		}

	}

	/** return number of chains  of model.
	 *
	 */
	@Override
	public int size(int modelnr) { return models.get(modelnr).size(); }

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
				// space group defined but no crystal cell: incomplete info, return false
				return  pdbHeader.getCrystallographicInfo().getCrystalCell() != null &&
						pdbHeader.getCrystallographicInfo().getCrystalCell().isCellReasonable();
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

	/** {@inheritDoc} */
	@Override
	public List<Chain> getChains(int modelIdx){
		return getModel(modelIdx);
	}

	/** {@inheritDoc} */
	@Override
	public List<Chain> getChains(){
		if (models.size()==0) {
			return new ArrayList<>(0);
		}
		return getChains(0);

	}

	@Override
	public List<Chain> getPolyChains() { 
		if (models.size()==0) {
			return new ArrayList<>(0);
		}
		return getPolyChains(0);
	}

	@Override
	public List<Chain> getPolyChains(int modelIdx) {
		return models.get(modelIdx).getPolyChains();
	}

	@Override
	public List<Chain> getNonPolyChains() { 
		if (models.size()==0) {
			return new ArrayList<>(0);
		}
		return  getNonPolyChains(0);
	}

	@Override
	public List<Chain> getNonPolyChains(int modelIdx) {
		return models.get(modelIdx).getNonPolyChains();
	}
	
	@Override
	public List<Chain> getWaterChains() {
		if (models.size()==0) {
			return new ArrayList<>(0);
		}
		return getWaterChains(0);
	}

	@Override
	public List<Chain> getWaterChains(int modelIdx) {
		return models.get(modelIdx).getWaterChains();
	}



	/** {@inheritDoc} */
	@Override
	public void setChains(int modelnr, List<Chain> chains){
		for (Chain c: chains){
			c.setStructure(this);
		}
		if (models.size()>modelnr) {
			models.remove(modelnr);
		}

		Model model = new Model();
		model.setChains(chains);
		models.add(modelnr, model);

	}

	/** Retrieve all Chains belonging to a model .
	 *
	 * @param modelnr  an int
	 * @return a List object
	 */
	@Override
	public List<Chain> getModel(int modelnr) {

		return models.get(modelnr).getChains();
	}

	/** {@inheritDoc} */
	@Override
	public Chain getChainByPDB(String authId, int modelnr)
			throws StructureException{

		Chain c = getPolyChainByPDB(authId, modelnr);
		
		if (c==null ) {
			throw new StructureException("Could not find chain with authId \"" + authId + "\"" + " for PDB id " + pdb_id + ", model "+modelnr);			
		}
		
		return c;
	}

	/** {@inheritDoc} */
	@Override
	public Chain getChain(String asymId, int modelnr) {

		List<Chain> chains = getChains(modelnr);
		for (Chain c : chains) {
			if (c.getId().equals(asymId)) {
				return c;
			}
		}
		return null;

	}

	/** {@inheritDoc} */
	@Override
	public Chain getChain(String asymId) {

		return getChain(asymId,0);

	}

	/** {@inheritDoc} */
	@Override
	public Chain getChainByPDB(String chainId)
			throws StructureException{
		if(nrModels() < 1 ) {
			throw new StructureException("No chains are present.");
		}
		return getChainByPDB(chainId,0);
	}

	@Override
	public Chain getPolyChain(String asymId) {
		return getPolyChain(asymId, 0);

	}
	
	@Override
	public Chain getPolyChain(String asymId, int modelIdx) {
		Model model = models.get(modelIdx);
		if (model==null) {
			return null;
		}
		List<Chain> polyChains = model.getPolyChains();
		for (Chain c : polyChains){
			if (c.getId().equals(asymId))
				return c;
		}
		return null;
	}


	@Override
	public Chain getNonPolyChain(String asymId) {
		return getNonPolyChain(asymId, 0);
	}
	
	@Override
	public Chain getNonPolyChain(String asymId, int modelIdx) {
		Model model = models.get(modelIdx);
		if (model==null) {
			return null;
		}
		
		List<Chain> nonpolyChains = model.getNonPolyChains();
		for (Chain c : nonpolyChains){
			if (c.getId().equals(asymId))
				return c;
		}

		return null;
	}

	@Override
	public Chain getPolyChainByPDB(String authId) {
		return getPolyChainByPDB(authId, 0);
	}

	@Override
	public Chain getPolyChainByPDB(String authId, int modelIdx) {
		Model model = models.get(modelIdx);
		if (model==null) {
			return null;
		}

		List<Chain> polyChains = model.getPolyChains();
		for (Chain c : polyChains){
			if (c.getName().equals(authId))
				return c;
		}

		return null;
	}

	@Override
	public List<Chain> getNonPolyChainsByPDB(String authId) {
		return getNonPolyChainsByPDB(authId, 0);
	}
	
	@Override
	public List<Chain> getNonPolyChainsByPDB(String authId, int modelIdx) {
		List<Chain> chains = new ArrayList<>();
		Model model = models.get(modelIdx);
		if (model==null) {
			return chains;
		}


		List<Chain> nonpolyChains = model.getNonPolyChains();
		for (Chain c : nonpolyChains){
			if (c.getName().equals(authId))
				chains.add(c);
		}

		return chains;
	}

	@Override
	public Chain getWaterChain(String asymId) {
		return getWaterChain(asymId, 0);
	}


	@Override
	public Chain getWaterChain(String asymId, int modelIdx) {
		Model model = models.get(modelIdx);
		if (model==null) {
			return null;
		}
		List<Chain> waterChains = model.getWaterChains();
		for (Chain c : waterChains){
			if (c.getId().equals(asymId))
				return c;
		}
		return null;
	}


	@Override
	public Chain getWaterChainByPDB(String authId) {
		return getWaterChainByPDB(authId, 0);
	}


	@Override
	public Chain getWaterChainByPDB(String authId, int modelIdx) {
		Model model = models.get(modelIdx);
		if (model==null) {
			return null;
		}
		List<Chain> waterChains = model.getWaterChains();
		for (Chain c : waterChains){
			if (c.getName().equals(authId))
				return c;
		}

		return null;
	}



	/** {@inheritDoc} */
	@Override
	public String toPDB() {
		FileConvert f = new FileConvert(this) ;
		return f.toPDB();
	}

	/** {@inheritDoc} */
	@Override
	public String toMMCIF() {
		FileConvert f = new FileConvert(this);
		return f.toMMCIF();
	}

	/** {@inheritDoc} */
	@Override
	public boolean hasChain(String authId) {
		int modelnr = 0;

		List<Chain> chains = getChains(modelnr);
		for (Chain c : chains) {
			// we check here with equals because we might want to distinguish between upper and lower case chains!
			if (c.getId().equals(authId)) {
				return true;
			}
		}
		return false;
	}

	/** {@inheritDoc} */
	@Override
	public boolean hasNonPolyChain(String asymId){
		int modelnr = 0;

		List<Chain> chains = models.get(modelnr).getNonPolyChains();
		for (Chain c : chains) {
			// we check here with equals because we might want to distinguish between upper and lower case chains!
			if (c.getId().equals(asymId)) {
				return true;
			}
		}
		return false;
	}

	/** {@inheritDoc} */
	@Override
	public boolean hasPdbChain(String authId) {
		int modelnr = 0;

		List<Chain> chains = getChains(modelnr);
		for (Chain c : chains) {
			// we check here with equals because we might want to distinguish between upper and lower case chains!
			if (c.getName().equals(authId)) {
				return true;
			}
		}
		return false;
	}

	/** {@inheritDoc} */
	@Override
	public void setEntityInfos(List<EntityInfo> molList){
		this.entityInfos = molList;
	}

	/** {@inheritDoc} */
	@Override
	public void addEntityInfo(EntityInfo entityInfo) {
		this.entityInfos.add(entityInfo);
	}

	/** {@inheritDoc} */
	@Override
	public List<EntityInfo> getEntityInfos() {
		return entityInfos;
	}

	/** {@inheritDoc} */
	@Override
	public EntityInfo getCompoundById(int molId) {
		return getEntityById(molId);
	}

	/** {@inheritDoc} */
	@Override
	public EntityInfo getEntityById(int entityId) {
		for (EntityInfo mol : this.entityInfos){
			if (mol.getMolId()==entityId){
				return mol;
			}
		}
		return null;
	}


	/** {@inheritDoc} */
	@Override
	public List<DBRef> getDBRefs() {
		return dbrefs;
	}


	/** {@inheritDoc} */
	@Override
	public void setDBRefs(List<DBRef> dbrefs) {
		if ( dbrefs == null)
			throw new IllegalArgumentException("trying to set dbrefs to null!");

		for( DBRef ref : dbrefs){
			ref.setParent(this);
		}
		this.dbrefs = dbrefs;
	}


	/** {@inheritDoc} */
	@Override
	public PDBHeader getPDBHeader() {
		return pdbHeader;
	}

	/** {@inheritDoc} */
	@Override
	public void setPDBHeader(PDBHeader pdbHeader){
		this.pdbHeader = pdbHeader;
	}

	/** {@inheritDoc} */
	@Override
	public List<Bond> getSSBonds(){
		return ssbonds;

	}

	/** {@inheritDoc} */
	@Override
	public void setSSBonds(List<Bond> ssbonds){
		this.ssbonds = ssbonds;
	}

	/**
	 * Adds a single disulfide Bond to this structure
	 *
	 * @param ssbond the SSBond.
	 */
	@Override
	public void addSSBond(Bond ssbond){
		ssbonds.add(ssbond);
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

	/** {@inheritDoc} */
	@Override
	public String getIdentifier() {
		//1. StructureIdentifier
		if(getStructureIdentifier() != null) {
			return getStructureIdentifier().getIdentifier();
		}
		//2. Name
		if(getName() != null) {
			return getName();
		}
		//3. PDBCode + ranges
		return toCanonical().getIdentifier();
	}

	/** {@inheritDoc} */
	@Deprecated
	@Override
	public String getPdbId() {
		return pdb_id;
	}

	/** {@inheritDoc} */
	@Override
	public void resetModels() {
		models = new ArrayList<>();
	}
	/** {@inheritDoc} */
	@Deprecated
	@Override
	public List<ResidueRange> getResidueRanges() {
		return toCanonical().getResidueRanges();
	}
	/** {@inheritDoc} */
	@Deprecated
	@Override
	public List<String> getRanges() {
		return ResidueRange.toStrings(getResidueRanges());
	}

	/**
	 * Creates a SubstructureIdentifier based on the residues in this Structure.
	 *
	 * Only the first and last residues of each chain are considered, so chains
	 * with gaps
	 * @return A {@link SubstructureIdentifier} with residue ranges constructed from each chain
	 */
	private SubstructureIdentifier toCanonical() {
		StructureIdentifier real = getStructureIdentifier();
		if(real != null) {
			try {
				return real.toCanonical();
			} catch (StructureException e) {
				// generate fake one if needed
			}
		}

		// No identifier set, so generate based on residues present in the structure
		List<ResidueRange> range = new ArrayList<>();
		for (Chain chain : getChains()) {
			List<Group> groups = chain.getAtomGroups();
			ListIterator<Group> groupsIt = groups.listIterator();
			if(!groupsIt.hasNext()) {
				continue; // no groups in chain
			}
			Group g = groupsIt.next();
			ResidueNumber first = g.getResidueNumber();

			//TODO Detect missing intermediate residues -sbliven, 2015-01-28
			//Already better than previous whole-chain representation

			// get last residue
			while(groupsIt.hasNext()) {
				g = groupsIt.next();
			}
			ResidueNumber last = g.getResidueNumber();

			range.add(new ResidueRange(chain.getName(),first,last));
		}
		return new SubstructureIdentifier(getPDBCode(),range);
	}





}
