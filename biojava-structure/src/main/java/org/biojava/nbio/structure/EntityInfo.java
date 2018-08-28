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
 * Created on 22.01.2007
 *
 */
package org.biojava.nbio.structure;


import org.biojava.nbio.structure.io.FileParsingParameters;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Serializable;
import java.util.*;

/**
 * An object to contain the info from the PDB header for a Molecule.
 * In mmCIF dictionary, it is called an Entity. In the case of polymers it
 * is defined as each group of sequence identical NCS-related chains
 *
 * Now PDB file format 3.2 aware - contains the new TAX_ID fields for the
 * organism studied and the expression system.
 *
 * @author Jules Jacobsen
 * @author Jose Duarte
 * @author Anthony Bradley
 * @since 1.5
 */
public class EntityInfo implements Serializable {

	private final static Logger logger = LoggerFactory.getLogger(EntityInfo.class);


	// TODO We could drop a lot of the stuff here that is PDB-file related (actually many PDB files don't contain many of these fields) - JD 2016-03-25
	//     The only really essential part of a EntityInfo is the member chains and the entity_id/mol_id
	// See also issue https://github.com/biojava/biojava/issues/219

	private static final long serialVersionUID = 2991897825657586356L;

	/**
	 * The list of chains that are described by this EntityInfo
	 */
	private List<Chain> chains;

	/**
	 * The Molecule identifier, called entity_id in mmCIF dictionary
	 */
	private int molId;

	/**
	 * A map to cache residue number mapping, between ResidueNumbers and index (1-based) in aligned sequences (SEQRES).
	 * Initialised lazily upon call to {@link #getAlignedResIndex(Group, Chain)}
	 * Keys are asym_ids of chains, values maps of residue numbers to indices.
	 */
	private Map<String, Map<ResidueNumber,Integer>> chains2pdbResNums2ResSerials;

	private String refChainId;
	private String description = null;
	private String title = null;
	/**
	 * The type of entity (polymer, non-polymer, water)
	 */
	private EntityType type = null;
	private List<String> synonyms = null;
	private List<String> ecNums = null;
	private String engineered = null;
	private String mutation = null;
	private String biologicalUnit = null;
	private String details = null;
	private String numRes = null;
	private String resNames = null;
	private String headerVars = null;
	private String synthetic = null;
	private String fragment = null;
	private String organismScientific = null;
	private String organismTaxId = null;
	private String organismCommon = null;
	private String strain = null;
	private String variant = null;
	private String cellLine = null;
	private String atcc = null;
	private String organ = null;
	private String tissue = null;
	private String cell = null;
	private String organelle = null;
	private String secretion = null;
	private String gene = null;
	private String cellularLocation = null;
	private String expressionSystem = null;
	private String expressionSystemTaxId = null;
	private String expressionSystemStrain = null;
	private String expressionSystemVariant = null;
	private String expressionSystemCellLine = null;
	private String expressionSystemAtccNumber = null;
	private String expressionSystemOrgan = null;
	private String expressionSystemTissue = null;
	private String expressionSystemCell = null;
	private String expressionSystemOrganelle = null;
	private String expressionSystemCellularLocation = null;
	private String expressionSystemVectorType = null;
	private String expressionSystemVector = null;
	private String expressionSystemPlasmid = null;
	private String expressionSystemGene = null;
	private String expressionSystemOtherDetails = null;

	private Long id;

	public EntityInfo () {
		chains = new ArrayList<Chain>();
		chains2pdbResNums2ResSerials = new HashMap<String, Map<ResidueNumber,Integer>>();
		molId = -1;
	}

	/**
	 * Constructs a new EntityInfo copying all data from the given one
	 * but not setting the Chains
	 * @param c
	 */
	public EntityInfo (EntityInfo c) {

		this.chains = new ArrayList<Chain>();

		this.chains2pdbResNums2ResSerials = new HashMap<String, Map<ResidueNumber,Integer>>();

		this.molId = c.molId;
		
		this.type = c.type;

		this.refChainId = c.refChainId;

		this.description = c.description;
		this.title = c.title;

		if (c.synonyms!=null) {
			this.synonyms = new ArrayList<String>();
			synonyms.addAll(c.synonyms);
		}
		if (c.ecNums!=null) {
			this.ecNums = new ArrayList<String>();
			ecNums.addAll(c.ecNums);
		}

		this.engineered = c.engineered;
		this.mutation = c.mutation;
		this.biologicalUnit = c.biologicalUnit;
		this.details = c.details;

		this.numRes = c.numRes;
		this.resNames = c.resNames;

		this.headerVars = c.headerVars;

		this.synthetic = c.synthetic;
		this.fragment = c.fragment;
		this.organismScientific = c.organismScientific;
		this.organismTaxId = c.organismTaxId;
		this.organismCommon = c.organismCommon;
		this.strain = c.strain;
		this.variant = c.variant;
		this.cellLine = c.cellLine;
		this.atcc = c.atcc;
		this.organ = c.organ;
		this.tissue = c.tissue;
		this.cell = c.cell;
		this.organelle = c.organelle;
		this.secretion = c.secretion;
		this.gene = c.gene;
		this.cellularLocation = c.cellularLocation;
		this.expressionSystem = c.expressionSystem;
		this.expressionSystemTaxId = c.expressionSystemTaxId;
		this.expressionSystemStrain = c.expressionSystemStrain;
		this.expressionSystemVariant = c.expressionSystemVariant;
		this.expressionSystemCellLine = c.expressionSystemCellLine;
		this.expressionSystemAtccNumber = c.expressionSystemAtccNumber;
		this.expressionSystemOrgan = c.expressionSystemOrgan;
		this.expressionSystemTissue = c.expressionSystemTissue;
		this.expressionSystemCell = c.expressionSystemCell;
		this.expressionSystemOrganelle = c.expressionSystemOrganelle;
		this.expressionSystemCellularLocation = c.expressionSystemCellularLocation;
		this.expressionSystemVectorType = c.expressionSystemVectorType;
		this.expressionSystemVector = c.expressionSystemVector;
		this.expressionSystemPlasmid = c.expressionSystemPlasmid;
		this.expressionSystemGene = c.expressionSystemGene;
		this.expressionSystemOtherDetails = c.expressionSystemOtherDetails;


	}

	@Override
	public String toString(){
		StringBuilder buf = new StringBuilder();
		buf.append("EntityInfo: ").append(molId).append(" ");
		buf.append(description==null?"(no name)":"("+description+")");
		buf.append(" asymIds: ");
		if (chains!=null) {
			for (int i=0;i<chains.size();i++) {
				buf.append(chains.get(i).getId());
				if (i!=chains.size()-1) buf.append(",");
			}
		} else {
			buf.append("no chains");
		}
		return buf.toString();
	}

	/**
	 * Get the representative Chain for this EntityInfo.
	 * We choose the Chain with the first asym_id chain identifier after
	 * lexicographical sorting (case insensitive),
	 * e.g. chain A if EntityInfo is composed of chains A,B,C,D,E
	 * @return
	 */
	public Chain getRepresentative() {

		List<String> chainIds = new ArrayList<String>();
		for (Chain chain:chains) {
			chainIds.add(chain.getId());
		}

		Collections.sort(chainIds, String.CASE_INSENSITIVE_ORDER);

		for (Chain chain:chains) {
			if (chain.getId().equals(chainIds.get(0))) {
				return chain;
			}
		}

		logger.error("Could not find a representative chain for EntityInfo '{}'", this.toString());

		return null;
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

	/**
	 * Return the list of member chain ids (asym ids) that are described by this EntityInfo,
	 * only unique chain IDs are contained in the list.
	 * Note that in the case of multimodel structures this will return just the unique
	 * chain identifiers whilst {@link #getChains()} will return a corresponding chain
	 * per model.
	 * @return the list of unique ChainIDs that are described by this EnityInfo
	 * @see #setChains(List)
	 * @see #getChains()
	 */
	public List<String> getChainIds() {

		Set<String> uniqChainIds = new TreeSet<String>();
		for (int i=0;i<getChains().size();i++) {
			uniqChainIds.add(getChains().get(i).getId());
		}

		return new ArrayList<String>(uniqChainIds);
	}

	/**
	 * Given a Group g of Chain c (member of this EnityInfo) return the corresponding position in the
	 * alignment of all member sequences (1-based numbering), i.e. the index (1-based) in the SEQRES sequence.
	 * This allows for comparisons of residues belonging to different chains of the same EnityInfo (entity).
	 * <p>
	 * If {@link FileParsingParameters#setAlignSeqRes(boolean)} is not used or SEQRES not present, a mapping
	 * will not be available and this method will return {@link ResidueNumber#getSeqNum()} for all residues, which
	 * in some cases will be correctly aligned indices (when no insertion codes are
	 * used and when all chains within the entity are numbered in the same way), but
	 * in general they will be neither unique (because of insertion codes) nor aligned.
	 * </p>
	 * @param g
	 * @param c
	 * @return the aligned residue index (1 to n), if no SEQRES groups are available at all then {@link ResidueNumber#getSeqNum()}
	 * is returned as a fall-back, if the group is not found in the SEQRES groups then -1 is returned
	 * for the given group and chain
	 * @throws IllegalArgumentException if the given Chain is not a member of this EnityInfo
	 * @see {@link Chain#getSeqResGroup(int)}
	 */
	public int getAlignedResIndex(Group g, Chain c) {

		boolean contained = false;
		for (Chain member:getChains()) {
			if (c.getId().equals(member.getId())) {
				contained = true;
				break;
			}
		}
		if (!contained)
			throw new IllegalArgumentException("Given chain with asym_id "+c.getId()+" is not a member of this entity: "+getChainIds().toString());

		if (!chains2pdbResNums2ResSerials.containsKey(c.getId())) {
			// we do lazy initialisation of the map
			initResSerialsMap(c);
		}
		// if no seqres groups are available at all the map will be null
		Map<ResidueNumber,Integer> map = chains2pdbResNums2ResSerials.get(c.getId());
		int serial;
		if (map!=null) {

			ResidueNumber resNum = g.getResidueNumber();
			// the resNum will be null for groups that are SEQRES only and not in ATOM,
			// still it can happen that a group is in ATOM in one chain but not in other of the same entity.
			// This is what we try to find out here (analogously to what we do in initResSerialsMap() ):
			if (resNum==null && c.getSeqResGroups()!=null && !c.getSeqResGroups().isEmpty()) {
				
				int index = c.getSeqResGroups().indexOf(g);

				resNum = findResNumInOtherChains(index, c);

			}

			if (resNum == null) {
				// still null, we really can't map
				serial = -1;
			}
			else {

				Integer alignedSerial = map.get(resNum);

				if (alignedSerial==null) {
					// the map doesn't contain this group, something's wrong: return -1
					serial = -1;
				} else {
					serial = alignedSerial;
				}
			}

		} else {
			// no seqres groups available we resort to using the pdb residue numbers are given
			serial = g.getResidueNumber().getSeqNum();
		}
		return serial;
	}

	private void initResSerialsMap(Chain c) {
		if (c.getSeqResGroups()==null || c.getSeqResGroups().isEmpty()) {
			logger.warn("No SEQRES groups found in chain with asym_id {}, will use residue numbers as given (no insertion codes, not necessarily aligned). "
					+ "Make sure your structure has SEQRES records and that you use FileParsingParameters.setAlignSeqRes(true)",
					c.getId());
			// we add a explicit null to the map so that we flag it as unavailable for this chain
			chains2pdbResNums2ResSerials.put(c.getId(), null);
			return;
		}

		Map<ResidueNumber,Integer> resNums2ResSerials = new HashMap<ResidueNumber, Integer>();
		chains2pdbResNums2ResSerials.put(c.getId(), resNums2ResSerials);

		for (int i=0;i<c.getSeqResGroups().size();i++) {

			// The seqres group will have a null residue number whenever its corresponding atom group doesn't exist
			// because it is missing in the electron density.
			// However, it can be observed in the density in other chains of the same entity,
			// to be complete we go and look for the residue number in other chains, so that we have a
			// seqres to atom mapping as complete as possible (with all known atom groups of any chain of this entity)

			ResidueNumber resNum = c.getSeqResGroup(i).getResidueNumber();

			if (resNum==null) {
				resNum = findResNumInOtherChains(i,c);
			}

			// NOTE that resNum will still be null here for cases where the residue
			// is missing in atom groups (not observed in density) in all chains
			// Thus the mapping will not be possible for residues that are only in SEQRES groups
			resNums2ResSerials.put(resNum, i+1);
		}
	}

	private ResidueNumber findResNumInOtherChains(int i, Chain chain) {

		// note that getChains contains all chains from all models, we'll just use first model found and skip the others
		for (Chain c: getFirstModelChains()) {
			
			if (c == chain) continue;			
			
			Group seqResGroup = c.getSeqResGroup(i);

			if (seqResGroup==null) {
				logger.warn("The SEQRES group is null for index {} in chain with asym_id {}, whilst it wasn't null in chain with asym_id {}",
						 i, c.getId(), chain.getId());
				continue;
			}

			if (seqResGroup.getResidueNumber()!=null) return seqResGroup.getResidueNumber();

		}

		return null;
	}

	/**
	 * Return the ref chain id value.
	 * @return the RefChainID
	 * @see #setRefChainId(String)
	 */
	public String getRefChainId() {
		return refChainId;
	}

	/**
	 * Return the ref chain id value.
	 * @param refChainId the RefChainID
	 * @see #getRefChainId()
	 */
	public void setRefChainId(String refChainId) {
		this.refChainId = refChainId;
	}

	/**
	 * Return the molecule identifier, called entity_id in mmCIF dictionary.
	 * @return the molecule id
	 * @see #setMolId(int)
	 */
	public int getMolId() {
		return molId;
	}

	/**
	 * Set the molecule identifier, called entity_id in mmCIF dictionary.
	 * @param molId the molecule id
	 * @see #getMolId()
	 */
	public void setMolId(int molId) {
		this.molId = molId;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String molName) {
		this.description = molName;
	}

	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public List<String> getSynonyms() {
		return synonyms;
	}

	public void setSynonyms(List<String> synonyms) {
		this.synonyms = synonyms;
	}

	public List<String> getEcNums() {
		return ecNums;
	}

	public void setEcNums(List<String> ecNums) {
		this.ecNums = ecNums;
	}

	public String getEngineered() {
		return engineered;
	}

	public void setEngineered(String engineered) {
		this.engineered = engineered;
	}

	public String getMutation() {
		return mutation;
	}

	public void setMutation(String mutation) {
		this.mutation = mutation;
	}

	public String getBiologicalUnit() {
		return biologicalUnit;
	}

	public void setBiologicalUnit(String biologicalUnit) {
		this.biologicalUnit = biologicalUnit;
	}

	public String getDetails() {
		return details;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	public String getNumRes() {
		return numRes;
	}

	public void setNumRes(String numRes) {
		this.numRes = numRes;
	}

	public String getResNames() {
		return resNames;
	}

	public void setResNames(String resNames) {
		this.resNames = resNames;
	}

	public String getHeaderVars() {
		return headerVars;
	}

	public void setHeaderVars(String headerVars) {
		this.headerVars = headerVars;
	}

	public String getSynthetic() {
		return synthetic;
	}

	public void setSynthetic(String synthetic) {
		this.synthetic = synthetic;
	}

	public String getFragment() {
		return fragment;
	}

	public void setFragment(String fragment) {
		this.fragment = fragment;
	}

	public String getOrganismScientific() {
		return organismScientific;
	}

	public void setOrganismScientific(String organismScientific) {
		this.organismScientific = organismScientific;
	}

	public String getOrganismTaxId() {
		return organismTaxId;
	}

	public void setOrganismTaxId(String organismTaxId) {
		this.organismTaxId = organismTaxId;
	}

	public String getOrganismCommon() {
		return organismCommon;
	}

	public void setOrganismCommon(String organismCommon) {
		this.organismCommon = organismCommon;
	}

	public String getStrain() {
		return strain;
	}

	public void setStrain(String strain) {
		this.strain = strain;
	}

	public String getVariant() {
		return variant;
	}

	public void setVariant(String variant) {
		this.variant = variant;
	}

	public String getCellLine() {
		return cellLine;
	}

	public void setCellLine(String cellLine) {
		this.cellLine = cellLine;
	}

	public String getAtcc() {
		return atcc;
	}

	public void setAtcc(String atcc) {
		this.atcc = atcc;
	}

	public String getOrgan() {
		return organ;
	}

	public void setOrgan(String organ) {
		this.organ = organ;
	}

	public String getTissue() {
		return tissue;
	}

	public void setTissue(String tissue) {
		this.tissue = tissue;
	}

	public String getCell() {
		return cell;
	}

	public void setCell(String cell) {
		this.cell = cell;
	}

	public String getOrganelle() {
		return organelle;
	}

	public void setOrganelle(String organelle) {
		this.organelle = organelle;
	}

	public String getSecretion() {
		return secretion;
	}

	public void setSecretion(String secretion) {
		this.secretion = secretion;
	}

	public String getGene() {
		return gene;
	}

	public void setGene(String gene) {
		this.gene = gene;
	}

	public String getCellularLocation() {
		return cellularLocation;
	}

	public void setCellularLocation(String cellularLocation) {
		this.cellularLocation = cellularLocation;
	}

	public String getExpressionSystem() {
		return expressionSystem;
	}

	public String getExpressionSystemTaxId() {
		return expressionSystemTaxId;
	}

	public void setExpressionSystemTaxId(String expressionSystemTaxId) {
		this.expressionSystemTaxId = expressionSystemTaxId;
	}

	public void setExpressionSystem(String expressionSystem) {
		this.expressionSystem = expressionSystem;
	}

	public String getExpressionSystemStrain() {
		return expressionSystemStrain;
	}

	public void setExpressionSystemStrain(String expressionSystemStrain) {
		this.expressionSystemStrain = expressionSystemStrain;
	}

	public String getExpressionSystemVariant() {
		return expressionSystemVariant;
	}

	public void setExpressionSystemVariant(String expressionSystemVariant) {
		this.expressionSystemVariant = expressionSystemVariant;
	}

	public String getExpressionSystemCellLine() {
		return expressionSystemCellLine;
	}

	public void setExpressionSystemCellLine(String expressionSystemCellLine) {
		this.expressionSystemCellLine = expressionSystemCellLine;
	}

	public String getExpressionSystemAtccNumber() {
		return expressionSystemAtccNumber;
	}

	public void setExpressionSystemAtccNumber(String expressionSystemAtccNumber) {
		this.expressionSystemAtccNumber = expressionSystemAtccNumber;
	}

	public String getExpressionSystemOrgan() {
		return expressionSystemOrgan;
	}

	public void setExpressionSystemOrgan(String expressionSystemOrgan) {
		this.expressionSystemOrgan = expressionSystemOrgan;
	}

	public String getExpressionSystemTissue() {
		return expressionSystemTissue;
	}

	public void setExpressionSystemTissue(String expressionSystemTissue) {
		this.expressionSystemTissue = expressionSystemTissue;
	}

	public String getExpressionSystemCell() {
		return expressionSystemCell;
	}

	public void setExpressionSystemCell(String expressionSystemCell) {
		this.expressionSystemCell = expressionSystemCell;
	}

	public String getExpressionSystemOrganelle() {
		return expressionSystemOrganelle;
	}

	public void setExpressionSystemOrganelle(String expressionSystemOrganelle) {
		this.expressionSystemOrganelle = expressionSystemOrganelle;
	}

	public String getExpressionSystemCellularLocation() {
		return expressionSystemCellularLocation;
	}

	public void setExpressionSystemCellularLocation(String expressionSystemCellularLocation) {
		this.expressionSystemCellularLocation = expressionSystemCellularLocation;
	}

	public String getExpressionSystemVectorType() {
		return expressionSystemVectorType;
	}

	public void setExpressionSystemVectorType(String expressionSystemVectorType) {
		this.expressionSystemVectorType = expressionSystemVectorType;
	}

	public String getExpressionSystemVector() {
		return expressionSystemVector;
	}

	public void setExpressionSystemVector(String expressionSystemVector) {
		this.expressionSystemVector = expressionSystemVector;
	}

	public String getExpressionSystemPlasmid() {
		return expressionSystemPlasmid;
	}

	public void setExpressionSystemPlasmid(String expressionSystemPlasmid) {
		this.expressionSystemPlasmid = expressionSystemPlasmid;
	}

	public String getExpressionSystemGene() {
		return expressionSystemGene;
	}

	public void setExpressionSystemGene(String expressionSystemGene) {
		this.expressionSystemGene = expressionSystemGene;
	}

	public String getExpressionSystemOtherDetails() {
		return expressionSystemOtherDetails;
	}

	public void setExpressionSystemOtherDetails(String expressionSystemOtherDetails) {
		this.expressionSystemOtherDetails = expressionSystemOtherDetails;
	}

	/**
	 * Get the list of chains that are part of this EntityInfo. Note that for multi-model
	 * structures chains from all models are returned.
	 *
	 * @return a List of Chain objects
	 */
	public List<Chain> getChains(){
		return this.chains;
	}

	private List<Chain> getFirstModelChains() {

		Map<String, Chain> firstModelChains = new LinkedHashMap<>();
		Set<String> lookupChainIds = new HashSet<>(getChainIds());

		for (Chain chain : chains) {
			if (lookupChainIds.contains(chain.getId())) {
				if (!firstModelChains.containsKey(chain.getId())) {
					firstModelChains.put(chain.getId(), chain);
				}
			}
		}

		return new ArrayList<>(firstModelChains.values());
	}
	
	 /**
	  * Add new Chain to this EntityInfo
	  * @param chain
	  */
	public void addChain(Chain chain){
		this.chains.add(chain);
	}

	/**
	 * Set the chains for this EntityInfo
	 * @param chains
	 */
	public void setChains(List<Chain> chains){
		this.chains = chains;
	}

	/**
	 * Get the type of entity this EntityInfo describes.
	 * Options are polymer, non-polymer or water.
	 * @return a string describing the type of entity. (polymer, non-polymer or water).
	 */
	public EntityType getType() {
		return this.type;
	}

	/**
	 * Set the type of entity this EntityInfo describes.
	 * Options are polymer, non-polymer or water.
	 * @param type a string describing the type of entity. (polymer, non-polymer or water).
	 */
	public void setType(EntityType type) {
		this.type = type;
	}
}
 