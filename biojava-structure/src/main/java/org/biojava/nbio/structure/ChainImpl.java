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
 * Created on 12.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.nbio.structure;


import org.biojava.nbio.structure.io.FileConvert;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;


/**
 * A Chain in a PDB file. It contains several groups which can be of
 * one of the types defined in the {@link GroupType} constants.
 *
 * @author Andreas Prlic
 * @author Jules Jacobsen
 * @since 1.4
 */
public class ChainImpl implements Chain {

	private final static Logger logger = LoggerFactory.getLogger(ChainImpl.class);

	private static final long serialVersionUID = 1990171805277911840L;

	/**
	 * The default chain identifier used to be an empty space
	 */
	private static final String DEFAULT_CHAIN_ID = "A";

	private String swissprot_id ;
	private String authId; // the 'public' chain identifier as assigned by authors in PDB files

	private List <Group> groups;
	private List<Group> seqResGroups;

	private EntityInfo entity;
	private Structure parent;

	private Map<String, Integer> pdbResnumMap;
	private String asymId; // the 'internal' chain identifier as used in mmCIF files


	private List<SeqMisMatch> seqMisMatches = null;
	/**
	 *  Constructs a ChainImpl object.
	 */
	public ChainImpl() {
		super();

		authId = DEFAULT_CHAIN_ID;
		groups = new ArrayList<>() ;

		seqResGroups = new ArrayList<>();
		pdbResnumMap = new HashMap<>();
		asymId = null;

	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public String getId() {
		return asymId;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public void setId(String asymId) {
		this.asymId = asymId;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public String getName() { return authId; }

	/** {@inheritDoc}
	 *
	 */
	@Override
	public void setName(String authId) { this.authId = authId; }

	/** {@inheritDoc}
	 *
	 */
	@Override
	@Deprecated
	public void setParent(Structure parent) {
		setStructure(parent);
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public void setStructure(Structure parent){
		this.parent = parent;
	}

	/** Returns the parent Structure of this chain.
	 *
	 * @return the parent Structure object
	 */
	@Override
	public Structure getStructure() {

		return parent;
	}


	/** Returns the parent Structure of this chain.
	 *
	 * @return the parent Structure object
	 * @deprecated  use getStructure instead.
	 */
	@Override
	@Deprecated
	public Structure getParent() {

		return getStructure();
	}

	/** Returns an identical copy of this Chain .
	 * @return an identical copy of this Chain
	 */
	@Override
	public Object clone() {
		// go through all groups and add to new Chain.
		ChainImpl n = new ChainImpl();
		// copy chain data:

		n.setId(getId());
		n.setName(getName());
		n.setSwissprotId ( getSwissprotId());

		// NOTE the EntityInfo will be reset at the parent level (Structure) if cloning is happening from parent level
		// here we don't deep-copy it and just keep the same reference, in case the cloning is happening at the Chain level only
		n.setEntityInfo(this.entity);


		for (Group group : groups) {
			Group g = (Group) group.clone();
			n.addGroup(g);
			g.setChain(n);
		}

		if (seqResGroups!=null){

			List<Group> tmpSeqRes = new ArrayList<>();

			// cloning seqres and atom groups is ugly, due to their
			// nested relationship (some of the atoms can be in the seqres, but not all)

			for (Group seqResGroup : seqResGroups) {

				if (seqResGroup==null) {
					tmpSeqRes.add(null);
					continue;
				}

				int i = groups.indexOf(seqResGroup);

				Group g ;

				if (i!=-1) {
					// group found in atom groups, we get the equivalent reference from the newly cloned atom groups
					g = n.getAtomGroup(i);
				} else {
					// group not found in atom groups, we clone the seqres group
					g = (Group) seqResGroup.clone();
				}
				g.setChain(n);
				tmpSeqRes.add(g);
			}

			n.setSeqResGroups(tmpSeqRes);
		}

		return n ;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public void setEntityInfo(EntityInfo mol) {
		this.entity = mol;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public EntityInfo getEntityInfo() {
		return this.entity;
	}

	/** set the Swissprot id of this chains .
	 * @param sp_id  a String specifying the swissprot id value
	 * @see #getSwissprotId
	 */
	@Override
	public void setSwissprotId(String sp_id){
		swissprot_id = sp_id ;
	}

	/** get the Swissprot id of this chains .
	 * @return a String representing the swissprot id value
	 * @see #setSwissprotId
	 */
	@Override
	public String getSwissprotId() {
		return swissprot_id ;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public void addGroup(Group group) {

		group.setChain(this);

		// Set the altlocs chain as well
		for(Group g : group.getAltLocs()) {
			g.setChain(this);
		}

		groups.add(group);

		// store the position internally for quick access of this group

		String pdbResnum = null ;
		ResidueNumber resNum = group.getResidueNumber();
		if ( resNum != null)
			pdbResnum = resNum.toString();
		if ( pdbResnum != null) {
			Integer pos = groups.size() - 1;
			// ARGH sometimes numbering in PDB files is confusing.
			// e.g. PDB: 1sfe
			/*
			 * ATOM    620  N   GLY    93     -24.320  -6.591   4.210  1.00 46.82           N
			 * ATOM    621  CA  GLY    93     -24.960  -6.849   5.497  1.00 47.35           C
			 * ATOM    622  C   GLY    93     -26.076  -5.873   5.804  1.00 47.24           C
			 * ATOM    623  O   GLY    93     -26.382  -4.986   5.006  1.00 47.56           O
			 *    and ...
			 * HETATM 1348  O   HOH    92     -21.853 -16.886  19.138  1.00 66.92           O
			 * HETATM 1349  O   HOH    93     -26.126   1.226  29.069  1.00 71.69           O
			 * HETATM 1350  O   HOH    94     -22.250 -18.060  -6.401  1.00 61.97           O
			 */

			// this check is to give in this case the entry priority that is an AminoAcid / comes first...
			// a good example of same residue number for 2 residues is 3th3, chain T, residue 201 (a LYS and a sugar BGC covalently attached to it) - JD 2016-03-09
			if (  pdbResnumMap.containsKey(pdbResnum)) {

				logger.warn("Adding residue {}({}) to chain {} but a residue with same residue number is already present: {}({}). Will add only the aminoacid residue (if any) to the lookup, lookups for that residue number won't work properly.",
						pdbResnum, group.getPDBName(), getChainID(), groups.get(pdbResnumMap.get(pdbResnum)).getResidueNumber(), groups.get(pdbResnumMap.get(pdbResnum)).getPDBName());
				if ( group instanceof AminoAcid)
					pdbResnumMap.put(pdbResnum,pos);
			} else
				pdbResnumMap.put(pdbResnum,pos);
		}

	}


	/**
	 * {@inheritDoc}
	 */
	@Override
	public Group getAtomGroup(int position) {

		return groups.get(position);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public List<Group> getAtomGroups(GroupType type){

		List<Group> tmp = new ArrayList<>() ;
		for (Group g : groups) {
			if (g.getType().equals(type)) {
				tmp.add(g);
			}
		}

		return tmp ;
	}


	/** {@inheritDoc}
	 *
	 */
	@Override
	public List<Group> getAtomGroups(){
		return groups ;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public void setAtomGroups(List<Group> groups){
		for (Group g:groups){
			g.setChain(this);
		}
		this.groups = groups;
	}

	@Override
	public Group[] getGroupsByPDB(ResidueNumber start, ResidueNumber end, boolean ignoreMissing)
			throws StructureException {
		// Short-circut for include all groups
		if(start == null && end == null) {
			return groups.toArray(new Group[groups.size()]);
		}


		List<Group> retlst = new ArrayList<>();

		boolean adding, foundStart;
		if( start == null ) {
			// start with first group
			adding = true;
			foundStart = true;
		} else {
			adding = false;
			foundStart = false;
		}

		
		for (Group g: groups){

			// Check for start
			if (!adding && start.equalsPositional(g.getResidueNumber())) {
				adding = true;
				foundStart = true;
			}

			// Check if past start
			if ( ignoreMissing && ! (foundStart && adding) ) {
				ResidueNumber pos = g.getResidueNumber();

				if ( start != null && start.compareToPositional(pos) <= 0) {
					foundStart = true;
					adding = true;
				}
			}

			if ( adding)
				retlst.add(g);

			// check for end
			if ( end != null && end.equalsPositional(g.getResidueNumber())) {
				if ( ! adding)
					throw new StructureException("did not find start PDB residue number " + start + " in chain " + authId);
				adding = false;
				break;
			}
			// check if past end
			if ( ignoreMissing && adding && end != null){

				ResidueNumber pos = g.getResidueNumber();
				if ( end.compareToPositional(pos) <= 0) {
					adding = false;
					break;
				}

			}
		}

		if ( ! foundStart){
			throw new StructureException("did not find start PDB residue number " + start + " in chain " + authId);
		}
		if ( end != null && adding && !ignoreMissing) {
			throw new StructureException("did not find end PDB residue number " + end + " in chain " + authId);
		}


		//not checking if the end has been found in this case...

		return retlst.toArray(new Group[retlst.size()] );
	}


	/**
	 * {@inheritDoc}
	 *
	 */
	@Override
	public Group getGroupByPDB(ResidueNumber resNum) throws StructureException {
		String pdbresnum = resNum.toString();
		if ( pdbResnumMap.containsKey(pdbresnum)) {
			Integer pos = pdbResnumMap.get(pdbresnum);
			return groups.get(pos);
		} else {
			throw new StructureException("unknown PDB residue number " + pdbresnum + " in chain " + authId);
		}
	}

	/**
	 * {@inheritDoc}
	 *
	 */
	@Override
	public Group[] getGroupsByPDB(ResidueNumber start, ResidueNumber end)
			throws StructureException {
		return getGroupsByPDB(start, end, false);
	}



	/**
	 * {@inheritDoc}
	 */
	@Override
	public int getSeqResLength() {
		//new method returns the length of the sequence defined in the SEQRES records
		return seqResGroups.size();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void   setChainID(String asymId) { this.asymId = asymId;   }


	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getChainID()           {	return this.asymId;  }



	/** String representation.
	 * @return String representation of the Chain
	 */
	@Override
	public String toString(){
		String newline = System.getProperty("line.separator");
		StringBuilder str = new StringBuilder();
		str.append("Chain asymId:").append(getChainID()).append(" authId:").append(getName()).append(newline);
		if ( entity != null ){
			if ( entity.getDescription() != null){
				str.append(entity.getDescription()).append(newline);
			}
		}
		str.append("total SEQRES length: ").append(getSeqResGroups().size()).append(" total ATOM length:")
		.append(getAtomLength()).append(" residues ").append(newline);

		return str.toString() ;

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Sequence<?> getBJSequence()  {

		String seq = getSeqResSequence();

		Sequence<AminoAcidCompound> s = null;

		try {
			s = new ProteinSequence(seq);
		} catch (CompoundNotFoundException e) {
			logger.error("Could not create sequence object from seqres sequence. Some unknown compound: {}",e.getMessage());
		}

		//TODO: return a DNA sequence if the content is DNA...
		return s;

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getAtomSequence(){


		List<Group> groups = getAtomGroups();
		StringBuilder sequence = new StringBuilder() ;

		for ( Group g: groups){
			ChemComp cc = g.getChemComp();

			if ( PolymerType.PROTEIN_ONLY.contains(cc.getPolymerType()) ||
					PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())){
				// an amino acid residue.. use for alignment
				String oneLetter= ChemCompGroupFactory.getOneLetterCode(cc);
				if ( oneLetter == null)
					oneLetter = Character.toString(StructureTools.UNKNOWN_GROUP_LABEL);
				sequence.append(oneLetter);
			}

		}
		return sequence.toString();


	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getSeqResSequence(){

		StringBuilder str = new StringBuilder();
		for (Group g : seqResGroups) {
			ChemComp cc = g.getChemComp();
			if ( cc == null) {
				logger.warn("Could not load ChemComp for group: ", g);
				str.append(StructureTools.UNKNOWN_GROUP_LABEL);
			} else if ( PolymerType.PROTEIN_ONLY.contains(cc.getPolymerType()) ||
					PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())){
				// an amino acid residue.. use for alignment
				String oneLetter= ChemCompGroupFactory.getOneLetterCode(cc);
				// AB oneLetter.length() should be one. e.g. in 1EMA it is 3 and this makes mapping residue to sequence impossible.
				if ( oneLetter == null || oneLetter.isEmpty() || oneLetter.equals("?")) {
					oneLetter = Character.toString(StructureTools.UNKNOWN_GROUP_LABEL);
				}
				str.append(oneLetter);
			} else {
				str.append(StructureTools.UNKNOWN_GROUP_LABEL);
			}
		}
		return str.toString();
	}
	
	/**
	 * Get the one letter sequence so that Sequence is guaranteed to
	 * be the same length as seqResGroups.
	 * Method related to https://github.com/biojava/biojava/issues/457
	 * @return a string of the sequence guaranteed to be the same length
	 * as seqResGroups.
	 */
	public String getSeqResOneLetterSeq(){

		StringBuilder str = new StringBuilder();
		for (Group g : seqResGroups) {
			ChemComp cc = g.getChemComp();
			if ( cc == null) {
				logger.warn("Could not load ChemComp for group: ", g);
				str.append(StructureTools.UNKNOWN_GROUP_LABEL);
			} else if ( PolymerType.PROTEIN_ONLY.contains(cc.getPolymerType()) ||
					PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())){
				// an amino acid residue.. use for alignment
				String oneLetter= ChemCompGroupFactory.getOneLetterCode(cc);
				// AB oneLetter.length() should be one. e.g. in 1EMA it is 3 and this makes mapping residue to sequence impossible.
				if ( oneLetter == null || oneLetter.isEmpty() || oneLetter.equals("?") || oneLetter.length()!=1) {
					oneLetter = Character.toString(StructureTools.UNKNOWN_GROUP_LABEL);
				}
				str.append(oneLetter);
			} else {
				str.append(StructureTools.UNKNOWN_GROUP_LABEL);
			}
		}
		return str.toString();
	}


	/**
	 * {@inheritDoc}
	 */
	@Override
	public Group getSeqResGroup(int position) {

		return seqResGroups.get(position);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public List<Group> getSeqResGroups(GroupType type) {
		List<Group> tmp = new ArrayList<>() ;
		for (Group g : seqResGroups) {
			if (g.getType().equals(type)) {
				tmp.add(g);
			}
		}

		return tmp ;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public List<Group> getSeqResGroups() {
		return seqResGroups;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public void setSeqResGroups(List<Group> groups){
		for (Group g: groups){
			if (g != null) {
				g.setChain(this);
			}
		}
		this.seqResGroups = groups;
	}


	/** {@inheritDoc}
	 *
	 */
	@Override
	public int getAtomLength() {

		return groups.size();
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public List<Group> getAtomLigands(){
		List<Group> ligands = new ArrayList<>();

		for (Group g : groups)
			if (!seqResGroups.contains(g) && !g.isWater())
				ligands.add(g);

		return ligands;
	}

	@Override
	public String getInternalChainID() {
		return asymId;
	}

	@Override
	public void setInternalChainID(String internalChainID) {
		this.asymId = internalChainID;

	}

	@Override
	public String toPDB() {
		return FileConvert.toPDB(this);
	}

	@Override
	public String toMMCIF() {
		return FileConvert.toMMCIF(this, true);
	}

	@Override
	public void setSeqMisMatches(List<SeqMisMatch> seqMisMatches) {
		this.seqMisMatches = seqMisMatches;
	}

	@Override
	public List<SeqMisMatch> getSeqMisMatches() {
		return seqMisMatches;
	}
	
	@Override
	public EntityType getEntityType() {
		if (getEntityInfo()==null) return null;
		return getEntityInfo().getType();
	}

	@Override
	public boolean isWaterOnly() {
		for (Group g : getAtomGroups()) {
			if (!g.isWater())
				return false;
		}
		return true;
	}

	@Override
	public boolean isPureNonPolymer() {
		for (Group g : getAtomGroups()) {

			//ChemComp cc = g.getChemComp();

			if ( 	g.isPolymeric() &&
					!g.isHetAtomInFile() ) {

				// important: the aminoacid or nucleotide residue can be in Atom records

				return false;
			}

		}
		return true;
	}

	@Override
	public GroupType getPredominantGroupType(){

		double ratioResiduesToTotal = StructureTools.RATIO_RESIDUES_TO_TOTAL;

		int sizeAminos = getAtomGroups(GroupType.AMINOACID).size();
		int sizeNucleotides = getAtomGroups(GroupType.NUCLEOTIDE).size();
		List<Group> hetAtoms = getAtomGroups(GroupType.HETATM);
		int sizeHetatoms = hetAtoms.size();
		int sizeWaters = 0;
		for (Group g : hetAtoms) {
			if (g.isWater())
				sizeWaters++;
		}
		int sizeHetatomsWithoutWater = sizeHetatoms - sizeWaters;

		int fullSize = sizeAminos + sizeNucleotides + sizeHetatomsWithoutWater;

		if ((double) sizeAminos / (double) fullSize > ratioResiduesToTotal)
			return GroupType.AMINOACID;

		if ((double) sizeNucleotides / (double) fullSize > ratioResiduesToTotal)
			return GroupType.NUCLEOTIDE;

		if ((double) (sizeHetatomsWithoutWater) / (double) fullSize > ratioResiduesToTotal)
			return GroupType.HETATM;

		// finally if neither condition works, we try based on majority, but log
		// it
		GroupType max;
		if (sizeNucleotides > sizeAminos) {
			if (sizeNucleotides > sizeHetatomsWithoutWater) {
				max = GroupType.NUCLEOTIDE;
			} else {
				max = GroupType.HETATM;
			}
		} else {
			if (sizeAminos > sizeHetatomsWithoutWater) {
				max = GroupType.AMINOACID;
			} else {
				max = GroupType.HETATM;
			}
		}
		logger.debug(
				"Ratio of residues to total for chain with asym_id {} is below {}. Assuming it is a {} chain. "
						+ "Counts: # aa residues: {}, # nuc residues: {}, # non-water het residues: {}, # waters: {}, "
						+ "ratio aa/total: {}, ratio nuc/total: {}",
				getId(), ratioResiduesToTotal, max, sizeAminos,
				sizeNucleotides, sizeHetatomsWithoutWater, sizeWaters,
				(double) sizeAminos / (double) fullSize,
				(double) sizeNucleotides / (double) fullSize);

		return max;
	}

	@Override
	public  boolean isProtein() {
		return getPredominantGroupType() == GroupType.AMINOACID;
	}

	@Override
	public  boolean isNucleicAcid() {
		return getPredominantGroupType() == GroupType.NUCLEOTIDE;
	}


}

