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
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.SeqRes2AtomAligner;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Serializable;
import java.util.*;


/**
 * A Chain in a PDB file. It contains several groups which can be of
 * one of the types defined in the {@link GroupType} constants.
 *
 * @author Andreas Prlic
 * @author Jules Jacobsen
 * @since 1.4
 */
public class ChainImpl implements Chain, Serializable {

	private final static Logger logger = LoggerFactory.getLogger(ChainImpl.class);

	private static final long serialVersionUID = 1990171805277911840L;

	/**
	 * The default chain identifier used to be an empty space
	 */
	public static String DEFAULT_CHAIN_ID = "A";

	private String swissprot_id ;
	private String chainID ; // the chain identifier as in PDB files

	private List <Group> groups;
	private List<Group> seqResGroups;

	private Long id;
	private Compound mol;
	private Structure parent;

	private Map<String, Integer> pdbResnumMap;
	private String internalChainID; // the chain identifier used in mmCIF files


	private List<SeqMisMatch> seqMisMatches = null;
	/**
	 *  Constructs a ChainImpl object.
	 */
	public ChainImpl() {
		super();

		chainID = DEFAULT_CHAIN_ID;
		groups = new ArrayList<Group>() ;

		seqResGroups = new ArrayList<Group>();
		pdbResnumMap = new HashMap<String,Integer>();
		internalChainID = null;

	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public Long getId() {
		return id;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public void setId(Long id) {
		this.id = id;
	}

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

		n.setChainID( getChainID());
		n.setSwissprotId ( getSwissprotId());

		// NOTE the Compound will be reset at the parent level (Structure) if cloning is happening from parent level
		// here we don't deep-copy it and just keep the same reference, in case the cloning is happening at the Chain level only
		n.setCompound(this.mol);

		n.setInternalChainID(internalChainID);

		for (Group group : groups) {
			Group g = (Group) group.clone();
			n.addGroup(g);
			g.setChain(n);
		}

		if (!seqResGroups.isEmpty()){

			// cloning seqres and atom groups is ugly, due to their
			// nested relationship (some of the atoms can be in the seqres, but not all)

			List<Group> tmpSeqRes = new ArrayList<Group>();
			for (Group seqResGroup : seqResGroups) {
				Group g = (Group) seqResGroup.clone();
				g.setChain(n);
				tmpSeqRes.add(g);
			}

			Chain tmp = new ChainImpl();
			// that's a bit confusing, but that's how to set the seqres so that SeqRes2AtomAligner can use them
			tmp.setAtomGroups(tmpSeqRes);

			// now match them up..
			SeqRes2AtomAligner seqresaligner = new SeqRes2AtomAligner();

			seqresaligner.mapSeqresRecords(n, tmp);


		}


		return n ;
	}



	/** {@inheritDoc}
	 *
	 */
	@Override
	public void setCompound(Compound mol) {
		this.mol = mol;
	}

	/** {@inheritDoc}
	 *
	 */
	@Override
	public Compound getCompound() {
		return this.mol;
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
         and ...
		 * HETATM 1348  O   HOH    92     -21.853 -16.886  19.138  1.00 66.92           O
		 * HETATM 1349  O   HOH    93     -26.126   1.226  29.069  1.00 71.69           O
		 * HETATM 1350  O   HOH    94     -22.250 -18.060  -6.401  1.00 61.97           O
		 */

			// this check is to give in this case the entry priority that is an AminoAcid / comes first...
			if (  pdbResnumMap.containsKey(pdbResnum)) {
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

		List<Group> tmp = new ArrayList<Group>() ;
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
	@Deprecated // TODO dmyersturnbull: why is this deprecated if it's declared in Chain?
	public Group[] getGroupsByPDB(ResidueNumber start, ResidueNumber end, boolean ignoreMissing)
			throws StructureException {

		if (! ignoreMissing )
			return getGroupsByPDB(start, end);


		List<Group> retlst = new ArrayList<Group>();

		String pdbresnumStart = start.toString();
		String pdbresnumEnd   = end.toString();


		int startPos = Integer.MIN_VALUE;
		int endPos   = Integer.MAX_VALUE;


		startPos = start.getSeqNum();
		endPos   = end.getSeqNum();



		boolean adding = false;
		boolean foundStart = false;

		for (Group g: groups){

			if ( g.getResidueNumber().toString().equals(pdbresnumStart)) {
				adding = true;
				foundStart = true;
			}

			if ( ! (foundStart && adding) ) {


				int pos = g.getResidueNumber().getSeqNum();

				if ( pos >= startPos) {
					foundStart = true;
					adding = true;
				}


			}

			if ( adding)
				retlst.add(g);

			if ( g.getResidueNumber().toString().equals(pdbresnumEnd)) {
				if ( ! adding)
					throw new StructureException("did not find start PDB residue number " + pdbresnumStart + " in chain " + chainID);
				adding = false;
				break;
			}
			if (adding){

				int pos = g.getResidueNumber().getSeqNum();
				if (pos >= endPos) {
					adding = false;
					break;
				}

			}
		}

		if ( ! foundStart){
			throw new StructureException("did not find start PDB residue number " + pdbresnumStart + " in chain " + chainID);
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
			throw new StructureException("unknown PDB residue number " + pdbresnum + " in chain " + chainID);
		}
	}

	/**
	 * {@inheritDoc}
	 *
	 */
	@Override
	public Group[] getGroupsByPDB(ResidueNumber start, ResidueNumber end)
			throws StructureException {

		String pdbresnumStart = start.toString();
		String pdbresnumEnd   = end.toString();

		List<Group> retlst = new ArrayList<Group>();

		Iterator<Group> iter = groups.iterator();
		boolean adding = false;
		boolean foundStart = false;

		while ( iter.hasNext()){
			Group g = iter.next();
			if ( g.getResidueNumber().toString().equals(pdbresnumStart)) {
				adding = true;
				foundStart = true;
			}

			if ( adding)
				retlst.add(g);

			if ( g.getResidueNumber().toString().equals(pdbresnumEnd)) {
				if ( ! adding)
					throw new StructureException("did not find start PDB residue number " + pdbresnumStart + " in chain " + chainID);
				adding = false;
				break;
			}
		}

		if ( ! foundStart){
			throw new StructureException("did not find start PDB residue number " + pdbresnumStart + " in chain " + chainID);
		}
		if ( adding) {
			throw new StructureException("did not find end PDB residue number " + pdbresnumEnd + " in chain " + chainID);
		}

		return retlst.toArray(new Group[retlst.size()] );
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
	public void   setChainID(String nam) { chainID = nam;   }


	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getChainID()           {	return chainID;  }



	/** String representation.
	 * @return String representation of the Chain
	 */
	@Override
	public String toString(){
		String newline = System.getProperty("line.separator");
		StringBuilder str = new StringBuilder();
		str.append("Chain >").append(getChainID()).append("<").append(newline);
		if ( mol != null ){
			if ( mol.getMolName() != null){
				str.append(mol.getMolName()).append(newline);
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

		String prop = System.getProperty(PDBFileReader.LOAD_CHEM_COMP_PROPERTY);

		if ( prop != null && prop.equalsIgnoreCase("true")){

			logger.info("The property {} is true. Will use the chemical component dictionary files to produce the atom sequence",
					PDBFileReader.LOAD_CHEM_COMP_PROPERTY);

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

		logger.info("The property {} is false or not set. Will use only amino acids for the atom sequence",
				PDBFileReader.LOAD_CHEM_COMP_PROPERTY);

		// not using ChemComp records...
		List<Group> aminos = getAtomGroups(GroupType.AMINOACID);
		StringBuilder sequence = new StringBuilder() ;
		for (Group amino : aminos) {
			if (amino instanceof AminoAcid) {
				AminoAcid a = (AminoAcid) amino;
				sequence.append(a.getAminoType());
			} else {
				// I suppose this can't happen, but just in case...
				logger.warn("Group {} is tagged as GroupType.AMINOACID but its class is not AminoAcid", amino.toString());
			}
		}

		return sequence.toString();

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getSeqResSequence(){

		String prop = System.getProperty(PDBFileReader.LOAD_CHEM_COMP_PROPERTY);

		if ( prop != null && prop.equalsIgnoreCase("true")) {

			logger.info("The property {} is true. Will use the chemical component dictionary files to produce the seqres sequence",
					PDBFileReader.LOAD_CHEM_COMP_PROPERTY);

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
					if ( oneLetter == null || oneLetter.isEmpty() || oneLetter.equals("?"))
						oneLetter = Character.toString(StructureTools.UNKNOWN_GROUP_LABEL);
					str.append(oneLetter);
				} else {
					str.append(StructureTools.UNKNOWN_GROUP_LABEL);
				}
			}
			return str.toString();
		}

		logger.info("The property {} is false or not set. Will use only amino acids for the seqres sequence",
				PDBFileReader.LOAD_CHEM_COMP_PROPERTY);

		StringBuilder str = new StringBuilder();
		for (Group group : seqResGroups) {
			if (group instanceof AminoAcid) {
				AminoAcid aa = (AminoAcid)group;
				str.append(aa.getAminoType()) ;
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
		List<Group> tmp = new ArrayList<Group>() ;
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
			g.setChain(this);
		}
		this.seqResGroups = groups;
	}

	protected void addSeqResGroup(Group g){
		seqResGroups.add(g);
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
		List<Group> ligands = new ArrayList<Group>();

		for (Group g : groups)
			if (!seqResGroups.contains(g) && !g.isWater())
				ligands.add(g);

		return ligands;
	}

	@Override
	public String getInternalChainID() {
		return internalChainID;
	}

	@Override
	public void setInternalChainID(String internalChainID) {
		this.internalChainID = internalChainID;

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
}

