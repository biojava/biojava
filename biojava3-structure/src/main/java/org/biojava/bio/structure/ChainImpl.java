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
package org.biojava.bio.structure;


import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.SeqRes2AtomAligner;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava3.core.exceptions.CompoundNotFoundException;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.template.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


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
	public Long getId() {
		return id;
	}

	/** {@inheritDoc}
	 *
	 */
	public void setId(Long id) {
		this.id = id;
	}

	/** {@inheritDoc}
	 *
	 */
	public void setParent(Structure parent) {
		this.parent = parent;
	}

	/** Returns the parent Structure of this chain.
	 *
	 * @return the parent Structure object
	 */

	public Structure getParent() {


		return parent;
	}


	/** Returns an identical copy of this Chain .
	 * @return an identical copy of this Chain
	 */
	public Object clone() {
		// go through all groups and add to new Chain.
		ChainImpl n = new ChainImpl();
		// copy chain data:

		n.setChainID( getChainID());
		n.setSwissprotId ( getSwissprotId());
		// TODO should the Compound be deep copied too? I'd keep it like this (actually StructureInterface.getContactOverlapScore depends on it NOT being deep copied!!) - JD 2014.12.10
		n.setCompound(this.getCompound());
		n.setInternalChainID(internalChainID);

		for (int i=0;i<groups.size();i++){
			Group g = groups.get(i);
			n.addGroup((Group)g.clone());
		}
		
		if (seqResGroups.size() > 0 ){

			// cloning seqres and atom groups is ugly, due to their
			// nested relationship (some of the atoms can be in the seqres, but not all)

			List<Group> tmpSeqRes = new ArrayList<Group>();
			for (int i=0;i<seqResGroups.size();i++){
				Group g = (Group)seqResGroups.get(i).clone();
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
	public void setCompound(Compound mol) {
		this.mol = mol;
	}

	/** {@inheritDoc}
	 *
	 */
	public Compound getCompound() {
		return this.mol;
	}

	/** set the Swissprot id of this chains .
	 * @param sp_id  a String specifying the swissprot id value
	 * @see #getSwissprotId
	 */

	public void setSwissprotId(String sp_id){
		swissprot_id = sp_id ;
	}

	/** get the Swissprot id of this chains .
	 * @return a String representing the swissprot id value
	 * @see #setSwissprotId
	 */
	public String getSwissprotId() {
		return swissprot_id ;
	}

	/** {@inheritDoc}
	 *
	 */
	public void addGroup(Group group) {

		group.setChain(this);

		groups.add(group);

		// store the position internally for quick access of this group

		String pdbResnum = null ;
		ResidueNumber resNum = group.getResidueNumber();
		if ( resNum != null)
			pdbResnum = resNum.toString();
		if ( pdbResnum != null) {
			Integer pos = new Integer(groups.size()-1);
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
	public Group getAtomGroup(int position) {

		return (Group)groups.get(position);
	}

	/**  
	 * {@inheritDoc}
	 */
	public List<Group> getAtomGroups(String type){
		List<Group> tmp = new ArrayList<Group>() ;
		for (int i=0;i<groups.size();i++){
			Group g = (Group)groups.get(i);
			if (g.getType().equals(type)){
				tmp.add(g);
			}
		}

		return tmp ;
	}


	/** {@inheritDoc}
	 *
	 */
	public List<Group> getAtomGroups(){
		return groups ;
	}

	/** {@inheritDoc}
	 *
	 */
	public void setAtomGroups(List<Group> groups){
		for (Group g:groups){
			g.setChain(this);
		}
		this.groups = groups;
	}

	/** {@inheritDoc}
	 *
	 */
	public Group[] getGroupsByPDB(String pdbresnumStart, String pdbresnumEnd, boolean ignoreMissing)
			throws StructureException {

		ResidueNumber start = ResidueNumber.fromString(pdbresnumStart);
		ResidueNumber end = ResidueNumber.fromString(pdbresnumEnd);

		if (! ignoreMissing )
			return getGroupsByPDB(start, end);

		return getGroupsByPDB(start, end, ignoreMissing);

	}

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

		return (Group[]) retlst.toArray(new Group[retlst.size()] );
	}


	/**
	 * {@inheritDoc}
	 *
	 */
	public Group getGroupByPDB(String pdbresnum) throws StructureException {
		ResidueNumber resNum = ResidueNumber.fromString(pdbresnum);
		return getGroupByPDB(resNum);

	}

	/**
	 * {@inheritDoc}
	 *
	 */
	public Group getGroupByPDB(ResidueNumber resNum) throws StructureException {
		String pdbresnum = resNum.toString();
		if ( pdbResnumMap.containsKey(pdbresnum)) {
			Integer pos = (Integer) pdbResnumMap.get(pdbresnum);
			return (Group) groups.get(pos.intValue());
		} else {
			throw new StructureException("unknown PDB residue number " + pdbresnum + " in chain " + chainID);
		}
	}

	/**
	 * {@inheritDoc}
	 *
	 */
	public Group[] getGroupsByPDB(String pdbresnumStart, String pdbresnumEnd)
			throws StructureException {
		ResidueNumber start = ResidueNumber.fromString(pdbresnumStart);
		ResidueNumber end = ResidueNumber.fromString(pdbresnumEnd);

		return getGroupsByPDB(start,end);
	}


	/**
	 * {@inheritDoc}
	 *
	 */
	public Group[] getGroupsByPDB(ResidueNumber start, ResidueNumber end)
			throws StructureException {

		String pdbresnumStart = start.toString();
		String pdbresnumEnd   = end.toString();

		List<Group> retlst = new ArrayList<Group>();

		Iterator<Group> iter = groups.iterator();
		boolean adding = false;
		boolean foundStart = false;

		while ( iter.hasNext()){
			Group g = (Group) iter.next();
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

		return (Group[]) retlst.toArray(new Group[retlst.size()] );
	}



	/** {@inheritDoc}
	 *
	 */
	public int getLengthAminos() {

		List<Group> g = getAtomGroups(GroupType.AMINOACID);
		return g.size() ;
	}

	/**
	 * {@inheritDoc}
	 */
	public int getSeqResLength() {
		//new method returns the length of the sequence defined in the SEQRES records
		return seqResGroups.size();
	}

	/**
	 * {@inheritDoc}
	 */
	public void   setChainID(String nam) { chainID = nam;   }


	/**
	 * {@inheritDoc}
	 */
	public String getChainID()           {	return chainID;  }



	/** String representation.
	 * @return String representation of the Chain
	 *  */
	public String toString(){
		String newline = System.getProperty("line.separator");
		StringBuilder str = new StringBuilder();
		str.append("Chain >"+getChainID()+"<"+newline) ;
		if ( mol != null ){
			if ( mol.getMolName() != null){
				str.append(mol.getMolName()).append(newline);
			}
		}
		str.append("total SEQRES length: " + getSeqResGroups().size() +
				" total ATOM length:" + getAtomLength() + " residues " + newline);

		// commented out the looping over residues, I thought it didn't help much, especially in debugging - JD 2014-12-10
		// loop over the residues
		//for ( int i = 0 ; i < seqResGroups.size();i++){
		//	Group gr = (Group) seqResGroups.get(i);
		//	str.append(gr.toString()).append(newline);
		//}
		return str.toString() ;

	}

	/** Convert the SEQRES groups of a Chain to a Biojava Sequence object.
	 *
	 * @return the SEQRES groups of the Chain as a Sequence object.
	 */
	public Sequence<?> getBJSequence()  {

		//List<Group> groups = c.getSeqResGroups();
		String seq = getSeqResSequence();

		//		String name = "";
		//		if ( this.getParent() != null )
		//			name = getParent().getPDBCode();
		//		name += "." + getName();

		Sequence<AminoAcidCompound> s = null;

		try {
			s = new ProteinSequence(seq);
		} catch (CompoundNotFoundException e) {
			logger.error("Could not create sequence object from seqres sequence. Some unknown compound: {}",e.getMessage());
		}

		//TODO: return a DNA sequence if the content is DNA...
		return s;

	}

	/** {@inheritDoc}
	 *
	 */
	public String getAtomSequence(){

		String prop = System.getProperty(PDBFileReader.LOAD_CHEM_COMP_PROPERTY);

		if ( prop != null && prop.equalsIgnoreCase("true")){


			List<Group> groups = getAtomGroups();
			StringBuffer sequence = new StringBuffer() ;

			for ( Group g: groups){
				ChemComp cc = g.getChemComp();

				if ( PolymerType.PROTEIN_ONLY.contains(cc.getPolymerType()) ||
						PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())){
					// an amino acid residue.. use for alignment
					String oneLetter= ChemCompGroupFactory.getOneLetterCode(cc);
					if ( oneLetter == null)
						oneLetter = "X";
					sequence.append(oneLetter);
				}

			}
			return sequence.toString();
		}

		// not using ChemCOmp records...		
		List<Group> aminos = getAtomGroups(GroupType.AMINOACID);
		StringBuffer sequence = new StringBuffer() ;
		for ( int i=0 ; i< aminos.size(); i++){
			AminoAcid a = (AminoAcid)aminos.get(i);
			sequence.append( a.getAminoType());
		}

		return sequence.toString();

	}

	/**
	 * {@inheritDoc}	 
	 */
	public String getSeqResSequence(){

		String prop = System.getProperty(PDBFileReader.LOAD_CHEM_COMP_PROPERTY);

		if ( prop != null && prop.equalsIgnoreCase("true")){
			StringBuffer str = new StringBuffer();
			for (Group g : seqResGroups) {
				ChemComp cc = g.getChemComp();
				if ( cc == null) {
					logger.warn("Could not load ChemComp for group: ", g);
					str.append("X");
				} else if ( PolymerType.PROTEIN_ONLY.contains(cc.getPolymerType()) ||
						PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())){
					// an amino acid residue.. use for alignment
					String oneLetter= ChemCompGroupFactory.getOneLetterCode(cc);
					if ( oneLetter == null ||  oneLetter.length()==0  || oneLetter.equals("?"))
						oneLetter = "X";
					str.append(oneLetter);
				} else {
					str.append("X");
				}
			}
			return str.toString();
		}

		StringBuffer str = new StringBuffer();
		for (Group group : seqResGroups) {
			if (group instanceof AminoAcid) {
				AminoAcid aa = (AminoAcid)group;
				str.append(aa.getAminoType()) ;
			} else {
				str.append("X");
			}
		}
		return str.toString();

	}


	/** {@inheritDoc}
	 *
	 */
	public Group getSeqResGroup(int position) {

		return seqResGroups.get(position);
	}

	/** {@inheritDoc}
	 *
	 */
	public List<Group> getSeqResGroups(String type) {
		List<Group> tmp = new ArrayList<Group>() ;
		for (int i=0;i<seqResGroups.size();i++){
			Group g = (Group)seqResGroups.get(i);
			if (g.getType().equals(type)){
				tmp.add(g);
			}
		}

		return tmp ;
	}

	/** {@inheritDoc}
	 *
	 */
	public List<Group> getSeqResGroups() {
		return seqResGroups;
	}

	/** {@inheritDoc}
	 *
	 */
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
	public int getAtomLength() {

		return groups.size();
	}

	/** {@inheritDoc}
	 *
	 */
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
}

