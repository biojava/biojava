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
 * created at Apr 26, 2008
 */
package org.biojava.nbio.structure.io.mmcif;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.DBRef;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.HetatomImpl;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.SeqMisMatch;
import org.biojava.nbio.structure.SeqMisMatchImpl;
import org.biojava.nbio.structure.Site;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.BondMaker;
import org.biojava.nbio.structure.io.ChargeAdder;
import org.biojava.nbio.structure.io.EntityFinder;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.SeqRes2AtomAligner;
import org.biojava.nbio.structure.io.mmcif.model.AtomSite;
import org.biojava.nbio.structure.io.mmcif.model.AtomSites;
import org.biojava.nbio.structure.io.mmcif.model.AuditAuthor;
import org.biojava.nbio.structure.io.mmcif.model.Cell;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompAtom;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompBond;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompDescriptor;
import org.biojava.nbio.structure.io.mmcif.model.DatabasePDBremark;
import org.biojava.nbio.structure.io.mmcif.model.DatabasePDBrev;
import org.biojava.nbio.structure.io.mmcif.model.DatabasePdbrevRecord;
import org.biojava.nbio.structure.io.mmcif.model.Entity;
import org.biojava.nbio.structure.io.mmcif.model.EntityPoly;
import org.biojava.nbio.structure.io.mmcif.model.EntityPolySeq;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcGen;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcNat;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcSyn;
import org.biojava.nbio.structure.io.mmcif.model.Exptl;
import org.biojava.nbio.structure.io.mmcif.model.PdbxAuditRevisionHistory;
import org.biojava.nbio.structure.io.mmcif.model.PdbxChemCompDescriptor;
import org.biojava.nbio.structure.io.mmcif.model.PdbxChemCompIdentifier;
import org.biojava.nbio.structure.io.mmcif.model.PdbxDatabaseStatus;
import org.biojava.nbio.structure.io.mmcif.model.PdbxEntityNonPoly;
import org.biojava.nbio.structure.io.mmcif.model.PdbxNonPolyScheme;
import org.biojava.nbio.structure.io.mmcif.model.PdbxPolySeqScheme;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.nbio.structure.io.mmcif.model.Refine;
import org.biojava.nbio.structure.io.mmcif.model.Struct;
import org.biojava.nbio.structure.io.mmcif.model.StructAsym;
import org.biojava.nbio.structure.io.mmcif.model.StructConn;
import org.biojava.nbio.structure.io.mmcif.model.StructKeywords;
import org.biojava.nbio.structure.io.mmcif.model.StructNcsOper;
import org.biojava.nbio.structure.io.mmcif.model.StructRef;
import org.biojava.nbio.structure.io.mmcif.model.StructRefSeq;
import org.biojava.nbio.structure.io.mmcif.model.StructRefSeqDif;
import org.biojava.nbio.structure.io.mmcif.model.StructSite;
import org.biojava.nbio.structure.io.mmcif.model.StructSiteGen;
import org.biojava.nbio.structure.io.mmcif.model.Symmetry;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.biojava.nbio.structure.xtal.SymoplibParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A MMcifConsumer implementation that builds an in-memory representation of the
 * content of a mmcif file as a BioJava Structure object.
 *
 * @author Andreas Prlic
 * @since 1.7
 */

public class SimpleMMcifConsumer implements MMcifConsumer {

	private static final Logger logger = LoggerFactory.getLogger(SimpleMMcifConsumer.class);

	private Structure structure;
	private Chain currentChain;
	private Group currentGroup;

	/**
	 * A temporary data structure to hold all parsed chains
	 */
	private ArrayList<List<Chain>> allModels; 
	/**
	 * The current set of chains per model
	 */
	private List<Chain>      currentModel;
	private List<Entity>     entities;
	/**
	 * Needed in header only mode to get mapping between asym ids and author ids
	 */
	private List<EntityPoly> entityPolys;
	private List<StructRef>  strucRefs;
	private List<Chain>      seqResChains;
	private List<Chain>      entityChains; // needed to link entities, chains and compounds...
	private List<StructAsym> structAsyms;  // needed to link entities, chains and compounds...
	private List<PdbxStructOperList> structOpers ; //
	private List<PdbxStructAssembly> strucAssemblies;
	private List<PdbxStructAssemblyGen> strucAssemblyGens;
	private List<EntitySrcGen> entitySrcGens;
	private List<EntitySrcNat> entitySrcNats;
	private List<EntitySrcSyn> entitySrcSyns;
	private List<StructConn> structConn;
	private List<StructNcsOper> structNcsOper;
	private List<StructRefSeqDif> sequenceDifs;
	private List<StructSiteGen> structSiteGens;
	
	private Matrix4d parsedScaleMatrix;



	/**
	 * A map of asym ids (internal chain ids) to entity ids extracted from
	 * the _struct_asym category
	 */
	private Map<String,String> asymId2entityId;

	/**
	 * A map of asym ids (internal chain ids) to author ids extracted from 
	 * the _entity_poly category. Used in header only parsing.
	 */
	private Map<String,String> asymId2authorId;

	private String currentNmrModelNumber ;

	private FileParsingParameters params;

	public  SimpleMMcifConsumer(){
		params = new FileParsingParameters();
		documentStart();

	}

	@Override
	public void newEntity(Entity entity) {
		logger.debug("New entity: {}",entity.toString());
		entities.add(entity);
	}

	@Override
	public void newEntityPoly(EntityPoly entityPoly) {
		entityPolys.add(entityPoly);
	}

	@Override
	public void newPdbxStructOperList(PdbxStructOperList structOper){

		structOpers.add(structOper);
	}

	@Override
	public void newStructAsym(StructAsym sasym){

		structAsyms.add(sasym);
	}

	private Entity getEntity(int entity_id){
		try {
			for (Entity e: entities){
				int eId = Integer.parseInt(e.getId());
				if  (eId== entity_id){
					return e;
				}
			}
		} catch (NumberFormatException e) {
			logger.warn("Entity id does not look like a number:", e.getMessage());
		}
		return null;
	}

	@Override
	public void newStructKeywords(StructKeywords kw){
		PDBHeader header = structure.getPDBHeader();
		if ( header == null)
			header = new PDBHeader();
		header.setDescription(kw.getPdbx_keywords());
		header.setClassification(kw.getPdbx_keywords());
	}

	@Override
	public void setStruct(Struct struct) {

		PDBHeader header = structure.getPDBHeader();
		if ( header == null)
			header = new PDBHeader();

		header.setTitle(struct.getTitle());
		header.setIdCode(struct.getEntry_id());
		//header.setDescription(struct.getPdbx_descriptor());
		//header.setClassification(struct.getPdbx_descriptor());
		//header.setDescription(struct.getPdbx_descriptor());



		structure.setPDBHeader(header);
		structure.setPDBCode(struct.getEntry_id());
	}

	/** initiate new group, either Hetatom, Nucleotide, or AminoAcid */
	private Group getNewGroup(String recordName,Character aminoCode1, long seq_id,String groupCode3) {

		Group g = ChemCompGroupFactory.getGroupFromChemCompDictionary(groupCode3);
		if ( g != null && !g.getChemComp().isEmpty()) {
			if ( g instanceof AminoAcidImpl) {
				AminoAcidImpl aa = (AminoAcidImpl) g;
				aa.setId(seq_id);
			} else if ( g instanceof NucleotideImpl) {
				NucleotideImpl nuc =  (NucleotideImpl) g;
				nuc.setId(seq_id);
			} else if ( g instanceof HetatomImpl) {
				HetatomImpl het = (HetatomImpl)g;
				het.setId(seq_id);
			}
			return g;
		}



		Group group;
		if ( recordName.equals("ATOM") ) {
			if (StructureTools.isNucleotide(groupCode3))  {
				// it is a nucleotide
				NucleotideImpl nu = new NucleotideImpl();
				group = nu;
				nu.setId(seq_id);

			} else if (aminoCode1==null || aminoCode1 == StructureTools.UNKNOWN_GROUP_LABEL){
				HetatomImpl h = new HetatomImpl();
				h.setId(seq_id);
				group = h;

			} else {
				AminoAcidImpl aa = new AminoAcidImpl() ;
				aa.setAminoType(aminoCode1);
				aa.setId(seq_id);
				group = aa ;
			}
		}
		else {
			if (StructureTools.isNucleotide(groupCode3))  {
				// it is a nucleotide
				NucleotideImpl nu = new NucleotideImpl();
				group = nu;
				nu.setId(seq_id);
			}
			else if (aminoCode1 != null ) {
				AminoAcidImpl aa = new AminoAcidImpl() ;
				aa.setAminoType(aminoCode1);
				aa.setId(seq_id);
				group = aa ;
			} else {
				HetatomImpl h = new HetatomImpl();
				h.setId(seq_id);
				group = h;
			}
		}
		return  group ;
	}

	/**
	 * Test if the given asymId is already present in the list of chains given. If yes, returns the chain
	 * otherwise returns null.
	 */
	private static Chain isKnownChain(String asymId, List<Chain> chains){

		for (int i = 0; i< chains.size();i++){
			Chain testchain =  chains.get(i);
			//System.out.println("comparing chainID >"+chainID+"< against testchain " + i+" >" +testchain.getName()+"<");
			if (asymId.equals(testchain.getId())) {
				//System.out.println("chain "+ chainID+" already known ...");
				return testchain;
			}
		}

		return null;
	}

	@Override
	public void newAtomSite(AtomSite atom) {

		if (params.isHeaderOnly()) return;

		// Warning: getLabel_asym_id is not the "chain id" in the PDB file
		// it is the internally used chain id.
		// later on we will fix this...

		// later one needs to map the asym id to the pdb_strand_id

		//TODO: add support for FileParsingParams.getMaxAtoms()

		boolean startOfNewChain = false;

		String asymId = atom.getLabel_asym_id();
		String authId = atom.getAuth_asym_id();

		String recordName    = atom.getGroup_PDB();
		String residueNumberS = atom.getAuth_seq_id();
		Integer residueNrInt = Integer.parseInt(residueNumberS);

		// the 3-letter name of the group:
		String groupCode3    = atom.getLabel_comp_id();

		boolean isHetAtomInFile = false;

		Character aminoCode1 = null;
		if ( recordName.equals("ATOM") )
			aminoCode1 = StructureTools.get1LetterCodeAmino(groupCode3);
		else {
			aminoCode1 = StructureTools.get1LetterCodeAmino(groupCode3);

			// for nucleotides this will be null..
			if (aminoCode1 != null &&  aminoCode1.equals(StructureTools.UNKNOWN_GROUP_LABEL))
				aminoCode1 = null;

			isHetAtomInFile = true;
		}
		String insCodeS = atom.getPdbx_PDB_ins_code();
		Character insCode = null;
		if (!  insCodeS.equals("?")) {
			insCode = insCodeS.charAt(0);
		}
		// we store the internal seq id in the Atom._id field
		// this is not a PDB file field but we need this to internally assign the insertion codes later
		// from the pdbx_poly_seq entries..

		long seq_id = -1;
		try {
			seq_id = Long.parseLong(atom.getLabel_seq_id());
		} catch (NumberFormatException e){
			// non polymer chains (ligands and small molecules) will have a label_seq_id set to '.', thus it is ok to
			// silently ignore this
			//logger.debug("Could not parse number for _atom_site.label_seq_id: "+e.getMessage());
		}

		String nmrModelNumber = atom.getPdbx_PDB_model_num();

		if ( currentNmrModelNumber == null) {
			currentNmrModelNumber = nmrModelNumber;
		}

		if (! currentNmrModelNumber.equals(nmrModelNumber)){
			currentNmrModelNumber = nmrModelNumber;

			// add previous data
			if ( currentChain != null ) {
				currentChain.addGroup(currentGroup);
				currentGroup.trimToSize();
			}

			// we came to the beginning of a new NMR model
			allModels.add(currentModel);
			currentModel = new ArrayList<Chain>();
			currentChain = null;
			currentGroup = null;
		}


		if (currentChain == null) {

			currentChain = new ChainImpl();
			currentChain.setName(authId);
			currentChain.setId(asymId);
			currentModel.add(currentChain);
			startOfNewChain = true;
		}

		//System.out.println("BEFORE: " + chain_id + " " + current_chain.getName());
		if ( ! asymId.equals(currentChain.getId()) ) {
			//logger.info("unknown chain. creating new chain. authId:" + authId + " asymId: " + asymId);
			startOfNewChain = true;

			// end up old chain...
			currentChain.addGroup(currentGroup);

			// see if old chain is known ...
			Chain testchain = isKnownChain(asymId,currentModel);

			if ( testchain == null) {
				//logger.info("unknown chain. creating new chain. authId:" + authId + " asymId: " + asymId);

				currentChain = new ChainImpl();
				currentChain.setName(authId);
				currentChain.setId(asymId);

			}   else {
				currentChain = testchain;
			}

			if ( ! currentModel.contains(currentChain))
				currentModel.add(currentChain);

		}


		ResidueNumber residueNumber = new ResidueNumber(authId,residueNrInt, insCode);

		if (currentGroup == null) {


			currentGroup = getNewGroup(recordName,aminoCode1,seq_id, groupCode3);

			currentGroup.setResidueNumber(residueNumber);
			currentGroup.setPDBName(groupCode3);
			currentGroup.setHetAtomInFile(isHetAtomInFile);
		}

		// SET UP THE ALT LOC GROUP
		Group altGroup = null;
		String altLocS = atom.getLabel_alt_id();
		Character altLoc = ' ';
		if ( altLocS.length()>0) {
			altLoc = altLocS.charAt(0);
			if ( altLoc.equals('.') )
				altLoc = ' ';

		}
		// If it's the start of the new chain 
		if ( startOfNewChain){
			currentGroup = getNewGroup(recordName,aminoCode1,seq_id, groupCode3);
			currentGroup.setResidueNumber(residueNumber);
			currentGroup.setPDBName(groupCode3);
			currentGroup.setHetAtomInFile(isHetAtomInFile);
		}
		// ANTHONY BRADLEY ADDED THIS -> WE ONLY WAN'T TO CHECK FOR ALT LOCS WHEN IT's NOT THE FIRST GROUP IN CHAIN
		else{
			// check if residue number is the same ...
			// insertion code is part of residue number
			if ( ! residueNumber.equals(currentGroup.getResidueNumber())) {
				//System.out.println("end of residue: "+current_group.getPDBCode()+" "+residueNrInt);
				currentChain.addGroup(currentGroup);
				currentGroup.trimToSize();
				currentGroup = getNewGroup(recordName,aminoCode1,seq_id,groupCode3);
				currentGroup.setPDBName(groupCode3);
				currentGroup.setResidueNumber(residueNumber);
				currentGroup.setHetAtomInFile(isHetAtomInFile);


			} else {
				// same residueNumber, but altLocs...
				// test altLoc

				if ( ! altLoc.equals(' ') && ( ! altLoc.equals('.'))) {
					logger.debug("found altLoc! " + altLoc + " " + currentGroup + " " + altGroup);
					altGroup = getCorrectAltLocGroup( altLoc,recordName,aminoCode1,groupCode3, seq_id);
					if (altGroup.getChain()==null) {
						altGroup.setChain(currentChain);
					}
				}
			}
		}
		//atomCount++;
		//System.out.println("fixing atom name for  >" + atom.getLabel_atom_id() + "< >" + fullname + "<");


		if ( params.isParseCAOnly() ){
			// yes , user wants to get CA only
			// only parse CA atoms...
			if (! (atom.getLabel_atom_id().equals(StructureTools.CA_ATOM_NAME) && atom.getType_symbol().equals("C"))) {
				//System.out.println("ignoring " + line);
				//atomCount--;
				return;
			}
		}

		//see if chain_id is one of the previous chains ...

		Atom a = convertAtom(atom);

		//see if chain_id is one of the previous chains ...
		if ( altGroup != null) {
			altGroup.addAtom(a);
			altGroup = null;
		}
		else {
			currentGroup.addAtom(a);
		}


		String atomName = a.getName();
		// make sure that main group has all atoms 
		// GitHub issue: #76
		if ( ! currentGroup.hasAtom(atomName)) {
			// Unless it's microheterogenity https://github.com/rcsb/codec-devel/issues/81
			if (currentGroup.getPDBName().equals(a.getGroup().getPDBName())) {
				if(!StructureTools.hasNonDeuteratedEquiv(a,currentGroup)){
					currentGroup.addAtom(a);
				}
			}

		}
	}

	/** 
	 * Convert a mmCIF AtomSite object to a BioJava Atom object
	 *
	 * @param atom the mmmcif AtomSite record
	 * @return an Atom
	 */
	private Atom convertAtom(AtomSite atom){


		Atom a = new AtomImpl();

		a.setPDBserial(Integer.parseInt(atom.getId()));
		a.setName(atom.getLabel_atom_id());

		double x = Double.parseDouble (atom.getCartn_x());
		double y = Double.parseDouble (atom.getCartn_y());
		double z = Double.parseDouble (atom.getCartn_z());
		a.setX(x);
		a.setY(y);
		a.setZ(z);

		float occupancy = Float.parseFloat (atom.getOccupancy());
		a.setOccupancy(occupancy);

		float temp = Float.parseFloat (atom.getB_iso_or_equiv());
		a.setTempFactor(temp);

		String alt = atom.getLabel_alt_id();
		if (( alt != null ) && ( alt.length() > 0) && (! alt.equals("."))){
			a.setAltLoc(new Character(alt.charAt(0)));
		} else {
			a.setAltLoc(new Character(' '));
		}

		Element element = Element.R;
		try {
			element = Element.valueOfIgnoreCase(atom.getType_symbol());
		}  catch (IllegalArgumentException e) {
			logger.info("Element {} was not recognised as a BioJava-known element, the element will be represented as the generic element {}", atom.getType_symbol(), Element.R.name());
		}
		a.setElement(element);

		return a;

	}


	private Group getCorrectAltLocGroup( Character altLoc,
			String recordName,
			Character aminoCode1,
			String groupCode3,
			long seq_id) {

		// see if we know this altLoc already;
		List<Atom> atoms = currentGroup.getAtoms();
		if ( atoms.size() > 0) {
			Atom a1 = atoms.get(0);
			// we are just adding atoms to the current group
			// probably there is a second group following later...
			if (a1.getAltLoc().equals(altLoc)) {

				return currentGroup;
			}
		}

		List<Group> altLocs = currentGroup.getAltLocs();
		for ( Group altLocG : altLocs ){
			atoms = altLocG.getAtoms();
			if ( atoms.size() > 0) {
				for ( Atom a1 : atoms) {
					if (a1.getAltLoc().equals( altLoc)) {

						return altLocG;
					}
				}
			}
		}

		// no matching altLoc group found.
		// build it up.

		if ( groupCode3.equals(currentGroup.getPDBName())) {
			if ( currentGroup.getAtoms().size() == 0) {
				//System.out.println("current group is empty " + current_group + " " + altLoc);
				return currentGroup;
			}
			//System.out.println("cloning current group " + current_group + " " + current_group.getAtoms().get(0).getAltLoc() + " altLoc " + altLoc);
			Group altLocG = (Group) currentGroup.clone();
			// drop atoms from cloned group...
			// https://redmine.open-bio.org/issues/3307
			altLocG.setAtoms(new ArrayList<Atom>());
			altLocG.getAltLocs().clear();
			currentGroup.addAltLoc(altLocG);
			return altLocG;
		}

		//	System.out.println("new  group " + recordName + " " + aminoCode1 + " " +groupCode3);
		//String recordName,Character aminoCode1, long seq_id,String groupCode3) {
		Group altLocG = getNewGroup(recordName,aminoCode1,seq_id,groupCode3);

		altLocG.setPDBName(groupCode3);
		altLocG.setResidueNumber(currentGroup.getResidueNumber());
		currentGroup.addAltLoc(altLocG);
		return altLocG;
	}

	/** 
	 * Start the parsing
	 */
	@Override
	public void documentStart() {
		structure = new StructureImpl();

		currentChain        = null;
		currentGroup 		= null;
		currentNmrModelNumber 	= null;
		//atomCount     		= 0;

		allModels     = new ArrayList<List<Chain>>();
		currentModel  = new ArrayList<Chain>();
		entities      = new ArrayList<Entity>();
		entityPolys   = new ArrayList<>();
		strucRefs     = new ArrayList<StructRef>();
		seqResChains  = new ArrayList<Chain>();
		entityChains  = new ArrayList<Chain>();
		structAsyms   = new ArrayList<StructAsym>();

		asymId2entityId = new HashMap<String,String>();
		asymId2authorId = new HashMap<>();
		structOpers   = new ArrayList<PdbxStructOperList>();
		strucAssemblies = new ArrayList<PdbxStructAssembly>();
		strucAssemblyGens = new ArrayList<PdbxStructAssemblyGen>();
		entitySrcGens = new ArrayList<EntitySrcGen>();
		entitySrcNats = new ArrayList<EntitySrcNat>();
		entitySrcSyns = new ArrayList<EntitySrcSyn>();
		structConn = new ArrayList<StructConn>();
		structNcsOper = new ArrayList<StructNcsOper>();
		sequenceDifs = new ArrayList<StructRefSeqDif>();
		structSiteGens = new ArrayList<StructSiteGen>();
	}


	@Override
	public void documentEnd() {

		// Expected that there is one current_chain that needs to be added to the model
		// When in headerOnly mode, no Atoms are read, and there will not be an active
		// current_chain.
		if ( currentChain != null ) {

			currentChain.addGroup(currentGroup);
			if (isKnownChain(currentChain.getId(),currentModel) == null) {
				currentModel.add(currentChain);
			}
		} else if (!params.isHeaderOnly()){
			logger.warn("current chain is null at end of document.");
		}

		allModels.add(currentModel);

		// this populates the asymId2authorId and asymId2entityId maps, needed in header only mode to get the mapping 
		// between the 2 chain identifiers.
		initMaps();

		for (StructAsym asym : structAsyms) {

			logger.debug("Entity {} matches asym_id: {}", asym.getEntity_id(), asym.getId() );

			Chain s = getEntityChain(asym.getEntity_id());
			Chain seqres = (Chain)s.clone();
			// to solve issue #160 (e.g. 3u7t)
			seqres = removeSeqResHeterogeneity(seqres);
			seqres.setId(asym.getId());
			if (asymId2authorId.get(asym.getId()) !=null ){ 
				seqres.setName(asymId2authorId.get(asym.getId()));
			} else {
				seqres.setName(asym.getId());
			}

			EntityType type = null;
			try {
				Entity ent = getEntity(Integer.parseInt(asym.getEntity_id()));
				type = EntityType.entityTypeFromString(ent.getType());
			} catch (NumberFormatException e) {
				logger.debug("Could not parse integer from entity id field {}", asym.getEntity_id());
			}

			// we'll only add seqres chains that are polymeric or unknown
			if (type==null || type==EntityType.POLYMER ) {
				seqResChains.add(seqres);	
			}

			logger.debug(" seqres: " + asym.getId() + " " + seqres + "<") ;
			// adding the entities to structure
			addEntities(asym);

		}

		if (structAsyms.isEmpty()) {
			logger.warn("No _struct_asym category in file, no SEQRES groups will be added.");
		}

		// entities
		// In addEntities above we created the entities if they were present in the file
		// Now we need to make sure that they are linked to chains and also that if they are not present in the file we need to add them now
		linkEntities();

		// now that we know the entities, we can add all chains to structure so that they are stored
		// properly as polymer/nonpolymer/water chains inside structure
		for (List<Chain> model:allModels) {
			structure.addModel(model);
		}

		// Only align if requested (default) and not when headerOnly mode with no Atoms.
		// Otherwise, we store the empty SeqRes Groups unchanged in the right chains.
		if ( params.isAlignSeqRes() && !params.isHeaderOnly() ){
			logger.debug("Parsing mode align_seqres, will parse SEQRES and align to ATOM sequence");
			alignSeqRes();
		} else {
			logger.debug("Parsing mode unalign_seqres, will parse SEQRES but not align it to ATOM sequence");
			SeqRes2AtomAligner.storeUnAlignedSeqRes(structure, seqResChains, params.isHeaderOnly());
		}


		// Now make sure all altlocgroups have all the atoms in all the groups
		StructureTools.cleanUpAltLocs(structure);

		// NOTE bonds and charges can only be done at this point that the chain id mapping is properly sorted out
		if (!params.isHeaderOnly()) {
			if ( params.shouldCreateAtomBonds()) {
				addBonds();
			}

			if ( params.shouldCreateAtomCharges()) {
				addCharges();
			}
		}

		if (!params.isHeaderOnly()) {

			// Do structure.setSites(sites) after any chain renaming to be like PDB.
			addSites();
		}



		// set the oligomeric state info in the header...
		if (params.isParseBioAssembly()) {

			// the more detailed mapping of chains to rotation operations happens in StructureIO...

			Map<Integer,BioAssemblyInfo> bioAssemblies = new LinkedHashMap<Integer, BioAssemblyInfo>();

			for ( PdbxStructAssembly psa : strucAssemblies){

				List<PdbxStructAssemblyGen> psags = new ArrayList<PdbxStructAssemblyGen>(1);

				for ( PdbxStructAssemblyGen psag: strucAssemblyGens ) {
					if ( psag.getAssembly_id().equals(psa.getId())) {
						psags.add(psag);
					}
				}

				BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

				// these are the transformations that need to be applied to our model
				List<BiologicalAssemblyTransformation> transformations = builder.getBioUnitTransformationList(psa, psags, structOpers);

				int bioAssemblyId = -1;
				try {
					bioAssemblyId = Integer.parseInt(psa.getId());
				} catch (NumberFormatException e) {
					logger.info("Could not parse a numerical bio assembly id from '{}'",psa.getId());
				}

				// if bioassembly id is not numerical we throw it away
				// this happens usually for viral capsid entries, like 1ei7
				// see issue #230 in github
				if (bioAssemblyId!=-1) {
					int mmSize = 0;
					// note that the transforms contain asym ids of both polymers and non-polymers
					// For the mmsize, we are only interested in the polymers
					for (BiologicalAssemblyTransformation transf:transformations) {
						Chain c = structure.getChain(transf.getChainId());
						if (c==null) {
							logger.info("Could not find asym id {} specified in struct_assembly_gen", transf.getChainId());
							continue;
						}
						if (c.getEntityType() == EntityType.POLYMER &&
							// for entries like 4kro, sugars are annotated as polymers but we
							// don't want them in the macromolecularSize count
							!c.getEntityInfo().getDescription().contains("SUGAR") ) {
								
								mmSize++;
							}
					}
					
					BioAssemblyInfo bioAssembly = new BioAssemblyInfo();
					bioAssembly.setId(bioAssemblyId);
					bioAssembly.setMacromolecularSize(mmSize);
					bioAssembly.setTransforms(transformations);
					bioAssemblies.put(bioAssemblyId,bioAssembly);
				}

			}
			structure.getPDBHeader().setBioAssemblies(bioAssemblies);
		}

		setStructNcsOps();
		
		setCrystallographicInfoMetadata();


		Map<String,List<SeqMisMatch>> misMatchMap = new HashMap<String, List<SeqMisMatch>>();
		for (StructRefSeqDif sdif : sequenceDifs) {
			SeqMisMatch misMatch = new SeqMisMatchImpl();
			misMatch.setDetails(sdif.getDetails());

			String insCode = sdif.getPdbx_pdb_ins_code();
			if ( insCode != null && insCode.equals("?"))
				insCode = null;
			misMatch.setInsCode(insCode);
			misMatch.setOrigGroup(sdif.getDb_mon_id());
			misMatch.setPdbGroup(sdif.getMon_id());
			misMatch.setPdbResNum(sdif.getPdbx_auth_seq_num());
			misMatch.setUniProtId(sdif.getPdbx_seq_db_accession_code());
			misMatch.setSeqNum(sdif.getSeq_num());


			List<SeqMisMatch> mms = misMatchMap.get(sdif.getPdbx_pdb_strand_id());
			if ( mms == null) {
				mms = new ArrayList<SeqMisMatch>();
				misMatchMap.put(sdif.getPdbx_pdb_strand_id(),mms);
			}
			mms.add(misMatch);

		}

		for (String chainId : misMatchMap.keySet()){

			Chain chain = structure.getPolyChainByPDB(chainId);

			if ( chain == null) {
				logger.warn("Could not set mismatches for chain with author id" + chainId);
				continue;
			}

			chain.setSeqMisMatches(misMatchMap.get(chainId));


		}

	}

	/**
	 * Here we link entities to chains.
	 * Also if entities are not present in file, this initialises the entities with some heuristics, see {@link org.biojava.nbio.structure.io.EntityFinder}
	 */
	private void linkEntities() {

		for (int i =0; i< allModels.size() ; i++){
			for (Chain chain : allModels.get(i)) {
				//logger.info("linking entities for " + chain.getId() + " "  + chain.getName());
				String entityId = asymId2entityId.get(chain.getId());

				if (entityId==null) {
					// this can happen for instance if the cif file didn't have _struct_asym category at all
					// and thus we have no asymId2entityId mapping at all
					logger.info("No entity id could be found for chain {}", chain.getId());
					continue;
				}
				int eId = Integer.parseInt(entityId);

				// Entities are not added for non-polymeric entities, if a chain is non-polymeric its entity won't be found.
				// TODO: add all entities and unique compounds and add methods to directly get polymer or non-polymer
				// asyms (chains).  Either create a unique StructureImpl or modify existing for a better representation of the
				// mmCIF internal data structures but is compatible with Structure interface.
				// Some examples of PDB entries with this kind of problem:
				//   - 2uub: asym_id X, chainName Z, entity_id 24: fully non-polymeric but still with its own chainName
				//   - 3o6j: asym_id K, chainName Z, entity_id 6 : a single water molecule
				//   - 1dz9: asym_id K, chainName K, entity_id 6 : a potassium ion alone

				EntityInfo entityInfo = structure.getEntityById(eId);
				if (entityInfo==null) {
					// Supports the case where the only chain members were from non-polymeric entity that is missing.
					// Solved by creating a new Compound(entity) to which this chain will belong.
					logger.info("Could not find an Entity for entity_id {}, for chain id {}, creating a new Entity.",
							eId, chain.getId());
					entityInfo = new EntityInfo();
					entityInfo.setMolId(eId);
					entityInfo.addChain(chain);
					if (chain.isWaterOnly()) {
						entityInfo.setType(EntityType.WATER);
					} else {
						entityInfo.setType(EntityType.NONPOLYMER);
					}
					chain.setEntityInfo(entityInfo);
					structure.addEntityInfo(entityInfo);
				} else {
					logger.debug("Adding chain with chain id {} (auth id {}) to Entity with entity_id {}",
							chain.getId(), chain.getName(), eId);
					entityInfo.addChain(chain);
					chain.setEntityInfo(entityInfo);
				}

			}

		}

		// if no entity information was present in file we then go and find the entities heuristically with EntityFinder
		List<EntityInfo> entityInfos = structure.getEntityInfos();
		if (entityInfos==null || entityInfos.isEmpty()) {

			List<List<Chain>> polyModels = new ArrayList<>();
			List<List<Chain>> nonPolyModels = new ArrayList<>();
			List<List<Chain>> waterModels = new ArrayList<>();

			for (List<Chain> model:allModels) {

				List<Chain> polyChains = new ArrayList<>();
				List<Chain> nonPolyChains = new ArrayList<>();
				List<Chain> waterChains = new ArrayList<>();

				polyModels.add(polyChains);
				nonPolyModels.add(nonPolyChains);
				waterModels.add(waterChains);

				for (Chain c:model) {

					// we only have entities for polymeric chains, all others are ignored for assigning entities
					if (c.isWaterOnly()) {
						waterChains.add(c);

					} else if (c.isPureNonPolymer()) {
						nonPolyChains.add(c);

					} else {
						polyChains.add(c);
					}
				}
			}

			entityInfos = EntityFinder.findPolyEntities(polyModels);
			EntityFinder.createPurelyNonPolyEntities(nonPolyModels, waterModels, entityInfos);


			structure.setEntityInfos(entityInfos);
		}

		// final sanity check: it can happen that from the annotated entities some are not linked to any chains
		// e.g. 3s26: a sugar entity does not have any chains associated to it (it seems to be happening with many sugar compounds)
		// we simply log it, this can sign some other problems if the entities are used down the line
		for (EntityInfo e:entityInfos) {
			if (e.getChains().isEmpty()) {
				logger.info("Entity {} '{}' has no chains associated to it",
						e.getMolId()<0?"with no entity id":e.getMolId(), e.getDescription());
			}
		}

	}

	private void addCharges() {
		ChargeAdder.addCharges(structure);
	}

	/**
	 * The method will return a new reference to a Chain with any consecutive groups
	 * having same residue numbers removed.
	 * This is necessary to solve the microheterogeneity issue in entries like 3u7t (see github issue #160)
	 * @param c
	 * @return
	 */
	private static Chain removeSeqResHeterogeneity(Chain c) {

		Chain trimmedChain = new ChainImpl();

		ResidueNumber lastResNum = null;

		for (Group g:c.getAtomGroups()) {

			// note we have to deep copy this, otherwise they stay linked and would get altered in addGroup(g)
			ResidueNumber currentResNum = new ResidueNumber(
					g.getResidueNumber().getChainName(),
					g.getResidueNumber().getSeqNum(),
					g.getResidueNumber().getInsCode());

			if (lastResNum == null || !lastResNum.equals(currentResNum) ) {
				trimmedChain.addGroup(g);
			} else {
				logger.debug("Removing seqres group because it seems to be repeated in entity_poly_seq, most likely has hetero='y': "+g);
			}

			lastResNum = currentResNum;

		}
		return trimmedChain;
	}

	private void addBonds() {
		BondMaker maker = new BondMaker(structure, params);
		maker.makeBonds();
		maker.formBondsFromStructConn(structConn);
	}

	private void alignSeqRes() {

		logger.debug("Parsing mode align_seqres, will align to ATOM to SEQRES sequence");

		// fix SEQRES residue numbering for all models

		for (int model=0;model<structure.nrModels();model++) {

			List<Chain> atomList   = structure.getModel(model);

			for (Chain seqResChain: seqResChains){

				// this extracts the matching atom chain from atomList
				Chain atomChain = SeqRes2AtomAligner.getMatchingAtomRes(seqResChain, atomList, true);

				if (atomChain == null) {
					// most likely there's no observed residues at all for the seqres chain: can't map
					// e.g. 3zyb: chains with asym_id L,M,N,O,P have no observed residues
					logger.info("Could not map SEQRES chain with asym_id={} to any ATOM chain. Most likely there's no observed residues in the chain.",
							seqResChain.getId());
					continue;
				}

				//map the atoms to the seqres...

				// we need to first clone the seqres so that they stay independent for different models
				List<Group> seqResGroups = new ArrayList<Group>();
				for (int i=0;i<seqResChain.getAtomGroups().size();i++) {
					seqResGroups.add((Group)seqResChain.getAtomGroups().get(i).clone());
				}

				for ( int seqResPos = 0 ; seqResPos < seqResGroups.size(); seqResPos++) {
					Group seqresG = seqResGroups.get(seqResPos);
					boolean found = false;
					for ( Group atomG: atomChain.getAtomGroups()) {

						int internalNr = getInternalNr (atomG);

						if (seqresG.getResidueNumber().getSeqNum() == internalNr ) {
							seqResGroups.set(seqResPos, atomG);
							found = true;
							break;
						}


					}
					if ( ! found)
						// so far the residue number has tracked internal numbering.
						// however there are no atom records, as such this can't be a PDB residue number...
						seqresG.setResidueNumber(null);
				}
				atomChain.setSeqResGroups(seqResGroups);

			}
		}
	}

	private int getInternalNr(Group atomG) {
		if ( atomG.getType().equals(GroupType.AMINOACID)) {
			AminoAcidImpl aa = (AminoAcidImpl) atomG;
			return new Long(aa.getId()).intValue();
		} else if ( atomG.getType().equals(GroupType.NUCLEOTIDE)) {
			NucleotideImpl nu = (NucleotideImpl) atomG;
			return new Long(nu.getId()).intValue();
		} else {
			HetatomImpl he = (HetatomImpl) atomG;
			return new Long(he.getId()).intValue();
		}
	}

	private void addEntities(StructAsym asym) {
		int eId = 0;
		try {
			eId = Integer.parseInt(asym.getEntity_id());
		} catch (NumberFormatException e) {
			logger.warn("Could not parse mol_id from string {}. Will use 0 for creating Entity",asym.getEntity_id());
		}
		Entity e = getEntity(eId);

		// for some mmCIF files like 1yrm all 3 of _entity_src_gen, _entity_src_nat and _pdbx_entity_src_syn are missing
		// we need to fill the Compounds in some other way:

		EntityInfo entityInfo = structure.getEntityById(eId);

		if (entityInfo==null) {
			//logger.info("Creating new EntityInfo " + eId + " " + e.getId() + " " + e.getPdbx_description());
			entityInfo = new EntityInfo();
			entityInfo.setMolId(eId);
			// we only add the compound if a polymeric one (to match what the PDB parser does)
			if (e!=null) {
				entityInfo.setDescription(e.getPdbx_description());

				EntityType eType = EntityType.entityTypeFromString(e.getType());
				if (eType!=null) {
					entityInfo.setType(eType);
				} else {
					logger.warn("Type '{}' is not recognised as a valid entity type for entity {}", e.getType(), eId);
				}
				addAncilliaryEntityData(asym, eId, e, entityInfo);
				structure.addEntityInfo(entityInfo);
				logger.debug("Adding Entity with entity id {} from _entity, with name: {}",eId, entityInfo.getDescription());
			}
		}
	}


	/**
	 * Add any extra information to the entity information.
	 * @param asym 
	 * @param entityId 
	 * @param entity 
	 * @param entityInfo 
	 */
	private void addAncilliaryEntityData(StructAsym asym, int entityId, Entity entity, EntityInfo entityInfo) {
		// Loop through each of the entity types and add the corresponding data
		// We're assuming if data is duplicated between sources it is consistent
		// This is a potentially huge assumption...


		for (EntitySrcGen esg : entitySrcGens) {

			if (! esg.getEntity_id().equals(asym.getEntity_id()))
				continue;

			addInformationFromESG(esg, entityId, entityInfo);

		}

		for (EntitySrcNat esn : entitySrcNats) {
			if (! esn.getEntity_id().equals(asym.getEntity_id()))
				continue;
			addInformationFromESN(esn, entityId, entityInfo);

		}

		for (EntitySrcSyn ess : entitySrcSyns) {
			if (! ess.getEntity_id().equals(asym.getEntity_id()))
				continue;
			addInfoFromESS(ess, entityId, entityInfo);

		}		
	}

	/**
	 * Add the information from an ESG to a compound.
	 * @param entitySrcInfo
	 * @param entityId
	 * @param c
	 */
	private void addInformationFromESG(EntitySrcGen entitySrcInfo, int entityId, EntityInfo c) {
		c.setAtcc(entitySrcInfo.getPdbx_gene_src_atcc());
		c.setCell(entitySrcInfo.getPdbx_gene_src_cell());
		c.setOrganismCommon(entitySrcInfo.getGene_src_common_name());
		c.setOrganismScientific(entitySrcInfo.getPdbx_gene_src_scientific_name());
		c.setOrganismTaxId(entitySrcInfo.getPdbx_gene_src_ncbi_taxonomy_id());
		c.setExpressionSystemTaxId(entitySrcInfo.getPdbx_host_org_ncbi_taxonomy_id());
		c.setExpressionSystem(entitySrcInfo.getPdbx_host_org_scientific_name());
	}

	/**
	 * Add the information to entity info from ESN.
	 * @param esn
	 * @param eId
	 * @param c
	 */
	private void addInformationFromESN(EntitySrcNat esn, int eId, EntityInfo c) {

		c.setAtcc(esn.getPdbx_atcc());
		c.setCell(esn.getPdbx_cell());
		c.setOrganismCommon(esn.getCommon_name());
		c.setOrganismScientific(esn.getPdbx_organism_scientific());
		c.setOrganismTaxId(esn.getPdbx_ncbi_taxonomy_id());

	}
	/**
	 * Add the information from ESS to Entity info.
	 * @param ess
	 * @param eId
	 * @param c
	 */
	private void addInfoFromESS(EntitySrcSyn ess, int eId, EntityInfo c) {
		c.setOrganismCommon(ess.getOrganism_common_name());
		c.setOrganismScientific(ess.getOrganism_scientific());
		c.setOrganismTaxId(ess.getNcbi_taxonomy_id());

	}

	private void initMaps() {


		if (structAsyms == null || structAsyms.isEmpty()) {
			logger.info("No _struct_asym category found in file. No asym id to entity_id mapping will be available");
			return;
		}

		Map<String, List<String>> entityId2asymId = new HashMap<>();

		for (StructAsym asym : structAsyms) {

			logger.debug("Entity {} matches asym_id: {}", asym.getEntity_id(), asym.getId() );

			asymId2entityId.put(asym.getId(), asym.getEntity_id());

			if (entityId2asymId.containsKey(asym.getEntity_id())) {
				List<String> asymIds = entityId2asymId.get(asym.getEntity_id());
				asymIds.add(asym.getId());
			} else {
				List<String> asymIds = new ArrayList<>();
				asymIds.add(asym.getId());
				entityId2asymId.put(asym.getEntity_id(), asymIds);
			}
		}

		if (entityPolys==null || entityPolys.isEmpty()) {
			logger.info("No _entity_poly category found in file. No asym id to author id mapping will be available for header only parsing");
			return;
		}

		for (EntityPoly ep:entityPolys) {
			if (ep.getPdbx_strand_id()==null) {
				logger.info("_entity_poly.pdbx_strand_id is null for entity {}. Won't be able to map asym ids to author ids for this entity.", ep.getEntity_id());
				continue;
			}
			String[] chainNames = ep.getPdbx_strand_id().split(",");
			List<String> asymIds = entityId2asymId.get(ep.getEntity_id());
			if (chainNames.length!=asymIds.size()) {
				logger.warn("The list of asym ids (from _struct_asym) and the list of author ids (from _entity_poly) for entity {} have different lengths! Can't provide a mapping from asym ids to author chain ids", ep.getEntity_id());
				continue;
			}
			for (int i=0; i<chainNames.length; i++) {
				asymId2authorId.put(asymIds.get(i), chainNames[i]);
			}
		}
	}
	
	private void setStructNcsOps() {
		
		ArrayList<Matrix4d> ncsOperators = new ArrayList<Matrix4d>();
		
		for (StructNcsOper sNcsOper:structNcsOper) {
			
			if (!sNcsOper.getCode().equals("generate")) continue;
			
			try {
				Matrix4d op = new Matrix4d();
				op.setElement(3, 0, 0.0);
				op.setElement(3, 1, 0.0);
				op.setElement(3, 2, 0.0);
				op.setElement(3, 3, 1.0);


				op.setElement(0, 0, Double.parseDouble(sNcsOper.getMatrix11()));
				op.setElement(0, 1, Double.parseDouble(sNcsOper.getMatrix12()));
				op.setElement(0, 2, Double.parseDouble(sNcsOper.getMatrix13()));

				op.setElement(1, 0, Double.parseDouble(sNcsOper.getMatrix21()));
				op.setElement(1, 1, Double.parseDouble(sNcsOper.getMatrix22()));
				op.setElement(1, 2, Double.parseDouble(sNcsOper.getMatrix23()));

				op.setElement(2, 0, Double.parseDouble(sNcsOper.getMatrix31()));
				op.setElement(2, 1, Double.parseDouble(sNcsOper.getMatrix32()));
				op.setElement(2, 2, Double.parseDouble(sNcsOper.getMatrix33()));

				op.setElement(0, 3, Double.parseDouble(sNcsOper.getVector1()));
				op.setElement(1, 3, Double.parseDouble(sNcsOper.getVector2()));
				op.setElement(2, 3, Double.parseDouble(sNcsOper.getVector3()));

				ncsOperators.add(op);
				
			} catch (NumberFormatException e) {
				logger.warn("Error parsing doubles in NCS operator list, skipping operator {}", structNcsOper.indexOf(sNcsOper)+1); 
			}

		}
		
		// we only set it if not empty, otherwise remains null
		if (ncsOperators.size()>0) {
			structure.getCrystallographicInfo().setNcsOperators(
					ncsOperators.toArray(new Matrix4d[ncsOperators.size()]));
		}
	}
	
	private void setCrystallographicInfoMetadata() {
		if (parsedScaleMatrix!=null) {
			
			PDBCrystallographicInfo crystalInfo = structure.getCrystallographicInfo();
			
			boolean nonStd = false;
			if (crystalInfo.getCrystalCell()!=null && !crystalInfo.getCrystalCell().checkScaleMatrix(parsedScaleMatrix)) {
				nonStd = true;
			}
			
			crystalInfo.setNonStandardCoordFrameConvention(nonStd); 
		}
	}


	/** This method will return the parsed protein structure, once the parsing has been finished
	 *
	 * @return a BioJava protein structure object
	 */
	public Structure getStructure() {

		return structure;
	}

	@Override
	public void newDatabasePDBrevRecord(DatabasePdbrevRecord record) {

		PDBHeader header = structure.getPDBHeader();

		if ( header == null) {
			header = new PDBHeader();
			structure.setPDBHeader(header);
		}

		List<DatabasePdbrevRecord> revRecords = header.getRevisionRecords();
		if ( revRecords == null) {
			revRecords = new ArrayList<DatabasePdbrevRecord>();
			header.setRevisionRecords(revRecords);
		}
		revRecords.add(record);


	}


	@Override
	public void newDatabasePDBrev(DatabasePDBrev dbrev) {
		
		logger.debug("got a database revision:" + dbrev);
		
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
		PDBHeader header = structure.getPDBHeader();

		if ( header == null) {
			header = new PDBHeader();
		}

		if (dbrev.getNum().equals("1")){

			try {
				Date dep = dateFormat.parse(dbrev.getDate_original());
				header.setDepDate(dep);

			} catch (ParseException e){
				logger.warn("Could not parse date string '{}', deposition date will be unavailable", dbrev.getDate_original());
			}

			try {
				Date rel = dateFormat.parse(dbrev.getDate());
				header.setRelDate(rel);

			} catch (ParseException e){
				logger.warn("Could not parse date string '{}', modification date will be unavailable", dbrev.getDate());
			}


		} else {
			try {

				Date mod = dateFormat.parse(dbrev.getDate());
				header.setModDate(mod);

			} catch (ParseException e){
				logger.warn("Could not parse date string '{}', modification date will be unavailable", dbrev.getDate());
			}
		}

		structure.setPDBHeader(header);
	}
	
	@Override
	public void newPdbxAuditRevisionHistory(PdbxAuditRevisionHistory history) {
		
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
		PDBHeader header = structure.getPDBHeader();

		if ( header == null) {
			header = new PDBHeader();
		}

        // first entry in revision history is the release date
		if (history.getOrdinal().equals("1")){
			try {
				Date releaseDate = dateFormat.parse(history.getRevision_date());
				header.setRelDate(releaseDate);
				
			} catch (ParseException e){
				logger.warn("Could not parse date string '{}', release date will be unavailable", history.getRevision_date());
			}
		} else {
			// all other dates are revision dates;
			// since this method may be called multiple times,
			// the last revision date will "stick"
			try {
				Date revisionDate = dateFormat.parse(history.getRevision_date());
				header.setModDate(revisionDate);
			} catch (ParseException e){
				logger.warn("Could not parse date string '{}', revision date will be unavailable", history.getRevision_date());
			}
		}

		structure.setPDBHeader(header);
	}
	
	@Override
	public void newPdbxDatabaseStatus(PdbxDatabaseStatus status) {

		// the deposition date field is only available in mmCIF 5.0

		if (status.getRecvd_initial_deposition_date() == null) {
			// skip this method for older mmCIF versions
			return;
		}
		
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
		PDBHeader header = structure.getPDBHeader();

		if (header == null) {
			header = new PDBHeader();
		}

		try {
			Date depositionDate = dateFormat.parse(status.getRecvd_initial_deposition_date());
			header.setDepDate(depositionDate);
		} catch (ParseException e){
			logger.warn("Could not parse date string '{}', deposition date will be unavailable", status.getRecvd_initial_deposition_date());
		}

		structure.setPDBHeader(header);
	}

	@Override
	public void newDatabasePDBremark(DatabasePDBremark remark) {
		//System.out.println(remark);
		String id = remark.getId();
		if (id.equals("2")){

			//this remark field contains the resolution information:
			String line = remark.getText();

			int i = line.indexOf("ANGSTROM");
			if ( i > 5) {
				// line contains ANGSTROM info...
				String resolution = line.substring(i-5,i).trim();
				// convert string to float
				float res = 99 ;
				try {
					res = Float.parseFloat(resolution);

				} catch (NumberFormatException e) {
					logger.info("could not parse resolution from line and ignoring it " + line);
					return ;


				}
				// support for old style header

				PDBHeader pdbHeader = structure.getPDBHeader();
				pdbHeader.setResolution(res);

			}

		}
	}

	@Override
	public void newRefine(Refine r){

		PDBHeader pdbHeader = structure.getPDBHeader();
		// RESOLUTION
		// in very rare cases (for instance hybrid methods x-ray + neutron diffraction, e.g. 3ins, 4n9m)
		// there are 2 resolution values, one for each method
		// we take the last one found so that behaviour is like in PDB file parsing
		if (pdbHeader.getResolution()!=PDBHeader.DEFAULT_RESOLUTION) {
			logger.warn("More than 1 resolution value present, will use last one {} and discard previous {} "
					,r.getLs_d_res_high(), String.format("%4.2f",pdbHeader.getResolution()));
		}
		try {
			pdbHeader.setResolution(Float.parseFloat(r.getLs_d_res_high()));
		} catch (NumberFormatException e){
			logger.info("Could not parse resolution from " + r.getLs_d_res_high() + " " + e.getMessage());
		}


		// RFREE
		if (pdbHeader.getRfree()!=PDBHeader.DEFAULT_RFREE) {
			logger.warn("More than 1 Rfree value present, will use last one {} and discard previous {} ",
					r.getLs_R_factor_R_free(), String.format("%4.2f",pdbHeader.getRfree()));
		}
		if (r.getLs_R_factor_R_free()==null) {
			// some entries like 2ifo haven't got this field at all
			logger.info("_refine.ls_R_factor_R_free not present, not parsing Rfree value");
		} else {
			try {
				pdbHeader.setRfree(Float.parseFloat(r.getLs_R_factor_R_free()));
			} catch (NumberFormatException e){
				// no rfree present ('?') is very usual, that's why we set it to debug
				logger.debug("Could not parse Rfree from string '{}'", r.getLs_R_factor_R_free());
			}
		}

		// RWORK
		if(pdbHeader.getRwork()!=PDBHeader.DEFAULT_RFREE) {
			logger.warn("More than 1 R work value present, will use last one {} and discard previous {} ",
					r.getLs_R_factor_R_work(), String.format("%4.2f",pdbHeader.getRwork()));
		}
		if(r.getLs_R_factor_R_work()==null){
			logger.info("_refine.ls_R_factor_R_work not present, not parsing R-work value");
		}
		else{
			try{
				pdbHeader.setRwork(Float.parseFloat(r.getLs_R_factor_R_work()));
			}
			catch (NumberFormatException e){
				logger.debug("Could not parse R-work from string '{}'", r.getLs_R_factor_R_work());
			}

		}

	}


	@Override
	public void newAuditAuthor(AuditAuthor aa){

		String name =  aa.getName();

		StringBuffer famName = new StringBuffer();
		StringBuffer initials = new StringBuffer();
		boolean afterComma = false;
		for ( char c: name.toCharArray()) {
			if ( c == ' ')
				continue;
			if ( c == ','){
				afterComma = true;
				continue;
			}

			if ( afterComma)
				initials.append(c);
			else
				famName.append(c);
		}

		StringBuffer newaa = new StringBuffer();
		newaa.append(initials);
		newaa.append(famName);

		PDBHeader header = structure.getPDBHeader();
		String auth = header.getAuthors();
		if (auth == null) {
			header.setAuthors(newaa.toString());
		}else {
			auth += "," + newaa.toString();
			header.setAuthors(auth);

		}
	}

	@Override
	public void newExptl(Exptl exptl) {

		PDBHeader pdbHeader = structure.getPDBHeader();
		String method = exptl.getMethod();
		pdbHeader.setExperimentalTechnique(method);

	}

	@Override
	public void newCell(Cell cell) {

		try {
			float a = Float.parseFloat(cell.getLength_a());
			float b = Float.parseFloat(cell.getLength_b());
			float c = Float.parseFloat(cell.getLength_c());
			float alpha = Float.parseFloat(cell.getAngle_alpha());
			float beta = Float.parseFloat(cell.getAngle_beta());
			float gamma = Float.parseFloat(cell.getAngle_gamma());

			CrystalCell xtalCell = new CrystalCell();
			xtalCell.setA(a);
			xtalCell.setB(b);
			xtalCell.setC(c);
			xtalCell.setAlpha(alpha);
			xtalCell.setBeta(beta);
			xtalCell.setGamma(gamma);

			if (!xtalCell.isCellReasonable()) {
				// If the entry describes a structure determined by a technique other than X-ray crystallography,
				// cell is (sometimes!) a = b = c = 1.0, alpha = beta = gamma = 90 degrees
				// if so we don't add and CrystalCell will be null
				logger.debug("The crystal cell read from file does not have reasonable dimensions (at least one dimension is below {}), discarding it.",
						CrystalCell.MIN_VALID_CELL_SIZE);
				return;
			}

			structure.getPDBHeader().getCrystallographicInfo().setCrystalCell(xtalCell);

		} catch (NumberFormatException e){
			structure.getPDBHeader().getCrystallographicInfo().setCrystalCell(null);
			logger.info("could not parse some cell parameters ("+e.getMessage()+"), ignoring _cell ");
		}
	}

	@Override
	public void newSymmetry(Symmetry symmetry) {
		String spaceGroup = symmetry.getSpace_group_name_H_M();
		SpaceGroup sg = SymoplibParser.getSpaceGroup(spaceGroup);
		if (sg==null) {
			logger.warn("Space group '"+spaceGroup+"' not recognised as a standard space group");
			structure.getPDBHeader().getCrystallographicInfo().setNonStandardSg(true);
		} else {
			structure.getPDBHeader().getCrystallographicInfo().setSpaceGroup(sg);
			structure.getPDBHeader().getCrystallographicInfo().setNonStandardSg(false);
		}
	}

	@Override
	public void newStructNcsOper(StructNcsOper sNcsOper) {
		structNcsOper.add(sNcsOper);
	}
	
	public void newAtomSites(AtomSites atomSites) {
		
		try {
			Matrix4d m = new Matrix4d(
				Double.parseDouble(atomSites.getFract_transf_matrix11()), Double.parseDouble(atomSites.getFract_transf_matrix12()), Double.parseDouble(atomSites.getFract_transf_matrix13()), Double.parseDouble(atomSites.getFract_transf_vector1()),
				Double.parseDouble(atomSites.getFract_transf_matrix21()), Double.parseDouble(atomSites.getFract_transf_matrix22()), Double.parseDouble(atomSites.getFract_transf_matrix23()), Double.parseDouble(atomSites.getFract_transf_vector2()),
				Double.parseDouble(atomSites.getFract_transf_matrix31()), Double.parseDouble(atomSites.getFract_transf_matrix32()), Double.parseDouble(atomSites.getFract_transf_matrix33()), Double.parseDouble(atomSites.getFract_transf_vector3()),
				0,0,0,1);

			parsedScaleMatrix = m;
		
		} catch (NumberFormatException e) {
			logger.warn("Some values in _atom_sites.fract_transf_matrix or _atom_sites.fract_transf_vector could not be parsed as numbers. Can't check whether coordinate frame convention is correct! Error: {}", e.getMessage());
			structure.getPDBHeader().getCrystallographicInfo().setNonStandardCoordFrameConvention(false);
			
			// in this case parsedScaleMatrix stays null and can't be used in documentEnd()
		}
	}

	@Override
	public void newStructRef(StructRef sref) {
		logger.debug(sref.toString());
		strucRefs.add(sref);
	}

	private StructRef getStructRef(String ref_id){
		for (StructRef structRef : strucRefs) {

			if (structRef.getId().equals(ref_id)){
				return structRef;
			}

		}
		return null;

	}

	/**
	 * create a DBRef record from the StrucRefSeq record:
	 * <pre>
	 * PDB record                    DBREF
	 * Field Name                    mmCIF Data Item
	 * Section                       n.a.
	 * PDB_ID_Code                   _struct_ref_seq.pdbx_PDB_id_code
	 * Strand_ID                     _struct_ref_seq.pdbx_strand_id
	 * Begin_Residue_Number          _struct_ref_seq.pdbx_auth_seq_align_beg
	 * Begin_Ins_Code                _struct_ref_seq.pdbx_seq_align_beg_ins_code
	 * End_Residue_Number            _struct_ref_seq.pdbx_auth_seq_align_end
	 * End_Ins_Code                  _struct_ref_seq.pdbx_seq_align_end_ins_code
	 * Database                      _struct_ref.db_name
	 * Database_Accession_No         _struct_ref_seq.pdbx_db_accession
	 * Database_ID_Code              _struct_ref.db_code
	 * Database_Begin_Residue_Number _struct_ref_seq.db_align_beg
	 * Databaes_Begin_Ins_Code       _struct_ref_seq.pdbx_db_align_beg_ins_code
	 * Database_End_Residue_Number   _struct_ref_seq.db_align_end
	 * Databaes_End_Ins_Code         _struct_ref_seq.pdbx_db_align_end_ins_code
	 * </pre>
	 *
	 *
	 */
	@Override
	public void newStructRefSeq(StructRefSeq sref) {
		DBRef r = new DBRef();

		r.setIdCode(sref.getPdbx_PDB_id_code());
		r.setDbAccession(sref.getPdbx_db_accession());
		r.setDbIdCode(sref.getPdbx_db_accession());

		r.setChainName(sref.getPdbx_strand_id());
		StructRef structRef = getStructRef(sref.getRef_id());
		if (structRef == null){
			logger.info("could not find StructRef " + sref.getRef_id() + " for StructRefSeq " + sref);
		} else {
			r.setDatabase(structRef.getDb_name());
			r.setDbIdCode(structRef.getDb_code());
		}

		int seqbegin;
		int seqend;
		try{
			seqbegin = Integer.parseInt(sref.getPdbx_auth_seq_align_beg());
			seqend   = Integer.parseInt(sref.getPdbx_auth_seq_align_end());
		}
		catch(NumberFormatException e){
			// this happens in a few entries, annotation error? e.g. 6eoj
			logger.warn("Couldn't parse pdbx_auth_seq_align_beg/end in _struct_ref_seq. Will not store dbref alignment info for accession {}. Error: {}", r.getDbAccession(), e.getMessage());
			return;
		}
		
		Character begin_ins_code = ' ';
		if (sref.getPdbx_seq_align_beg_ins_code() != null ) {
		    begin_ins_code = new Character(sref.getPdbx_seq_align_beg_ins_code().charAt(0));
		}
		
		Character end_ins_code = ' ';
		if (sref.getPdbx_seq_align_end_ins_code() != null) {
		    end_ins_code = new Character(sref.getPdbx_seq_align_end_ins_code().charAt(0));
		}

		if (begin_ins_code == '?')
			begin_ins_code = ' ';

		if (end_ins_code == '?')
			end_ins_code = ' ';

		r.setSeqBegin(seqbegin);
		r.setInsertBegin(begin_ins_code);

		r.setSeqEnd(seqend);
		r.setInsertEnd(end_ins_code);

		int dbseqbegin = Integer.parseInt(sref.getDb_align_beg());
		int dbseqend   = Integer.parseInt(sref.getDb_align_end());
		
		Character db_begin_in_code = ' ';
		if (sref.getPdbx_db_align_beg_ins_code() != null) {
		    db_begin_in_code = new Character(sref.getPdbx_db_align_beg_ins_code().charAt(0));
		}
		
		Character db_end_in_code = ' ';
		if (sref.getPdbx_db_align_end_ins_code() != null) {
		    db_end_in_code = new Character(sref.getPdbx_db_align_end_ins_code().charAt(0));
		}

		if (db_begin_in_code == '?')
			db_begin_in_code = ' ';

		if (db_end_in_code == '?')
			db_end_in_code = ' ';


		r.setDbSeqBegin(dbseqbegin);
		r.setIdbnsBegin(db_begin_in_code);

		r.setDbSeqEnd(dbseqend);
		r.setIdbnsEnd(db_end_in_code);

		List<DBRef> dbrefs = structure.getDBRefs();
		if ( dbrefs == null)
			dbrefs = new ArrayList<DBRef>();
		dbrefs.add(r);

		logger.debug(r.toPDB());

		structure.setDBRefs(dbrefs);

	}

	@Override
	public void newStructRefSeqDif(StructRefSeqDif sref) {
		sequenceDifs.add(sref);
	}

	private Chain getEntityChain(String entity_id){

		for (Chain chain : entityChains) {
			if ( chain.getId().equals(entity_id)){

				return chain;
			}
		}
		// does not exist yet, so create...

		Chain	chain = new ChainImpl();
		chain.setId(entity_id);
		entityChains.add(chain);

		return chain;

	}

	//private Chain getSeqResChain(String chainID){
	//	return getChainFromList(seqResChains, chainID);
	//}


	/**
	 * Data items in the ENTITY_SRC_GEN category record details of
	 * the source from which the entity was obtained in cases
	 * where the source was genetically manipulated.  The
	 * following are treated separately:  items pertaining to the tissue
	 * from which the gene was obtained, items pertaining to the host
	 * organism for gene expression and items pertaining to the actual
	 * producing organism (plasmid).
	 */
	@Override
	public void newEntitySrcGen(EntitySrcGen entitySrcGen){

		// add to internal list. Map to Compound object later on...
		entitySrcGens.add(entitySrcGen);
	}

	@Override
	public void newEntitySrcNat(EntitySrcNat entitySrcNat){

		// add to internal list. Map to Compound object later on...
		entitySrcNats.add(entitySrcNat);
	}

	@Override
	public void newEntitySrcSyn(EntitySrcSyn entitySrcSyn){

		// add to internal list. Map to Compound object later on...
		entitySrcSyns.add(entitySrcSyn);
	}

	/**
	 * The EntityPolySeq object provide the amino acid sequence objects for the Entities.
	 * Later on the entities are mapped to the BioJava {@link Chain} and {@link EntityInfo} objects.
	 * @param epolseq the EntityPolySeq record for one amino acid
	 */
	@Override
	public void newEntityPolySeq(EntityPolySeq epolseq) {

		logger.debug("NEW entity poly seq " + epolseq);

		int eId = -1;
		try {
			eId = Integer.parseInt(epolseq.getEntity_id());
		} catch (NumberFormatException e) {
			logger.warn("Could not parse entity id from EntityPolySeq: "+e.getMessage());
		}
		Entity e = getEntity(eId);

		if (e == null){
			logger.info("Could not find entity "+ epolseq.getEntity_id()+". Can not match sequence to it.");
			return;
		}

		Chain entityChain = getEntityChain(epolseq.getEntity_id());

		// first we check through the chemcomp provider, if it fails we do some heuristics to guess the type of group
		// TODO some of this code is analogous to getNewGroup() and we should try to unify them - JD 2016-03-08

		Group g = ChemCompGroupFactory.getGroupFromChemCompDictionary(epolseq.getMon_id());
		//int seqId = Integer.parseInt(epolseq.getNum());
		if ( g != null && !g.getChemComp().isEmpty()) {
			if ( g instanceof AminoAcidImpl) {
				AminoAcidImpl aa = (AminoAcidImpl) g;
				aa.setRecordType(AminoAcid.SEQRESRECORD);
				//aa.setId(seqId);
			}
		} else {

			if (epolseq.getMon_id().length()==3 && StructureTools.get1LetterCodeAmino(epolseq.getMon_id())!=null){
				AminoAcidImpl a = new AminoAcidImpl();
				a.setRecordType(AminoAcid.SEQRESRECORD);
				Character code1 = StructureTools.get1LetterCodeAmino(epolseq.getMon_id());
				a.setAminoType(code1);
				g = a;

			} else if ( StructureTools.isNucleotide(epolseq.getMon_id())) {
				// the group is actually a nucleotide group...
				NucleotideImpl n = new NucleotideImpl();
				g = n;

			} else {
				logger.debug("Residue {} {} is not a standard aminoacid or nucleotide, will create a het group for it", epolseq.getNum(),epolseq.getMon_id());
				HetatomImpl h = new HetatomImpl();
				g = h;

			}


		}
		// at this stage we don't know about author residue numbers (insertion codes)
		// we abuse now the ResidueNumber field setting the internal residue numbers (label_seq_id, strictly sequential and follow the seqres sequence 1 to n)
		// later the actual ResidueNumbers (author residue numbers) have to be corrected in alignSeqRes()
		g.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));

		g.setPDBName(epolseq.getMon_id());

		entityChain.addGroup(g);

	}

	@Override
	public void newPdbxPolySeqScheme(PdbxPolySeqScheme ppss) {

		//if ( headerOnly)
		//	return;

		// replace the group asym ids with the real PDB ids!
		// replaceGroupSeqPos(ppss);  // This might be incorrect in some pdb, to use auth_seq_id of the pdbx_poly_seq_scheme.


	}


	@Override
	public void newPdbxNonPolyScheme(PdbxNonPolyScheme ppss) {

		//if (headerOnly)
		//	return;

		// merge the EntityPolySeq info and the AtomSite chains into one...
		//already known ignore:

	}

	@Override
	public void newPdbxEntityNonPoly(PdbxEntityNonPoly pen){
		// TODO: do something with them...
		// not implemented yet...
		logger.debug(pen.getEntity_id() + " " + pen.getName() + " " + pen.getComp_id());

	}

	@Override
	public void newChemComp(ChemComp c) {
		// TODO: do something with them...

	}

	@Override
	public void newGenericData(String category, List<String> loopFields,
			List<String> lineData) {

		//logger.debug("unhandled category so far: " + category);
	}

	@Override
	public FileParsingParameters getFileParsingParameters()
	{
		return params;
	}

	@Override
	public void setFileParsingParameters(FileParsingParameters params)
	{
		this.params = params;

	}

	@Override
	public void newChemCompDescriptor(ChemCompDescriptor ccd) {

		// TODO nothing happening here yet.

	}



	public List<PdbxStructOperList> getStructOpers() {
		return structOpers;
	}

	@Override
	public void newPdbxStrucAssembly(PdbxStructAssembly strucAssembly) {
		strucAssemblies.add(strucAssembly);

	}

	public List<PdbxStructAssembly> getStructAssemblies(){
		return strucAssemblies;
	}

	@Override
	public void newPdbxStrucAssemblyGen(PdbxStructAssemblyGen strucAssembly) {
		strucAssemblyGens.add(strucAssembly);

	}

	public List<PdbxStructAssemblyGen> getStructAssemblyGens(){
		return strucAssemblyGens;
	}

	@Override
	public void newChemCompAtom(ChemCompAtom atom) {

	}

	@Override
	public void newPdbxChemCompIndentifier(PdbxChemCompIdentifier id) {

	}

	@Override
	public void newChemCompBond(ChemCompBond bond) {

	}

	@Override
	public void newPdbxChemCompDescriptor(PdbxChemCompDescriptor desc) {

	}

	@Override
	public void newStructConn(StructConn structConn) {
		this.structConn.add(structConn);
	}

	@Override
	public void newStructSiteGen(StructSiteGen siteGen) { this.structSiteGens.add(siteGen);	}

	@Override
	public void newStructSite(StructSite structSite) {

		if (params.isHeaderOnly()) {
			return;
		}

		// Simply implement the method.
		List<Site> sites = structure.getSites();
		if (sites == null) sites = new ArrayList<Site>();

		Site site = null;
		for (Site asite : sites) {
			if (asite.getSiteID().equals(structSite.getId())) {
				site = asite; 		// Prevent duplicate siteIds
			}
		}
		boolean addSite = false;
		if (site == null) { site = new Site(); addSite = true; }
		site.setSiteID(structSite.getId());
		site.setDescription(structSite.getDetails());
		// site.setPdbxEvidenceCode(structSite.getPdbxEvidenceCode()); // TODO - add addition fields in Sites
		if (addSite) sites.add(site);

		structure.setSites(sites);
	}

	/**
	 * Build sites in a BioJava Structure using the original author chain id & residue numbers.
	 * Sites are built from struct_site_gen records that have been parsed.
	 */
	private void addSites() {
		List<Site> sites = structure.getSites();
		if (sites == null) sites = new ArrayList<Site>();

		for (StructSiteGen siteGen : structSiteGens) {
			// For each StructSiteGen, find the residues involved, if they exist then
			String site_id = siteGen.getSite_id(); // multiple could be in same site.
			if (site_id == null) site_id = "";
			String comp_id = siteGen.getLabel_comp_id();  // PDBName

			// Assumption: the author chain ID and residue number for the site is consistent with the original
			// author chain id and residue numbers.

			String asymId = siteGen.getLabel_asym_id(); // chain name
			String authId = siteGen.getAuth_asym_id(); // chain Id
			String auth_seq_id = siteGen.getAuth_seq_id(); // Res num

			String insCode = siteGen.getPdbx_auth_ins_code();
			if ( insCode != null && insCode.equals("?"))
				insCode = null;

			// Look for asymID = chainID and seqID = seq_ID.  Check that comp_id matches the resname.
			Group g = null;
			try {
				Chain chain = structure.getChain(asymId);

				if (null != chain) {
					try {
						Character insChar = null;
						if (null != insCode && insCode.length() > 0) insChar = insCode.charAt(0);
						g = chain.getGroupByPDB(new ResidueNumber(null, Integer.parseInt(auth_seq_id), insChar));
					} catch (NumberFormatException e) {
						logger.warn("Could not lookup residue : " + authId + auth_seq_id);
					}
				}
			} catch (StructureException e) {
				logger.warn("Problem finding residue in site entry " + siteGen.getSite_id() + " - " + e.getMessage(), e.getMessage());
			}

			if (g != null) {
				// 2. find the site_id, if not existing, create anew.
				Site site = null;
				for (Site asite: sites) {
					if (site_id.equals(asite.getSiteID())) site = asite;
				}

				boolean addSite = false;

				// 3. add this residue to the site.
				if (site == null) {
					addSite = true;
					site = new Site();
					site.setSiteID(site_id);
				}

				List<Group> groups = site.getGroups();
				if (groups == null) groups = new ArrayList<Group>();

				// Check the self-consistency of the residue reference from auth_seq_id and chain_id
				if (!comp_id.equals(g.getPDBName())) {
					logger.warn("comp_id doesn't match the residue at " + authId + " " + auth_seq_id + " - skipping");
				} else {
					groups.add(g);
					site.setGroups(groups);
				}
				if (addSite) sites.add(site);
			}
		}
		structure.setSites(sites);
	}
}
