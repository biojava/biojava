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
import org.biojava.nbio.structure.Compound;
import org.biojava.nbio.structure.DBRef;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.HetatomImpl;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.SSBond;
import org.biojava.nbio.structure.SSBondImpl;
import org.biojava.nbio.structure.SeqMisMatch;
import org.biojava.nbio.structure.SeqMisMatchImpl;
import org.biojava.nbio.structure.Site;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.BondMaker;
import org.biojava.nbio.structure.io.ChargeAdder;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.LigandConnectMaker;
import org.biojava.nbio.structure.io.SeqRes2AtomAligner;
import org.biojava.nbio.structure.io.mmcif.model.AtomSite;
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
import org.biojava.nbio.structure.io.mmcif.model.EntityPolySeq;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcGen;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcNat;
import org.biojava.nbio.structure.io.mmcif.model.EntitySrcSyn;
import org.biojava.nbio.structure.io.mmcif.model.Exptl;
import org.biojava.nbio.structure.io.mmcif.model.PdbxChemCompDescriptor;
import org.biojava.nbio.structure.io.mmcif.model.PdbxChemCompIdentifier;
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

/** A MMcifConsumer implementation that build a in-memory representation of the
 * content of a mmcif file as a BioJava Structure object.
 *  @author Andreas Prlic
 *  @since 1.7
 */

public class SimpleMMcifConsumer implements MMcifConsumer {

	private static final Logger logger = LoggerFactory.getLogger(SimpleMMcifConsumer.class);


	private Structure structure;
	private Chain current_chain;
	private Group current_group;


	private List<Chain>      current_model;
	private List<Entity>     entities;
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

	/**
	 * A map of asym ids (internal chain ids) to strand ids (author chain ids) 
	 * extracted from pdbx_poly_seq_scheme/pdbx_non_poly_seq_scheme categories
	 */
	private Map<String,String> asymStrandId;

	/**
	 * A map of asym ids (internal chain ids) to strand ids (author chain ids) 
	 * extracted from the information in _atom_sites category. Will be used
	 * if no mapping is found in pdbx_poly_seq_scheme/pdbx_non_poly_seq_scheme
	 */
	private Map<String,String> asymId2StrandIdFromAtomSites;

	/**
	 * A map of asym ids (internal chain ids) to entity ids extracted from 
	 * the _struct_asym category
	 */
	private Map<String,String> asymId2entityId;

	private String current_nmr_model ;

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

		if ( params.isLoadChemCompInfo() ){
			Group g = ChemCompGroupFactory.getGroupFromChemCompDictionary(groupCode3);
			if ( g != null) {
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

	/** test if the chain is already known (is in current_model
	 * ArrayList) and if yes, returns the chain
	 * if no -> returns null
	 */
	private Chain isKnownChain(String chainID, List<Chain> chains){

		for (int i = 0; i< chains.size();i++){
			Chain testchain =  chains.get(i);
			//System.out.println("comparing chainID >"+chainID+"< against testchain " + i+" >" +testchain.getName()+"<");
			if (chainID.equals(testchain.getChainID())) {
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


		// Here we can either use Internal or public facing chain ids]
		String chain_id;
		if(params.isUseInternalChainIds()==true){
			chain_id = atom.getLabel_asym_id();	
		}
		else{
			chain_id = atom.getAuth_asym_id();
		}


		String recordName    = atom.getGroup_PDB();
		String residueNumberS = atom.getAuth_seq_id();
		Integer residueNrInt = Integer.parseInt(residueNumberS);

		// the 3-letter name of the group:
		String groupCode3    = atom.getLabel_comp_id();

		Character aminoCode1 = null;
		if ( recordName.equals("ATOM") )
			aminoCode1 = StructureTools.get1LetterCodeAmino(groupCode3);
		else {
			aminoCode1 = StructureTools.get1LetterCodeAmino(groupCode3);

			// for nucleotides this will be null..
			if (aminoCode1 != null &&  aminoCode1.equals(StructureTools.UNKNOWN_GROUP_LABEL)) 
				aminoCode1 = null;
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

		String nmrModel = atom.getPdbx_PDB_model_num();

		if ( current_nmr_model == null) {
			current_nmr_model = nmrModel;
		}

		if (! current_nmr_model.equals(nmrModel)){
			current_nmr_model = nmrModel;

			// add previous data
			if ( current_chain != null ) {
				current_chain.addGroup(current_group);
				current_group.trimToSize();
			}

			// we came to the beginning of a new NMR model
			structure.addModel(current_model);
			current_model = new ArrayList<Chain>();
			current_chain = null;
			current_group = null;
		}


		if (current_chain == null) {
			current_chain = new ChainImpl();
			current_chain.setChainID(chain_id);
			current_model.add(current_chain);
			startOfNewChain = true;
		}

		//System.out.println("BEFORE: " + chain_id + " " + current_chain.getName());
		if ( ! chain_id.equals(current_chain.getChainID()) ) {
			startOfNewChain = true;

			// end up old chain...
			current_chain.addGroup(current_group);

			// see if old chain is known ...
			Chain testchain ;
			testchain = isKnownChain(current_chain.getChainID(),current_model);

			//System.out.println("trying to re-using known chain " + current_chain.getName() + " " + chain_id);		
			if ( testchain != null && testchain.getChainID().equals(chain_id)){
				//System.out.println("re-using known chain " + current_chain.getName() + " " + chain_id);				

			} else {

				testchain = isKnownChain(chain_id,current_model);
			}

			if ( testchain == null) {
				//System.out.println("unknown chain. creating new chain.");

				current_chain = new ChainImpl();
				current_chain.setChainID(chain_id);

			}   else {
				current_chain = testchain;
			}

			if ( ! current_model.contains(current_chain))
				current_model.add(current_chain);

		} 


		ResidueNumber residueNumber = new ResidueNumber(chain_id,residueNrInt, insCode);

		if (current_group == null) {

			current_group = getNewGroup(recordName,aminoCode1,seq_id, groupCode3);

			current_group.setResidueNumber(residueNumber);
			current_group.setPDBName(groupCode3);
		}

		if ( startOfNewChain){
			current_group = getNewGroup(recordName,aminoCode1,seq_id, groupCode3);

			current_group.setResidueNumber(residueNumber);
			current_group.setPDBName(groupCode3);
		}

		Group altGroup = null;
		String altLocS = atom.getLabel_alt_id();
		Character altLoc = ' ';
		if ( altLocS.length()>0) {
			altLoc = altLocS.charAt(0);
			if ( altLoc.equals('.') )
				altLoc = ' ';

		}

		// check if residue number is the same ...
		// insertion code is part of residue number
		if ( ! residueNumber.equals(current_group.getResidueNumber())) {
			//System.out.println("end of residue: "+current_group.getPDBCode()+" "+residueNrInt);
			current_chain.addGroup(current_group);
			current_group.trimToSize();
			current_group = getNewGroup(recordName,aminoCode1,seq_id,groupCode3);
			current_group.setPDBName(groupCode3);
			current_group.setResidueNumber(residueNumber);


			//                        System.out.println("Made new group:  " + groupCode3 + " " + resNum + " " + iCode);

		} else {
			// same residueNumber, but altLocs...

			// test altLoc
			if ( ! altLoc.equals(' ') && ( ! altLoc.equals('.'))) {												
				logger.debug("found altLoc! " + altLoc + " " + current_group + " " + altGroup);
				altGroup = getCorrectAltLocGroup( altLoc,recordName,aminoCode1,groupCode3, seq_id);
				if (altGroup.getChain()==null) {
					altGroup.setChain(current_chain);
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

		// filling the map in case there's no pdbx_poly_seq_scheme/pdbx_non_poly_seq_scheme in the file
		asymId2StrandIdFromAtomSites.put(atom.getLabel_asym_id(), atom.getAuth_asym_id());

		//see if chain_id is one of the previous chains ...

		Atom a = convertAtom(atom);

		//see if chain_id is one of the previous chains ...
		if ( altGroup != null) {
			altGroup.addAtom(a);
			altGroup = null;
		}
		else {
			current_group.addAtom(a);
		}


		// make sure that main group has all atoms
		// GitHub issue: #76
		if ( ! current_group.hasAtom(a.getName())) {
			current_group.addAtom(a);
		}


		//System.out.println(">" + atom.getLabel_atom_id()+"< " + a.getGroup().getPDBName() + " " + a.getGroup().getChemComp()  );

		//System.out.println(current_group);

	}

	/** convert a MMCif AtomSite object to a BioJava Atom object
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
			String recordName, Character aminoCode1, String groupCode3, long seq_id) {

		// see if we know this altLoc already;
		List<Atom> atoms = current_group.getAtoms();
		if ( atoms.size() > 0) {
			Atom a1 = atoms.get(0);
			// we are just adding atoms to the current group
			// probably there is a second group following later...
			if (a1.getAltLoc().equals(altLoc)) {

				return current_group;
			}
		}

		List<Group> altLocs = current_group.getAltLocs();
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

		if ( groupCode3.equals(current_group.getPDBName())) {
			if ( current_group.getAtoms().size() == 0) {
				//System.out.println("current group is empty " + current_group + " " + altLoc);
				return current_group;
			}
			//System.out.println("cloning current group " + current_group + " " + current_group.getAtoms().get(0).getAltLoc() + " altLoc " + altLoc);
			Group altLocG = (Group) current_group.clone();
			// drop atoms from cloned group...
			// https://redmine.open-bio.org/issues/3307
			altLocG.setAtoms(new ArrayList<Atom>());
			altLocG.getAltLocs().clear();
			current_group.addAltLoc(altLocG);
			return altLocG;	
		}

		//	System.out.println("new  group " + recordName + " " + aminoCode1 + " " +groupCode3);
		//String recordName,Character aminoCode1, long seq_id,String groupCode3) {
		Group altLocG = getNewGroup(recordName,aminoCode1,seq_id,groupCode3);

		altLocG.setPDBName(groupCode3);
		altLocG.setResidueNumber(current_group.getResidueNumber());
		current_group.addAltLoc(altLocG);
		return altLocG;
	}

	/** Start the parsing
	 *
	 */
	@Override
	public void documentStart() {
		structure = new StructureImpl();

		current_chain 		= null;
		current_group 		= null;
		current_nmr_model 	= null;
		//atomCount     		= 0;

		current_model = new ArrayList<Chain>();
		entities      = new ArrayList<Entity>();
		strucRefs     = new ArrayList<StructRef>();
		seqResChains  = new ArrayList<Chain>();
		entityChains  = new ArrayList<Chain>();
		structAsyms   = new ArrayList<StructAsym>();
		asymStrandId  = new HashMap<String, String>();
		asymId2StrandIdFromAtomSites = new HashMap<String, String>();
		asymId2entityId = new HashMap<String,String>();		
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
		if ( current_chain != null ) {

			current_chain.addGroup(current_group);
			if (isKnownChain(current_chain.getChainID(),current_model) == null) {
				current_model.add(current_chain);
			}
		} else if (!params.isHeaderOnly()){
			logger.warn("current chain is null at end of document.");			
		}

		structure.addModel(current_model);

		// Goal is to reproduce the PDB files exactly:
		// What has to be done is to use the auth_mon_id for the assignment. For this

		// map entities to Chains and Compound objects...


		for (StructAsym asym : structAsyms) {
			logger.debug("Entity {} matches asym_id: {}", asym.getEntity_id(), asym.getId() );

			asymId2entityId.put(asym.getId(), asym.getEntity_id());

			Chain s = getEntityChain(asym.getEntity_id());
			Chain seqres = (Chain)s.clone();
			// to solve issue #160 (e.g. 3u7t)
			seqres = removeSeqResHeterogeneity(seqres);
			seqres.setChainID(asym.getId());

			seqResChains.add(seqres);
			logger.debug(" seqres: " + asym.getId() + " " + seqres + "<") ;

			// adding the compounds (entities)
			addCompounds(asym);

		}

		if (structAsyms.isEmpty()) {
			logger.warn("No _struct_asym category in file, no SEQRES groups will be added."); 
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

		// Various tasks that require Atoms can be skipped.
		if (!params.isHeaderOnly()) {
			if ( params.shouldCreateAtomBonds()) {
				addBonds();
			}

			// Adds the equivalent CONECT-style records as a PDB has for ligands.
			if ( params.isCreateLigandConects()) {
				addLigandConnections();
			}

			if ( params.shouldCreateAtomCharges()) {
				addCharges();
			}
		}

		if (asymStrandId.isEmpty()) {
			logger.warn("No pdbx_poly_seq_scheme/pdbx_non_poly_seq_scheme categories present. Will use chain id mapping from _atom_sites category");

			asymStrandId = asymId2StrandIdFromAtomSites;
		}

		// mismatching Author assigned chain IDS and PDB internal chain ids:
		// fix the chain IDS in the current model:
		if(params.isUseInternalChainIds()==false){
			for (int i =0; i< structure.nrModels() ; i++){
				List<Chain> model = structure.getModel(i);

				List<Chain> pdbChains = new ArrayList<Chain>();

				for (Chain chain : model) {
					for (String asym : asymStrandId.keySet()) {
						if ( chain.getChainID().equals(asym)){
							String newChainId = asymStrandId.get(asym);

							logger.debug("Renaming chain with asym_id {} ({} atom groups) to author_asym_id/strand_id  {}", 
									asym, chain.getAtomGroups().size(), newChainId);

							chain.setChainID(newChainId);
							chain.setInternalChainID(asym);
							// set chain of all groups
							for(Group g : chain.getAtomGroups()) {
								ResidueNumber resNum = g.getResidueNumber();
								if(resNum != null)
									resNum.setChainId(newChainId);
							}
							for(Group g : chain.getSeqResGroups()) {
								ResidueNumber resNum = g.getResidueNumber();
								if(resNum != null)
									resNum.setChainId(newChainId);
							}
							Chain known =  isKnownChain(chain.getChainID(), pdbChains);
							if ( known == null ){
								pdbChains.add(chain);
							} else {
								// and now we join the 2 chains together again, because in cif files the data can be split up...
								for ( Group g : chain.getAtomGroups()){
									known.addGroup(g);
								}
							}

							break;
						}
					}
				}

				structure.setModel(i,pdbChains);
			}
		}

		// compounds (entities)
		// In addCompounds above we created the compounds if they were present in the file
		// Now we need to make sure that they are linked to chains and also that if they are not present in the file we need to add them now
		linkCompounds();

		if (!params.isHeaderOnly()) {
			createSSBonds();

			// Do structure.setSites(sites) after any chain renaming to be like PDB.
			addSites();
		}



		// set the oligomeric state info in the header...
		if (params.isParseBioAssembly()) {

			// the more detailed mapping of chains to rotation operations happens in StructureIO...

			Map<Integer,BioAssemblyInfo> bioAssemblies = new HashMap<Integer, BioAssemblyInfo>();

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

				int mmSize = 0;
				int bioAssemblyId = -1;
				try {
					bioAssemblyId = Integer.parseInt(psa.getId());
				} catch (NumberFormatException e) {
					logger.info("Could not parse a numerical bio assembly id from '{}'",psa.getId());
				}
				try {
					mmSize = Integer.parseInt(psa.getOligomeric_count());					
				} catch (NumberFormatException e) {
					if (bioAssemblyId!=-1)
						// if we have a numerical id, then it's unusual to have no oligomeric size: we warn about it
						logger.warn("Could not parse oligomeric count from '{}' for biological assembly id {}",
								psa.getOligomeric_count(),psa.getId());
					else 
						// no numerical id (PAU,XAU in virus entries), it's normal to have no oligomeric size
						logger.info("Could not parse oligomeric count from '{}' for biological assembly id {}",
								psa.getOligomeric_count(),psa.getId());
				}

				// if bioassembly id is not numerical we throw it away
				// this happens usually for viral capsid entries, like 1ei7
				// see issue #230 in github
				if (bioAssemblyId!=-1) {
					BioAssemblyInfo bioAssembly = new BioAssemblyInfo();
					bioAssembly.setId(bioAssemblyId);
					bioAssembly.setMacromolecularSize(mmSize);
					bioAssembly.setTransforms(transformations);
					bioAssemblies.put(bioAssemblyId,bioAssembly);
				}

			}
			structure.getPDBHeader().setBioAssemblies(bioAssemblies);
		}

		ArrayList<Matrix4d> ncsOperators = new ArrayList<Matrix4d>();
		for (StructNcsOper sNcsOper:structNcsOper) {
			if (sNcsOper.getCode().equals("generate")) {
				ncsOperators.add(sNcsOper.getOperator());
			}
		}
		// we only set it if not empty, otherwise remains null
		if (ncsOperators.size()>0) {
			structure.getCrystallographicInfo().setNcsOperators(
					ncsOperators.toArray(new Matrix4d[ncsOperators.size()]));
		}


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
			try {
				Chain c = structure.getChainByPDB(chainId);
				c.setSeqMisMatches(misMatchMap.get(chainId));
			} catch (Exception e){
				logger.warn("could not set mismatches for chain " + chainId);

			}
		}

	}

	/**
	 * Here we link compounds (entities) to chains.
	 * Also if compounds are not present in file, this initialises the compounds with some heuristics, see {@link CompoundFinder}
	 */
	private void linkCompounds() {


		for (int i =0; i< structure.nrModels() ; i++){
			for (Chain chain : structure.getModel(i)) {

				String entityId = asymId2entityId.get(chain.getInternalChainID());
				if (entityId==null) {
					// this can happen for instance if the cif file didn't have _struct_asym category at all
					// and thus we have no asymId2entityId mapping at all
					logger.warn("No entity id could be found for chain {}", chain.getInternalChainID());					
					continue;
				}
				int eId = Integer.parseInt(entityId);

				// Compounds are not added for non-polymeric entities, if a chain is non-polymeric its compound won't be found. 
				// TODO: add all entities and unique compounds and add methods to directly get polymer or non-polymer
				// asyms (chains).  Either create a unique StructureImpl or modify existing for a better representation of the
				// mmCIF internal data structures but is compatible with Structure interface.
				// Some examples of PDB entries with this kind of problem:
				//   - 2uub: asym_id X, chainId Z, entity_id 24: fully non-polymeric but still with its own chainId
				//   - 3o6j: asym_id K, chainId Z, entity_id 6 : a single water molecule
				//   - 1dz9: asym_id K, chainId K, entity_id 6 : a potassium ion alone 

				Compound compound = structure.getCompoundById(eId);
				if (compound==null) {
					// Supports the case where the only chain members were from non-polymeric entity that is missing.
					// Solved by creating a new Compound(entity) to which this chain will belong.
					logger.warn("Could not find a compound for entity_id {}, for chain id {}, creating a new compound.",
							eId, chain.getChainID());
					compound = new Compound();
					compound.setMolId(eId);
					compound.addChain(chain);
					chain.setCompound(compound);
					structure.addCompound(compound);
				} else {
					logger.debug("Adding chain with chain id {} (asym id {}) to compound with entity_id {}",
							chain.getChainID(), chain.getInternalChainID(), eId);
					compound.addChain(chain);
					chain.setCompound(compound);
				}

			}			

		}

		// to make sure we have Compounds linked to chains, we call getCompounds() which will lazily initialise the
		// compounds using heuristics (see CompoundFinder) in the case that they were not explicitly present in the file
		List<Compound> compounds = structure.getCompounds();

		// final sanity check: it can happen that from the annotated compounds some are not linked to any chains
		// e.g. 3s26: a sugar entity does not have any chains associated to it (it seems to be happening with many sugar compounds)
		// we simply log it, this can sign some other problems if the compounds are used down the line
		for (Compound compound:compounds) {
			if (compound.getChains().isEmpty()) {
				logger.info("Compound {} '{}' has no chains associated to it",
						compound.getId()==null?"with no entity id":compound.getId(), compound.getMolName());
			}
		}

	}

	private void addCharges() {
		ChargeAdder adder = new ChargeAdder(structure);
		adder.addCharges();
	}

	/**
	 * The method will return a new reference to a Chain with any consecutive groups 
	 * having same residue numbers removed.
	 * This is necessary to solve the microheterogeneity issue in entries like 3u7t (see github issue #160)
	 * @param c
	 * @return
	 */
	private Chain removeSeqResHeterogeneity(Chain c) {

		Chain trimmedChain = new ChainImpl();

		ResidueNumber lastResNum = null;

		for (Group g:c.getAtomGroups()) {

			// note we have to deep copy this, otherwise they stay linked and would get altered in addGroup(g) 
			ResidueNumber currentResNum = new ResidueNumber(
					g.getResidueNumber().getChainId(),
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
		BondMaker maker = new BondMaker(structure);
		maker.makeBonds();	
	}

	private void addLigandConnections(){
		LigandConnectMaker maker = new LigandConnectMaker(structure);
		maker.addLigandConnections();
	}

	private void alignSeqRes() {

		logger.debug("Parsing mode align_seqres, will align to ATOM to SEQRES sequence");

		// fix SEQRES residue numbering for all models

		for (int model=0;model<structure.nrModels();model++) {

			List<Chain> atomList   = structure.getModel(model);

			for (Chain seqResChain: seqResChains){

				// this extracts the matching atom chain from atomList
				Chain atomChain = SeqRes2AtomAligner.getMatchingAtomRes(seqResChain, atomList);

				if (atomChain == null) {
					// most likely there's no observed residues at all for the seqres chain: can't map
					// e.g. 3zyb: chains with asym_id L,M,N,O,P have no observed residues
					logger.warn("Could not map SEQRES chain with asym_id={} to any ATOM chain. Most likely there's no observed residues in the chain.",
							seqResChain.getChainID());
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

	private void addCompounds(StructAsym asym) {
		int eId = 0;
		try {
			eId = Integer.parseInt(asym.getEntity_id());
		} catch (NumberFormatException e) {
			logger.warn("Could not parse mol_id from string {}. Will use 0 for creating Compound",asym.getEntity_id());
		}
		Entity e = getEntity(eId);

		for (EntitySrcGen esg : entitySrcGens) {

			if (! esg.getEntity_id().equals(asym.getEntity_id()))
				continue;

			// found the matching EntitySrcGen
			// get the corresponding Entity
			Compound c = structure.getCompoundById(eId);
			if ( c == null){
				if (e!=null) { 
					if (e.getType().equals("polymer")) {
						c = createNewCompoundFromESG(esg, eId);
						c.setMolName(e.getPdbx_description());
						structure.addCompound(c);
						logger.debug("Adding Compound with entity id {} from _entity_src_syn, with name: {}",eId,c.getMolName());
					} else if (e.getType().equals("non-solvent")) {
						// TODO handle non-polymer compounds.
					} else if (e.getType().equals("water")) {
						// TODO handle solvent entity.
					} else {
						logger.warn("Could not add entity id " + esg.getEntity_id() + " that has unknown _entity.type");
					}
				} 
			}

		}

		for (EntitySrcNat esn : entitySrcNats) {
			if (! esn.getEntity_id().equals(asym.getEntity_id()))
				continue;

			// found the matching EntitySrcGen
			// get the corresponding Entity
			Compound c = structure.getCompoundById(eId);
			if ( c == null){		
				if (e!=null && e.getType().equals("polymer")) {
					c = createNewCompoundFromESN(esn, eId);
					c.setMolName(e.getPdbx_description());
					structure.addCompound(c);
					logger.debug("Adding Compound with entity id {} from _entity_src_syn, with name: {}",eId,c.getMolName());
				}
			}

		}

		for (EntitySrcSyn ess : entitySrcSyns) {
			if (! ess.getEntity_id().equals(asym.getEntity_id()))
				continue;

			// found the matching EntitySrcGen
			// get the corresponding Entity
			Compound c = structure.getCompoundById(eId);
			if ( c == null){	
				if (e!=null && e.getType().equals("polymer")) {
					c = createNewCompoundFromESS(ess, eId);
					c.setMolName(e.getPdbx_description());
					structure.addCompound(c);
					logger.debug("Adding Compound with entity id {} from _entity_src_syn, with name: {}",eId,c.getMolName());
				}
			}
		}

		// for some mmCIF files like 1yrm all 3 of _entity_src_gen, _entity_src_nat and _pdbx_entity_src_syn are missing
		// we need to fill the Compounds in some other way:

		Compound c = structure.getCompoundById(eId);

		if (c==null) {
			c = new Compound();
			c.setMolId(eId);

			// we only add the compound if a polymeric one (to match what the PDB parser does)
			if (e!=null && e.getType().equals("polymer")) {
				c.setMolName(e.getPdbx_description());
				structure.addCompound(c);
				logger.debug("Adding Compound with entity id {} from _entity, with name: {}",eId, c.getMolName());
			}
		}
	}

	private Compound createNewCompoundFromESG(EntitySrcGen esg, int eId) {

		Compound c = new Compound();
		c.setMolId(eId);
		c.setAtcc(esg.getPdbx_gene_src_atcc());
		c.setCell(esg.getPdbx_gene_src_cell());
		c.setOrganismCommon(esg.getGene_src_common_name());
		c.setOrganismScientific(esg.getPdbx_gene_src_scientific_name());
		c.setOrganismTaxId(esg.getPdbx_gene_src_ncbi_taxonomy_id());
		c.setExpressionSystemTaxId(esg.getPdbx_host_org_ncbi_taxonomy_id());
		c.setExpressionSystem(esg.getPdbx_host_org_scientific_name());
		return c;

	}

	private Compound createNewCompoundFromESN(EntitySrcNat esn, int eId) {

		Compound c = new Compound();

		c.setMolId(eId);
		c.setAtcc(esn.getPdbx_atcc());
		c.setCell(esn.getPdbx_cell());
		c.setOrganismCommon(esn.getCommon_name());
		c.setOrganismScientific(esn.getPdbx_organism_scientific());
		c.setOrganismTaxId(esn.getPdbx_ncbi_taxonomy_id());

		return c;

	}

	private Compound createNewCompoundFromESS(EntitySrcSyn ess, int eId) {

		Compound c = new Compound();

		c.setMolId(eId);
		c.setOrganismCommon(ess.getOrganism_common_name());
		c.setOrganismScientific(ess.getOrganism_scientific());
		c.setOrganismTaxId(ess.getNcbi_taxonomy_id());


		return c;

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
		//System.out.println("got a database revision:" + dbrev);
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
				Date mod = dateFormat.parse(dbrev.getDate());
				header.setModDate(mod);

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
		if (sg==null) logger.warn("Space group '"+spaceGroup+"' not recognised as a standard space group"); 

		structure.getPDBHeader().getCrystallographicInfo().setSpaceGroup(sg); 
	}

	@Override
	public void newStructNcsOper(StructNcsOper sNcsOper) {
		structNcsOper.add(sNcsOper);
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

	/** create a DBRef record from the StrucRefSeq record:
	 *  <pre>
  PDB record 					DBREF
  Field Name 					mmCIF Data Item
  Section   	  				n.a.
  PDB_ID_Code   	  			_struct_ref_seq.pdbx_PDB_id_code
  Strand_ID   	 			 	_struct_ref_seq.pdbx_strand_id
  Begin_Residue_Number   	  	_struct_ref_seq.pdbx_auth_seq_align_beg
  Begin_Ins_Code   	  			_struct_ref_seq.pdbx_seq_align_beg_ins_code
  End_Residue_Number   	  		_struct_ref_seq.pdbx_auth_seq_align_end
  End_Ins_Code   	  			_struct_ref_seq.pdbx_seq_align_end_ins_code
  Database   	  				_struct_ref.db_name
  Database_Accession_No   	  	_struct_ref_seq.pdbx_db_accession
  Database_ID_Code   	  		_struct_ref.db_code
  Database_Begin_Residue_Number	_struct_ref_seq.db_align_beg
  Databaes_Begin_Ins_Code   	_struct_ref_seq.pdbx_db_align_beg_ins_code
  Database_End_Residue_Number  	_struct_ref_seq.db_align_end
  Databaes_End_Ins_Code   	  	_struct_ref_seq.pdbx_db_align_end_ins_code
  </pre>
	 *
	 *
	 */
	@Override
	public void newStructRefSeq(StructRefSeq sref) {
		//if (DEBUG)
		//	System.out.println(sref);
		DBRef r = new DBRef();


		//if (DEBUG)
		//	System.out.println( " " + sref.getPdbx_PDB_id_code() + " " + sref.getPdbx_db_accession());
		r.setIdCode(sref.getPdbx_PDB_id_code());
		r.setDbAccession(sref.getPdbx_db_accession());
		r.setDbIdCode(sref.getPdbx_db_accession());

		r.setChainId(sref.getPdbx_strand_id());
		StructRef structRef = getStructRef(sref.getRef_id());
		if (structRef == null){
			logger.warn("could not find StructRef " + sref.getRef_id() + " for StructRefSeq " + sref);
		} else {
			r.setDatabase(structRef.getDb_name());
			r.setDbIdCode(structRef.getDb_code());
		}


		int seqbegin = Integer.parseInt(sref.getPdbx_auth_seq_align_beg());
		int seqend   = Integer.parseInt(sref.getPdbx_auth_seq_align_end());
		Character begin_ins_code = new Character(sref.getPdbx_seq_align_beg_ins_code().charAt(0));
		Character end_ins_code   = new Character(sref.getPdbx_seq_align_end_ins_code().charAt(0));

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
		Character db_begin_in_code = new Character(sref.getPdbx_db_align_beg_ins_code().charAt(0));
		Character db_end_in_code   = new Character(sref.getPdbx_db_align_end_ins_code().charAt(0));

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

	private static Chain getChainFromList(List<Chain> chains, String name){
		for (Chain chain : chains) {
			if ( chain.getChainID().equals(name)){

				return chain;
			}
		}
		// does not exist yet, so create...

		Chain	chain = new ChainImpl();
		chain.setChainID(name);
		chains.add(chain);

		return chain;
	}

	private Chain getEntityChain(String entity_id){

		return getChainFromList(entityChains,entity_id);
	}

	//private Chain getSeqResChain(String chainID){
	//	return getChainFromList(seqResChains, chainID);
	//}


	/** Data items in the ENTITY_SRC_GEN category record details of
               the source from which the entity was obtained in cases
               where the source was genetically manipulated.  The
               following are treated separately:  items pertaining to the tissue
               from which the gene was obtained, items pertaining to the host
               organism for gene expression and items pertaining to the actual
               producing organism (plasmid).
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

	/** The EntityPolySeq object provide the amino acid sequence objects for the Entities.
	 * Later on the entities are mapped to the BioJava Chain and Compound objects.
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


		// create group from epolseq;
		// by default this are the SEQRES records...


		if (epolseq.getMon_id().length()==3 && StructureTools.get1LetterCodeAmino(epolseq.getMon_id())!=null){
			AminoAcid g = new AminoAcidImpl();

			g.setRecordType(AminoAcid.SEQRESRECORD);

			g.setPDBName(epolseq.getMon_id());

			Character code1 = StructureTools.get1LetterCodeAmino(epolseq.getMon_id());
			g.setAminoType(code1);

			g.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
			// ARGH at this stage we don't know about insertion codes
			// this has to be obtained from _pdbx_poly_seq_scheme
			entityChain.addGroup(g);

		} else if ( StructureTools.isNucleotide(epolseq.getMon_id())) {
			// the group is actually a nucleotide group...
			NucleotideImpl n = new NucleotideImpl();

			n.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
			n.setPDBName(epolseq.getMon_id());
			entityChain.addGroup(n);				
		} else {				
			logger.debug("Residue {} {} is not a standard aminoacid or nucleotide, will create a het group for it", epolseq.getNum(),epolseq.getMon_id());
			HetatomImpl h = new HetatomImpl();				
			h.setPDBName(epolseq.getMon_id());
			h.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
			entityChain.addGroup(h);

		}

	}

	@Override
	public void newPdbxPolySeqScheme(PdbxPolySeqScheme ppss) {

		//if ( headerOnly)
		//	return;

		// replace the group asym ids with the real PDB ids!
		// replaceGroupSeqPos(ppss);  // This might be incorrect in some pdb, to use auth_seq_id of the pdbx_poly_seq_scheme.

		// merge the EntityPolySeq info and the AtomSite chains into one...
		//already known ignore:
		if (asymStrandId.containsKey(ppss.getAsym_id()))
			return;

		// this is one of the internal mmcif rules it seems...
		if ( ppss.getPdb_strand_id() == null) {
			asymStrandId.put(ppss.getAsym_id(), ppss.getAuth_mon_id());
			return;
		}

		//System.out.println(ppss.getAsym_id() + " = " + ppss.getPdb_strand_id());

		asymStrandId.put(ppss.getAsym_id(), ppss.getPdb_strand_id());

	}


	@Override
	public void newPdbxNonPolyScheme(PdbxNonPolyScheme ppss) {

		//if (headerOnly)
		//	return;

		// merge the EntityPolySeq info and the AtomSite chains into one...
		//already known ignore:
		if (asymStrandId.containsKey(ppss.getAsym_id()))
			return;

		// this is one of the interal mmcif rules it seems...
		if ( ppss.getPdb_strand_id() == null) {
			asymStrandId.put(ppss.getAsym_id(), ppss.getAsym_id());
			return;
		}

		asymStrandId.put(ppss.getAsym_id(), ppss.getPdb_strand_id());

	}

	@Override
	public void newPdbxEntityNonPoly(PdbxEntityNonPoly pen){
		// TODO: do something with them...
		// not implemented yet...
		//System.out.println(pen.getEntity_id() + " " + pen.getName() + " " + pen.getComp_id());
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

	private void createSSBonds() {
		List<SSBond> bonds = structure.getSSBonds();
		if (bonds == null) bonds = new ArrayList<SSBond>();

		final String symop = "1_555"; // For now - accept bonds within origin asymmetric unit.

		// For SSBond equivalent, parse through the struct_conn records
		int internalId = 0;
		for (StructConn conn : structConn) {

			if (!conn.getConn_type_id().equals("disulf")) continue;


			String ptnr1_chainId = conn.getPtnr1_auth_asym_id();

			String insCode1 = "";
			if (!conn.getPdbx_ptnr1_PDB_ins_code().equals("?")) insCode1 = conn.getPdbx_ptnr1_PDB_ins_code();
			String insCode2 = "";
			if (!conn.getPdbx_ptnr2_PDB_ins_code().equals("?")) insCode2 = conn.getPdbx_ptnr2_PDB_ins_code();

			String ptnr1_seqId = conn.getPtnr1_auth_seq_id() + insCode1 ;
			String ptnr2_chainId = conn.getPtnr2_auth_asym_id();
			String ptnr2_seqId = conn.getPtnr2_auth_seq_id() + insCode2;
			// conn.getId() would equal disulf#.

			// if we can find both of these residues - 
			Group s1 = lookupResidue(ptnr1_chainId, ptnr1_seqId);
			Group s2 = lookupResidue(ptnr2_chainId, ptnr2_seqId);

			// and is SS - then we should create a new disulfide bond.
			if (null != s1 && null != s2) {
				if ("CYS".equals(s1.getPDBName()) && "CYS".equals(s2.getPDBName()) ) {

					// TODO: when issue 220 is implemented, add robust symmetry handling 
					// to allow disulfide bonds between symmetry-related molecules.
					if (!conn.getPtnr1_symmetry().equals(symop) || !conn.getPtnr2_symmetry().equals(symop) ) {
						logger.info("Skipping ss bond between groups {} and {} belonging to different symmetry partners, because it is not supported yet", s1.getResidueNumber(), s2.getResidueNumber());
						continue;
					}


					SSBondImpl bond = new SSBondImpl();

					bond.setSerNum(internalId++); // An internal label what bond # 
					bond.setChainID1(ptnr1_chainId);
					bond.setResnum1(ptnr1_seqId);
					if ("?".equals(conn.getPdbx_ptnr1_PDB_ins_code())) {
						conn.setPdbx_ptnr1_PDB_ins_code(null);
					}
					bond.setInsCode1(conn.getPdbx_ptnr1_PDB_ins_code());
					bond.setChainID2(ptnr2_chainId);
					bond.setResnum2(ptnr2_seqId);
					if ("?".equals(conn.getPdbx_ptnr2_PDB_ins_code())) {
						conn.setPdbx_ptnr2_PDB_ins_code(null);
					}
					bond.setInsCode2(conn.getPdbx_ptnr2_PDB_ins_code());

					if (bonds.contains(bond)) {
						if (s1.getAltLocs().isEmpty() && s1.getAltLocs().isEmpty()) {
							// no alt locs in either group: there's something weird with the file having to ss bonds for same pair
							logger.warn("SS bond between residues {} and {} is repeated in file and neither group has alt locations, won't add it to list of SS bonds. ", s1.getResidueNumber(), s2.getResidueNumber());
						} else {
							// there is alt locs in either group. This is normal: an ssbond per alt loc is in file, e.g. 3dvf
							logger.info("SS bond between residues {} and {} is repeated in file and one of the 2 groups has alt locations. Adding only 1 SS bond for it. ", s1.getResidueNumber(), s2.getResidueNumber());
						}
					} else {					
						bonds.add(bond);
					}
				}
			}
		}

		structure.setSSBonds(bonds);
	}

	/**
	 * Lookup a residue - not found exceptions are handled with logger warnings.
	 * @param chainId
	 * @param seqId
	 * @return Successful = Group, Failure = null
	 */
	private Group lookupResidue(String chainId, String seqId) {
		try {
			Chain chain = structure.getChainByPDB(chainId);
			if (null != chain) {
				try {
					ResidueNumber resNum = null;
					if (Character.isAlphabetic( seqId.charAt(seqId.length()-1 )) ) {
						resNum = new ResidueNumber(chainId, Integer.parseInt(seqId.substring(0,seqId.length()-1)), seqId.charAt(seqId.length()-1));
					} else {
						resNum = new ResidueNumber(chainId, Integer.parseInt(seqId), ' ');
					}

					return chain.getGroupByPDB(resNum);
				} catch (NumberFormatException e) {
					logger.warn("Could not parse number for residue number {} specified in _struct_conn record", chainId + "-" + seqId);
				}
			}
		} catch (StructureException e) {
			logger.warn("Problem finding residue " + chainId + "-" + seqId + " specified in _struct_conn record." );
		}
		// Could not find.
		return null;
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
			String chain_id = siteGen.getAuth_asym_id(); // ChainID
			String auth_seq_id = siteGen.getAuth_seq_id(); // Res num

			String insCode = siteGen.getPdbx_auth_ins_code();
			if ( insCode != null && insCode.equals("?"))
				insCode = null;

			// Look for asymID = chainID and seqID = seq_ID.  Check that comp_id matches the resname.
			Group g = null;
			try {
				Chain chain = structure.getChainByPDB(chain_id);
				if (null != chain) {
					try {
						Character insChar = null;
						if (null != insCode && insCode.length() > 0) insChar = insCode.charAt(0);
						g = chain.getGroupByPDB(new ResidueNumber(chain_id, Integer.parseInt(auth_seq_id), insChar));
					} catch (NumberFormatException e) {
						logger.warn("Could not lookup residue : " + chain_id + auth_seq_id);
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
					logger.warn("comp_id doesn't match the residue at " + chain_id + auth_seq_id + " - skipping");
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
