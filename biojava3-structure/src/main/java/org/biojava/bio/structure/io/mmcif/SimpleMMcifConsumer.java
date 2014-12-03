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
package org.biojava.bio.structure.io.mmcif;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Matrix4d;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Compound;
import org.biojava.bio.structure.DBRef;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.NucleotideImpl;
import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.UnknownPdbAminoAcidException;
import org.biojava.bio.structure.io.BondMaker;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.SeqRes2AtomAligner;
import org.biojava.bio.structure.io.mmcif.model.AtomSite;
import org.biojava.bio.structure.io.mmcif.model.AuditAuthor;
import org.biojava.bio.structure.io.mmcif.model.Cell;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava.bio.structure.io.mmcif.model.ChemCompAtom;
import org.biojava.bio.structure.io.mmcif.model.ChemCompBond;
import org.biojava.bio.structure.io.mmcif.model.ChemCompDescriptor;
import org.biojava.bio.structure.io.mmcif.model.DatabasePDBremark;
import org.biojava.bio.structure.io.mmcif.model.DatabasePDBrev;
import org.biojava.bio.structure.io.mmcif.model.Entity;
import org.biojava.bio.structure.io.mmcif.model.EntityPolySeq;
import org.biojava.bio.structure.io.mmcif.model.EntitySrcGen;
import org.biojava.bio.structure.io.mmcif.model.EntitySrcNat;
import org.biojava.bio.structure.io.mmcif.model.EntitySrcSyn;
import org.biojava.bio.structure.io.mmcif.model.Exptl;
import org.biojava.bio.structure.io.mmcif.model.PdbxChemCompDescriptor;
import org.biojava.bio.structure.io.mmcif.model.PdbxChemCompIdentifier;
import org.biojava.bio.structure.io.mmcif.model.PdbxEntityNonPoly;
import org.biojava.bio.structure.io.mmcif.model.PdbxNonPolyScheme;
import org.biojava.bio.structure.io.mmcif.model.PdbxPolySeqScheme;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.bio.structure.io.mmcif.model.Refine;
import org.biojava.bio.structure.io.mmcif.model.Struct;
import org.biojava.bio.structure.io.mmcif.model.StructAsym;
import org.biojava.bio.structure.io.mmcif.model.StructConn;
import org.biojava.bio.structure.io.mmcif.model.StructKeywords;
import org.biojava.bio.structure.io.mmcif.model.StructNcsOper;
import org.biojava.bio.structure.io.mmcif.model.StructRef;
import org.biojava.bio.structure.io.mmcif.model.StructRefSeq;
import org.biojava.bio.structure.io.mmcif.model.Symmetry;
import org.biojava.bio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.bio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.bio.structure.xtal.CrystalCell;
import org.biojava.bio.structure.xtal.SpaceGroup;
import org.biojava.bio.structure.xtal.SymoplibParser;
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

	private Map<String,String> asymStrandId;
	
	private Map<String,String> asymId2entityId;

	private String current_nmr_model ;

	private FileParsingParameters params;

	public  SimpleMMcifConsumer(){
		params = new FileParsingParameters();
		documentStart();

	}

	public void newEntity(Entity entity) {
		logger.debug("New entity: {}",entity.toString());
		entities.add(entity);
	}

	public void newPdbxStructOperList(PdbxStructOperList structOper){

		structOpers.add(structOper);
	}

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

	public void newStructKeywords(StructKeywords kw){
		PDBHeader header = structure.getPDBHeader();
		if ( header == null)
			header = new PDBHeader();
		header.setDescription(kw.getPdbx_keywords());
		header.setClassification(kw.getPdbx_keywords());
	}

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
			if (aminoCode1 == null)  {
				// it is a nucleotide
				NucleotideImpl nu = new NucleotideImpl();
				group = nu;
				nu.setId(seq_id);

			} else if (aminoCode1 == StructureTools.UNKNOWN_GROUP_LABEL){
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
			if (aminoCode1 != null ) {
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

	public void newAtomSite(AtomSite atom) {

		// Warning: getLabel_asym_id is not the "chain id" in the PDB file
		// it is the internally used chain id.
		// later on we will fix this...

		// later one needs to map the asym id to the pdb_strand_id

		//TODO: add support for MAX_ATOMS

		boolean startOfNewChain = false;

		//String chain_id      = atom.getAuth_asym_id();
		String chain_id      = atom.getLabel_asym_id();		
				
		String recordName    = atom.getGroup_PDB();
		String residueNumberS = atom.getAuth_seq_id();
		Integer residueNrInt = Integer.parseInt(residueNumberS);

		//String residueNumberSeqres = atom.getLabel_seq_id();
		// the 3-letter name of the group:
		String groupCode3    = atom.getLabel_comp_id();
		if ( groupCode3.length() == 1){
			groupCode3 = "  " + groupCode3;
		}
		if ( groupCode3.length() == 2){
			groupCode3 = " " + groupCode3;
		}
		Character aminoCode1 = null;
		if ( recordName.equals("ATOM") )
			aminoCode1 = StructureTools.get1LetterCode(groupCode3);
		else {
			aminoCode1 = StructureTools.get1LetterCode(groupCode3);
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
				altGroup = getCorrectAltLocGroup( altLoc,recordName,aminoCode1,groupCode3, seq_id);
				//System.out.println("found altLoc! " + altLoc + " " + current_group + " " + altGroup);
			}
		}

		if ( params.isHeaderOnly())
			return;

		//atomCount++;
		//System.out.println("fixing atom name for  >" + atom.getLabel_atom_id() + "< >" + fullname + "<");

		
		if ( params.isParseCAOnly() ){
			// yes , user wants to get CA only
			// only parse CA atoms...
			if (! (atom.getLabel_atom_id().equals("CA") && atom.getType_symbol().equals("C"))) {
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

		double occupancy = Double.parseDouble(atom.getOccupancy());
		a.setOccupancy(occupancy);

		double temp = Double.parseDouble(atom.getB_iso_or_equiv());
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
		}  catch (IllegalArgumentException e){}
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
		asymId2entityId = new HashMap<String,String>();
		structOpers   = new ArrayList<PdbxStructOperList>();
		strucAssemblies = new ArrayList<PdbxStructAssembly>();
		strucAssemblyGens = new ArrayList<PdbxStructAssemblyGen>();
		entitySrcGens = new ArrayList<EntitySrcGen>();
		entitySrcNats = new ArrayList<EntitySrcNat>();
		entitySrcSyns = new ArrayList<EntitySrcSyn>();
		structConn = new ArrayList<StructConn>();
		structNcsOper = new ArrayList<StructNcsOper>();
	}


	public void documentEnd() {

		// a problem occurred earlier so current_chain = null ...
		// most likely the buffered reader did not provide data ...
		if ( current_chain != null ) {

			current_chain.addGroup(current_group);
			if (isKnownChain(current_chain.getChainID(),current_model) == null) {
				current_model.add(current_chain);
			}
		} else {
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

			
			for (EntitySrcGen esg : entitySrcGens) {
				String eId = esg.getEntity_id();
				if (! eId.equals(asym.getEntity_id()))
					continue;

				int eIdInt = 0;
				try {
					eIdInt = Integer.parseInt(eId);
				} catch (NumberFormatException e) {
					logger.warn("Could not parse mol_id from string {}. Will use 0 for matching EntitySrcGen",eId);
				}
				
				// found the matching EntitySrcGen
				// get the corresponding Entity
				Compound c = structure.getCompoundById(eIdInt);
				if ( c == null){					
					c = createNewCompoundFromESG(esg, eIdInt);
					// add to chain
					List<Compound> compounds  = structure.getCompounds();
					compounds.add(c);
					logger.debug("Adding Compound with entity id {} from _entity_src_syn, with name: {}",eIdInt,c.getMolName());
				}

			}

			for (EntitySrcNat esn : entitySrcNats) {
				String eId = esn.getEntity_id();
				if (! eId.equals(asym.getEntity_id()))
					continue;

				int eIdInt = 0;
				try {
					eIdInt = Integer.parseInt(eId);
				} catch (NumberFormatException e) {
					logger.warn("Could not parse mol_id from string {}. Will use 0 for matching EntitySrcGen",eId);
				}
				// found the matching EntitySrcGen
				// get the corresponding Entity
				Compound c = structure.getCompoundById(eIdInt);
				if ( c == null){					
					c = createNewCompoundFromESN(esn, eIdInt);
					// add to chain
					List<Compound> compounds  = structure.getCompounds();
					compounds.add(c);
					logger.debug("Adding Compound with entity id {} from _entity_src_syn, with name: {}",eIdInt,c.getMolName());
				}

			}

			for (EntitySrcSyn ess : entitySrcSyns) {
				String eId = ess.getEntity_id();
				if (! eId.equals(asym.getEntity_id()))
					continue;

				int eIdInt = 0;
				try {
					eIdInt = Integer.parseInt(eId);
				} catch (NumberFormatException e) {
					logger.warn("Could not parse mol_id from string {}. Will use 0 for matching EntitySrcGen",eId);
				}
				// found the matching EntitySrcGen
				// get the corresponding Entity
				Compound c = structure.getCompoundById(eIdInt);
				if ( c == null){					
					c = createNewCompoundFromESS(ess, eIdInt);
					// add to chain
					List<Compound> compounds  = structure.getCompounds();
					compounds.add(c);
					logger.debug("Adding Compound with entity id {} from _entity_src_syn, with name: {}",eIdInt,c.getMolName());
				}
			}
			
			// for some mmCIF files like 1yrm all 3 of _entity_src_gen, _entity_src_nat and _pdbx_entity_src_syn are missing
			// we need to fill the Compounds in some other way:

			int eId = 0;
			try {
				eId = Integer.parseInt(asym.getEntity_id());
			} catch (NumberFormatException e) {
				logger.warn("Could not parse mol_id from string {}. Will use 0 for creating Compound",asym.getEntity_id());
			}

			Compound c = structure.getCompoundById(eId);

			if (c==null) {
				c = new Compound();
				c.setMolId(eId);
				Entity e = getEntity(eId);
				// we only add the compound if a polymeric one (to match what the PDB parser does)
				if (e!=null && e.getType().equals("polymer")) {
					if ( e != null)
						c.setMolName(e.getPdbx_description());
					List<Compound> compounds = structure.getCompounds();
					compounds.add(c);
					logger.debug("Adding Compound with entity id {} from _entity, with name: {}",eId, c.getMolName());
				}
			}

		}

		if ( params.isAlignSeqRes() ){
		
			SeqRes2AtomAligner aligner = new SeqRes2AtomAligner();			
			//aligner.align(structure,seqResChains);

			// fix SEQRES residue numbering
			List<Chain> atomList   = structure.getModel(0);
			for (Chain seqResChain: seqResChains){
		
					Chain atomChain = aligner.getMatchingAtomRes(seqResChain, atomList);

					//map the atoms to the seqres...

					List<Group> seqResGroups = seqResChain.getAtomGroups(); 

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


		addBonds();

		//TODO: add support for these:
		//structure.setConnections(connects);
		

		// mismatching Author assigned chain IDS and PDB internal chain ids:
		// fix the chain IDS in the current model:

		Set<String> asymIds = asymStrandId.keySet();

		for (int i =0; i< structure.nrModels() ; i++){
			List<Chain>model = structure.getModel(i);

			List<Chain> pdbChains = new ArrayList<Chain>();

			for (Chain chain : model) {
				for (String asym : asymIds) {
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
			
			// finally setting chains to compounds and compounds to chains now that we have the final chains
			for (Chain chain:pdbChains) {
				String entityId = asymId2entityId.get(chain.getInternalChainID());
				int eId = Integer.parseInt(entityId);
				Compound compound = structure.getCompoundById(eId);
				logger.debug("Adding chain with chain id {} (internal chain id {}) to compound with entity_id {}",
						chain.getChainID(), chain.getInternalChainID(), eId);
				compound.addChain(chain);
				chain.setCompound(compound);

			}			
			
		}
		

		// set the oligomeric state info in the header...

		PDBHeader header = structure.getPDBHeader();
		header.setNrBioAssemblies(strucAssemblies.size());

		// the more detailed mapping of chains to rotation operations happens in StructureIO...
		// TODO clean this up and move it here...
		//header.setBioUnitTranformationMap(tranformationMap);
		Map<String,List<BiologicalAssemblyTransformation>> transformationMap = new HashMap<String, List<BiologicalAssemblyTransformation>>();
		//int total = strucAssemblies.size();

		//for ( int defaultBioAssembly = 1 ; defaultBioAssembly <= total; defaultBioAssembly++){
		
		for ( PdbxStructAssembly psa : strucAssemblies){
		//List<ModelTransformationMatrix>tmp = getBioUnitTransformationList(pdbId, i +1);

			//PdbxStructAssembly psa = strucAssemblies.get(asmbl.getId());
			List<PdbxStructAssemblyGen> psags = new ArrayList<PdbxStructAssemblyGen>(1);

			for ( PdbxStructAssemblyGen psag: strucAssemblyGens ) {
				if ( psag.getAssembly_id().equals(psa.getId())) {
					psags.add(psag);
				}
			}

			//System.out.println("psags: " + psags.size());
			BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

			// these are the transformations that need to be applied to our model
			List<BiologicalAssemblyTransformation> transformations = builder.getBioUnitTransformationList(psa, psags, structOpers);

			transformationMap.put(psa.getId(),transformations);
			//System.out.println("mmcif header: " + (defaultBioAssembly+1) + " " + transformations.size() +" " +  transformations);

		}
		structure.getPDBHeader().setBioUnitTranformationMap(transformationMap);

		ArrayList<Matrix4d> ncsOperators = new ArrayList<Matrix4d>();
		for (StructNcsOper sNcsOper:structNcsOper) {
			if (sNcsOper.getCode().equals("generate")) {
				ncsOperators.add(sNcsOper.getOperator());
			}
		}
		// we only set it if not empty, otherwise remains null
		if (ncsOperators.size()>0) {
			structure.getCrystallographicInfo().setNcsOperators(
					(Matrix4d[]) ncsOperators.toArray(new Matrix4d[ncsOperators.size()]));
		}
		
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

	private Compound createNewCompoundFromESG(EntitySrcGen esg, int eId) {

		Entity e = getEntity(eId);
		Compound c = new Compound();
		c.setMolId(eId);
		
		if ( e != null)
			c.setMolName(e.getPdbx_description());
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

		Entity e = getEntity(eId);
		Compound c = new Compound();
		
		c.setMolId(eId);
		
		if ( e != null)
			c.setMolName(e.getPdbx_description());
		c.setAtcc(esn.getPdbx_atcc());
		c.setCell(esn.getPdbx_cell());
		c.setOrganismCommon(esn.getCommon_name());
		c.setOrganismScientific(esn.getPdbx_organism_scientific());
		c.setOrganismTaxId(esn.getPdbx_ncbi_taxonomy_id());

		return c;

	}

	private Compound createNewCompoundFromESS(EntitySrcSyn ess, int eId) {

		Entity e = getEntity(eId);
		Compound c = new Compound();
		
		c.setMolId(eId);
		
		if ( e != null)
			c.setMolName(e.getPdbx_description());


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

	public void newDatabasePDBrev(DatabasePDBrev dbrev) {
		//System.out.println("got a database revision:" + dbrev);
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
		PDBHeader header = structure.getPDBHeader();

		if ( header == null) {
			header = new PDBHeader();
		}


		if (dbrev.getNum().equals("1")){

			try {

				String date = dbrev.getDate_original();
				//System.out.println(date);
				Date dep = dateFormat.parse(date);
				//System.out.println(dep);
				header.setDepDate(dep);
				Date mod = dateFormat.parse(dbrev.getDate());

				header.setModDate(mod);

			} catch (ParseException e){
				e.printStackTrace();
			}
		} else {
			try {

				Date mod = dateFormat.parse(dbrev.getDate());
				header.setModDate(mod);

			} catch (ParseException e){
				e.printStackTrace();
			}
		}

		structure.setPDBHeader(header);
	}

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
		try {
			pdbHeader.setRfree(Float.parseFloat(r.getLs_R_factor_R_free()));
		} catch (NumberFormatException e){
			// no rfree present ('?') is very usual, that's why we set it to debug
			logger.debug("Could not parse Rfree from string '{}'", r.getLs_R_factor_R_free());
		}

		
	}


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

	public void newExptl(Exptl exptl) {

		PDBHeader pdbHeader = structure.getPDBHeader();
		String method = exptl.getMethod();
		pdbHeader.setExperimentalTechnique(method);

	}
	
	public void newCell(Cell cell) {
		
		try {
			float a = Float.parseFloat(cell.getLength_a());
			float b = Float.parseFloat(cell.getLength_b());
			float c = Float.parseFloat(cell.getLength_c());
			float alpha = Float.parseFloat(cell.getAngle_alpha());
			float beta = Float.parseFloat(cell.getAngle_beta());
			float gamma = Float.parseFloat(cell.getAngle_gamma());
			// If the entry describes a structure determined by a technique other than X-ray crystallography,
		    // cell is (sometimes!) a = b = c = 1.0, alpha = beta = gamma = 90 degrees
			// if so we don't add and CrystalCell will be null
			if (a == 1.0f && b == 1.0f && c == 1.0f && 
	        		alpha == 90.0f && beta == 90.0f && gamma == 90.0f ) {
	        	return;
	        } 
		
			CrystalCell xtalCell = new CrystalCell(); 
			structure.getPDBHeader().getCrystallographicInfo().setCrystalCell(xtalCell);
			xtalCell.setA(a);
			xtalCell.setB(b);
			xtalCell.setC(c);
			xtalCell.setAlpha(alpha);
			xtalCell.setBeta(beta);
			xtalCell.setGamma(gamma);
			
			
			
		} catch (NumberFormatException e){
			structure.getPDBHeader().getCrystallographicInfo().setCrystalCell(null);
			logger.info("could not parse some cell parameters ("+e.getMessage()+"), ignoring _cell ");
		}
		try {
			// if Z parsing fails it is not so important
			structure.getPDBHeader().getCrystallographicInfo().setZ(Integer.parseInt(cell.getZ_PDB()));
		} catch (NumberFormatException e) {
			logger.info("could not parse some the Z parameter from _cell ");
		}
	}
	
	public void newSymmetry(Symmetry symmetry) {
        String spaceGroup = symmetry.getSpace_group_name_H_M();
		SpaceGroup sg = SymoplibParser.getSpaceGroup(spaceGroup);
        if (sg==null) logger.warn("Space group '"+spaceGroup+"' not recognised as a standard space group"); 

		structure.getPDBHeader().getCrystallographicInfo().setSpaceGroup(sg); 
	}

	public void newStructNcsOper(StructNcsOper sNcsOper) {
		structNcsOper.add(sNcsOper);
	}
	
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

		AminoAcid g = new AminoAcidImpl();

		g.setRecordType(AminoAcid.SEQRESRECORD);

		try {
			
			if (epolseq.getMon_id().length()==3){
				g.setPDBName(epolseq.getMon_id());
				
				Character code1 = StructureTools.convert_3code_1code(epolseq.getMon_id());
				g.setAminoType(code1);

				g.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
				// ARGH at this stage we don;t know about insertion codes
				// this has to be obtained from _pdbx_poly_seq_scheme
				entityChain.addGroup(g);

			} else if ( StructureTools.isNucleotide(epolseq.getMon_id())) {
				// the group is actually a nucleotide group...
				NucleotideImpl n = new NucleotideImpl();

				n.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
				n.setPDBName(epolseq.getMon_id());
				entityChain.addGroup(n);				
			} else {				
				HetatomImpl h = new HetatomImpl();
				
				h.setPDBName(epolseq.getMon_id());
				//h.setAminoType('X');
				h.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
				entityChain.addGroup(h);

			}
		} catch (UnknownPdbAminoAcidException ex){
			logger.debug("Residue {} {} is not a standard aminoacid, will create a het group for it", epolseq.getNum(),epolseq.getMon_id());
			HetatomImpl h = new HetatomImpl();
			h.setPDBName(epolseq.getMon_id());
			//h.setAminoType('X');
			h.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
			entityChain.addGroup(h);

		}
	}


	/* returns the chains from all models that have the provided chainId
	 *
	 */
	private List<Chain> getChainsFromAllModels(String chainId){
		List<Chain> chains = new ArrayList<Chain>();


		for (int i=0 ; i < structure.nrModels();i++){
			List<Chain> model = structure.getModel(i);
			for (Chain c: model){
				if (c.getChainID().equals(chainId)) {
					chains.add(c);
				}
			}
		}

		return chains;
	}

	/** finds the residue in the internal representation and fixes the residue number and insertion code
	 *
	 * @param ppss
	 */
	private void replaceGroupSeqPos(PdbxPolySeqScheme ppss){

		if (ppss.getAuth_seq_num().equals("?"))
			return;

		//logger.info("replacegroupSeqPos " + ppss);
		// at this stage we are still using the internal asym ids...
		List<Chain> matchinChains = getChainsFromAllModels(ppss.getAsym_id());

		long sid = Long.parseLong(ppss.getSeq_id());
		for (Chain c: matchinChains){
			Group target = null;
			for (Group g: c.getAtomGroups()){

				if ( g instanceof AminoAcidImpl){
					AminoAcidImpl aa = (AminoAcidImpl)g;
					if (aa.getId() == sid ) {
						// found it:
						target = g;
						break;
					}
				}
				else if ( g instanceof NucleotideImpl) {
					NucleotideImpl n = (NucleotideImpl)g;
					if ( n.getId() == sid) {
						target = g;
						break;
					}
				} else if ( g instanceof HetatomImpl){
					HetatomImpl h = (HetatomImpl)g;
					if ( h.getId() == sid){
						target =h;
						break;
					}
				}
			}
			if (target == null){
				logger.info("could not find group at seq. position " + 
						ppss.getSeq_id() + " in internal chain " + c.getChainID() + ". " + ppss);
				continue;
			}

			if (! target.getPDBName().trim().equals(ppss.getMon_id())){
				logger.info("could not match PdbxPolySeqScheme to chain:" + target.getPDBName() + " " + ppss);
				continue;
			}

			// fix the residue number to the one used in the PDB files...
			Integer pdbResNum = Integer.parseInt(ppss.getAuth_seq_num());
			// check the insertion code...
			String insCodeS = ppss.getPdb_ins_code();
			Character insCode = null;
			if ( ( insCodeS != null) && (! insCodeS.equals(".")) && insCodeS.length()>0){
				//pdbResNum += insCode
				insCode = insCodeS.charAt(0);
			}

			ResidueNumber residueNumber = new ResidueNumber(null, pdbResNum,insCode);
			//logger.info("setting residue number for " + target +" to: " + residueNumber);
			target.setResidueNumber(residueNumber);
		}
	}
	public void newPdbxPolySeqScheme(PdbxPolySeqScheme ppss) {

		//if ( headerOnly)
		//	return;

		// replace the group asym ids with the real PDB ids!
		replaceGroupSeqPos(ppss);

		// merge the EntityPolySeq info and the AtomSite chains into one...
		//already known ignore:
		if (asymStrandId.containsKey(ppss.getAsym_id()))
			return;

		// this is one of the interal mmcif rules it seems...
		if ( ppss.getPdb_strand_id() == null) {
			asymStrandId.put(ppss.getAsym_id(), ppss.getAuth_mon_id());
			return;
		}

		//System.out.println(ppss.getAsym_id() + " = " + ppss.getPdb_strand_id());

		asymStrandId.put(ppss.getAsym_id(), ppss.getPdb_strand_id());

	}


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

	public void newPdbxEntityNonPoly(PdbxEntityNonPoly pen){
		// TODO: do something with them...
		// not implemented yet...
		//System.out.println(pen.getEntity_id() + " " + pen.getName() + " " + pen.getComp_id());
	}

	public void newChemComp(ChemComp c) {
		// TODO: do something with them...

	}

	public void newGenericData(String category, List<String> loopFields,
			List<String> lineData) {

		//logger.debug("unhandled category so far: " + category);		
	}

	public FileParsingParameters getFileParsingParameters()
	{
		return params;
	}

	public void setFileParsingParameters(FileParsingParameters params)
	{
		this.params = params;

	}

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
		// TODO Auto-generated method stub

	}

	@Override
	public void newPdbxChemCompIndentifier(PdbxChemCompIdentifier id) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newChemCompBond(ChemCompBond bond) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newPdbxChemCompDescriptor(PdbxChemCompDescriptor desc) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newStructConn(StructConn structConn) {
		this.structConn.add(structConn);
	}

}


