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
import java.util.logging.Logger;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.DBRef;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.NucleotideImpl;
import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.UnknownPdbAminoAcidException;
import org.biojava.bio.structure.io.FileParsingParameters;

import org.biojava.bio.structure.io.PDBParseException;
import org.biojava.bio.structure.io.SeqRes2AtomAligner;
import org.biojava.bio.structure.io.mmcif.model.AtomSite;
import org.biojava.bio.structure.io.mmcif.model.AuditAuthor;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;

import org.biojava.bio.structure.io.mmcif.model.ChemCompAtom;
import org.biojava.bio.structure.io.mmcif.model.ChemCompBond;
import org.biojava.bio.structure.io.mmcif.model.ChemCompDescriptor;
import org.biojava.bio.structure.io.mmcif.model.DatabasePDBremark;
import org.biojava.bio.structure.io.mmcif.model.DatabasePDBrev;
import org.biojava.bio.structure.io.mmcif.model.Entity;
import org.biojava.bio.structure.io.mmcif.model.EntityPolySeq;
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
import org.biojava.bio.structure.io.mmcif.model.StructKeywords;
import org.biojava.bio.structure.io.mmcif.model.StructRef;
import org.biojava.bio.structure.io.mmcif.model.StructRefSeq;
import org.biojava.bio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;

/** A MMcifConsumer implementation that build a in-memory representation of the
 * content of a mmcif file as a BioJava Structure object.
 *  @author Andreas Prlic
 *  @since 1.7
 */

public class SimpleMMcifConsumer implements MMcifConsumer {

	boolean DEBUG = false;

	Structure structure;
	Chain current_chain;
	Group current_group;
	int atomCount;


	List<Chain>      current_model;
	List<Entity>     entities;
	List<StructRef>  strucRefs;
	List<Chain>      seqResChains;
	List<Chain>      entityChains; // needed to link entities, chains and compounds...
	List<StructAsym> structAsyms;  // needed to link entities, chains and compounds...
	List<PdbxStructOperList> structOpers ; //
	List<PdbxStructAssembly> strucAssemblies;
	List<PdbxStructAssemblyGen> strucAssemblyGens;


	Map<String,String> asymStrandId;

	String current_nmr_model ;

	FileParsingParameters params;

	public static Logger logger =  Logger.getLogger("org.biojava.bio.structure");

	public  SimpleMMcifConsumer(){
		params = new FileParsingParameters();
		documentStart();
	}

	public void newEntity(Entity entity) {
		if (DEBUG)
			System.out.println(entity);
		entities.add(entity);
	}

	public void newPdbxStructOperList(PdbxStructOperList structOper){

		structOpers.add(structOper);
	}

	public void newStructAsym(StructAsym sasym){

		structAsyms.add(sasym);
	}

	private Entity getEntity(String entity_id){
		for (Entity e: entities){
			if  (e.getId().equals(entity_id)){
				return e;
			}
		}
		return null;
	}

	@SuppressWarnings("deprecation")
	public void newStructKeywords(StructKeywords kw){
		PDBHeader header = structure.getPDBHeader();
		if ( header == null)
			header = new PDBHeader();
		header.setDescription(kw.getPdbx_keywords());
		header.setClassification(kw.getPdbx_keywords());
		Map<String, Object> h = structure.getHeader();
		h.put("classification", kw.getPdbx_keywords());      
	}

	@SuppressWarnings("deprecation")
	public void setStruct(Struct struct) {
		//System.out.println(struct);

		PDBHeader header = structure.getPDBHeader();
		if ( header == null)
			header = new PDBHeader();

		header.setTitle(struct.getTitle());
		header.setIdCode(struct.getEntry_id());
		//header.setDescription(struct.getPdbx_descriptor());
		//header.setClassification(struct.getPdbx_descriptor());
		//header.setDescription(struct.getPdbx_descriptor());
		//System.out.println(struct.getPdbx_model_details());

		Map<String, Object> h = structure.getHeader();
		h.put("title", struct.getTitle());
		//h.put("classification", struct.getPdbx_descriptor());


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
			HetatomImpl h = new HetatomImpl();
			h.setId(seq_id);
			group = h;
		}		
		//System.out.println("new group type: "+ group.getType() );
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


	/** during mmcif parsing the full atom name string gets truncated, fix this...
	 *
	 * @param name
	 * @return
	 */
	private String fixFullAtomName(String name){

		if (name.equals("N")){
			return " N  ";
		}
		if (name.equals("CA")){
			return " CA ";
		}
		if (name.equals("C")){
			return " C  ";
		}
		if (name.equals("O")){
			return " O  ";
		}
		if (name.equals("CB")){
			return " CB ";
		}
		if (name.equals("CG"))
			return " CG ";

		if (name.length() == 2)
			return " " + name + " ";

		if (name.length() == 1)
			return " " + name + "  ";

		if (name.length() == 3)
			return " " + name ;

		return name;
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
		String fullname      = fixFullAtomName(atom.getLabel_atom_id());		
		String recordName    = atom.getGroup_PDB();
		String residueNumberS = atom.getAuth_seq_id();
		Integer residueNrInt = Integer.parseInt(residueNumberS);
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
			}

			// we came to the beginning of a new NMR model
			structure.setNmr(true);
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

			//current_group.setPDBCode(residueNumber);

			current_group.setResidueNumber(residueNumber);
			try { 
				current_group.setPDBName(groupCode3);
			} catch (PDBParseException e){
				System.err.println(e.getMessage());
			}
		}

		if ( startOfNewChain){
			current_group = getNewGroup(recordName,aminoCode1,seq_id, groupCode3);

			// current_group.setPDBCode(residueNumber);
			current_group.setResidueNumber(residueNumber);
			try {
				current_group.setPDBName(groupCode3);
			} catch (PDBParseException e){
				e.printStackTrace();
			}
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

			current_group = getNewGroup(recordName,aminoCode1,seq_id,groupCode3);
			//current_group.setPDBCode(pdbCode);
			try {
				current_group.setPDBName(groupCode3);
			} catch (PDBParseException e){
				e.printStackTrace();
			}
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

		atomCount++;
		//System.out.println("fixing atom name for  >" + atom.getLabel_atom_id() + "< >" + fullname + "<");

		if ( params.isParseCAOnly() ){
			// yes , user wants to get CA only
			// only parse CA atoms...
			if (! fullname.equals(" CA ")){
				//System.out.println("ignoring " + line);
				atomCount--;
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

		a.setFullName(fixFullAtomName(atom.getLabel_atom_id()));


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
			if ( current_group.getAtoms().size() == 0)
				return current_group;
			//System.out.println("cloning current group");
			Group altLocG = (Group) current_group.clone();
			current_group.addAltLoc(altLocG);
			return altLocG;	
		}


		Group altLocG = getNewGroup(recordName,aminoCode1,seq_id,groupCode3);

		try {
			altLocG.setPDBName(groupCode3);
		} catch (PDBParseException e) {
			e.printStackTrace();
		}
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
		atomCount     		= 0;

		current_model = new ArrayList<Chain>();
		entities      = new ArrayList<Entity>();
		strucRefs     = new ArrayList<StructRef>();
		seqResChains  = new ArrayList<Chain>();
		entityChains  = new ArrayList<Chain>();
		structAsyms   = new ArrayList<StructAsym>();
		asymStrandId  = new HashMap<String, String>();
		structOpers   = new ArrayList<PdbxStructOperList>();
		strucAssemblies = new ArrayList<PdbxStructAssembly>();
		strucAssemblyGens = new ArrayList<PdbxStructAssemblyGen>();
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
			if ( DEBUG){
				System.err.println("current chain is null at end of document.");
			}
		}

		structure.addModel(current_model);

		// Goal is to reproduce the PDB files exactly:
		// What has to be done is to use the auth_mon_id for the assignment. For this

		// map entities to Chains and Compound objects...


		for (StructAsym asym : structAsyms) {
			if ( DEBUG )
				System.out.println("entity " + asym.getEntity_id() + " matches asym id:" + asym.getId() );

			Chain s = getEntityChain(asym.getEntity_id());
			Chain seqres = (Chain)s.clone();
			seqres.setChainID(asym.getId());


			seqResChains.add(seqres);
			if ( DEBUG )
				System.out.println(" seqres: " + asym.getId() + " " + seqres + "<") ;

		}


		if ( params.isAlignSeqRes() ){

			SeqRes2AtomAligner aligner = new SeqRes2AtomAligner();			
			aligner.align(structure,seqResChains);
		}


		//TODO: add support for these:

		//structure.setConnections(connects);
		//structure.setCompounds(compounds);
		//linkChains2Compound(structure);

		// mismatching Author assigned chain IDS and PDB internal chain ids:
		// fix the chain IDS in the current model:

		Set<String> asymIds = asymStrandId.keySet();

		for (int i =0; i< structure.nrModels() ; i++){
			List<Chain>model = structure.getModel(i);

			List<Chain> pdbChains = new ArrayList<Chain>();

			for (Chain chain : model) {
				for (String asym : asymIds) {
					if ( chain.getChainID().equals(asym)){

						if (DEBUG)
							System.out.println("renaming " + asym  + " to : " + asymStrandId.get(asym));

						chain.setChainID(asymStrandId.get(asym));
						chain.setInternalChainID(asym);
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


		// set the oligomeric state info in the header...

		PDBHeader header = structure.getPDBHeader();
		header.setNrBioAssemblies(strucAssemblies.size());

		// the more detailed mapping of chains to rotation operations happens in StructureIO...
		// TODO clean this up and move it here...
		//header.setBioUnitTranformationMap(tranformationMap);
		Map<Integer,List<ModelTransformationMatrix>> transformationMap = new HashMap<Integer, List<ModelTransformationMatrix>>();
		int total = strucAssemblies.size();

		for ( int defaultBioAssembly = 1 ; defaultBioAssembly <= total; defaultBioAssembly++){
			//List<ModelTransformationMatrix>tmp = getBioUnitTransformationList(pdbId, i +1);

			PdbxStructAssembly psa = strucAssemblies.get(defaultBioAssembly-1);
			List<PdbxStructAssemblyGen> psags = new ArrayList<PdbxStructAssemblyGen>(1);

			for ( PdbxStructAssemblyGen psag: strucAssemblyGens ) {
				if ( psag.getAssembly_id().equals(defaultBioAssembly+"")) {
					psags.add(psag);
				}
			}

			//System.out.println("psags: " + psags.size());
			BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

			// these are the transformations that need to be applied to our model
			List<ModelTransformationMatrix> transformations = builder.getBioUnitTransformationList(psa, psags, structOpers);

			transformationMap.put(defaultBioAssembly,transformations);
			//System.out.println("mmcif header: " + (defaultBioAssembly+1) + " " + transformations.size() +" " +  transformations);

		}
		structure.getPDBHeader().setBioUnitTranformationMap(transformationMap);


	}


	/** This method will return the parsed protein structure, once the parsing has been finished
	 *
	 * @return a BioJava protein structure object
	 */
	public Structure getStructure() {

		return structure;
	}

	@SuppressWarnings("deprecation")
	public void newDatabasePDBrev(DatabasePDBrev dbrev) {
		//System.out.println("got a database revision:" + dbrev);
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
		PDBHeader header = structure.getPDBHeader();
		Map<String, Object> h = structure.getHeader();

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
				h.put("depDate", date);
				Date mod = dateFormat.parse(dbrev.getDate());

				header.setModDate(mod);
				h.put("revDate",dbrev.getDate());

			} catch (ParseException e){
				e.printStackTrace();
			}
		} else {
			try {

				Date mod = dateFormat.parse(dbrev.getDate());
				header.setModDate(mod);
				h.put("revDate",dbrev.getDate());

			} catch (ParseException e){
				e.printStackTrace();
			}
		}

		structure.setPDBHeader(header);
	}

	@SuppressWarnings("deprecation")
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
					System.err.println(e.getMessage());
					System.err.println("could not parse resolution from line and ignoring it " + line);
					return ;


				}
				// support for old style header

				Map<String,Object> header = structure.getHeader();
				header.put("resolution",new Float(res));
				structure.setHeader(header);

				PDBHeader pdbHeader = structure.getPDBHeader();
				pdbHeader.setResolution(res);

			}

		}
	}

	@SuppressWarnings("deprecation")
	public void newRefine(Refine r){
		// copy the resolution to header
		PDBHeader pdbHeader = structure.getPDBHeader();
		try {
			pdbHeader.setResolution(Float.parseFloat(r.getLs_d_res_high()));
		} catch (NumberFormatException e){
			logger.warning("could not parse resolution from " + r.getLs_d_res_high() + " " + e.getMessage());
		}
		Map<String,Object> header = structure.getHeader();
		header.put("resolution",pdbHeader.getResolution());
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

	@SuppressWarnings("deprecation")
	public void newExptl(Exptl exptl) {

		PDBHeader pdbHeader = structure.getPDBHeader();
		String method = exptl.getMethod();
		String old = pdbHeader.getTechnique();
		if ( (old != null) && (! old.equals(""))){
			method = old+"; " + method;
		}
		pdbHeader.setTechnique(method);
		Map<String,Object> header = structure.getHeader();
		header.put("technique",method);

	}

	public void newStructRef(StructRef sref) {
		if (DEBUG)
			System.out.println(sref);
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


		//TODO: make DBRef chain IDs a string for chainIDs that are longer than one char...
		r.setChainId(new Character(sref.getPdbx_strand_id().charAt(0)));
		StructRef structRef = getStructRef(sref.getRef_id());
		if (structRef == null){
			logger.warning("could not find StructRef " + sref.getRef_id() + " for StructRefSeq " + sref);
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

		if ( DEBUG)
			System.out.println(r.toPDB());

		structure.setDBRefs(dbrefs);

	}

	private Chain getChainFromList(List<Chain> chains, String name){
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

	/** The EntityPolySeq object provide the amino acid sequence objects for the Entities.
	 * Later on the entities are mapped to the BioJava Chain and Compound objects.
	 * @param epolseq the EntityPolySeq record for one amino acid
	 */
	public void newEntityPolySeq(EntityPolySeq epolseq) {

		if (DEBUG)
			System.out.println("NEW entity poly seq " + epolseq);

		Entity e = getEntity(epolseq.getEntity_id());

		if (e == null){
			System.err.println("could not find entity "+ epolseq.getEntity_id()+". Can not match sequence to it.");
			return;
		}

		Chain entityChain = getEntityChain(epolseq.getEntity_id());


		// create group from epolseq;
		// by default this are the SEQRES records...

		AminoAcid g = new AminoAcidImpl();

		g.setRecordType(AminoAcid.SEQRESRECORD);

		try {
			g.setPDBName(epolseq.getMon_id());

			Character code1 = StructureTools.convert_3code_1code(epolseq.getMon_id());
			g.setAminoType(code1);
			g.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
			// ARGH at this stage we don;t know about insertion codes
			// this has to be obtained from _pdbx_poly_seq_scheme
			entityChain.addGroup(g);

		}  catch (PDBParseException ex) {
			if ( StructureTools.isNucleotide(epolseq.getMon_id())) {
				// the group is actually a nucleotide group...
				NucleotideImpl n = new NucleotideImpl();
				n.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
				try {
					n.setPDBName(epolseq.getMon_id());
				} catch (PDBParseException pex){
					pex.printStackTrace();
				}
				entityChain.addGroup(n);
			}
			else {
				logger.warning(ex.getMessage() + " creating a hetatom called XXX ");
				HetatomImpl h = new HetatomImpl();
				try {
					h.setPDBName(epolseq.getMon_id());
					//h.setAminoType('X');
					h.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
					entityChain.addGroup(h);

				} catch (PDBParseException exc) {
					System.err.println("this is a helpless case and I am dropping group " + epolseq.getMon_id());
				}
				//ex.printStackTrace();
			}
		} catch (UnknownPdbAminoAcidException ex){
			//logger.warning("no sure what to do with:" + epolseq.getMon_id()+ " " + ex.getMessage());
			HetatomImpl h = new HetatomImpl();
			try {
				h.setPDBName(epolseq.getMon_id());
				//h.setAminoType('X');
				h.setResidueNumber(ResidueNumber.fromString(epolseq.getNum()));
				entityChain.addGroup(h);

			} catch (PDBParseException exc) {
				System.err.println("this is a helpless case and I am dropping group " + epolseq.getMon_id() + " " + ex.getMessage());
			}
			//System.err.println(ex.getMessage());
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

			if (! target.getPDBName().equals(ppss.getMon_id())){
				logger.info("could not match PdbxPolySeqScheme to chain:" + ppss);
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
	}

	public void newChemComp(ChemComp c) {
		// TODO: do something with them...

	}

	public void newGenericData(String category, List<String> loopFields,
			List<String> lineData) {

		if (DEBUG) {
			System.err.println("unhandled category so far: " + category);
		}
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

		// todo: nothing happening here yet.

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





}


