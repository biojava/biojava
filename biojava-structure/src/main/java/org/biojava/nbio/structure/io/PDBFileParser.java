/*
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
 * Created on 16.03.2004
 *
 */
package org.biojava.nbio.structure.io;

import static java.lang.Math.min;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Author;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.DBRef;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupIterator;
import org.biojava.nbio.structure.HetatomImpl;
import org.biojava.nbio.structure.JournalArticle;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Site;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompAtom;
import org.biojava.nbio.structure.io.util.PDBTemporaryStorageUtils.LinkRecord;
import org.biojava.nbio.structure.secstruc.SecStrucInfo;
import org.biojava.nbio.structure.secstruc.SecStrucType;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.biojava.nbio.structure.xtal.SymoplibParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class implements the actual PDB file parsing. Do not access it directly, but
 * via the PDBFileReader class.
 *
 * <h2>Parsing</h2>
 *
 * During the PDBfile parsing several Flags can be set. See the {@link #setFileParsingParameters(FileParsingParameters)} methods.
 *
 *
 * <p>
 * To provide excessive memory usage for large PDB files, there is the ATOM_CA_THRESHOLD.
 * If more Atoms than this threshold are being parsed in a PDB file, the parser will automatically
 * switch to a C-alpha only representation.
 *
 * <p>
 * The result of the parsing of the PDB file is a new {@link Structure} object.
 *
 * <p>
 * For more documentation on how to work with the Structure API please
 * see <a href="http://biojava.org/wiki/BioJava:CookBook#Protein_Structure" target="_top">
 * http://biojava.org/wiki/BioJava:CookBook#Protein_Structure</a>
 *
 *
 *
 *
 * <h2>Example</h2>
 * <p>
 * Q: How can I get a Structure object from a PDB file?
 * <p>
 * A:
 * <pre>
 * public {@link Structure} loadStructure(String pathToPDBFile){
 * 	// The PDBFileParser is wrapped by the PDBFileReader
 * 	{@link PDBFileReader} pdbreader = new {@link PDBFileReader}();
 *
 * 	{@link Structure} structure = null;
 * 	try{
 * 		structure = pdbreader.getStructure(pathToPDBFile);
 * 		System.out.println(structure);
 * 	} catch (IOException e) {
 * 		e.printStackTrace();
 * 	}
 * 	return structure;
 * }
 * </pre>
 *
 *
 * @author Andreas Prlic
 * @author Jules Jacobsen
 * @author Jose Duarte
 * @since 1.4
 */
public class PDBFileParser  {



	private static final Logger logger = LoggerFactory.getLogger(PDBFileParser.class);

	// for printing
	private static final String NEWLINE = System.getProperty("line.separator");


	// required for parsing:
	private String pdbId; //the actual id of the entry
	private Structure     structure;
	private List<List<Chain>> allModels; // a temp data structure to keep all models
	private List<Chain>   currentModel; // contains the ATOM records for each model
	private Chain         currentChain;
	private Group         currentGroup;

 	private List<Chain>   seqResChains; // contains all the chains for the SEQRES records
	//we're going to work on the assumption that the files are current -
	//if the pdb_HEADER_Handler detects a legacy format, this will be changed to true.
	//if true then lines will be truncated at 72 characters in certain cases
	//(pdb_COMPOUND_handler for example)
	private boolean isLegacyFormat = false;

	private boolean blankChainIdsPresent = false;
	
	// for re-creating the biological assembly
	private PDBBioAssemblyParser bioAssemblyParser = null;

	private PDBHeader pdbHeader;
	private PDBCrystallographicInfo crystallographicInfo;
	private JournalArticle journalArticle;
	private List<Map<String, Integer>> connects ;
	private List<Map<String,String>> helixList;
	private List<Map<String,String>> strandList;
	private List<Map<String,String>> turnList;

	private int lengthCheck ;

	private boolean isLastCompndLine = false;
	private boolean isLastSourceLine = false;
	private EntityInfo current_compound;
	private List<EntityInfo> entities = new ArrayList<EntityInfo>();
	private HashMap<Integer,List<String>> compoundMolIds2chainIds = new HashMap<Integer, List<String>>();
	private List<String> compndLines = new ArrayList<String>();
	private List<String> sourceLines = new ArrayList<String>();
	private List<String> journalLines = new ArrayList<String>();
	private List<DBRef> dbrefs;
	private Map<String, Site> siteMap = new LinkedHashMap<String, Site>();
	private Map<String, List<ResidueNumber>> siteToResidueMap = new LinkedHashMap<String, List<ResidueNumber>>();

	private List<SSBondImpl> ssbonds = new ArrayList<>();
	
    // for storing LINK until we have all the atoms parsed
    private List<LinkRecord> linkRecords;

	private Matrix4d currentNcsOp;
	private List<Matrix4d> ncsOperators;

	// for parsing COMPOUND and SOURCE Header lines
	private int prevMolId;
	private String previousContinuationField;
	private String continuationField;
	private String continuationString;

	private DateFormat dateFormat;

	// for rfree parsing
	private float rfreeStandardLine = -1;
	private float rfreeNoCutoffLine = -1;

	private static  final List<String> compndFieldValues = new ArrayList<String>(
			Arrays.asList(
					"MOL_ID:", "MOLECULE:", "CHAIN:", "SYNONYM:",
					"EC:", "FRAGMENT:", "ENGINEERED:", "MUTATION:",
					"BIOLOGICAL_UNIT:", "OTHER_DETAILS:"
					));


	private static final List<String> ignoreCompndFieldValues = new ArrayList<String>(
			Arrays.asList(
					"HETEROGEN:","ENGINEEREED:","FRAGMENT,",
					"MUTANT:","SYNTHETIC:"
					));
	// ENGINEEREED in pdb219d

	private static final List<String> sourceFieldValues = new ArrayList<String>(
			Arrays.asList("ENGINEERED:", "MOL_ID:", "SYNTHETIC:", "FRAGMENT:",
					"ORGANISM_SCIENTIFIC:", "ORGANISM_COMMON:",
					"ORGANISM_TAXID:","STRAIN:",
					"VARIANT:", "CELL_LINE:", "ATCC:", "ORGAN:", "TISSUE:",
					"CELL:", "ORGANELLE:", "SECRETION:", "GENE:",
					"CELLULAR_LOCATION:", "EXPRESSION_SYSTEM:",
					"EXPRESSION_SYSTEM_TAXID:",
					"EXPRESSION_SYSTEM_STRAIN:", "EXPRESSION_SYSTEM_VARIANT:",
					"EXPRESSION_SYSTEM_CELL_LINE:",
					"EXPRESSION_SYSTEM_ATCC_NUMBER:",
					"EXPRESSION_SYSTEM_ORGAN:", "EXPRESSION_SYSTEM_TISSUE:",
					"EXPRESSION_SYSTEM_CELL:", "EXPRESSION_SYSTEM_ORGANELLE:",
					"EXPRESSION_SYSTEM_CELLULAR_LOCATION:",
					"EXPRESSION_SYSTEM_VECTOR_TYPE:",
					"EXPRESSION_SYSTEM_VECTOR:", "EXPRESSION_SYSTEM_PLASMID:",
					"EXPRESSION_SYSTEM_GENE:", "OTHER_DETAILS:"));

	private int atomCount;

	// parsing options:

	private int atomCAThreshold ;

	private int loadMaxAtoms;

	private boolean atomOverflow;

	/** flag to tell parser to only read Calpha coordinates **/
	private boolean parseCAonly;


	private FileParsingParameters params;
	
	private boolean startOfMolecule;
	private boolean startOfModel;

	public PDBFileParser() {
		params = new FileParsingParameters();

		allModels = new ArrayList<>();
		structure     = null;
		currentModel  = null;
		currentChain  = null;
		currentGroup  = null;
		// we initialise to true since at the beginning of the file we are always starting a new molecule 
		startOfMolecule = true;
		startOfModel = true;

		
		pdbHeader 	  = new PDBHeader();
		crystallographicInfo = new PDBCrystallographicInfo();
		connects      = new ArrayList<Map<String,Integer>>() ;


		helixList     = new ArrayList<Map<String,String>>();
		strandList    = new ArrayList<Map<String,String>>();
		turnList      = new ArrayList<Map<String,String>>();
		current_compound = null;
		dbrefs        = new ArrayList<DBRef>();
		siteMap = null;
		dateFormat = new SimpleDateFormat("dd-MMM-yy", Locale.US);
		atomCount = 0;
		atomOverflow = false;
		parseCAonly = false;
		
		// this SHOULD not be done
		// DONOT:setFileParsingParameters(params);
		// set the correct max values for parsing...
		loadMaxAtoms = params.getMaxAtoms();
		atomCAThreshold = params.getAtomCaThreshold();
		
        linkRecords = new ArrayList<LinkRecord>();

		blankChainIdsPresent = false;
		
	}

	/** initiate new resNum, either Hetatom, Nucleotide, or AminoAcid */
	private Group getNewGroup(String recordName,Character aminoCode1, String aminoCode3) {

		Group g =  ChemCompGroupFactory.getGroupFromChemCompDictionary(aminoCode3);
		if ( g != null && !g.getChemComp().isEmpty())
			return g;


		Group group;
		if (aminoCode1 == null || StructureTools.UNKNOWN_GROUP_LABEL == aminoCode1 ){
			group = new HetatomImpl();

		} else if(StructureTools.isNucleotide(aminoCode3))  {
			// it is a nucleotide
			NucleotideImpl nu = new NucleotideImpl();
			group = nu;

		} else {
			AminoAcidImpl aa = new AminoAcidImpl() ;
			aa.setAminoType(aminoCode1);
			group = aa ;
		}

		//		System.out.println("new resNum type: "+ resNum.getType() );
		return  group ;
	}



	// Handler methods to deal with PDB file records properly.
	/**
	 Handler for
	 HEADER Record Format
	 <pre>
	 COLUMNS        DATA TYPE       FIELD           DEFINITION
	 ----------------------------------------------------------------------------------
	 1 -  6        Record name     "HEADER"
	 11 - 50        String(40)      classification  Classifies the molecule(s)
	 51 - 59        Date            depDate         Deposition date.  This is the date
	 the coordinates were received by
	 the PDB
	 63 - 66        IDcode          idCode          This identifier is unique within PDB
	</pre>
	 */
	private void pdb_HEADER_Handler(String line) {

		String classification  = null;
		String deposition_date = null;
		String pdbCode         = null;

		int len = line.trim().length();
		if(len > 10) {
			classification  = line.substring (10, min(len,50)).trim() ;
			pdbHeader.setClassification(classification);
		}
		if(len > 50) {
			deposition_date = line.substring (50, min(len,59)).trim() ;
			try {
				Date dep = dateFormat.parse(deposition_date);
				pdbHeader.setDepDate(dep);

			} catch (ParseException e){
				logger.info("Could not parse deposition date string '"+deposition_date+"'. Will continue without deposition date");
			}
		}
		if(len > 62) {
			pdbCode         = line.substring (62, min(len,66)).trim() ;
			pdbId = pdbCode;

			logger.debug("Parsing entry " + pdbId);


			structure.setPDBCode(pdbCode);
			pdbHeader.setIdCode(pdbCode);
		}

		//*really* old files (you'll need to hunt to find these as they
		//should have been remediated) have headers like below. Plus the
		//pdbId at positions 72-76 is present in every line

		//HEADER    PROTEINASE INHIBITOR (TRYPSIN)          05-OCT-84   5PTI      5PTI   3
		//HEADER    TRANSFERASE (ACYLTRANSFERASE)           02-SEP-92   1LAC      1LAC   2
		if (len > 66) {
			if (pdbId.equals(line.substring (72, 76))){
				isLegacyFormat = true;
				logger.warn(pdbId + " is a LEGACY entry - this will most likely not parse correctly.");
			}
		}

	}


	/** 
	 * Parses the following record:
	 * <pre>
	 *  COLUMNS      DATA  TYPE      FIELD         DEFINITION
	 * ------------------------------------------------------------------------------------
	 *  1 -  6      Record name     "AUTHOR"
	 *  9 - 10      Continuation    continuation  Allows concatenation of multiple records.
	 * 11 - 79      List            authorList    List of the author names, separated
	 *                                            by commas.
	 *
	 * </pre>
	 * @param line
	 */
	private void pdb_AUTHOR_Handler(String line) {

		String authors = line.substring(10).trim();

		String auth = pdbHeader.getAuthors();
		if (auth == null){
			pdbHeader.setAuthors(authors);
		} else {
			auth +=  authors;
			pdbHeader.setAuthors(auth);
		}

	}



	/** 
	 * Parses the following record:
	 *
	 * <pre>
	 * COLUMNS       DATA TYPE        FIELD        DEFINITION
	 * --------------------------------------------------------------------
	 *  1 -  6       Record name      "HELIX "
	 *  8 - 10       Integer          serNum       Serial number of the helix.
	 *                                             This starts at 1 and increases
	 *                                             incrementally.
	 * 12 - 14       LString(3)       helixID      Helix identifier. In addition
	 *                                             to a serial number, each helix is
	 *                                             given an alphanumeric character
	 *                                             helix identifier.
	 * 16 - 18       Residue name     initResName  Name of the initial residue.
	 * 20            Character        initChainID  Chain identifier for the chain
	 *                                             containing this helix.
	 * 22 - 25       Integer          initSeqNum   Sequence number of the initial
	 *                                             residue.
	 * 26            AChar            initICode    Insertion code of the initial
	 *                                             residue.
	 * 28 - 30       Residue name     endResName   Name of the terminal residue of
	 *                                             the helix.
	 * 32            Character        endChainID   Chain identifier for the chain
	 *                                             containing this helix.
	 * 34 - 37       Integer          endSeqNum    Sequence number of the terminal
	 *                                             residue.
	 * 38            AChar            endICode     Insertion code of the terminal
	 *                                             residue.
	 * 39 - 40       Integer          helixClass   Helix class (see below).
	 * 41 - 70       String           comment      Comment about this helix.
	 * 72 - 76       Integer          length       Length of this helix.
	 * </pre>
	 */
	private void pdb_HELIX_Handler(String line){

		if (params.isHeaderOnly()) return;

		if (line.length()<38) {
			logger.info("HELIX line has length under 38. Ignoring it.");
			return;
		}

		String initResName = line.substring(15,18).trim();
		String initChainId = line.substring(19,20);
		String initSeqNum  = line.substring(21,25).trim();
		String initICode   = line.substring(25,26);
		String endResName  = line.substring(27,30).trim();
		String endChainId  = line.substring(31,32);
		String endSeqNum   = line.substring(33,37).trim();
		String endICode    = line.substring(37,38);

		//System.out.println(initResName + " " + initChainId + " " + initSeqNum + " " + initICode + " " +
		//        endResName + " " + endChainId + " " + endSeqNum + " " + endICode);

		Map<String,String> m = new HashMap<String,String>();

		m.put("initResName",initResName);
		m.put("initChainId", initChainId);
		m.put("initSeqNum", initSeqNum);
		m.put("initICode", initICode);
		m.put("endResName", endResName);
		m.put("endChainId", endChainId);
		m.put("endSeqNum",endSeqNum);
		m.put("endICode",endICode);

		helixList.add(m);

	}

	/**
	 * Handler for
	 * <pre>
	 *       COLUMNS     DATA TYPE        FIELD           DEFINITION
	 * --------------------------------------------------------------
	 *  1 -  6     Record name      "SHEET "
	 *  8 - 10     Integer          strand       Strand number which starts at 1
	 *                                           for each strand within a sheet
	 *                                           and increases by one.
	 * 12 - 14     LString(3)       sheetID      Sheet identifier.
	 * 15 - 16     Integer          numStrands   Number of strands in sheet.
	 * 18 - 20     Residue name     initResName  Residue name of initial residue.
	 * 22          Character        initChainID  Chain identifier of initial
	 *                                           residue in strand.
	 * 23 - 26     Integer          initSeqNum   Sequence number of initial
	 *                                           residue in strand.
	 * 27          AChar            initICode    Insertion code of initial residue
	 *                                           in strand.
	 * 29 - 31     Residue name     endResName   Residue name of terminal residue.
	 * 33          Character        endChainID   Chain identifier of terminal
	 *                                           residue.
	 * 34 - 37     Integer          endSeqNum    Sequence number of terminal
	 *                                           residue.
	 * 38          AChar            endICode     Insertion code of terminal
	 *                                           residue.
	 * 39 - 40     Integer          sense        Sense of strand with respect to
	 *                                           previous strand in the sheet. 0
	 *                                           if first strand, 1 if parallel,
	 *                                           -1 if anti-parallel.
	 * 42 - 45     Atom             curAtom      Registration. Atom name in
	 *                                           current strand.
	 * 46 - 48     Residue name     curResName   Registration. Residue name in
	 *                                           current strand.
	 * 50          Character        curChainId   Registration. Chain identifier in
	 *                                           current strand.
	 * 51 - 54     Integer          curResSeq    Registration. Residue sequence
	 *                                           number in current strand.
	 * 55          AChar            curICode     Registration. Insertion code in
	 *                                           current strand.
	 * 57 - 60     Atom             prevAtom     Registration. Atom name in
	 *                                           previous strand.
	 * 61 - 63     Residue name     prevResName  Registration. Residue name in
	 *                                           previous strand.
	 * 65          Character        prevChainId  Registration. Chain identifier in
	 *                                           previous strand.
	 * 66 - 69     Integer          prevResSeq   Registration. Residue sequence
	 *                                           number in previous strand.
	 * 70          AChar            prevICode    Registration. Insertion code in
	 *                                               previous strand.
	 * </pre>
	 */
	private void pdb_SHEET_Handler( String line){

		if (params.isHeaderOnly()) return;

		if (line.length()<38) {
			logger.info("SHEET line has length under 38. Ignoring it.");
			return;
		}

		String initResName = line.substring(17,20).trim();
		String initChainId = line.substring(21,22);
		String initSeqNum  = line.substring(22,26).trim();
		String initICode   = line.substring(26,27);
		String endResName  = line.substring(28,31).trim();
		String endChainId  = line.substring(32,33);
		String endSeqNum   = line.substring(33,37).trim();
		String endICode    = line.substring(37,38);

		//System.out.println(initResName + " " + initChainId + " " + initSeqNum + " " + initICode + " " +
		//        endResName + " " + endChainId + " " + endSeqNum + " " + endICode);

		Map<String,String> m = new HashMap<String,String>();

		m.put("initResName",initResName);
		m.put("initChainId", initChainId);
		m.put("initSeqNum", initSeqNum);
		m.put("initICode", initICode);
		m.put("endResName", endResName);
		m.put("endChainId", endChainId);
		m.put("endSeqNum",endSeqNum);
		m.put("endICode",endICode);

		strandList.add(m);
	}


	/**
	 * Handler for TURN lines
	 * <pre>
	 * COLUMNS      DATA TYPE        FIELD         DEFINITION
	 * --------------------------------------------------------------------
	 *  1 -  6      Record name      "TURN "
	 *  8 - 10      Integer          seq           Turn number; starts with 1 and
	 *                                             increments by one.
	 * 12 - 14      LString(3)       turnId        Turn identifier
	 * 16 - 18      Residue name     initResName   Residue name of initial residue in
	 *                                             turn.
	 * 20           Character        initChainId   Chain identifier for the chain
	 *                                             containing this turn.
	 * 21 - 24      Integer          initSeqNum    Sequence number of initial residue
	 *                                             in turn.
	 * 25           AChar            initICode     Insertion code of initial residue
	 *                                             in turn.
	 * 27 - 29      Residue name     endResName    Residue name of terminal residue
	 *                                             of turn.
	 * 31           Character        endChainId    Chain identifier for the chain
	 *                                             containing this turn.
	 * 32 - 35      Integer          endSeqNum     Sequence number of terminal
	 *                                             residue of turn.
	 * 36           AChar            endICode      Insertion code of terminal residue
	 *                                             of turn.
	 * 41 - 70      String           comment       Associated comment.
	 * </pre>
	 * @param line
	 */
	private void pdb_TURN_Handler( String line){

		if (params.isHeaderOnly()) return;

		if (line.length()<36) {
			logger.info("TURN line has length under 36. Ignoring it.");
			return;
		}

		String initResName = line.substring(15,18).trim();
		String initChainId = line.substring(19,20);
		String initSeqNum  = line.substring(20,24).trim();
		String initICode   = line.substring(24,25);
		String endResName  = line.substring(26,29).trim();
		String endChainId  = line.substring(30,31);
		String endSeqNum   = line.substring(31,35).trim();
		String endICode    = line.substring(35,36);

		//System.out.println(initResName + " " + initChainId + " " + initSeqNum + " " + initICode + " " +
		//        endResName + " " + endChainId + " " + endSeqNum + " " + endICode);

		Map<String,String> m = new HashMap<String,String>();

		m.put("initResName",initResName);
		m.put("initChainId", initChainId);
		m.put("initSeqNum", initSeqNum);
		m.put("initICode", initICode);
		m.put("endResName", endResName);
		m.put("endChainId", endChainId);
		m.put("endSeqNum",endSeqNum);
		m.put("endICode",endICode);

		turnList.add(m);
	}

	/**
	 * Handler for
	 * REVDAT Record format:
	 * <pre>
	 *
	 * COLUMNS       DATA TYPE      FIELD         DEFINITION
	 * ----------------------------------------------------------------------------------
	 * 1 -  6       Record name    "REVDAT"
	 * 8 - 10       Integer        modNum        Modification number.
	 * 11 - 12       Continuation   continuation  Allows concatenation of multiple
	 * records.
	 * 14 - 22       Date           modDate       Date of modification (or release for
	 * new entries).  This is not repeated
	 * on continuation lines.
	 * 24 - 28       String(5)      modId         Identifies this particular
	 * modification.  It links to the
	 * archive used internally by PDB.
	 * This is not repeated on continuation
	 * lines.
	 * 32            Integer        modType       An integer identifying the type of
	 * modification.  In case of revisions
	 * with more than one possible modType,
	 * the highest value applicable will be
	 * assigned.
	 * 40 - 45       LString(6)     record        Name of the modified record.
	 * 47 - 52       LString(6)     record        Name of the modified record.
	 * 54 - 59       LString(6)     record        Name of the modified record.
	 * 61 - 66       LString(6)     record        Name of the modified record.
	 * </pre>
	 */
	private void pdb_REVDAT_Handler(String line) {

		// keep the first as latest modified date and the last as release date
		Date modDate = pdbHeader.getModDate();

		if ( modDate == null || modDate.equals(new Date(0)) ) {
			
			// modified date is still uninitialized
			String modificationDate = line.substring (13, 22).trim() ;

			try {
				Date dep = dateFormat.parse(modificationDate);
				pdbHeader.setModDate(dep);
				pdbHeader.setRelDate(dep);
			} catch (ParseException e){
				logger.info("Could not parse revision date string '"+modificationDate+"'. ");
			}

		} else {
			
			// set as the release date
			String releaseDate = line.substring (13, 22).trim() ;

			try {
				Date dep = dateFormat.parse(releaseDate);
				pdbHeader.setRelDate(dep);
			} catch (ParseException e){
				logger.info("Could not parse revision date string '"+releaseDate+"'. ");
			}
		}
	}

	/** 
	 * Handler for
	 * SEQRES record format
	 * SEQRES records contain the amino acid or nucleic acid sequence of residues in each chain of the macromolecule that was studied.
	 * <p>
	 * Record Format:
	 * <p>
	 * <pre>
	 * COLUMNS        DATA TYPE       FIELD         DEFINITION
	 * ---------------------------------------------------------------------------------
	 * 1 -  6        Record name     "SEQRES"
	 * 9 - 10        Integer         serNum        Serial number of the SEQRES record
	 * for the current chain.  Starts at 1
	 * and increments by one each line.
	 * Reset to 1 for each chain.
	 * 12             Character       chainID       Chain identifier.  This may be any
	 * single legal character, including a
	 * blank which is used if there is
	 * only one chain.
	 * 14 - 17        Integer         numRes        Number of residues in the chain.
	 * This value is repeated on every
	 * record.
	 * 20 - 22        Residue name    resName       Residue name.
	 * 24 - 26        Residue name    resName       Residue name.
	 * 28 - 30        Residue name    resName       Residue name.
	 * 32 - 34        Residue name    resName       Residue name.
	 * 36 - 38        Residue name    resName       Residue name.
	 * 40 - 42        Residue name    resName       Residue name.
	 * 44 - 46        Residue name    resName       Residue name.
	 * 48 - 50        Residue name    resName       Residue name.
	 * 52 - 54        Residue name    resName       Residue name.
	 * 56 - 58        Residue name    resName       Residue name.
	 * 60 - 62        Residue name    resName       Residue name.
	 * 64 - 66        Residue name    resName       Residue name.
	 * 68 - 70        Residue name    resName       Residue name.
	 * </pre>
	 * @author Jules Jacobsen
	 */
	private void pdb_SEQRES_Handler(String line) {

		/*
		 *          1         2         3         4         5         6         7
		 * 1234567890123456789012345678901234567890123456789012345678901234567890
		 * SEQRES   1 A  376  LYS PRO VAL THR VAL LYS LEU VAL ASP SER GLN ALA THR
		 * SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU
		 * SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR CYS ASN
		 * SEQRES   1 B   30  PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU
		 * SEQRES   2 B   30  ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR
		 * SEQRES   3 B   30  THR PRO LYS ALA
		 * SEQRES   1 C   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU
		 * SEQRES   2 C   21  TYR GLN LEU GLU ASN TYR CYS ASN
		 * SEQRES   1 D   30  PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU
		 * SEQRES   2 D   30  ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR
		 * SEQRES   3 D   30  THR PRO LYS ALA
		 */

		String recordName = line.substring(0, 6).trim();
		String chainID    = line.substring(11, 12);
		String newLength   = line.substring(13,17).trim();
		String subSequence = line.substring(18);

		if ( lengthCheck == -1 ){
			lengthCheck = Integer.parseInt(newLength);
		}

		StringTokenizer subSequenceResidues = new StringTokenizer(subSequence);

		Character aminoCode1 = null;
		if (! recordName.equals(AminoAcid.SEQRESRECORD)) {
			// should not have been called
			return;
		}

		currentChain = isKnownChain(chainID, seqResChains);
		if ( currentChain == null) {

			currentChain = new ChainImpl();
			currentChain.setId(chainID);
			currentChain.setName(chainID);

		}

		while (subSequenceResidues.hasMoreTokens()) {

			String threeLetter = subSequenceResidues.nextToken();

			aminoCode1 = StructureTools.get1LetterCode(threeLetter);

			//if (aminoCode1 == null) {
			// could be a nucleotide...
			// but getNewGroup takes care of that and converts ATOM records with aminoCode1 == nnull to nucleotide...
			//}
			currentGroup = getNewGroup("ATOM", aminoCode1, threeLetter);

			currentGroup.setPDBName(threeLetter);

			if ( currentGroup instanceof AminoAcid){
				AminoAcid aa = (AminoAcid)currentGroup;
				aa.setRecordType(AminoAcid.SEQRESRECORD);
			}
			// add the current resNum to the new chain.
			currentChain.addGroup(currentGroup);

		}
		Chain test = isKnownChain(chainID, seqResChains);

		if ( test == null)
			seqResChains.add(currentChain);

		if (currentGroup != null)
			currentGroup.trimToSize();

		currentGroup = null;
		currentChain = null;

		//		 the current chain is finished!
		//if ( current_chain.getLength() != lengthCheck ){
		//	System.err.println("the length of chain " + current_chain.getName() + "(" +
		//			current_chain.getLength() + ") does not match the expected " + lengthCheck);
		//}

		lengthCheck = Integer.parseInt(newLength);

	}



	/** 
	 * Handler for
	 * TITLE Record Format
	 * <pre>
	 COLUMNS        DATA TYPE       FIELD          DEFINITION
	 ----------------------------------------------------------------------------------
	 1 -  6        Record name     "TITLE "
	 9 - 10        Continuation    continuation   Allows concatenation of multiple
	 records.
	 11 - 70        String          title          Title of the experiment.
	 * </pre>
     *
	 */
	private void pdb_TITLE_Handler(String line) {
		String title;
		if ( line.length() > 79)
			title = line.substring(10,80).trim();
		else
			title = line.substring(10,line.length()).trim();

		String t = pdbHeader.getTitle();
		if ( (t != null) && (! t.equals("")) ){
			if (t.endsWith("-"))
				t += ""; // if last line ends with a hyphen then we don't add space
			else
				t += " ";
		}
		else t = "";

		t += title;

		pdbHeader.setTitle(t);
	}

	/**
	 * JRNL handler.
	 * The JRNL record contains the primary literature citation that describes the experiment which resulted
	 * in the deposited coordinate set. There is at most one JRNL reference per entry. If there is no primary
	 * reference, then there is no JRNL reference. Other references are given in REMARK 1.
	 *
	 * Record Format
	 * <pre>
	 * COLUMNS       DATA TYPE     FIELD         DEFINITION
	 * -----------------------------------------------------------------------
	 * 1 -  6       Record name   "JRNL  "
	 *
	 * 13 - 70       LString        text         See Details below.
	 * </pre>
	 */
	private void pdb_JRNL_Handler(String line) {
		//add the strings to the journalLines
		//the actual JournalArticle is then built when the whole entry is being
		//finalized with triggerEndFileChecks()
		//JRNL        TITL   NMR SOLUTION STRUCTURE OF RECOMBINANT TICK           1TAP  10
		if (line.substring(line.length() - 8, line.length() - 4).equals(pdbId)) {
			//trim off the trailing PDB id from legacy files.
			//are we really trying to still cater for these museum pieces?

			logger.debug("trimming legacy PDB id from end of JRNL section line");

			line = line.substring(0, line.length() - 8);
			journalLines.add(line);
		} else {
			journalLines.add(line);
		}
	}

	/**
	 * This should not be accessed directly, other than by </code>makeCompounds</code>. It still deals with the same
	 * lines in a similar manner but if not accessed from </code>makeCompounds</code> the last element will be
	 * missing. Don't say I didn't warn you.
	 *
	 * @param line
	 */
	private void pdb_COMPND_Handler(String line) {

		logger.debug("previousContinuationField  is "
					+ previousContinuationField);
		logger.debug("current continuationField  is "
					+ continuationField);
		logger.debug("current continuationString is "
					+ continuationString);
		logger.debug("current compound           is "
					+ current_compound);


		// In legacy PDB files the line ends with the PDB code and a serial number, chop those off!
		//format version 3.0 onwards will have 80 characters in a line
		//		if (line.length() > 72) {
		if (isLegacyFormat) {
			//                    if (DEBUG) {
			//                        System.out.println("We have a legacy file - truncating line length to 71 characters:");
			//                        System.out.println(line);
			//                    }
			line = line.substring(0, 72);
		}

		line = line.substring(10, line.length());


		String[] fieldList = line.trim().split("\\s+");
		int fl = fieldList.length;
		if ((fl >0 ) && compndFieldValues.contains(fieldList[0])) {

			continuationField = fieldList[0];
			if (previousContinuationField.equals("")) {
				previousContinuationField = continuationField;
			}

		} else if (fl>0) {
			// the ':' character indicates the end of a field name and should be invalid as part the first data token
			// e.g. obsolete file 1hhb has a malformed COMPND line that can only be caught with this kind of check
			if (fieldList[0].contains(":") ) {
				logger.info("COMPND line does not follow the PDB 3.0 format. Note that COMPND parsing is not supported any longer in format 2.3 or earlier");
				return;
			}

		} else {

			// the line will be added as data to the previous field
		}

		line = line.replace(continuationField, "").trim();

		StringTokenizer compndTokens = new StringTokenizer(line);

		//		System.out.println("PDBFileParser.pdb_COMPND_Handler: Tokenizing '" + line + "'");

		while (compndTokens.hasMoreTokens()) {
			String token = compndTokens.nextToken();

			if (previousContinuationField.equals("")) {
				previousContinuationField = continuationField;
			}

			if (previousContinuationField.equals(continuationField)
					&& compndFieldValues.contains(continuationField)) {

				logger.debug("Still in field " + continuationField);
				logger.debug("token = " + token);

				continuationString = continuationString.concat(token + " ");

				logger.debug("continuationString = "
							+ continuationString);

			}
			if (!continuationField.equals(previousContinuationField)) {

				if (continuationString.equals("")) {
					continuationString = token;

				} else {

					compndValueSetter(previousContinuationField,
							continuationString);
					previousContinuationField = continuationField;
					continuationString = token + " ";
				}
			} else if (ignoreCompndFieldValues.contains(token)) {
				// this field shall be ignored
				//continuationField = token;
			}
		}
		if (isLastCompndLine) {
			// final line in the section - finish off the compound
			//			System.out.println("[pdb_COMPND_Handler] Final COMPND line - Finishing off final MolID header.");
			compndValueSetter(continuationField, continuationString);
			continuationString = "";
			if (current_compound!=null) entities.add(current_compound);
		}
	}

	/**
	 * Set the value in the current molId object
	 * @param field
	 * @param value
	 */
	private void compndValueSetter(String field, String value) {

		value = value.trim().replace(";", "");
		if (field.equals("MOL_ID:")) {

			int i = -1;
			try {
				i = Integer.valueOf(value);
			} catch (NumberFormatException e){
				logger.warn("Value '{}' does not look like a number, while trying to parse COMPND MOL_ID line.",value);
			}
			if (i>0 && prevMolId!=i) {

				if (current_compound!=null) entities.add(current_compound);

				logger.debug("Initialising new Compound with mol_id {}", i);

				current_compound = new EntityInfo();

				current_compound.setMolId(i);
				
				// we will set polymer for all defined compounds in PDB file (non-polymer compounds are not defined in header) - JD 2016-03-25
				current_compound.setType(EntityType.POLYMER);

				prevMolId = i;
			}

		}

		// if for some reason (e.g. missing mol_id line) the current_compound is null we can't add anything to it, return
		if (current_compound==null) {
			return;
		}

		if (field.equals("MOLECULE:")) {
			current_compound.setDescription(value);

		}
		if (field.equals("CHAIN:")) {
			//System.out.println(value);
			StringTokenizer chainTokens = new StringTokenizer(value, ",");
			List<String> chains = new ArrayList<String>();

			while (chainTokens.hasMoreTokens()) {
				String chainID = chainTokens.nextToken().trim();
				// NULL is used in old PDB files to represent empty chain DI
				if (chainID.equals("NULL"))
					chainID = " ";
				chains.add(chainID);
			}
			compoundMolIds2chainIds.put(current_compound.getMolId(),chains);

		}
		if (field.equals("SYNONYM:")) {

			StringTokenizer synonyms = new StringTokenizer(value, ",");
			List<String> names = new ArrayList<String>();

			while (synonyms.hasMoreTokens()) {
				names.add(synonyms.nextToken());

				current_compound.setSynonyms(names);
			}

		}

		if (field.equals("EC:")) {

			StringTokenizer ecNumTokens = new StringTokenizer(value, ",");
			List<String> ecNums = new ArrayList<String>();

			while (ecNumTokens.hasMoreTokens()) {
				ecNums.add(ecNumTokens.nextToken());

				current_compound.setEcNums(ecNums);
			}

		}
		if (field.equals("FRAGMENT:")) {

			current_compound.setFragment(value);

		}
		if (field.equals("ENGINEERED:")) {

			current_compound.setEngineered(value);

		}
		if (field.equals("MUTATION:")) {

			current_compound.setMutation(value);

		}
		if (field.equals("BIOLOGICAL_UNIT:")) {

			current_compound.setBiologicalUnit(value);

		}
		if (field.equals("OTHER_DETAILS:")) {

			current_compound.setDetails(value);

		}

	}


	/** 
	 * Handler for
	 * SOURCE Record format
	 *
	 * The SOURCE record specifies the biological and/or chemical source of each biological molecule in the entry. Sources are described by both the common name and the scientific name, e.g., genus and species. Strain and/or cell-line for immortalized cells are given when they help to uniquely identify the biological entity studied.
	 * Record Format
	 * <pre>
	 * COLUMNS   DATA TYPE         FIELD          DEFINITION
	 * -------------------------------------------------------------------------------
	 *  1 -  6   Record name       "SOURCE"
	 *  9 - 10   Continuation      continuation   Allows concatenation of multiple records.
	 * 11 - 70   Specification     srcName        Identifies the source of the macromolecule in
	 *            list                            a token: value format.
	 * </pre>
	 * @param line the line to be parsed
	 */
	private void pdb_SOURCE_Handler(String line) {
		// works in the same way as the pdb_COMPND_Handler.
		String continuationNr = line.substring(9, 10).trim();



		logger.debug("current continuationNo     is "
				+ continuationNr);
		logger.debug("previousContinuationField  is "
				+ previousContinuationField);
		logger.debug("current continuationField  is "
				+ continuationField);
		logger.debug("current continuationString is "
				+ continuationString);
		logger.debug("current compound           is "
				+ current_compound);


		// following the docs, the last valid character should be 79, chop off the rest
		if (line.length() > 79) {
			line = line.substring(0, 79);
		}

		line = line.substring(10, line.length());

		logger.debug("LINE: >" + line + "<");

		String[] fieldList = line.split("\\s+");

		if (!fieldList[0].equals("")
				&& sourceFieldValues.contains(fieldList[0])) {
			//			System.out.println("[PDBFileParser.pdb_COMPND_Handler] Setting continuationField to '" + fieldList[0] + "'");
			continuationField = fieldList[0];
			if (previousContinuationField.equals("")) {
				previousContinuationField = continuationField;
			}

		} else if ((fieldList.length > 1) && ( sourceFieldValues.contains(fieldList[1]))) {
			//			System.out.println("[PDBFileParser.pdb_COMPND_Handler] Setting continuationField to '" + fieldList[1] + "'");
			continuationField = fieldList[1];
			if (previousContinuationField.equals("")) {
				previousContinuationField = continuationField;
			}

		} else {
			if (continuationNr.equals("")) {

				logger.debug("looks like an old PDB file");

				continuationField = "MOLECULE:";
				if (previousContinuationField.equals("")) {
					previousContinuationField = continuationField;
				}
			}

		}

		line = line.replace(continuationField, "").trim();

		StringTokenizer compndTokens = new StringTokenizer(line);

		//		System.out.println("PDBFileParser.pdb_COMPND_Handler: Tokenizing '" + line + "'");

		while (compndTokens.hasMoreTokens()) {
			String token = compndTokens.nextToken();

			if (previousContinuationField.equals("")) {
				//				System.out.println("previousContinuationField is empty. Setting to : " + continuationField);
				previousContinuationField = continuationField;
			}

			if (previousContinuationField.equals(continuationField)
					&& sourceFieldValues.contains(continuationField)) {

				logger.debug("Still in field " + continuationField);

				continuationString = continuationString.concat(token + " ");

				logger.debug("continuationString = "
							+ continuationString);
			}
			if (!continuationField.equals(previousContinuationField)) {

				if (continuationString.equals("")) {
					continuationString = token;

				} else {

					sourceValueSetter(previousContinuationField,
							continuationString);
					previousContinuationField = continuationField;
					continuationString = token + " ";
				}
			} else if (ignoreCompndFieldValues.contains(token)) {
				// this field shall be ignored
				//continuationField = token;
			}
		}
		if (isLastSourceLine) {
			// final line in the section - finish off the compound
			//			System.out.println("[pdb_SOURCE_Handler] Final SOURCE line - Finishing off final MolID header.");
			sourceValueSetter(continuationField, continuationString);
			continuationString = "";
			//compounds.add(current_compound);
		}

	}


	/** 
	 * Set the value in the current molId object
	 *
	 * @param field
	 * @param value
	 */
	private void sourceValueSetter(String field, String value) {

		value = value.trim().replace(";", "");
		//		System.out.println("[sourceValueSetter] " + field);
		if (field.equals("MOL_ID:")) {

			try {
				current_compound = entities.get(Integer.valueOf(value) - 1);
			} catch (NumberFormatException e){
				logger.info("could not process SOURCE MOL_ID record correctly:" + e.getMessage());
				return;
			}


			//			System.out.println("[sourceValueSetter] Fetching compound " + value + " " + current_compound.getMolId());

		}
		if (field.equals("SYNTHETIC:")) {
			current_compound.setSynthetic(value);
		} else if (field.equals("FRAGMENT:")) {
			current_compound.setFragment(value);
		} else if (field.equals("ORGANISM_SCIENTIFIC:")) {
			current_compound.setOrganismScientific(value);
		} else if (field.equals("ORGANISM_TAXID:")) {
			current_compound.setOrganismTaxId(value);
		} else if (field.equals("ORGANISM_COMMON:")) {
			current_compound.setOrganismCommon(value);
		} else if (field.equals("STRAIN:")) {
			current_compound.setStrain(value);
		} else if (field.equals("VARIANT:")) {
			current_compound.setVariant(value);
		} else if (field.equals("CELL_LINE:")) {
			current_compound.setCellLine(value);
		} else if (field.equals("ATCC:")) {
			current_compound.setAtcc(value);
		} else if (field.equals("ORGAN:")) {
			current_compound.setOrgan(value);
		} else if (field.equals("TISSUE:")) {
			current_compound.setTissue(value);
		} else if (field.equals("CELL:")) {
			current_compound.setCell(value);
		} else if (field.equals("ORGANELLE:")) {
			current_compound.setOrganelle(value);
		} else if (field.equals("SECRETION:")) {
			current_compound.setSecretion(value);
		} else if (field.equals("GENE:")) {
			current_compound.setGene(value);
		} else if (field.equals("CELLULAR_LOCATION:")) {
			current_compound.setCellularLocation(value);
		} else if (field.equals("EXPRESSION_SYSTEM:")) {
			current_compound.setExpressionSystem(value);
		} else if (field.equals("EXPRESSION_SYSTEM_TAXID:")) {
			current_compound.setExpressionSystemTaxId(value);
		} else if (field.equals("EXPRESSION_SYSTEM_STRAIN:")) {
			current_compound.setExpressionSystemStrain(value);
		} else if (field.equals("EXPRESSION_SYSTEM_VARIANT:")) {
			current_compound.setExpressionSystemVariant(value);
		} else if (field.equals("EXPRESSION_SYSTEM_CELL_LINE:")) {
			current_compound.setExpressionSystemCellLine(value);
		} else if (field.equals("EXPRESSION_SYSTEM_ATCC_NUMBER:")) {
			current_compound.setExpressionSystemAtccNumber(value);
		} else if (field.equals("EXPRESSION_SYSTEM_ORGAN:")) {
			current_compound.setExpressionSystemOrgan(value);
		} else if (field.equals("EXPRESSION_SYSTEM_TISSUE:")) {
			current_compound.setExpressionSystemTissue(value);
		} else if (field.equals("EXPRESSION_SYSTEM_CELL:")) {
			current_compound.setExpressionSystemCell(value);
		} else if (field.equals("EXPRESSION_SYSTEM_ORGANELLE:")) {
			current_compound.setExpressionSystemOrganelle(value);
		} else if (field.equals("EXPRESSION_SYSTEM_CELLULAR_LOCATION:")) {
			current_compound.setExpressionSystemCellularLocation(value);
		} else if (field.equals("EXPRESSION_SYSTEM_VECTOR_TYPE:")) {
			current_compound.setExpressionSystemVectorType(value);
		} else if (field.equals("EXPRESSION_SYSTEM_VECTOR:")) {
			current_compound.setExpressionSystemVector(value);
		} else if (field.equals("EXPRESSION_SYSTEM_PLASMID:")) {
			current_compound.setExpressionSystemPlasmid(value);
		} else if (field.equals("EXPRESSION_SYSTEM_GENE:")) {
			current_compound.setExpressionSystemGene(value);
		} else if (field.equals("OTHER_DETAILS:")) {
			current_compound.setExpressionSystemOtherDetails(value);
		}

	}

	/**
	 * Handler for REMARK lines
	 */
	private void pdb_REMARK_Handler(String line) {

		if ( line == null || line.length() < 11)
			return;


		if (line.startsWith("REMARK 800")) {
			pdb_REMARK_800_Handler(line);

		}  else if ( line.startsWith("REMARK 350")){

			if ( params.isParseBioAssembly()) {

				if (bioAssemblyParser == null){
					bioAssemblyParser = new PDBBioAssemblyParser();
				}

				bioAssemblyParser.pdb_REMARK_350_Handler(line);
			}

		// REMARK 3 (for R free)
		// note: if more than 1 value present (occurring in hybrid experimental technique entries, e.g. 3ins, 4n9m)
		// then last one encountered will be taken
		} else if (line.startsWith("REMARK   3   FREE R VALUE")) {

			// Rfree annotation is not very consistent in PDB format, it varies depending on the software
			// Here we follow this strategy:
			// a) take the '(NO CUTOFF)' value if the only one available (shelx software, e.g. 1x7q)
			// b) don't take it if also a line without '(NO CUTOFF)' is present (CNX software, e.g. 3lak)

			Pattern pR = Pattern.compile("^REMARK   3   FREE R VALUE\\s+(?:\\(NO CUTOFF\\))?\\s+:\\s+(\\d?\\.\\d+).*");
			Matcher mR = pR.matcher(line);
			if (mR.matches()) {
				try {
					rfreeNoCutoffLine = Float.parseFloat(mR.group(1));
				} catch (NumberFormatException e) {
					logger.info("Rfree value "+mR.group(1)+" does not look like a number, will ignore it");
				}
			}
			pR = Pattern.compile("^REMARK   3   FREE R VALUE\\s+:\\s+(\\d?\\.\\d+).*");
			mR = pR.matcher(line);
			if (mR.matches()) {
				try {
					rfreeStandardLine = Float.parseFloat(mR.group(1));
				} catch (NumberFormatException e) {
					logger.info("Rfree value '{}' does not look like a number, will ignore it", mR.group(1));
				}
			}

		// REMARK 3 RESOLUTION (contains more info than REMARK 2, for instance multiple resolutions in hybrid experimental technique entries)
		// note: if more than 1 value present (occurring in hybrid experimental technique entries, e.g. 3ins, 4n9m)
		// then last one encountered will be taken
		} else if (line.startsWith("REMARK   3   RESOLUTION RANGE HIGH")){
			Pattern pR = Pattern.compile("^REMARK   3   RESOLUTION RANGE HIGH \\(ANGSTROMS\\) :\\s+(\\d+\\.\\d+).*");
			Matcher mR = pR.matcher(line);
			if (mR.matches()) {
				try {
					float res = Float.parseFloat(mR.group(1));
					if (pdbHeader.getResolution()!=PDBHeader.DEFAULT_RESOLUTION) {
						logger.warn("More than 1 resolution value present, will use last one {} and discard previous {} "
								,mR.group(1), String.format("%4.2f",pdbHeader.getResolution()));
					}
					pdbHeader.setResolution(res);
				} catch (NumberFormatException e) {
					logger.info("Could not parse resolution '{}', ignoring it",mR.group(1));
				}
			}
		}

	}






	/** 
	 * Handler for
	 * EXPDTA Record Format
	<pre>
	 COLUMNS       DATA TYPE      FIELD         DEFINITION
	 -------------------------------------------------------------------------------
	 1 -  6       Record name    "EXPDTA"
	 9 - 10       Continuation   continuation  Allows concatenation of multiple
	 records.
	 11 - 70       SList          technique     The experimental technique(s) with
	 optional comment describing the
	 sample or experiment.

	 allowed techniques are:
	 ELECTRON DIFFRACTION
	 FIBER DIFFRACTION
	 FLUORESCENCE TRANSFER
	 NEUTRON DIFFRACTION
	 NMR
	 THEORETICAL MODEL
	 X-RAY DIFFRACTION
	</pre>
	 */
	private void pdb_EXPDTA_Handler(String line) {

		String technique  ;
		if (line.length() > 69)
			technique = line.substring (10, 70).trim() ;
		else
			technique = line.substring(10).trim();

		for (String singleTechnique: technique.split(";\\s+")) {
			pdbHeader.setExperimentalTechnique(singleTechnique);
		}


	}

	/** 
	 * Handler for
	 * CRYST1 Record Format
	 * The CRYST1 record presents the unit cell parameters, space group, and Z value.
	 * If the entry describes a structure determined by a technique other than X-ray crystallography,
	 * CRYST1 contains a = b = c = 1.0, alpha = beta = gamma = 90 degrees, space group = P 1, and Z =1.
	 * <pre>
	 * COLUMNS DATA TYPE    FIELD          DEFINITION
	 * -------------------------------------------------------------
	 *  1 - 6  Record name  "CRYST1"
	 *  7 - 15 Real(9.3)    a              a (Angstroms).
	 * 16 - 24 Real(9.3)    b              b (Angstroms).
	 * 25 - 33 Real(9.3)    c              c (Angstroms).
	 * 34 - 40 Real(7.2)    alpha          alpha (degrees).
	 * 41 - 47 Real(7.2)    beta           beta (degrees).
	 * 48 - 54 Real(7.2)    gamma          gamma (degrees).
	 * 56 - 66 LString      sGroup         Space group.
	 * 67 - 70 Integer      z              Z value.
	 * </pre>
	 */
	private void pdb_CRYST1_Handler(String line) {
		// for badly formatted files (e.g. phenix-produced ones), there's no z and the min length is 58 (e.g. for SG 'P 1')
		if (line.length() < 58) {
			logger.warn("CRYST1 record has fewer than 58 columns: will ignore it");
			return;
		}

		float a;
		float b;
		float c;
		float alpha;
		float beta;
		float gamma;
		String spaceGroup = "";

		try {
			a = Float.parseFloat(line.substring(6,15).trim());
			b = Float.parseFloat(line.substring(15,24).trim());
			c = Float.parseFloat(line.substring(24,33).trim());
			alpha = Float.parseFloat(line.substring(33,40).trim());
			beta = Float.parseFloat(line.substring(40,47).trim());
			gamma = Float.parseFloat(line.substring(47,54).trim());
		} catch (NumberFormatException e) {
			logger.info("could not parse CRYST1 record ("+e.getMessage()+") from line and ignoring it " + line);
			return ;
		}
		if (line.length()>=66) {
			// for well formatted files
			spaceGroup = line.substring(55,66).trim();
		} else {
			// for not-so-well formatted files, e.g. phenix-produced ones: they lack a Z value
			spaceGroup = line.substring(55,line.length()).trim();
		}

		CrystalCell xtalCell = new CrystalCell();
		xtalCell.setA(a);
		xtalCell.setB(b);
		xtalCell.setC(c);
		xtalCell.setAlpha(alpha);
		xtalCell.setBeta(beta);
		xtalCell.setGamma(gamma);

		if (!xtalCell.isCellReasonable()) {
			// If the entry describes a structure determined by a technique other than X-ray crystallography,
		    // CRYST1 contains a = b = c = 1.0, alpha = beta = gamma = 90 degrees, space group = P 1, and Z =1.
			// if so we don't add the crystal cell and it remains null
			logger.debug("The crystal cell read from file does not have reasonable dimensions (at least one dimension is below {}), discarding it.",
					CrystalCell.MIN_VALID_CELL_SIZE);
		} else {
			crystallographicInfo.setCrystalCell(xtalCell);
		}

		SpaceGroup sg = SymoplibParser.getSpaceGroup(spaceGroup);
		if (sg==null) {
			logger.warn("Space group '"+spaceGroup+"' not recognised as a standard space group");
			crystallographicInfo.setNonStandardSg(true);
		} else {
			crystallographicInfo.setSpaceGroup(sg);
			crystallographicInfo.setNonStandardSg(false);
		}
	}

	/**
	 * Handler for MTRIXn records. They specify extra NCS operators (usually in virus entries)
	 *
	 * See http://www.wwpdb.org/documentation/format33/sect8.html#MTRIXn
	 * <pre>
	 * COLUMNS        DATA TYPE     FIELD         DEFINITION
	 * -------------------------------------------------------------
	 *
	 *  1 -  6        Record name   "MTRIXn"      n=1, 2, or 3
	 *  8 - 10        Integer       serial        Serial number.
	 * 11 - 20        Real(10.6)    m[n][1]       Mn1
	 * 21 - 30        Real(10.6)    m[n][2]       Mn2
	 * 31 - 40        Real(10.6)    m[n][3]       Mn3
	 * 46 - 55        Real(10.5)    v[n]          Vn
	 * 60             Integer       iGiven        1
	 *
	 * </pre>
	 * Note that we ignore operators with iGiven==1
	 * 
	 * @param line
	 */
	private void pdb_MTRIXn_Handler(String line) {

		// don't process incomplete records
		if (line.length() < 60) {
			logger.info("MTRIXn record has fewer than 60 columns: will ignore it");
			return;
		}


		try {

			int rowIndex = Integer.parseInt(line.substring(5,6));
			double col1Value = Double.parseDouble(line.substring(10,20));
			double col2Value = Double.parseDouble(line.substring(20,30));
			double col3Value = Double.parseDouble(line.substring(30,40));
			double translValue = Double.parseDouble(line.substring(45,55));
			int iGiven = 0;
			if (!line.substring(59,60).trim().equals("")) {
				iGiven = Integer.parseInt(line.substring(59,60));
			}

			if (iGiven == 1) return;

			if (ncsOperators==null) {
				// we initialise on first pass
				ncsOperators = new ArrayList<Matrix4d>();
			}

			if (currentNcsOp==null) {
				currentNcsOp = new Matrix4d(1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1); // initialised to identity
			}

			currentNcsOp.setElement(rowIndex-1, 0, col1Value);
			currentNcsOp.setElement(rowIndex-1, 1, col2Value);
			currentNcsOp.setElement(rowIndex-1, 2, col3Value);
			currentNcsOp.setElement(rowIndex-1, 3, translValue);


			if (rowIndex==3) {
				ncsOperators.add(currentNcsOp);
				// we initialise for next matrix to come
				currentNcsOp = new Matrix4d(1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1); // initialised to identity
			}

		} catch (NumberFormatException e) {
			logger.info("Could not parse a number in MTRIXn record ("+e.getMessage()+") from line: >" + line+"<");
		}
	}

	/**
	 * Handler for ATOM.
	 * Record Format:
	 *
	 * <pre>
	 * ATOM      1  N   ASP A  15     110.964  24.941  59.191  1.00 83.44           N
	 *
	 * COLUMNS        DATA TYPE       FIELD         DEFINITION
	 * ---------------------------------------------------------------------------------
	 * 1 -  6        Record name     "ATOM  "
	 * 7 - 11        Integer         serial        Atom serial number.
	 * 13 - 16        Atom            name          Atom name.
	 * 17             Character       altLoc        Alternate location indicator.
	 * 18 - 20        Residue name    resName       Residue name.
	 * 22             Character       chainID       Chain identifier.
	 * 23 - 26        Integer         resSeq        Residue sequence number.
	 * 27             AChar           iCode         Code for insertion of residues.
	 * 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in Angstroms.
	 * 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in Angstroms.
	 * 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in Angstroms.
	 * 55 - 60        Real(6.2)       occupancy     Occupancy.
	 * 61 - 66        Real(6.2)       tempFactor    Temperature factor.
	 * 73 - 76        LString(4)      segID         Segment identifier, left-justified.
	 * 77 - 78        LString(2)      element       Element symbol, right-justified.
	 * 79 - 80        LString(2)      charge        Charge on the atom.
	 * </pre>
	 */
	private void  pdb_ATOM_Handler(String line)	{

		if ( params.isHeaderOnly())
			return;

		// let's first get the chain name which will serve to identify if we are starting a new molecule
		String chainName      = line.substring(21,22);
		
		if (chainName.equals(" ")) {
			blankChainIdsPresent = true;
		}
		
		if (currentChain!=null && !currentChain.getName().equals(chainName)) {
			// new chain name: another molecule coming
			startOfMolecule = true;
		}
		
		if (startOfMolecule) {
			// we add last chain if there was one
			if (currentChain!=null) {
				currentModel.add(currentChain);
				// let's not forget adding the last group to the finishing chain
				if (currentGroup!=null) {
					currentChain.addGroup(currentGroup);
				}
			}
			// we initialise the new molecule to come
			currentChain = new ChainImpl();
			// note that the chainId (asym id) is set properly later in assignAsymIds
			currentChain.setId(chainName);
			currentChain.setName(chainName);
			
		}

		if (startOfModel) {
			// we add last model if there was one
			if (currentModel!=null) {
				allModels.add(currentModel);
			}
			// we initialise the model to come
			currentModel = new ArrayList<>();
		}
		
		
		// let's get the residue number and see if we need to start a new group

		String groupCode3     = line.substring(17,20).trim();
		String resNum  = line.substring(22,26).trim();
		Character iCode = line.substring(26,27).charAt(0);
		if ( iCode == ' ')
			iCode = null;
		ResidueNumber residueNumber = new ResidueNumber(chainName, Integer.valueOf(resNum), iCode);

		//recordName      groupCode3
		//|                |    resNum
		//|                |    |   iCode
		//|     |          | |  |   ||
		//ATOM      1  N   ASP A  15     110.964  24.941  59.191  1.00 83.44           N
		//ATOM   1964  N   ARG H 221A      5.963 -16.715  27.669  1.00 28.59           N

		Character aminoCode1 = StructureTools.get1LetterCode(groupCode3);

		String recordName     = line.substring (0, 6).trim ();

		boolean isHetAtomInFile = false;
		
		if (recordName.equals("HETATM") ){
			// HETATOM RECORDS are treated slightly differently
			// some modified amino acids that we want to treat as amino acids
			// can be found as HETATOM records
			if ( aminoCode1 != null && aminoCode1.equals(StructureTools.UNKNOWN_GROUP_LABEL))
					aminoCode1 = null;
			
			isHetAtomInFile = true;
		}

		if ( startOfMolecule) {

			currentGroup = getNewGroup(recordName, aminoCode1, groupCode3);

			currentGroup.setPDBName(groupCode3);
			currentGroup.setResidueNumber(residueNumber);
			currentGroup.setHetAtomInFile(isHetAtomInFile);

		}
		
		// resetting states
		startOfModel = false;
		startOfMolecule = false;


		Character altLoc   = new Character(line.substring (16, 17).charAt(0));
		Group altGroup = null;


		// check if residue number is the same ...
		if ( ! residueNumber.equals(currentGroup.getResidueNumber())) {

			currentChain.addGroup(currentGroup);
			currentGroup.trimToSize();

			currentGroup = getNewGroup(recordName, aminoCode1, groupCode3);

			currentGroup.setPDBName(groupCode3);
			currentGroup.setResidueNumber(residueNumber);
			currentGroup.setHetAtomInFile(isHetAtomInFile);

		} else {
			// same residueNumber, but altLocs...

			// test altLoc
			if ( ! altLoc.equals(' ')) {
				logger.debug("found altLoc! " + currentGroup + " " + altGroup);
				altGroup = getCorrectAltLocGroup( altLoc,recordName,aminoCode1,groupCode3);
				if ( altGroup.getChain() == null) {
					// need to set current chain
					altGroup.setChain(currentChain);
				}

			}
		}

		atomCount++;

		if ( atomCount == atomCAThreshold ) {
			// throw away the SEQRES lines - too much to deal with...
			logger.warn("more than " + atomCAThreshold + " atoms in this structure, ignoring the SEQRES lines");
			seqResChains.clear();

			switchCAOnly();

		}



		if ( atomCount == loadMaxAtoms){
			logger.warn("File has more atoms than max specified in parsing parameters ({}). Ignoring atoms after line: {}", loadMaxAtoms, line);
			return;
		}
		if ( atomCount > loadMaxAtoms){
			return;
		}


		//          1         2         3         4         5         6
		//012345678901234567890123456789012345678901234567890123456789
		//ATOM      1  N   MET     1      20.154  29.699   5.276   1.0
		//ATOM    112  CA  ASP   112      41.017  33.527  28.371  1.00  0.00
		//ATOM     53  CA  MET     7      23.772  33.989 -21.600  1.00  0.00           C
		//ATOM    112  CA  ASP   112      37.613  26.621  33.571     0     0


		String fullname = line.substring (12, 16);

		// check for CA only if requested
		if ( parseCAonly ){
			// yes , user wants to get CA only
			// only parse CA atoms...
			if (! fullname.equals(" CA ")){
				//System.out.println("ignoring " + line);
				atomCount--;
				return;
			}
		}

		if ( params.getAcceptedAtomNames() != null) {

			boolean found = false;
			for (String ok : params.getAcceptedAtomNames()){
				//System.out.println(ok + "< >" + fullname +"<");

				if ( ok.equals(fullname.trim())) {
					found = true;
					break;
				}
			}
			if ( ! found) {
				atomCount--;
				return;
			}
		}
		// create new atom

		int pdbnumber = Integer.parseInt (line.substring (6, 11).trim ());
		AtomImpl atom = new AtomImpl() ;
		atom.setPDBserial(pdbnumber) ;

		atom.setAltLoc(altLoc);
		atom.setName(fullname.trim());

		double x = Double.parseDouble (line.substring (30, 38).trim());
		double y = Double.parseDouble (line.substring (38, 46).trim());
		double z = Double.parseDouble (line.substring (46, 54).trim());

		double[] coords = new double[3];
		coords[0] = x ;
		coords[1] = y ;
		coords[2] = z ;
		atom.setCoords(coords);

		float occu  = 1.0f;
		if ( line.length() > 59 ) {
			try {
				// occu and tempf are sometimes not used :-/
				occu = Float.parseFloat (line.substring (54, 60).trim());
			}  catch (NumberFormatException e){}
		}

		float tempf = 0.0f;
		if ( line.length() > 65) {
			try {
				tempf = Float.parseFloat (line.substring (60, 66).trim());
			}  catch (NumberFormatException e){}
		}

		atom.setOccupancy(  occu  );
		atom.setTempFactor( tempf );




		// Parse element from the element field. If this field is
		// missing (i.e. misformatted PDB file), then parse the
		// element from the chemical component.
		Element element = Element.R;
		boolean guessElement = true;
		if ( line.length() > 77 ) {
			// parse element from element field
			String elementSymbol = line.substring(76, 78).trim();
			if (elementSymbol.isEmpty()) {
				logger.info("Element column was empty for atom {} {}. Assigning atom element "
						+ "from Chemical Component Dictionary information", fullname.trim(), pdbnumber);
			} else {
			
				try {
					element = Element.valueOfIgnoreCase(elementSymbol);
					guessElement = false;
				}  catch (IllegalArgumentException e){
					logger.info("Element {} of atom {} {} was not recognised. Assigning atom element "
							+ "from Chemical Component Dictionary information", elementSymbol, 
							fullname.trim(), pdbnumber);
				}
			}
		} else {
			logger.info("Missformatted PDB file: element column of atom {} {} is not present. "
					+ "Assigning atom element from Chemical Component Dictionary information",
					fullname.trim(), pdbnumber);
		}
		if (guessElement) {
			String elementSymbol = null;
			if (currentGroup.getChemComp() != null) {
				for (ChemCompAtom a : currentGroup.getChemComp().getAtoms()) {
					if (a.getAtom_id().equals(fullname.trim())) {
						elementSymbol = a.getType_symbol();
						break;
					}
				}
				if (elementSymbol == null) {
					logger.info("Atom name {} was not found in the Chemical Component Dictionary information of {}. "
							+ "Assigning generic element R to it", fullname.trim(), currentGroup.getPDBName());
				} else {
					try {
						element = Element.valueOfIgnoreCase(elementSymbol);
					} catch (IllegalArgumentException e) {
						// this can still happen for cases like UNK
						logger.info("Element symbol {} found in chemical component dictionary for Atom {} {} could not be recognised as a known element. "
								+ "Assigning generic element R to it", elementSymbol, fullname.trim(), pdbnumber);
					}
				}			
			} else {
				logger.warn("Chemical Component Dictionary information was not found for Atom name {}. "
						+ "Assigning generic element R to it", fullname.trim());
			}
			
		}
		atom.setElement(element);


		//see if chain_id is one of the previous chains ...
		if ( altGroup != null) {
			altGroup.addAtom(atom);
			altGroup = null;
		}
		else {
			currentGroup.addAtom(atom);
		}


		// make sure that main group has all atoms
		// GitHub issue: #76
		if ( ! currentGroup.hasAtom(atom.getName())) {
			currentGroup.addAtom(atom);
		}



	}


	private Group getCorrectAltLocGroup( Character altLoc,
			String recordName, Character aminoCode1, String groupCode3) {

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
		Group altLocG = getNewGroup(recordName,aminoCode1,groupCode3);


		altLocG.setPDBName(groupCode3);

		altLocG.setResidueNumber(currentGroup.getResidueNumber());
		currentGroup.addAltLoc(altLocG);
		return altLocG;
	}

	private void switchCAOnly(){
		parseCAonly = true;


		currentModel = CAConverter.getRepresentativeAtomsOnly(currentModel);

		for ( int i =0; i< structure.nrModels() ; i++){
			//  iterate over all known models ...
			List<Chain> model = structure.getModel(i);
			model = CAConverter.getRepresentativeAtomsOnly(model);
			structure.setModel(i,model);
		}

		currentChain = CAConverter.getRepresentativeAtomsOnly(currentChain);

	}


	/** safes repeating a few lines ... */
	private Integer conect_helper (String line,int start,int end) {
		if (line.length() < end) return null;
		
		String sbond = line.substring(start,end).trim();
		int bond  = -1 ;
		Integer b = null ;

		if ( ! sbond.equals("")) {
			bond = Integer.parseInt(sbond);
			b = new Integer(bond);
		}

		return b ;
	}

	/**
	 * Handler for CONECT Record Format
	<pre>
	 COLUMNS         DATA TYPE        FIELD           DEFINITION
	 ---------------------------------------------------------------------------------
	 1 -  6         Record name      "CONECT"
	 7 - 11         Integer          serial          Atom serial number
	 12 - 16         Integer          serial          Serial number of bonded atom
	 17 - 21         Integer          serial          Serial number of bonded atom
	 22 - 26         Integer          serial          Serial number of bonded atom
	 27 - 31         Integer          serial          Serial number of bonded atom
	 32 - 36         Integer          serial          Serial number of hydrogen bonded
	 atom
	 37 - 41         Integer          serial          Serial number of hydrogen bonded
	 atom
	 42 - 46         Integer          serial          Serial number of salt bridged
	 atom
	 47 - 51         Integer          serial          Serial number of hydrogen bonded
	 atom
	 52 - 56         Integer          serial          Serial number of hydrogen bonded
	 atom
	 57 - 61         Integer          serial          Serial number of salt bridged
	 atom
	 </pre>
	 */
	private void pdb_CONECT_Handler(String line) {

		if ( atomOverflow) {
			return ;
		}
		if (params.isHeaderOnly()) {
			return;
		}
		
		// this try .. catch is e.g. to catch 1gte which has wrongly formatted lines...
		try {
			int atomserial = Integer.parseInt (line.substring(6 ,11).trim());
			Integer bond1      = conect_helper(line,11,16);
			Integer bond2      = conect_helper(line,16,21);
			Integer bond3      = conect_helper(line,21,26);
			Integer bond4      = conect_helper(line,26,31);
			Integer hyd1       = conect_helper(line,31,36);
			Integer hyd2       = conect_helper(line,36,41);
			Integer salt1      = conect_helper(line,41,46);
			Integer hyd3       = conect_helper(line,46,51);
			Integer hyd4       = conect_helper(line,51,56);
			Integer salt2      = conect_helper(line,56,61);

			//System.out.println(atomserial+ " "+ bond1 +" "+bond2+ " " +bond3+" "+bond4+" "+
			//		   hyd1+" "+hyd2 +" "+salt1+" "+hyd3+" "+hyd4+" "+salt2);
			HashMap<String, Integer> cons = new HashMap<String, Integer>();
			cons.put("atomserial",new Integer(atomserial));

			if ( bond1 != null) cons.put("bond1",bond1);
			if ( bond2 != null) cons.put("bond2",bond2);
			if ( bond3 != null) cons.put("bond3",bond3);
			if ( bond4 != null) cons.put("bond4",bond4);
			if ( hyd1  != null) cons.put("hydrogen1",hyd1);
			if ( hyd2  != null) cons.put("hydrogen2",hyd2);
			if ( salt1 != null) cons.put("salt1",salt1);
			if ( hyd3  != null) cons.put("hydrogen3",hyd3);
			if ( hyd4  != null) cons.put("hydrogen4",hyd4);
			if ( salt2 != null) cons.put("salt2",salt2);

			connects.add(cons);
		} catch (NumberFormatException e){
			logger.info("could not parse CONECT line correctly ("+e.getMessage()+"), at line : " + line);
			return;
		}
	}

	/**
	 * Handler for MODEL Record Format
	 * <pre>
	 * COLUMNS       DATA TYPE      FIELD         DEFINITION
	 * ----------------------------------------------------------------------
	 * 1 -  6       Record name    "MODEL "
	 * 11 - 14       Integer        serial        Model serial number.
	 * </pre>
	 */
	private void pdb_MODEL_Handler(String line) {

		if (params.isHeaderOnly()) return;
		
		// new model: we start a new molecule
		startOfMolecule = true;
		startOfModel = true;

	}
	
	/**
	 * Handler for TER record. The record is used in deposited PDB files and many others,
	 * but it's often forgotten by some softwares. In any case it helps identifying the 
	 * start of ligand molecules so we use it for that.
	 */
	private void pdb_TER_Handler() {
		startOfMolecule = true;		
	}


	/**
	 * DBREF handler
	 * <pre>
	 * COLUMNS       DATA TYPE          FIELD          DEFINITION
	 * ----------------------------------------------------------------
	 *  1 - 6        Record name        "DBREF "
	 *  8 - 11       IDcode             idCode         ID code of this entry.
	 * 13            Character          chainID        Chain identifier.
	 * 15 - 18       Integer            seqBegin       Initial sequence number
	 *                                                 of the PDB sequence segment.
	 * 19            AChar              insertBegin    Initial insertion code
	 *                                                 of the PDB sequence segment.
	 * 21 - 24       Integer            seqEnd         Ending sequence number
	 *                                                 of the PDB sequence segment.
	 * 25            AChar              insertEnd      Ending insertion code
	 *                                                 of the PDB sequence segment.
	 * 27 - 32       LString            database       Sequence database name.
	 * 34 - 41       LString            dbAccession    Sequence database accession code.
	 * 43 - 54      LString            dbIdCode        Sequence database
	 *                                                 identification code.
	 * 56 - 60      Integer            dbseqBegin      Initial sequence number of the
	 *                                                 database seqment.
	 * 61           AChar              idbnsBeg        Insertion code of initial residue
	 *                                                 of the segment, if PDB is the
	 *                                                 reference.
	 * 63 - 67      Integer            dbseqEnd        Ending sequence number of the
	 *                                                 database segment.
	 * 68           AChar              dbinsEnd        Insertion code of the ending
	 *                                                 residue of the segment, if PDB is
	 *                                                 the reference.
	 * </pre>
	 */
	private void pdb_DBREF_Handler(String line){

		logger.debug("Parsing DBREF " + line);

		DBRef dbref = new DBRef();
		String idCode      = line.substring(7,11);
		String chainName     = line.substring(12,13);
		String seqBegin    = line.substring(14,18);
		String insertBegin = line.substring(18,19);
		String seqEnd      = line.substring(20,24);
		String insertEnd   = line.substring(24,25);
		String database    = line.substring(26,32);
		String dbAccession = line.substring(33,41);
		String dbIdCode    = line.substring(42,54);
		String dbseqBegin  = line.substring(55,60);
		String idbnsBeg    = line.substring(60,61);
		String dbseqEnd    = line.substring(62,67);
		// Support implicit space character at end
		String dbinsEnd;
		if(line.length() >= 68)
			dbinsEnd       = line.substring(67,68);
		else
			dbinsEnd       = " ";

		dbref.setIdCode(idCode);
		dbref.setChainName(chainName);
		dbref.setSeqBegin(intFromString(seqBegin));
		dbref.setInsertBegin(insertBegin.charAt(0));
		dbref.setSeqEnd(intFromString(seqEnd));
		dbref.setInsertEnd(insertEnd.charAt(0));
		dbref.setDatabase(database.trim());
		dbref.setDbAccession(dbAccession.trim());
		dbref.setDbIdCode(dbIdCode.trim());
		dbref.setDbSeqBegin(intFromString(dbseqBegin));
		dbref.setIdbnsBegin(idbnsBeg.charAt(0));
		dbref.setDbSeqEnd(intFromString(dbseqEnd));
		dbref.setIdbnsEnd(dbinsEnd.charAt(0));

		//System.out.println(dbref.toPDB());
		dbrefs.add(dbref);
	}


	/**
	 * Process the disulfide bond info provided by an SSBOND record
	 *
	 * <pre>
	COLUMNS        DATA TYPE       FIELD         DEFINITION
	-------------------------------------------------------------------
	 1 -  6        Record name     "SSBOND"
	 8 - 10        Integer         serNum       Serial number.
	12 - 14        LString(3)      "CYS"        Residue name.
	16             Character       chainID1     Chain identifier.
	18 - 21        Integer         seqNum1      Residue sequence number.
	22             AChar           icode1       Insertion code.
	26 - 28        LString(3)      "CYS"        Residue name.
	30             Character       chainID2     Chain identifier.
	32 - 35        Integer         seqNum2      Residue sequence number.
	36             AChar           icode2       Insertion code.
	60 - 65        SymOP           sym1         Symmetry oper for 1st resid
	67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
	 * </pre>
	 */
	private void pdb_SSBOND_Handler(String line){

		if (params.isHeaderOnly()) return;

		if (line.length()<36) {
			logger.info("SSBOND line has length under 36. Ignoring it.");
			return;
		}

		String chain1      = line.substring(15,16);
		String seqNum1     = line.substring(17,21).trim();
		String icode1      = line.substring(21,22);
		String chain2      = line.substring(29,30);
		String seqNum2     = line.substring(31,35).trim();
		String icode2      = line.substring(35,36);

		if (line.length()>=72) {
			String symop1 = line.substring(59, 65).trim();
			String symop2 = line.substring(66, 72).trim();

			// until we implement proper treatment of symmetry in biojava #220, we can't deal with sym-related parteners properly, skipping them
			if (!symop1.equals("") && !symop2.equals("") && // in case the field is missing
					(!symop1.equals("1555") || !symop2.equals("1555")) ) {
				logger.info("Skipping ss bond between groups {} and {} belonging to different symmetry partners, because it is not supported yet", seqNum1+icode1, seqNum2+icode2);
				return;
			}
		}

		if (icode1.equals(" "))
			icode1 = "";
		if (icode2.equals(" "))
			icode2 = "";

		SSBondImpl ssbond = new SSBondImpl();

		ssbond.setChainID1(chain1);
		ssbond.setResnum1(seqNum1);
		ssbond.setChainID2(chain2);
		ssbond.setResnum2(seqNum2);
		ssbond.setInsCode1(icode1);
		ssbond.setInsCode2(icode2);
		ssbonds.add(ssbond);
	}


	/**
	 * Takes care of LINK records. These take the format of:
	 *
	 * <pre>
	 * COLUMNS        DATA TYPE       FIELD       DEFINITION
	 * --------------------------------------------------------------------------------
	 *  1 -  6        Record name     "LINK  "
	 * 13 - 16        Atom            name1       Atom name.
	 * 17             Character       altLoc1     Alternate location indicator.
	 * 18 - 20        Residue name    resName1    Residue name.
	 * 22             Character       chainID1    Chain identifier.
	 * 23 - 26        Integer         resSeq1     Residue sequence number.
	 * 27             AChar           iCode1      Insertion code.
	 * 43 - 46        Atom            name2       Atom name.
	 * 47             Character       altLoc2     Alternate location indicator.
	 * 48 - 50        Residue name    resName2    Residue name.
	 * 52             Character       chainID2    Chain identifier.
	 * 53 - 56        Integer         resSeq2     Residue sequence number.
	 * 57             AChar           iCode2      Insertion code.
	 * 60 - 65        SymOP           sym1        Symmetry operator for 1st atom.
	 * 67 - 72        SymOP           sym2        Symmetry operator for 2nd atom.
	 * </pre>
	 *
	 * (From http://www.wwpdb.org/documentation/format32/sect6.html#LINK)
	 *
	 * @param line the LINK record line to parse.
	 */
	private void pdb_LINK_Handler(String line) {

		if (params.isHeaderOnly()) return;
		
		// Check for the minimal set of fields.
		if (line.length()<56) {
			logger.info("LINK line has length under 56. Ignoring it.");
			return;
		}

		int len = line.length();
		
		String name1 = line.substring(12, 16).trim();
		String altLoc1 = line.substring(16, 17).trim();
		String resName1 = line.substring(17, 20).trim();
		String chainID1 = line.substring(21, 22).trim();
		String resSeq1 = line.substring(22, 26).trim();
		String iCode1 = line.substring(26, 27).trim();

		String name2 = line.substring(42, 46).trim();
		String altLoc2 = line.substring(46, 47).trim();
		String resName2 = line.substring(47, 50).trim();
		String chainID2 = line.substring(51, 52).trim();
		String resSeq2 = line.substring(52, 56).trim();
		String iCode2 = null;  // Might get trimmed if blank.
		if (len > 56) iCode2 = line.substring(56, 57).trim();

		String sym1 = null;
		if (len > 64) sym1 = line.substring(59, 65).trim();
		String sym2 = null;
		if (len > 71) sym2 = line.substring(66, 72).trim();

		linkRecords.add(new LinkRecord(
				name1, altLoc1, resName1, chainID1, resSeq1, iCode1,
				name2, altLoc2, resName2, chainID2, resSeq2, iCode2,
				sym1, sym2));
	}

	/**
	 * Handler for the SITE records. <br>
	 *
	 * <pre>
	 *
	 * COLUMNS	DATA TYPE 		FIELD 		DEFINITION
	 * ---------------------------------------------------------------------------------
	 * 1 - 6	Record name 	"SITE "
	 * 8 - 10 	Integer 		seqNum 		Sequence number.
	 * 12 - 14 	LString(3)		siteID 		Site name.
	 * 16 - 17 	Integer 		numRes 		Number of residues that compose the siteResidues.
	 * 19 - 21 	Residue name	resName1	Residue name for first residue that
	 * 										creates the siteResidues.
	 * 23 		Character 		chainID1 	Chain identifier for first residue of siteResidues.
	 * 24 - 27 	Integer 		seq1 		Residue sequence number for first residue
	 * 										of the siteResidues.
	 * 28 		AChar 			iCode1 		Insertion code for first residue of the siteResidues.
	 *
	 * example:
	 *          1         2         3         4         5         6         7         8
	 * 12345678901234567890123456789012345678901234567890123456789012345678901234567890
	 * SITE     1 AC1  3 HIS A  94 HIS A   96  HIS A 119
	 * SITE     1 AC2  5 ASN A  62 GLY A   63  HIS A  64  HOH A 328
	 * SITE     2 AC2  5 HOH A 634
	 * SITE     1 AC3  5 GLN A 136 GLN A  137  PRO A 138  GLU A 205
	 * SITE     2 AC3  5 CYS A 206
	 * SITE     1 AC4 11 HIS A  64 HIS A   94  HIS A  96  HIS A 119
	 * SITE     2 AC4 11 LEU A 198 THR A  199  THR A 200  TRP A 209
	 * SITE     3 AC4 11 HOH A 572 HOH A  582  HOH A 635
	 * </pre>
	 * @param line the SITE line record being currently read
	 * @author Amr AL-Hossary
	 * @author Jules Jacobsen
	 */
	private void pdb_SITE_Handler(String line){

		if (params.isHeaderOnly()) return;

		//  make a map of: SiteId to List<ResidueNumber>

		logger.debug("Site Line:"+line);


		String siteID = line.substring(11, 14);
		//fetch the siteResidues from the map
		List<ResidueNumber> siteResidues = siteToResidueMap.get(siteID);

		//if the siteResidues doesn't yet exist, make a new one.
		if (siteResidues == null || ! siteToResidueMap.containsKey(siteID.trim())){
			siteResidues = new ArrayList<ResidueNumber>();
			siteToResidueMap.put(siteID.trim(), siteResidues);

			logger.debug(String.format("New Site made: %s %s", siteID,  siteResidues));
			logger.debug("Now made " + siteMap.size() + " sites");

		}

		logger.debug(String.format("SiteId: %s", siteID));


		//line = 'SITE     1 AC1  6 ARG H 221A LYS H 224  HOH H 403  HOH H 460'
		//line.substring(18) = 'ARG H 221A LYS H 224  HOH H 403  HOH H 460'
		line = line.substring(18);
		String groupString = null;
		//groupString = 'ARG H 221A'
		//keep iterating through chunks of 10 characters - these are the groups in the siteResidues
		while (!(groupString = line.substring(0, 10)).equals("          ")) {
			//groupstring: 'ARG H 221A'

			logger.debug("groupString: '" + groupString + "'");

			//set the residue name
			//residueName = 'ARG'
			String residueName = groupString.substring(0, 3);
			Character aminoCode1 = StructureTools.get1LetterCode(residueName);
			if (aminoCode1 != null) {
				if (aminoCode1.equals(StructureTools.UNKNOWN_GROUP_LABEL)) {
					aminoCode1 = null;
				}
			}

			//this is already in the right format, so no need to fiddle with it...
			//pdbCode = 'H 221A'
			//                    String pdbCode = groupString.substring(4, 10).trim();
			String chainId = groupString.substring(4, 5);
			Integer resNum = Integer.valueOf(groupString.substring(5, 9).trim());
			Character insCode = groupString.substring(9, 10).charAt(0);
			//set insCode to null as a measure to prevent storing thousands of empty Strings
			//- the empty value is returned using Group.getInsCode()
			//                    if (insCode.equals(" ")) {
			//                        insCode = null;
			//                    }

			logger.debug(String.format("Site: %s: 'resName:%s resNum:%s insCode:%s'", siteID, residueName, resNum, insCode));

			//make a new resNum with the data - this will be linked up with a site later
			ResidueNumber residueNumber = new ResidueNumber();


			logger.debug("pdbCode: '" + resNum + insCode + "'");

			residueNumber.setChainName(chainId);
			residueNumber.setSeqNum(resNum);
			residueNumber.setInsCode(insCode);
			//add the resNum to the groups
			siteResidues.add(residueNumber);

			logger.debug("Adding residueNumber " + residueNumber + " to site " + siteID);

			line = line.substring(11);
		}

		logger.debug("Current SiteMap (contains "+ siteToResidueMap.keySet().size() + " sites):");
		for (String key : siteToResidueMap.keySet()) {
			logger.debug(key + " : " + siteToResidueMap.get(key));
		}

	}

	//Site variable related to parsing the REMARK 800 records.
	Site site;
	private void pdb_REMARK_800_Handler(String line){

		if (params.isHeaderOnly()) return;

		// 'REMARK 800 SITE_IDENTIFIER: CAT                                                 '
		line = line.substring(11);
		String[] fields = line.split(": ");

		if (fields.length == 2) {
			if (fields[0].equals("SITE_IDENTIFIER")) {
				//                    remark800Counter++;
				String siteID = fields[1].trim();

				logger.debug("siteID: '" + siteID +"'");

				//fetch the siteResidues from the map
				site = siteMap.get(siteID);

				//if the siteResidues doesn't yet exist, make a new one.
				if (site == null || !siteID.equals(site.getSiteID())) {
					site = new Site(siteID, new ArrayList<Group>());
					siteMap.put(site.getSiteID(), site);

					logger.debug("New Site made: " + site);
					logger.debug("Now made " + siteMap.size() + " sites");

				}
			}
			if (fields[0].equals("EVIDENCE_CODE")) {
				//                    remark800Counter++;
				String evCode = fields[1].trim();

				logger.debug("evCode: '" + evCode +"'");

				//fetch the siteResidues from the map
				site.setEvCode(evCode);
			}
			if (fields[0].equals("SITE_DESCRIPTION")) {
				//                    remark800Counter++;
				String desc = fields[1].trim();

				logger.debug("desc: '" + desc +"'");

				//fetch the siteResidues from the map
				site.setDescription(desc);

				logger.debug("Finished making REMARK 800 for site " + site.getSiteID());
				logger.debug(site.remark800toPDB());

			}
		}
	}

	private int intFromString(String intString){
		int val = Integer.MIN_VALUE;
		try {
			val = Integer.parseInt(intString.trim());
		} catch (NumberFormatException ex){
			logger.info("Could not parse a number: " + ex.getMessage());
		}
		return val;
	}



	/** 
	 * Finds in the given list of chains the first one that has as name the given chainID.
	 * If no such Chain can be found it returns null.
	 */
	private static Chain isKnownChain(String chainID, List<Chain> chains){

		for (int i = 0; i< chains.size();i++){
			Chain testchain =  chains.get(i);
			if (chainID.equals(testchain.getName())) {
				return testchain;
			}
		}

		return null;
	}



	private BufferedReader getBufferedReader(InputStream inStream)
			throws IOException {

		BufferedReader buf ;
		if (inStream == null) {
			throw new IOException ("input stream is null!");
		}

		buf = new BufferedReader (new InputStreamReader (inStream));
		return buf ;

	}



	/**
	 * Parse a PDB file and return a datastructure implementing
	 * PDBStructure interface.
	 *
	 * @param inStream  an InputStream object
	 * @return a Structure object
	 * @throws IOException
	 */
	public Structure parsePDBFile(InputStream inStream)
			throws IOException
	{

		BufferedReader buf = getBufferedReader(inStream);

		return parsePDBFile(buf);

	}

	/**
	 * Parse a PDB file and return a datastructure implementing
	 * PDBStructure interface.
	 *
	 * @param buf  a BufferedReader object
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public  Structure parsePDBFile(BufferedReader buf)
			throws IOException
	{
		// set the correct max values for parsing...
		loadMaxAtoms = params.getMaxAtoms();
		atomCAThreshold = params.getAtomCaThreshold();


		// (re)set structure

		allModels = new ArrayList<>();
		structure     = new StructureImpl() ;
		currentModel  = null;
		currentChain  = null;
		currentGroup  = null;
		// we initialise to true since at the beginning of the file we are always starting a new molecule 
		startOfMolecule = true;
		startOfModel = true;

		seqResChains  = new ArrayList<Chain>();
		siteMap = new LinkedHashMap<String, Site>();
		pdbHeader     = new PDBHeader();
		connects      = new ArrayList<Map<String,Integer>>();
		previousContinuationField = "";
		continuationField = "";
		continuationString = "";
		current_compound = null;
		sourceLines.clear();
		compndLines.clear();
		isLastCompndLine = false;
		isLastSourceLine = false;
		prevMolId = -1;
		entities.clear();
		helixList.clear();
		strandList.clear();
		turnList.clear();
		lengthCheck = -1;
		atomCount = 0;
		atomOverflow = false;
		linkRecords = new ArrayList<LinkRecord>();
		siteToResidueMap.clear();
		
		blankChainIdsPresent = false;

		parseCAonly = params.isParseCAOnly();

		String line = null;

		while ((line = buf.readLine()) != null) {

			// ignore empty lines
			if ( line.equals("") ||
					(line.equals(NEWLINE))){
				continue;
			}


			// ignore short TER and END lines
			if ( line.startsWith("END")) {
				continue;
			}

			if ( line.length() < 6 && !line.startsWith("TER")) {
				logger.info("Found line length below 6. Ignoring it, line: >" + line +"<" );
				continue;
			}

			String recordName = null;
			if (line.length()<6)
				recordName = line.trim();
			else
				recordName = line.substring (0, 6).trim ();

			try {
				if (recordName.equals("ATOM"))
					pdb_ATOM_Handler(line);
				else if (recordName.equals("SEQRES"))
					pdb_SEQRES_Handler(line);
				else if (recordName.equals("HETATM"))
					pdb_ATOM_Handler(line);
				else if (recordName.equals("MODEL"))
					pdb_MODEL_Handler(line);
				else if (recordName.equals("TER"))
					pdb_TER_Handler();
				else if (recordName.equals("HEADER"))
					pdb_HEADER_Handler(line);
				else if (recordName.equals("AUTHOR"))
					pdb_AUTHOR_Handler(line);
				else if (recordName.equals("TITLE"))
					pdb_TITLE_Handler(line);
				else if (recordName.equals("SOURCE"))
					sourceLines.add(line); //pdb_SOURCE_Handler
				else if (recordName.equals("COMPND"))
					compndLines.add(line); //pdb_COMPND_Handler
				else if (recordName.equals("JRNL"))
					pdb_JRNL_Handler(line);
				else if (recordName.equals("EXPDTA"))
					pdb_EXPDTA_Handler(line);
				else if (recordName.equals("CRYST1"))
					pdb_CRYST1_Handler(line);
				else if (recordName.startsWith("MTRIX"))
					pdb_MTRIXn_Handler(line);
				else if (recordName.equals("REMARK"))
					pdb_REMARK_Handler(line);
				else if (recordName.equals("CONECT"))
					pdb_CONECT_Handler(line);
				else if (recordName.equals("REVDAT"))
					pdb_REVDAT_Handler(line);
				else if (recordName.equals("DBREF"))
					pdb_DBREF_Handler(line);
				else if (recordName.equals("SITE"))
					pdb_SITE_Handler(line);
				else if (recordName.equals("SSBOND"))
					pdb_SSBOND_Handler(line);
				else if (recordName.equals("LINK"))
					pdb_LINK_Handler(line);
				else if ( params.isParseSecStruc()) {
					if ( recordName.equals("HELIX") ) pdb_HELIX_Handler (  line ) ;
					else if (recordName.equals("SHEET")) pdb_SHEET_Handler(line ) ;
					else if (recordName.equals("TURN")) pdb_TURN_Handler(   line ) ;
				}
			} catch (StringIndexOutOfBoundsException | NullPointerException ex) {
				logger.info("Unable to parse [" + line + "]");
			} 
		}

		makeCompounds(compndLines, sourceLines);

		triggerEndFileChecks();

		if (params.shouldCreateAtomBonds()) {
			formBonds();
		}

		if ( params.shouldCreateAtomCharges()) {
			addCharges();
		}

		if ( params.isParseSecStruc() && !params.isHeaderOnly())
			setSecStruc();

		// Now correct the alternate location group
		StructureTools.cleanUpAltLocs(structure);

		return structure;

	}


	/**
	 * Add the charges to the Structure
	 */
	private void addCharges() {
		ChargeAdder.addCharges(structure);
	}

	/**
	 * This is the new method for building the COMPND and SOURCE records. Now each method is self-contained.
	 * @author Jules Jacobsen
	 * @param  compoundList
	 * @param  sourceList
	 */
	private void makeCompounds(List<String> compoundList,
			List<String> sourceList) {
		//		System.out.println("[makeCompounds] making compounds from compoundLines");

		for (String line : compoundList) {
			if (compoundList.indexOf(line) + 1 == compoundList.size()) {
				//				System.out.println("[makeCompounds] Final line in compoundLines.");
				isLastCompndLine = true;
			}
			pdb_COMPND_Handler(line);

		}
		//		System.out.println("[makeCompounds] adding sources to compounds from sourceLines");
		// since we're starting again from the first compound, reset it here
		if ( entities.size() == 0){
			current_compound = new EntityInfo();
		} else {
			current_compound = entities.get(0);
		}
		for (String line : sourceList) {
			if (sourceList.indexOf(line) + 1 == sourceList.size()) {
				//				System.out.println("[makeCompounds] Final line in sourceLines.");
				isLastSourceLine = true;
			}
			pdb_SOURCE_Handler(line);
		}

	}

	/**
	 * Handles creation of all bonds. Looks at LINK records, SSBOND (Disulfide
	 * bonds), peptide bonds, and intra-residue bonds.
	 * <p>
	 * Note: the current implementation only looks at the first model of each
	 * structure. This may need to be fixed in the future.
	 */
	private void formBonds() {

		BondMaker maker = new BondMaker(structure, params);
		
		// LINK records should be preserved, they are the way that
		// inter-residue bonds are created for ligands such as trisaccharides, unusual polymers.
        // The analogy in mmCIF is the _struct_conn record.
		for (LinkRecord linkRecord : linkRecords) {
            maker.formLinkRecordBond(linkRecord);
        }

		maker.formDisulfideBonds(ssbonds);

		maker.makeBonds();
	}



	private void triggerEndFileChecks(){

		// we need to add the last chain and model, checking for nulls (e.g. the file could be completely empty of ATOM lines)
		if (currentChain!=null && currentGroup!=null) {
			currentChain.addGroup(currentGroup);
		}
		if (currentModel!=null && currentChain!=null) {
			currentModel.add(currentChain);
		}
		if (currentModel!=null) {
			allModels.add(currentModel);
		}
		
		if (blankChainIdsPresent) {
			// from biojava 5.0 there's limited support for old pdb files with blank chain ids
			logger.warn("Found some blank chain ids in PDB file. Please note that support for them has been discontinued and things might not work properly.");
		}

		// reordering chains following the mmcif model and assigning entities
		assignChainsAndEntities();
		structure.setEntityInfos(entities);
		

		
		// header data
		
		Date modDate = pdbHeader.getModDate();
		if ( modDate.equals(new Date(0)) ) {
			// modification date = deposition date
			Date depositionDate = pdbHeader.getDepDate();

			if (! depositionDate.equals(modDate)){
				// depDate is 0000-00-00
				pdbHeader.setDepDate(depositionDate);
			}

		}
		
		structure.setPDBHeader(pdbHeader);
		structure.setCrystallographicInfo(crystallographicInfo);

		//set the JournalArticle, if there is one
		if (!journalLines.isEmpty()) {
			buildjournalArticle();
			pdbHeader.setJournalArticle(journalArticle);
		}
		
		structure.setDBRefs(dbrefs);

		// Only align if requested (default) and not when headerOnly mode with no Atoms.
		// Otherwise, we store the empty SeqRes Groups unchanged in the right chains.
		if ( params.isAlignSeqRes() && !params.isHeaderOnly() && !seqResChains.isEmpty()){
			logger.debug("Parsing mode align_seqres, will parse SEQRES and align to ATOM sequence");
			SeqRes2AtomAligner aligner = new SeqRes2AtomAligner();
			aligner.align(structure,seqResChains);

		} else {
			logger.debug("Parsing mode unalign_seqres, will parse SEQRES but not align it to ATOM sequence");
			SeqRes2AtomAligner.storeUnAlignedSeqRes(structure, seqResChains, params.isHeaderOnly());
		}


		
		//associate the temporary Groups in the siteMap to the ones
		if (!params.isHeaderOnly()) {
			// Only can link SITES if Atom Groups were parsed.
			linkSitesToGroups(); // will work now that setSites is called
		}

		if ( bioAssemblyParser != null){
			bioAssemblyParser.setMacromolecularSizes();
			pdbHeader.setBioAssemblies(bioAssemblyParser.getTransformationMap());
		}

		if (ncsOperators !=null && ncsOperators.size()>0) {
			crystallographicInfo.setNcsOperators(
				ncsOperators.toArray(new Matrix4d[ncsOperators.size()]));
		}


		// rfree end file check
		// Rfree annotation is not very consistent in PDB format, it varies depending on the software
		// Here we follow this strategy:
		// a) take the '(NO CUTOFF)' value if the only one available (shelx software, e.g. 1x7q)
		// b) don't take it if also a line without '(NO CUTOFF)' is present (CNX software, e.g. 3lak)

		if (rfreeNoCutoffLine>0 && rfreeStandardLine<0) {
			pdbHeader.setRfree(rfreeNoCutoffLine);
		} else if (rfreeNoCutoffLine>0 && rfreeStandardLine>0) {
			pdbHeader.setRfree(rfreeStandardLine);
		} else if (rfreeNoCutoffLine<0 && rfreeStandardLine>0) {
			pdbHeader.setRfree(rfreeStandardLine);
		} // otherwise it remains default value: PDBHeader.DEFAULT_RFREE


		
	}

	private void setSecStruc(){

		setSecElement(helixList, SecStrucInfo.PDB_AUTHOR_ASSIGNMENT,
				SecStrucType.helix4);
		setSecElement(strandList, SecStrucInfo.PDB_AUTHOR_ASSIGNMENT,
				SecStrucType.extended);
		setSecElement(turnList, SecStrucInfo.PDB_AUTHOR_ASSIGNMENT,
				SecStrucType.turn);

		//Now insert random coil to the Groups that did not have SS information
		GroupIterator gi = new GroupIterator(structure);
		while (gi.hasNext()){
			Group g = gi.next();
			if (g.hasAminoAtoms()){
				if (g.getProperty(Group.SEC_STRUC) == null){
					SecStrucInfo ss = new SecStrucInfo(g,
							SecStrucInfo.PDB_AUTHOR_ASSIGNMENT,
							SecStrucType.coil);
					g.setProperty(Group.SEC_STRUC, ss);
				}
			}
		}

	}

	private void setSecElement(List<Map<String,String>> secList, String assignment, SecStrucType type){


		Iterator<Map<String,String>> iter = secList.iterator();
		nextElement:
			while (iter.hasNext()){
				Map<String,String> m = iter.next();

				// assign all residues in this range to this secondary structure type
				// String initResName = (String)m.get("initResName");
				String initChainId = m.get("initChainId");
				String initSeqNum  = m.get("initSeqNum" );
				String initICode   = m.get("initICode" );
				// String endResName  = (String)m.get("endResName" );
				String endChainId  = m.get("endChainId" );
				String endSeqNum   = m.get("endSeqNum");
				String endICode    = m.get("endICode");

				if (initICode.equals(" "))
					initICode = "";
				if (endICode.equals(" "))
					endICode = "";

				GroupIterator gi = new GroupIterator(structure);
				boolean inRange = false;
				while (gi.hasNext()){
					Group g = gi.next();
					Chain c = g.getChain();

					if (c.getName().equals(initChainId)){

						String pdbCode = initSeqNum + initICode;
						if ( g.getResidueNumber().toString().equals(pdbCode)  ) {
							inRange = true;
						}
					}
					if ( inRange){
						if (g.hasAminoAtoms()) {
							SecStrucInfo ss = new SecStrucInfo(g, assignment, type);
							g.setProperty(Group.SEC_STRUC, ss);
						}

					}
					if ( c.getName().equals(endChainId)){
						String pdbCode = endSeqNum + endICode;
						if (pdbCode.equals(g.getResidueNumber().toString())){
							inRange = false;
							continue nextElement;
						}
					}
				}
			}
	}

	/**
	 * Gets all chains with given chainName from given models list
	 * @param chainName
	 * @param polyModels
	 * @return
	 */
	private static List<List<Chain>> findChains(String chainName, List<List<Chain>> polyModels) {
		List<List<Chain>> models = new ArrayList<>();

		for (List<Chain> chains:polyModels) {
			List<Chain> matchingChains = new ArrayList<>();
			models.add(matchingChains);
			for (Chain c:chains) {
				if (c.getName().equals(chainName)) {
					matchingChains.add(c);
				}
			}
		}
		return models;
	}
	
	/**
	 * Split the given chain (containing non-polymer groups and water groups only) 
	 * into individual chains per non-polymer group and individual chains per contiguous sets of water groups. 
	 * @param chain
	 * @return a list of lists of size 2: first list is the split non-poly chains, second list is the split water chains 
	 */
	private static List<List<Chain>> splitNonPolyChain(Chain chain) {
		List<Chain> splitNonPolys = new ArrayList<>();
		List<Chain> waterChains = new ArrayList<>();
		
		Chain split = null;
		boolean previousGroupIsWater = false;
		
		for (Group g:chain.getAtomGroups()){
			
			if (!previousGroupIsWater) {
				// add last one if there's one
				if (split!=null) {
					splitNonPolys.add(split);
				}
				split = new ChainImpl();
				split.setName(chain.getName());
			} else if (!g.isWater()) { 
				// previous group is water and this group is not water: we change from a water chain to a non-poly
				// we'll need to add now the water chain to the list of water chains
				waterChains.add(split);
				split = new ChainImpl();
				split.setName(chain.getName());
			}
			
			if (g.isWater()) {
				previousGroupIsWater = true;
			} else {
				previousGroupIsWater = false;
				
			}
						
			// this should include alt locs (referenced from the main group)
			split.addGroup(g);
			
		}
		
		// adding the last split chain: either to water or non-poly depending on what was the last seen group
		if (split!=null) {
			if (previousGroupIsWater)
				waterChains.add(split);
			else
				splitNonPolys.add(split);
		}

		
		List<List<Chain>> all = new ArrayList<>(2);
		all.add(splitNonPolys);
		all.add(waterChains);

		return all;
	}
	
	/**
	 * Assign asym ids following the rules used by the PDB to assign asym ids in mmCIF files
	 * @param polys
	 * @param nonPolys
	 * @param waters
	 */
	private void assignAsymIds(List<List<Chain>> polys, List<List<Chain>> nonPolys, List<List<Chain>> waters) {
		
		for (int i=0; i<polys.size(); i++) {
			String asymId = "A";		

			for (Chain poly:polys.get(i)) {
				poly.setId(asymId);
				asymId = getNextAsymId(asymId);
			}
			for (Chain nonPoly:nonPolys.get(i)) {
				nonPoly.setId(asymId);
				asymId = getNextAsymId(asymId);			
			}
			for (Chain water:waters.get(i)) {
				water.setId(asymId);
				asymId = getNextAsymId(asymId);			
			}
		}
	}
	
	/**
	 * Gets the next asym id given an asymId, according to the convention followed by 
	 * mmCIF files produced by the PDB
	 * i.e.: A,B,...,Z,AA,BA,CA,...,ZA,AB,BB,CB,...,ZB,.......,ZZ,AAA,BAA,CAA,...
	 * @param asymId
	 * @return
	 */
	private String getNextAsymId(String asymId) {
		if (asymId.length()==1) {
			if (!asymId.equals("Z")) {
				return Character.toString(getNextChar(asymId.charAt(0)));
			} else {
				return "AA";
			}
		} else if (asymId.length()==2) {
			if (asymId.equals("ZZ")) {
				return "AAA";
			}
			char[] c = new char[2];
			asymId.getChars(0, 2, c, 0);
			c[0] = getNextChar(c[0]);
			if (c[0]=='A') {
				c[1] = getNextChar(c[1]);
			} 
			return new String(c);
		} else if (asymId.length()==3) {
			char[] c = new char[3];
			asymId.getChars(0, 3, c, 0);
			c[0] = getNextChar(c[0]);
			if (c[0]=='A') {
				c[1] = getNextChar(c[1]);
				if (c[1]=='A') {
					c[2] = getNextChar(c[2]);
				}
			}
			return new String(c);
		}
		return null;
	}
	
	private char getNextChar(char c) {
		if (c!='Z') {
			return ((char)(c+1));
		} else {
			return 'A';
		}
	}
	
	/** 
	 * Here we assign chains following the mmCIF data model:
	 * one chain per polymer, one chain per non-polymer group and 
	 * several water chains.
	 * <p>
	 * Subsequently we assign entities for them: either from those read from 
	 * COMPOUND records or from those found heuristically through {@link EntityFinder} 
	 *
	 */
	private void assignChainsAndEntities(){
		
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
		
		List<List<Chain>> splitNonPolyModels = new ArrayList<>();
		for (int i=0; i<nonPolyModels.size(); i++) {
			List<Chain> nonPolyModel = nonPolyModels.get(i);
			List<Chain> waterModel = waterModels.get(i);
			
			List<Chain> splitNonPolys = new ArrayList<>();
			splitNonPolyModels.add(splitNonPolys);
			
			for (Chain nonPoly:nonPolyModel) {
				List<List<Chain>> splits = splitNonPolyChain(nonPoly);
				splitNonPolys.addAll(splits.get(0));
				waterModel.addAll(splits.get(1));
			}
		}
		
		
		// now we have all chains as in mmcif, let's assign ids following the mmcif rules
		assignAsymIds(polyModels, splitNonPolyModels, waterModels);
		

		if (!entities.isEmpty()) {
			// if the file contained COMPOUND records then we can assign entities to the poly chains
			for (EntityInfo comp : entities){
				List<String> chainIds = compoundMolIds2chainIds.get(comp.getMolId());
				if ( chainIds == null)
					continue;
				for ( String chainId : chainIds) {

					List<List<Chain>> models = findChains(chainId, polyModels);

					for (List<Chain> matchingChains:models) {
						for (Chain chain:matchingChains) {
							comp.addChain(chain);
							chain.setEntityInfo(comp);
						}

						if (matchingChains.isEmpty()) {
							// usually if this happens something is wrong with the PDB header
							// e.g. 2brd - there is no Chain A, although it is specified in the header
							// Some bona-fide cases exist, e.g. 2ja5, chain N is described in SEQRES
							// but the authors didn't observe in the density so it's completely missing
							// from the ATOM lines
							logger.warn("Could not find polymeric chain {} to link to entity {}. The chain will be missing in the entity.", chainId, comp.getMolId());
						}
					}
				}
			}
			
		} else {

			logger.info("Entity information (COMPOUND record) not found in file. Will assign entities heuristically");
			// if no entity information was present in file we then go and find the entities heuristically with EntityFinder
			entities = EntityFinder.findPolyEntities(polyModels);

		}
		
		// now we assign entities to the nonpoly and water chains
		EntityFinder.createPurelyNonPolyEntities(splitNonPolyModels, waterModels, entities);


		// in some rare cases purely non-polymer or purely water chain are present in pdb files
		// see https://github.com/biojava/biojava/pull/394
		// these case should be covered by the above

		
		// now that we have entities in chains we add the chains to the structure
		
		for (int i=0;i<allModels.size();i++) {
			List<Chain> model = new ArrayList<>();
			model.addAll(polyModels.get(i));
			model.addAll(splitNonPolyModels.get(i));
			model.addAll(waterModels.get(i));
			structure.addModel(model);
		}


	}
	
	/**
	 * Links the Sites in the siteMap to the Groups in the Structure via the
	 * siteToResidueMap ResidueNumber.
	 * @author Jules Jacobsen
	 * @return
	 */
	private void linkSitesToGroups() {

		//System.out.println("LINK SITES TO GROUPS:" + siteToResidueMap.keySet().size());

		//link the map of siteIds : <ResidueNumber> with the sites by using ResidueNumber to get the correct group back.
		//the return list

		if ( siteMap == null || siteToResidueMap == null){
			logger.info("Sites can not be linked to residues!");

			return;
		}

		List<Site> sites = null;
		//check that there are chains with which to associate the groups
		if (structure.getChains().isEmpty()) {
			sites = new ArrayList<Site>(siteMap.values());
			logger.info("No chains to link Site Groups with - Sites will not be present in the Structure");
			return;
		}

		//check that the keys in the siteMap and SiteToResidueMap are equal
		if (! siteMap.keySet().equals(siteToResidueMap.keySet())) {
			logger.info("Not all sites have been properly described in the PDB " + pdbId + " header - some Sites will not be present in the Structure");
			logger.debug(siteMap.keySet() + " | " + siteToResidueMap.keySet());
			//return;
		}

		//so we have chains - associate the siteResidues-related groups with the ones
		//already in in the chains
		for (String key : siteMap.keySet()) {
			Site currentSite = siteMap.get(key);
			List<ResidueNumber> linkedGroups = siteToResidueMap.get(key);
			if ( linkedGroups == null)
				continue;
			for (ResidueNumber residueNumber : linkedGroups) {

				String pdbCode = residueNumber.toString();
				String chain = residueNumber.getChainName();
				//                    System.out.println("chain: '" + chain + "'");
				//                    String resNum = resNum.getSeqNum().toString();
				//                    System.out.println("resNum: '" + resNum + "'");

				Group linkedGroup = null;
				try {
					//TODO: implement findGroup(ResidueNumber resNum)
					linkedGroup = structure.findGroup(chain, pdbCode);
				} catch (StructureException ex) {
					logger.info("Can't find group " + pdbCode + " in chain " + chain + " in order to link up SITE records (PDB ID " + pdbId +")");
					continue;
				}

				//                    System.out.println("Adding group: " + linkedGroup.getSeqNum() + " to site " + site.getSiteID());
				currentSite.getGroups().add(linkedGroup);
			}
		}

		//System.out.println("SITEMAP: " + siteMap);

		sites = new ArrayList<Site>(siteMap.values());
		structure.setSites(sites);
		//System.out.println("STRUCTURE SITES: " + structure.getSites().size());
		//            for (Site site : structure.getSites()) {
		//                System.out.println(site);
		//            }
		//            System.out.println("Linked Site Groups with Chains");

	}

	private void buildjournalArticle() {

		logger.debug("building new JournalArticle");
		//            for (String line : journalLines) {
		//                System.out.println(line);
		//            }

		this.journalArticle = new JournalArticle();
		//        JRNL        AUTH   M.HAMMEL,G.SFYROERA,D.RICKLIN,P.MAGOTTI,
		//        JRNL        AUTH 2 J.D.LAMBRIS,B.V.GEISBRECHT
		//        JRNL        TITL   A STRUCTURAL BASIS FOR COMPLEMENT INHIBITION BY
		//        JRNL        TITL 2 STAPHYLOCOCCUS AUREUS.
		//        JRNL        REF    NAT.IMMUNOL.                  V.   8   430 2007
		//        JRNL        REFN                   ISSN 1529-2908
		//        JRNL        PMID   17351618
		//        JRNL        DOI    10.1038/NI1450
		StringBuffer auth = new StringBuffer();
		StringBuffer titl = new StringBuffer();
		StringBuffer edit = new StringBuffer();
		StringBuffer ref = new StringBuffer();
		StringBuffer publ = new StringBuffer();
		StringBuffer refn = new StringBuffer();
		StringBuffer pmid = new StringBuffer();
		StringBuffer doi = new StringBuffer();

		for (String line : journalLines) {
			if ( line.length() < 19 ) {
				logger.info("can not process Journal line: " + line);
				continue;
			}
			//            System.out.println("'" + line + "'");
			String subField = line.substring(12, 16);
			//            System.out.println("'" + subField + "'");
			if (subField.equals("AUTH")) {
				auth.append(line.substring(19, line.length()).trim());

				logger.debug("AUTH '" + auth.toString() + "'");

			}
			if (subField.equals("TITL")) {
				//add a space to the end of a line so that when wrapped the
				//words on the join won't be concatenated
				titl.append(line.substring(19, line.length()).trim()).append(" ");

				logger.debug("TITL '" + titl.toString() + "'");

			}
			if (subField.equals("EDIT")) {
				edit.append(line.substring(19, line.length()).trim());

				logger.debug("EDIT '" + edit.toString() + "'");

			}
			//        JRNL        REF    NAT.IMMUNOL.                  V.   8   430 2007
			if (subField.equals("REF ")) {
				ref.append(line.substring(19, line.length()).trim()).append(" ");

				logger.debug("REF '" + ref.toString() + "'");

			}
			if (subField.equals("PUBL")) {
				publ.append(line.substring(19, line.length()).trim()).append(" ");

				logger.debug("PUBL '" + publ.toString() + "'");

			}
			//        JRNL        REFN                   ISSN 1529-2908
			if (subField.equals("REFN")) {
				if ( line.length() < 35 ) {
					logger.info("can not process Journal REFN line: " + line);
					continue;
				}
				refn.append(line.substring(35, line.length()).trim());

				logger.debug("REFN '" + refn.toString() + "'");

			}
			//        JRNL        PMID   17351618
			if (subField.equals("PMID")) {
				pmid.append(line.substring(19, line.length()).trim());

				logger.debug("PMID '" + pmid.toString() + "'");

			}
			//        JRNL        DOI    10.1038/NI1450
			if (subField.equals("DOI ")) {
				doi.append(line.substring(19, line.length()).trim());

				logger.debug("DOI '" + doi.toString() + "'");

			}
		}

		//now set the parts of the JournalArticle
		journalArticle.setAuthorList(authorBuilder(auth.toString()));
		journalArticle.setEditorList(authorBuilder(edit.toString()));
		journalArticle.setRef(ref.toString());
		JournalParser journalParser = new JournalParser(ref.toString());
		journalArticle.setJournalName(journalParser.getJournalName());
		if (!journalArticle.getJournalName().equals("TO BE PUBLISHED")) {
			journalArticle.setIsPublished(true);
		}
		journalArticle.setVolume(journalParser.getVolume());
		journalArticle.setStartPage(journalParser.getStartPage());
		journalArticle.setPublicationDate(journalParser.getPublicationDate());
		journalArticle.setPublisher(publ.toString().trim());
		journalArticle.setTitle(titl.toString().trim());
		journalArticle.setRefn(refn.toString().trim());
		journalArticle.setPmid(pmid.toString().trim());
		journalArticle.setDoi(doi.toString().trim());


		logger.debug("Made JournalArticle:");
		logger.debug(journalArticle.toString());

	}

	//inner class to deal with all the journal info
	private class JournalParser {

		private String journalName;
		private String volume;
		private String startPage;
		private int publicationDate;


		public JournalParser(String ref) {

			logger.debug("JournalParser init '" + ref + "'");


			if (ref.equals("TO BE PUBLISHED ")) {
				journalName = ref.trim();

				logger.debug(String.format("JournalParser found journalString '%s'", journalName));

				return;
			}

			if (ref.length() < 48) {
				logger.info("REF line too short - must be at least 48 characters to be valid for parsing.");
				journalName = "";
				volume = "";
				startPage = "";
				publicationDate = 0;
				return;
			}
			//can be multi line:
			//REF    PHILOS.TRANS.R.SOC.LONDON,    V. 293    53 1981
			//REF  2 SER.B

			//or

			//REF    GLYCOGEN PHOSPHORYLASE B:                1 1991
			//REF  2 DESCRIPTION OF THE PROTEIN
			//REF  3 STRUCTURE

			//but usually single line
			//REF    NUCLEIC ACIDS RES.                         2009
			//REF    MOL.CELL                                   2009
			//REF    NAT.STRUCT.MOL.BIOL.          V.  16   238 2009
			//REF    ACTA CRYSTALLOGR.,SECT.F      V.  65   199 2009
			//check if the date is present at the end of the line.
			//                             09876543210987654321
			//'J.BIOL.CHEM.                  V. 280 23000 2005 '
			//'J.AM.CHEM.SOC.                V. 130 16011 2008 '
			//'NAT.STRUCT.MOL.BIOL.          V.  16   238 2009'
			String volumeInformation = ref.substring(30, 48);

			logger.debug(String.format("Parsing volumeInformation: '%s'", volumeInformation));

			//volumeInformation: 'V. 293    53 1981 '
			//                      String dateString = ref.substring(ref.length() - 5 , ref.length() - 1).trim();
			//			String startPageString = ref.substring(ref.length() - 11 , ref.length() - 6).trim();
			//			String volumeString = ref.substring(ref.length() - 16 , ref.length() - 12).trim();
			//			String journalString = ref.substring(0 , ref.length() - 18).trim();
			String dateString = volumeInformation.substring(volumeInformation.length() - 5 , volumeInformation.length() - 1).trim();
			String startPageString = volumeInformation.substring(volumeInformation.length() - 11 , volumeInformation.length() - 6).trim();
			String volumeString = volumeInformation.substring(volumeInformation.length() - 16 , volumeInformation.length() - 12).trim();
			//for the journal string we need to remove the volume information which might be in the middle of the string (e.g. 1gpb, 3pfk)
			String journalString = ref.substring(0 , 29).trim() + " " + ref.substring(30, ref.length() - 1).replace(volumeInformation.trim(), "").trim();
			journalString = journalString.trim();
			//                        System.out.println("journalString: " + journalString);

			logger.debug(String.format("JournalParser found volumeString '%s'", volumeString));
			logger.debug(String.format("JournalParser found startPageString '%s'", startPageString));
			logger.debug(String.format("JournalParser found dateString '%s'", dateString));
			logger.debug(String.format("JournalParser found journalString '%s'", journalString));


			if (!dateString.equals("    ")) {
				try {
					publicationDate = Integer.valueOf(dateString);
				} catch (NumberFormatException nfe) {
					logger.info(dateString + " is not a valid integer for a date in JRNL sub-section REF line 1");
				}
				//				if (DEBUG) {
				//					System.out.println("JournalParser set date " + publicationDate);
				//				}
			}

			if (!startPageString.equals("    ")) {
				startPage = startPageString;
				//				if (DEBUG) {
				//					System.out.println("JournalParser set startPage " + startPage);
				//				}
			}

			if (!volumeString.equals("    ")) {
				volume = volumeString;
				//				if (DEBUG) {
				//					System.out.println("JournalParser set volume " + volume);
				//				}
			}

			if (!journalString.equals("    ")) {
				journalName = journalString;

				logger.debug("JournalParser set journalName " + journalName);

			}
		}

		private String getJournalName() {
			return journalName;
		}

		private int getPublicationDate() {
			return publicationDate;
		}

		private String getStartPage() {
			return startPage;
		}

		private String getVolume() {
			return volume;
		}
	}

	private List<Author> authorBuilder(String authorString) {
		ArrayList<Author> authorList = new ArrayList<Author>();

		if (authorString.equals("")) {
			return authorList;
		}

		String[] authors = authorString.split(",");
		//        if (DEBUG) {
		//            for (int i = 0; i < authors.length; i++) {
		//                String string = authors[i];
		//                System.out.println("authorBuilder author: '" + string + "'");
		//            }
		//        }
		//        AUTH   SEATTLE STRUCTURAL GENOMICS CENTER FOR INFECTIOUS
		//        AUTH 2 DISEASE (SSGCID)
		//        or
		//        AUTH   E.DOBROVETSKY,A.DONG,A.SEITOVA,B.DUNCAN,L.CROMBET,
		//        AUTH 2 M.SUNDSTROM,C.H.ARROWSMITH,A.M.EDWARDS,C.BOUNTRA,
		//        AUTH 3 A.BOCHKAREV,D.COSSAR,
		//        AUTH 4 STRUCTURAL GENOMICS CONSORTIUM (SGC)
		//        or
		//        AUTH   T.-C.MOU,S.R.SPRANG,N.MASADA,D.M.F.COOPER
		if (authors.length == 1) {
			//only one element means it's a consortium only
			Author author = new Author();
			author.setSurname(authors[0]);

			logger.debug("Set consortium author name " + author.getSurname());

			authorList.add(author);
		} else {
			for (int i = 0; i < authors.length; i++) {
				String authorFullName = authors[i];

				logger.debug("Building author " + authorFullName);

				Author author = new Author();
				String regex = "\\.";
				String[] authorNames = authorFullName.split(regex);
				//                if (DEBUG) {
				//                    System.out.println("authorNames size " + authorNames.length);
				//                    for (int j = 0; j < authorNames.length; j++) {
				//                        String name = authorNames[j];
				//                        System.out.println("split authName '" + name + "'");
				//
				//                    }
				//                }
				if (authorNames.length == 0) {
					author.setSurname(authorFullName);

					logger.debug("Unable to split using '" + regex + "' Setting whole name " + author.getSurname());

				}
				//again there might be a consortium name so there may be no elements
				else if (authorNames.length == 1) {
					author.setSurname(authorNames[0]);

					logger.debug("Set consortium author name in multiple author block " + author.getSurname
								());

				} else {
					String initials = "";
					for (int j = 0; j < authorNames.length - 1; j++) {
						String initial = authorNames[j];
						//                        if (DEBUG) {
						//                            System.out.println("adding initial '" + initial + "'");
						//                        }
						//build the initials back up again
						initials += initial + ".";
					}

					logger.debug("built initials '" + initials + "'");

					author.setInitials(initials);
					//surname is always last
					int lastName = authorNames.length - 1;
					String surname = authorNames[lastName];

					logger.debug("built author surname " + surname);

					author.setSurname(surname);

				}
				authorList.add(author);
			}
		}
		return authorList;
	}

	public void setFileParsingParameters(FileParsingParameters params)
	{
		this.params= params;

		// set the correct max values for parsing...
		loadMaxAtoms = params.getMaxAtoms();
		atomCAThreshold = params.getAtomCaThreshold();

	}

	public FileParsingParameters getFileParsingParameters(){
		return params;
	}


}
