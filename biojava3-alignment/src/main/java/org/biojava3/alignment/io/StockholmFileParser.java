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
 * Created on August 13, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.io;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import org.biojava3.alignment.io.StockholmFileAnnotation.StockholmFileAnnotationReference;
import org.biojava3.core.exceptions.ParserException;
import org.biojava3.core.util.InputStreamProvider;

/**
 * Stockholm file parser.<br>
 * for more information about the format refer to 
 * <ul>
 * <li><a href="ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/userman.txt">ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/userman.txt</a>.</li>
 * <li><a href="ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/USERMAN">ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/USERMAN</a>.</li>
 * <li><a href="http://sonnhammer.sbc.su.se/Stockholm.html">http://sonnhammer.sbc.su.se/Stockholm.html</a>.</li>
 * </ul>
 * <pre>
 * Pfam DESCRIPTION OF FIELDS

   Compulsory fields:
   ------------------

   AC   Accession number:           Accession number in form PFxxxxx.version or PBxxxxxx.
   ID   Identification:             One word name for family.
   DE   Definition:                 Short description of family.
   AU   Author:                     Authors of the entry.
   SE   Source of seed:             The source suggesting the seed members belong to one family.
   GA   Gathering method:           Search threshold to build the full alignment.
   TC   Trusted Cutoff:             Lowest sequence score and domain score of match in the full alignment.
   NC   Noise Cutoff:               Highest sequence score and domain score of match not in full alignment.
   TP   Type:                       Type of family -- presently Family, Domain, Motif or Repeat.
   SQ   Sequence:                   Number of sequences in alignment.
   //                               End of alignment.

   Optional fields:
   ----------------

   DC   Database Comment:           Comment about database reference.
   DR   Database Reference:         Reference to external database.
   RC   Reference Comment:          Comment about literature reference.
   RN   Reference Number:           Reference Number.
   RM   Reference Medline:          Eight digit medline UI number.
   RT   Reference Title:            Reference Title.
   RA   Reference Author:           Reference Author
   RL   Reference Location:         Journal location.
   PI   Previous identifier:        Record of all previous ID lines.
   KW   Keywords:                   Keywords.
   CC   Comment:                    Comments.
   NE   Pfam accession:             Indicates a nested domain.
   NL   Location:                   Location of nested domains - sequence ID, start and end of insert.
   WK   Wikipedia Reference:        Reference to wikipedia.

   Obsolete fields:
   -----------
   AL   Alignment method of seed:   The method used to align the seed members.
   AM   Alignment Method:	    The order ls and fs hits are aligned to the model to build the full align.

 * </pre>
 * 
 * @since 3.0.5
 * @author Amr AL-Hossary
 * @author Marko Vaz
 * 
 */
public class StockholmFileParser {

	/**indicates reading as much as possible, without limits */
	public static final int INFINITY = -1;
	/** #=GF &lt;feature&gt; &lt;Generic per-File annotation, free text&gt; */
	private static final String GENERIC_PER_FILE_ANNOTATION = "GF";
	/** #=GC &lt;feature&gt; &lt;Generic per-Column annotation, exactly 1 char per column&gt;  */
	private static final String GENERIC_PER_CONSENSUS_ANNOTATION = "GC";
	/** #=GS &lt;seqname&gt; &lt;feature&gt; &lt;Generic per-Sequence annotation, free text&gt; */
	private static final String GENERIC_PER_SEQUENCE_ANNOTATION = "GS";
	/** #=GR &lt;seqname&gt; &lt;feature&gt; &lt;Generic per-Residue annotation, exactly 1 char per residue&gt; */
	private static final String GENERIC_PER_RESIDUE_ANNOTATION = "GR";

	// COMPULSORY FIELDS
	/** Accession number in form PFxxxxx (Pfam) or RFxxxxx (Rfam).*/
	private static final String GF_ACCESSION_NUMBER = "AC";
	/** One word name for family.*/
	private static final String GF_IDENTIFICATION = "ID";
	/** Short description of family. */
	private static final String GF_DEFINITION = "DE";
	/** Authors of the entry.*/
	private static final String GF_AUTHOR = "AU";
	/**Indicates the order that ls and fs matches are aligned to 
	 * the model to give the full alignment. (OBSOLETE IN HMMER3)*/
	private static final String GF_ALIGNMENT_METHOD = "AM";
	/**	Command line used to generate the model*/
	private static final String GF_BUILD_METHOD = "BM";
	/** Command line used to perform the search*/
	private static final String GF_SEARCH_METHOD = "SM";
	/**	The source suggesting the seed members belong to one family.*/
	private static final String GF_SOURCE_SEED = "SE";
	/**	The source (prediction or publication) of the consensus RNA secondary structure used by Rfam. */
	private static final String GF_SOURCE_STRUCTURE = "SS";
	/** Search threshold to build the full alignment.*/
	private static final String GF_GATHERING_THRESHOLD = "GA";
	/** Lowest sequence score (and domain score for Pfam) of match in the full alignment.*/
	private static final String GF_TRUSTED_CUTOFF = "TC";
	/**	Highest sequence score (and domain score for Pfam) of match not in full alignment.*/
	private static final String GF_NOISE_CUTOFF = "NC";
	/**	Type of family -- presently Family, Domain, Motif or Repeat for Pfam.
	 * -- a tree with roots Gene, Intron or Cis-reg for Rfam. */
	private static final String GF_TYPE_FIELD = "TP";
	/** Number of sequences in alignment, and start of MSA.*/
	private static final String GF_SEQUENCE = "SQ";


	// OPTIONAL FIELDS
	
	/** Comment about database reference.*/
	private static final String GF_DB_COMMENT = "DC";
	/** Reference to external database.*/
	private static final String GF_DB_REFERENCE = "DR";
	/** Comment about literature reference. */
	private static final String GF_REFERENCE_COMMENT = "RC";
	/** Reference Number.*/
	private static final String GF_REFERENCE_NUMBER = "RN";
	/** Eight digit medline UI number.*/
	private static final String GF_REFERENCE_MEDLINE = "RM";
	/**	Reference Title.*/
	private static final String GF_REFERENCE_TITLE = "RT";
	/**	Reference Author.*/
	private static final String GF_REFERENCE_AUTHOR = "RA";
	/** Journal Location.*/
	private static final String GF_REFERENCE_LOCALTION = "RL";
	/** Record of all previous ID lines.*/
	private static final String GF_PREVIOUS_IDS = "PI";
	/** Keywords*/
	private static final String GF_KEYWORDS = "KW";
	/** Comments*/
	private static final String GF_COMMENT = "CC";
	/** Indicates a nested domain*/
	private static final String GF_PFAM_ACCESSION = "NE";
	/** Location of nested domains - sequence ID, start and end of insert.*/
	private static final String GF_LOCATION = "NL";
	/** Wikipedia page*/
	private static final String GF_WIKIPEDIA_LINK = "WK";
	/** Clan accession*/
	private static final String GF_CLAN = "CL";
	/** Used for listing Clan membership*/
	private static final String GF_MEMBERSHIP = "MB";

	/** FOR EMBEDDING TREES **/
	 
	/** A tree in New Hampshire eXtended format.*/
	private static final String GF_NEW_HAMPSHIRE = "NH";
	/** A unique identifier for the next tree.*/
	private static final String GF_TREE_ID = "TN";

	
	// OTHER
	
	/** A method used to set the bit score threshold based on the ratio of
	 * expected false positives to true positives. Floating point number between 0 and 1. */
	private static final String GF_FALSE_DISCOVERY_RATE = "FR";

	//	#=GS <seqname> <feature> <Generic per-Sequence annotation, free text>
	
	private static final String GS_ACCESSION_NUMBER = "AC";
	private static final String GS_DESCRIPTION = "DE";
	private static final String GS_DATABASE_REFERENCE = "DR";
	private static final String GS_ORGANISM_SPECIES = "OS";
	private static final String GS_ORGANISM_CLASSIFICATION = "OC";
	private static final String GS_LOOK = "LO";

//	#=GR <seqname> <feature> <Generic per-Residue annotation, exactly 1 char per residue>
	
	/** For RNA [.,;<>(){}[]AaBb...],<br>
	 * For protein [HGIEBTSCX]	 */
	private static final String GR_SECONDARY_STRUCTURE = "SS";
	/**[0-9X]<br>
	 * (0=0%-10%; ...; 9=90%-100%)*/
	private static final String GR_SURFACE_ACCESSIBILITY = "SA";
	
	/**[Mio] */
	private static final String GR_TRANS_MEMBRANE = "TM";
	/** [0-9*]<br>
	 * (0=0.00-0.05; 1=0.05-0.15; *=0.95-1.00)	 */
	private static final String GR_POSTERIOR_PROBABILITY = "PP";
	/** [*] */
	private static final String GR_LIGAND_BINDING = "LI";
	/** [*] */
	private static final String GR_ACTIVE_SITE = "AS";
	/** [*] */
	private static final String GR_AS_PFAM_PREDICTED = "pAS";
	/** [*] */
	private static final String GR_AS_SWISSPROT = "sAS";
	/** [0-2] */
	private static final String GR_INTRON = "IN";

//	#=GC <feature> <Generic per-Column annotation, exactly 1 char per column>
	
	private static final String GC_SEQUENSE_CONSENSUS = "seq_cons";
	private static final String GC_SECONDARY_STRUCTURE = "SS_cons";
	private static final String GC_SURFACE_ACCESSIBILITY = "SA_cons";
	private static final String GC_TRANS_MEMBRANE = "TM_cons";
	private static final String GC_POSTERIOR_PROBABILITY = "PP_cons";
	private static final String GC_LIGAND_BINDING = "LI_cons";
	private static final String GC_ACTIVE_SITE = "AS_cons";
	private static final String GC_AS_PFAM_PREDICTED = "pAS_cons";
	private static final String GC_AS_SWISSPROT = "sAS_cons";
	private static final String GC_INTRON = "IN_cons";
	/** Often the consensus RNA or protein sequence is used as a reference
	 * Any non-gap character (eg. x's) can indicate consensus/conserved/match
	 * columns .'s or -'s indicate insert columns ~'s indicate unaligned
	 * insertions Upper and lower case can be used to discriminate strong and
	 * weakly conserved residues respectively
	 */
	private static final String GC_REFERENCE_ANNOTATION = "RF";
	/**Indicates which columns in an alignment should be masked, such that the
	 * emission probabilities for match states corresponding to those columns
	 * will be the background distribution. */
	private static final String GC_MODEL_MASK = "MM";

	private StockholmStructure stockholmStructure;
//	private boolean endFile = false;
	private static final int STATUS_OUTSIDE_FILE = 0;
	private static final int STATUS_INSIDE_FILE  = 10;
	private static final int STATUS_IN_SEQUENCE  = 20;
	
	private int status=STATUS_OUTSIDE_FILE;
	Scanner internalScanner= null;
	private InputStream cashedInputStream;


	/**
	 * Parses a Stockholm file and returns a {@link StockholmStructure} object with its content.<br>
	 * This function is meant to be used for single access to specific 
	 * file and it closes the file after doing its assigned job. Any subsequent call 
	 * to {@link #parseNext(int)} will throw an exception or will function with unpredicted behavior.
	 * 
	 * @param filename complete(?) path to the file from where to read the content
	 * @return stockholm file content
	 * @throws IOException when an exception occurred while opening/reading/closing the file+
	 * @throws ParserException if unexpected format is encountered
	 */
	public StockholmStructure parse(String filename) throws IOException,ParserException{
		InputStream inStream = new InputStreamProvider().getInputStream(filename);
		StockholmStructure structure = parse(inStream);
		inStream.close();
		return structure;
	}
	/**
	 * Parses a Stockholm file and returns a {@link StockholmStructure} object with its content.<br>
	 * This function doesn't close the file after doing its assigned job; to allow for further calls of {@link #parseNext(int)}.
	 * @see #parseNext(int)
	 * 
	 * @param filename file from where to read the content. see {@link InputStreamProvider} for more details.
	 * @param max maximum number of files to read, {@link #INFINITY} for all.
	 * @return a vector of {@link StockholmStructure} containing parsed structures.
	 * @throws IOException when an exception occurred while opening/reading/closing the file.
	 * @throws ParserException if unexpected format is encountered
	 */
	public List<StockholmStructure> parse(String filename, int max) throws IOException,ParserException{
		InputStreamProvider isp = new InputStreamProvider();
		InputStream inStream = isp.getInputStream(filename);
		return parse(inStream, max);
	}

	/**parses {@link InputStream} and returns a the first contained alignment in a {@link StockholmStructure} object.
	 * Used mainly for multiple files within the same input stream, (e.g. when 
	 * reading from Pfam flat files. <br>
	 * This method leaves the stream open for further calls of {@link #parseNext(int)}.
	 * @see #parseNext(int)
	 * @param inStream the {@link InputStream} containing the file to read.
	 * @return a {@link StockholmStructure} object representing file contents.
	 * @throws IOException 
	 * @throws ParserException 
	 */
	public StockholmStructure parse(InputStream inStream) throws ParserException, IOException {
		return parse(inStream,1).get(0);
	}

	/**parses an {@link InputStream} and returns at maximum <code>max</code> objects contained in
	 * that file.<br>
	 * This method leaves the stream open for further calls of {@link #parse(InputStream, int)} (same function) or {@link #parseNext(int)}.
	 * 
	 * @see #parseNext(int)
	 * @param inStream the stream to parse
	 * @param max maximum number of structures to try to 
	 * parse, {@link #INFINITY} to try to obtain as much as possible. 
	 * @return a {@link List} of {@link StockholmStructure} objects. If there are no more
	 * structures, an empty list is returned.
	 * @throws IOException in case an I/O Exception occurred.
	 */
	public List<StockholmStructure> parse(InputStream inStream, int max) throws IOException {
		if (max < INFINITY) {
			throw new IllegalArgumentException("max can't be -ve value "+max);
		}
		if (inStream != this.cashedInputStream) {
			this.cashedInputStream=inStream;
			this.internalScanner=null;
		}
		
		if (internalScanner == null) {
			internalScanner= new Scanner(inStream);
		}
		ArrayList<StockholmStructure> structures= new ArrayList<StockholmStructure>();
		while (max != INFINITY && max-- >0) {
			StockholmStructure structure = parse(internalScanner);
			if(structure != null){
				structures.add(structure);
			}else {
				break;
			}
		}
		return structures;
	}
	
	/**Tries to parse and return as maximum as <code>max</code> structures in the last used file or input stream.<br>
	 * Please consider calling either {@link #parse(InputStream)}, 
	 * {@link #parse(InputStream, int)}, or {@link #parse(String, int)} before calling this function.
	 * @param max
	 * @return
	 * @throws IOException
	 */
	public List<StockholmStructure> parseNext(int max) throws IOException {
		return parse(this.cashedInputStream, max);
	}

	/**
	 * Parses a Stockholm file and returns a {@link StockholmStructure} object with its content.
	 * This method returns just after reaching the end of structure delimiter line ("//"), leaving any remaining empty lines unconsumed.
	 * 
	 * @param scanner from where to read the file content
	 * @return Stockholm file content, <code>null</code> if couldn't or no more structures.
	 * @throws IOException 
	 * @throws Exception
	 */
	StockholmStructure parse(Scanner scanner) throws ParserException, IOException {
		if (scanner == null) {
			if(internalScanner != null){
				scanner = internalScanner;
			}else {
				throw new IllegalArgumentException("No Scanner defined");
			}
		}
		String line = null;
		int linesCount = 0;
		try {
			while(scanner.hasNextLine()){
				line = scanner.nextLine();
				// if the file is empty
				//this condition will not happen, just left in case we decided to go for buffereedReader again for performance purpose.
				if (linesCount == 0 && line == null) {
					throw new IOException("Could not parse Stockholm file, BufferedReader returns null!");
				}

				// ignore empty lines
				if ((/*status==STATUS_INSIDE_FILE &&*/ line == null) || line.trim().length()==0) {
					continue;
				}

				if (line.startsWith("#=G")) {
//					// comment line or metadata
//					line = line.substring(1).trim();
//					line = line.substring(1).trim();
					if (line.startsWith(GENERIC_PER_FILE_ANNOTATION,2)) {
						// #=GF <featurename> <generic per-file annotation, free text>
						int firstSpaceIndex=line.indexOf(' ', 5);
						String featureName=line.substring(5, firstSpaceIndex);
						String value= line.substring(firstSpaceIndex).trim();
						handleFileAnnotation(featureName,value);
					} else if (line.startsWith(GENERIC_PER_CONSENSUS_ANNOTATION,2)) {
						//Being in a consensus means we are no longer in a sequence.
						this.status=STATUS_INSIDE_FILE;
						// #=GC <featurename> <generic per-column annotation, exactly 1 char per column>
						int firstSpaceIndex=line.indexOf(' ', 5);
						String featureName=line.substring(5, firstSpaceIndex);
						String value= line.substring(firstSpaceIndex).trim();
						handleConsensusAnnotation(featureName, value);
					} else if (line.startsWith(GENERIC_PER_SEQUENCE_ANNOTATION,2)) {
						// #=GS <seqname> <featurename> <generic per-sequence annotation, free text>
						int index1=line.indexOf(' ', 5);
						String seqName=line.substring(5, index1);
						while (line.charAt(++index1)<= ' ')//i.e. white space
							;//keep advancing
						int index2=line.indexOf(' ', index1);
						String featureName=line.substring(index1, index2);
						String value= line.substring(index2).trim();
						handleSequenceAnnotation(seqName,featureName,value);
					} else if (line.startsWith(GENERIC_PER_RESIDUE_ANNOTATION,2)) {
						// #=GR <seqname> <featurename> <generic per-sequence AND per-column mark-up, exactly 1 character per column>
						int index1=line.indexOf(' ', 5);
						String seqName=line.substring(5, index1);
						while (line.charAt(++index1)== ' ') 
							;//keep advancing
						int index2=line.indexOf(' ', index1);
						String featureName=line.substring(index1, index2);
						String value= line.substring(index2).trim();
						handleResidueAnnotation(seqName,featureName,value);
					}
				}else if(line.startsWith("# STOCKHOLM")) { // it is the header line
					if (status== STATUS_OUTSIDE_FILE) {
						status= STATUS_INSIDE_FILE;
						String[] header = line.split("\\s+");
						this.stockholmStructure = new StockholmStructure();
						this.stockholmStructure.getFileAnnotation().setFormat(header[1]);
						this.stockholmStructure.getFileAnnotation().setVersion(header[2]);
					} else {
						throw new ParserException("Uexpected Format line: ["+line+"]");
					}
				} else if (line.trim().equals("//")) {
					status=STATUS_OUTSIDE_FILE;
					break;//should we just break immediately or jump next empty lines?
				} else /*if (!line.startsWith("#")) */{
					if (status == STATUS_IN_SEQUENCE) {
						// This line corresponds to a sequence. Something like:
						// O83071/192-246 MTCRAQLIAVPRASSLAEAIACAQKMRVSRVPVYERS
						handleSequenceLine(line);
//					}else if (status==STATUS_OUTSIDE_FILE) {
//						throw new ParserException("The end of file character was allready reached but there are still sequence lines");
					}else {
						System.err.println("Error: Unknown or unexpected line [" +line+"].\nPlease contact the Biojava team.");
						throw new ParserException("Error: Unknown or unexpected line [" +line+"].");
					}
				}
				linesCount++;
			}
		} catch (IOException e) {
			e.printStackTrace();
			throw new IOException("Error parsing Stockholm file");
		}
		StockholmStructure structure = this.stockholmStructure;
		this.stockholmStructure=null;
		if (structure != null) {
			int length = -1;
			Map<String, StringBuffer> sequences = structure.getSequences();
			for (String sequencename : sequences.keySet()) {
				StringBuffer sequence = sequences.get(sequencename);
				if (length == -1) {
					length = sequence.length();
				} else if (length != sequence.length()) {
					throw new RuntimeException(
							"Sequences have different lengths");
				}
			}
		}
		return structure;
	}

	/**
	 * Handles a line that corresponds to a sequence. <br>
	 * e.g.: COATB_BPIKE/30-81
	 * AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA<br>
	 * N.B.: This function can't tolerate sequences with intrinsic white space.
	 * @param line
	 *            the line to be parsed
	 * @throws Exception
	 */
	private void handleSequenceLine(String line) throws ParserException {
		String[] lineContent = line.split("\\s+");
		if (lineContent.length != 2) {
			throw new ParserException("Could not split sequence line into sequence name and sequence:\n"
							+ line);
		}
		stockholmStructure.addSequence(lineContent[0], lineContent[1]);
	}

	/**
	 * #=GF &lt;feature&gt; &lt;Generic per-File annotation, free text&gt;
	 * @param featureName 
	 * @param value the line to be parsed
	 */
	private void handleFileAnnotation(String featureName, String value) {
		if (featureName.equals(GF_ACCESSION_NUMBER)) {
			stockholmStructure.getFileAnnotation().setGFAccessionNumber(value);
		} else if (featureName.equals(GF_IDENTIFICATION)) {
			stockholmStructure.getFileAnnotation().setGFIdentification(value);
		} else if (featureName.equals(GF_DB_REFERENCE)) { 
			stockholmStructure.getFileAnnotation().addDBReference(value);
		} else if (featureName.equals(GF_DEFINITION)) {
			stockholmStructure.getFileAnnotation().setGFDefinition(value);
		} else if (featureName.equals(GF_AUTHOR)) {
			stockholmStructure.getFileAnnotation().setGFAuthors(value);
		} else if (featureName.equals(GF_ALIGNMENT_METHOD)) {
			stockholmStructure.getFileAnnotation().setAlignmentMethod(value);
		} else if (featureName.equals(GF_BUILD_METHOD)) {
			stockholmStructure.getFileAnnotation().addGFBuildMethod(value);
		} else if (featureName.equals(GF_SEARCH_METHOD)) {
			stockholmStructure.getFileAnnotation().setGFSearchMethod(value);
		} else if (featureName.equals(GF_SOURCE_SEED)) {
			stockholmStructure.getFileAnnotation().setGFSourceSeed(value);
		} else if (featureName.equals(GF_SOURCE_STRUCTURE)) {
			stockholmStructure.getFileAnnotation().setGFSourceStructure(value);
		} else if (featureName.equals(GF_GATHERING_THRESHOLD)) {
			stockholmStructure.getFileAnnotation().setGFGatheringThreshs(value);
		} else if (featureName.equals(GF_TRUSTED_CUTOFF)) {
			stockholmStructure.getFileAnnotation().setGFTrustedCutoffs(value);
		} else if (featureName.equals(GF_NOISE_CUTOFF)) {
			stockholmStructure.getFileAnnotation().setGFNoiseCutoffs(value);
		} else if (featureName.equals(GF_TYPE_FIELD)) {
			stockholmStructure.getFileAnnotation().setGFTypeField(value);
		} else if (featureName.equals(GF_PREVIOUS_IDS)) {
			stockholmStructure.getFileAnnotation().setGFPreviousIDs(value);
		} else if (featureName.equals(GF_SEQUENCE)) {
			status=STATUS_IN_SEQUENCE;
			stockholmStructure.getFileAnnotation().setGFNumSequences(value);
		} else if (featureName.equals(GF_DB_COMMENT)) {
			stockholmStructure.getFileAnnotation().setGFDBComment(value);
//		} else if (featureName.equals(GF_DB_REFERENCE)) {
//			stockholmStructure.getFileAnnotation().addDBReference(value);
		} else if (featureName.equals(GF_REFERENCE_COMMENT)) {
			stockholmStructure.getFileAnnotation().setGFRefComment(value);
		} else if (featureName.equals(GF_REFERENCE_NUMBER)) {
			StockholmFileAnnotationReference reference = new StockholmFileAnnotationReference();
			stockholmStructure.getFileAnnotation().getReferences().add(reference);
		} else if (featureName.equals(GF_REFERENCE_MEDLINE)) {
			stockholmStructure.getFileAnnotation().getReferences().lastElement().setRefMedline(value);
		} else if (featureName.equals(GF_REFERENCE_TITLE)) {
			stockholmStructure.getFileAnnotation().getReferences().lastElement().addToRefTitle(value);
		} else if (featureName.equals(GF_REFERENCE_AUTHOR)) {
			stockholmStructure.getFileAnnotation().getReferences().lastElement().addToRefAuthor(value);
		} else if (featureName.equals(GF_REFERENCE_LOCALTION)) {
			stockholmStructure.getFileAnnotation().getReferences().lastElement().setRefLocation(value);
		} else if (featureName.equals(GF_KEYWORDS)) {
			stockholmStructure.getFileAnnotation().setGFKeywords(value);
		} else if (featureName.equals(GF_COMMENT)) {
			stockholmStructure.getFileAnnotation().addToGFComment(value);
		} else if (featureName.equals(GF_PFAM_ACCESSION)) {
			stockholmStructure.getFileAnnotation().setGFPfamAccession(value);
		} else if (featureName.equals(GF_LOCATION)) {
			stockholmStructure.getFileAnnotation().setGFLocation(value);
		} else if (featureName.equals(GF_WIKIPEDIA_LINK)) {
			stockholmStructure.getFileAnnotation().setGFWikipediaLink(value);
		} else if (featureName.equals(GF_CLAN)) {
			stockholmStructure.getFileAnnotation().setGFClan(value);
		} else if (featureName.equals(GF_MEMBERSHIP)) {
			stockholmStructure.getFileAnnotation().setGFMembership(value);
		} else if (featureName.equals(GF_NEW_HAMPSHIRE)) {
			stockholmStructure.getFileAnnotation().addGFNewHampshire(value);
		} else if (featureName.equals(GF_TREE_ID)) {
			stockholmStructure.getFileAnnotation().addGFTreeID(value);
		} else if (featureName.equals(GF_FALSE_DISCOVERY_RATE)) {
			stockholmStructure.getFileAnnotation().addGFFalseDiscoveryRate(value);
		} else {
			// unknown feature
			System.err.println("Warning: Unknown File Feature [" +featureName+"].\nPlease contact the Biojava team.");
		}
	}

	/**usually a single line of:<br>
	 * #=GC &lt;feature&gt; &lt;Generic per-Column annotation, exactly 1 char per column&gt;
	 * @param featureName the feature name :)
	 * @param value the line to be parsed.
	 */
	private void handleConsensusAnnotation(String featureName, String value) {
		if (featureName.equals(GC_SECONDARY_STRUCTURE)) {
			stockholmStructure.getConsAnnotation().setSecondaryStructure(value);
		} else if (featureName.equals(GC_SEQUENSE_CONSENSUS)) {
			stockholmStructure.getConsAnnotation().setSequenceConsensus(value);
		} else if (featureName.equals(GC_SURFACE_ACCESSIBILITY)) {
			stockholmStructure.getConsAnnotation().setSurfaceAccessibility(value);
		} else if (featureName.equals(GC_TRANS_MEMBRANE)) {
			stockholmStructure.getConsAnnotation().setTransMembrane(value);
		} else if (featureName.equals(GC_POSTERIOR_PROBABILITY)) {
			stockholmStructure.getConsAnnotation().setPosteriorProbability(value);
		} else if (featureName.equals(GC_LIGAND_BINDING)) {
			stockholmStructure.getConsAnnotation().setLigandBinding(value);
		} else if (featureName.equals(GC_ACTIVE_SITE)) {
			stockholmStructure.getConsAnnotation().setActiveSite(value);
		} else if (featureName.equals(GC_AS_PFAM_PREDICTED)) {
			stockholmStructure.getConsAnnotation().setAsPFamPredicted(value);
		} else if (featureName.equals(GC_AS_SWISSPROT)) {
			stockholmStructure.getConsAnnotation().setAsSwissProt(value);
		} else if (featureName.equals(GC_INTRON)) {
			stockholmStructure.getConsAnnotation().setIntron(value);
		} else if (featureName.equals(GC_REFERENCE_ANNOTATION)) {
			stockholmStructure.getConsAnnotation().setReferenceAnnotation(value);
		} else if (featureName.equals(GC_MODEL_MASK)) {
			stockholmStructure.getConsAnnotation().setModelMask(value);
		} else {
			// unknown feature
			System.err.println("Warning: Unknown Consensus Feature [" +featureName+"].\nPlease contact the Biojava team.");
		}
	}

	/**
	 * #=GS &lt;seqname&gt; &lt;feature&gt; &lt;Generic per-Sequence annotation, free text&gt;
	 * 
	 * @param line the line to be parsed
	 */
	private void handleSequenceAnnotation(String seqName, String featureName,String value) {
		if (featureName.equals(GS_ACCESSION_NUMBER)) {
			stockholmStructure.addGSAccessionNumber(seqName, value);
		} else if (featureName.equals(GS_DESCRIPTION)) {
			stockholmStructure.addGSDescription(seqName, value);
		} else if (featureName.equals(GS_DATABASE_REFERENCE)) {
			stockholmStructure.addGSdbReference(seqName, value);
		} else if (featureName.equals(GS_ORGANISM_SPECIES)) {
			stockholmStructure.addGSOrganismSpecies(seqName, value);
		} else if (featureName.equals(GS_ORGANISM_CLASSIFICATION)) {
			stockholmStructure.addGSOrganismClassification(seqName, value);
		} else if (featureName.equals(GS_LOOK)) {
			stockholmStructure.addGSLook(seqName, value);
		} else {
			// unknown feature
			System.err.println("Warning: Unknown Sequence Feature [" +featureName+"].\nPlease contact the Biojava team.");
		}
	}

	/**
	 * #=GR &lt;seqname&gt; &lt;feature&gt; &lt;Generic per-Residue annotation, exactly 1 char per residue&gt;
	 * 
	 * @param line
	 *            the line to be parsed
	 */
	private void handleResidueAnnotation(String seqName, String featureName,String value) {

		if (featureName.equals(GR_SURFACE_ACCESSIBILITY)) {
			stockholmStructure.addSurfaceAccessibility(seqName, value);
		} else if (featureName.equals(GR_TRANS_MEMBRANE)) {
			stockholmStructure.addTransMembrane(seqName, value);
		} else if (featureName.equals(GR_POSTERIOR_PROBABILITY)) {
			stockholmStructure.addPosteriorProbability(seqName, value);
		} else if (featureName.equals(GR_LIGAND_BINDING)) {
			stockholmStructure.addLigandBinding(seqName, value);
		} else if (featureName.equals(GR_ACTIVE_SITE)) {
			stockholmStructure.addActiveSite(seqName, value);
		} else if (featureName.equals(GR_AS_PFAM_PREDICTED)) {
			stockholmStructure.addASPFamPredicted(seqName, value);
		} else if (featureName.equals(GR_AS_SWISSPROT)) {
			stockholmStructure.addASSwissProt(seqName, value);
		} else if (featureName.equals(GR_INTRON)) {
			stockholmStructure.addIntron(seqName, value);
		} else if (featureName.equals(GR_SECONDARY_STRUCTURE)) {
			stockholmStructure.addSecondaryStructure(seqName, value);
		} else {
			// unknown feature
			System.err.println("Warning: Unknown Residue Feature [" +featureName+"].\nPlease contact the Biojava team.");
		}
	}
}
