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
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.DataSource;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.template.GenbankHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.AbstractSequence.AnnotationType;
import org.biojava3.core.exceptions.ParserException;

/**
 * Use GenbankReaderHelper as an example of how to use this class where GenbankReaderHelper should be the
 * primary class used to read Genbank files
 * -- copied from original FastReader by Scooter Willis ;lt;willishf at gmail dot com&gt;
 * @author Karl Nicholas 

 */
public class GenbankReader<S extends AbstractSequence<C>, C extends Compound> {

    /**
     * The name of this format
     */
    public static final String GENBANK_FORMAT = "GENBANK";
    
    protected static final String LOCUS_TAG =           "LOCUS";
    protected static final String DEFINITION_TAG =      "DEFINITION";
    protected static final String ACCESSION_TAG =       "ACCESSION";
    protected static final String VERSION_TAG =         "VERSION";
    protected static final String KEYWORDS_TAG =        "KEYWORDS";
    //                                                  "SEGMENT"
    protected static final String SOURCE_TAG =          "SOURCE";
    protected static final String ORGANISM_TAG =        "ORGANISM";
    protected static final String REFERENCE_TAG =       "REFERENCE";
    protected static final String AUTHORS_TAG =         "AUTHORS";
    protected static final String CONSORTIUM_TAG =      "CONSRTM";
    protected static final String TITLE_TAG =           "TITLE";
    protected static final String JOURNAL_TAG =         "JOURNAL";
    protected static final String PUBMED_TAG =          "PUBMED";
    protected static final String MEDLINE_TAG =         "MEDLINE"; //deprecated
    protected static final String REMARK_TAG =          "REMARK";
    protected static final String COMMENT_TAG =         "COMMENT";
    protected static final String FEATURE_TAG =         "FEATURES";
    protected static final String BASE_COUNT_TAG_FULL = "BASE COUNT"; //deprecated
    protected static final String BASE_COUNT_TAG =      "BASE";
    //                                                  "CONTIG"
    protected static final String START_SEQUENCE_TAG =  "ORIGIN";
    protected static final String END_SEQUENCE_TAG =    "//";
    // locus line
    protected static final Pattern lp = Pattern.compile("^(\\S+)\\s+\\d+\\s+(bp|aa)\\s{1,4}([dms]s-)?(\\S+)?\\s+(circular|linear)?\\s*(\\S+)?\\s*(\\S+)?$");
    // version line
    protected static final Pattern vp = Pattern.compile("^(\\S*?)(\\.(\\d+))?(\\s+GI:(\\S+))?$");
    // reference line
    protected static final Pattern refRange = Pattern.compile("^\\s*(\\d+)\\s+to\\s+(\\d+)$");
    protected static final Pattern refp = Pattern.compile("^(\\d+)\\s*(?:(\\((?:bases|residues)\\s+(\\d+\\s+to\\s+\\d+(\\s*;\\s*\\d+\\s+to\\s+\\d+)*)\\))|\\(sites\\))?");
    // dbxref line
    protected static final Pattern dbxp = Pattern.compile("^([^:]+):(\\S+)$");
    //sections start at a line and continue till the first line afterwards with a
    //non-whitespace first character
    //we want to match any of the following as a new section within a section
    //  \s{0,8} word \s{0,7} value
    //  \s{21} /word = value
    //  \s{21} /word
    protected static final Pattern sectp = Pattern.compile("^(\\s{0,8}(\\S+)\\s{0,7}(.*)|\\s{21}(/\\S+?)=(.*)|\\s{21}(/\\S+))$");
    
    protected static final Pattern readableFiles = Pattern.compile(".*(g[bp]k*$|\\u002eg[bp].*)");
    protected static final Pattern headerLine = Pattern.compile("^LOCUS.*");

    SequenceCreatorInterface<C> sequenceCreator;
    GenbankHeaderParserInterface<S,C> headerParser;
    BufferedReaderBytesRead br;
    InputStreamReader isr;
    FileInputStream fi = null;
    long fileIndex = 0;
    long sequenceIndex = 0;
    String line = "";
    String header= "";
    
    /**
     * If you are going to use FileProxyProteinSequenceCreator then do not use this constructor because we need details about
     * local file offsets for quick reads. InputStreams does not give you the name of the stream to access quickly via file seek. A seek in
     * an inputstream is forced to read all the data so you don't gain anything.
     * @param br
     * @param headerParser
     * @param sequenceCreator
     */
    public GenbankReader(InputStream is, GenbankHeaderParserInterface<S,C> headerParser,
    		SequenceCreatorInterface<C> sequenceCreator) {
        this.headerParser = headerParser;
        isr = new InputStreamReader(is);
        this.br = new BufferedReaderBytesRead(isr);
        this.sequenceCreator = sequenceCreator;
    }

    /**
     * If you are going to use the FileProxyProteinSequenceCreator then you
     * need to use this constructor because we need details about
     * the location of the file.
     * @param file
     * @param headerParser
     * @param sequenceCreator
     * @throws FileNotFoundException if the file does not exist, is a directory 
     * 	rather than a regular file, or for some other reason cannot be opened
     * 	for reading.
     * @throws SecurityException if a security manager exists and its checkRead
     * 	method denies read access to the file.
     */
    public GenbankReader(File file, GenbankHeaderParserInterface<S,C> headerParser, 
    		SequenceCreatorInterface<C> sequenceCreator) throws FileNotFoundException {
        this.headerParser = headerParser;
        fi = new FileInputStream(file);
        isr = new InputStreamReader(fi);
        this.br = new BufferedReaderBytesRead(isr);
        this.sequenceCreator = sequenceCreator;
    }

    /**
     * The parsing is done in this method.<br>
     * This method tries to process all the available Genbank records 
     * in the File or InputStream, closes the underlying resource, 
     * and return the results in {@link LinkedHashMap}.<br>
     * You don't need to call {@link #close()} after calling this method.
     * @see #process(int)
     * @return {@link HashMap} containing all the parsed Genbank records 
     * present, starting current fileIndex onwards.
     * @throws IOException if an error occurs reading the input file
     */
    public LinkedHashMap<String,S> process() throws IOException {
    	LinkedHashMap<String,S> sequences = process(-1);
    	close();
    	return sequences;
    }

    private String sectionKey = null;
//    private NCBITaxon tax = null;
    private String accession = null;
    private String identifier = null;
    private String seqName = null;


    /**
     * This method tries to parse maximum <code>max</code> records from
     * the open File or InputStream, and leaves the underlying resource open.<br>
     * Subsequent calls to the same method continue parsing the rest of the file.<br>
     * This is particularly useful when dealing with very big data files,
     * (e.g. NCBI nr database), which can't fit into memory and will take long
     * time before the first result is available.<br>
     * <b>N.B.</b>
     * <ul>
     * <li>This method ca't be called after calling its NO-ARGUMENT twin.</li> 
     * <li>remember to close the underlying resource when you are done.</li> 
     * </ul>
     * @see #process()
     * @author Amr AL-Hossary
     * @since 3.0.6
     * @param max maximum number of records to return, <code>-1</code> for infinity.
     * @return {@link HashMap} containing maximum <code>max</code> parsed Genbank records 
     * present, starting current fileIndex onwards.
     * @throws IOException if an error occurs reading the input file
     */
    public LinkedHashMap<String,S> process(int max) throws IOException {
        LinkedHashMap<String,S> sequences = new LinkedHashMap<String,S>();

        // Get an ordered list of key->value pairs in array-tuples
        List section = null;
        try{
            do {
                section = this.readSection();
                sectionKey = ((String[])section.get(0))[0];
                if(sectionKey == null){
                    throw new ParserException("Section key was null");
                }
                // process section-by-section
                if (sectionKey.equals(LOCUS_TAG)) {
                    String loc = ((String[])section.get(0))[1];
                    header = loc;
                    Matcher m = lp.matcher(loc);
                    if (m.matches()) {
                    	seqName = m.group(1);
                    }
                } else if (sectionKey.equals(DEFINITION_TAG)) {
                } else if (sectionKey.equals(ACCESSION_TAG)) {
                    // if multiple accessions, store only first as accession,
                    // and store rest in annotation
                    String[] accs = ((String[])section.get(0))[1].split("\\s+");
                    accession = accs[0].trim();
                } else if (sectionKey.equals(VERSION_TAG)) {
                } else if (sectionKey.equals(KEYWORDS_TAG)) {
                } else if (sectionKey.equals(SOURCE_TAG)) {
                    // ignore - can get all this from the first feature
                } else if (sectionKey.equals(REFERENCE_TAG) ) {
                } else if (sectionKey.equals(COMMENT_TAG) ) {
                } else if (sectionKey.equals(FEATURE_TAG) ) {
                } else if (sectionKey.equals(BASE_COUNT_TAG)) {
                    // ignore - can calculate from sequence content later if needed
                } else if (sectionKey.equals(START_SEQUENCE_TAG) ) {

                	// our first line is ignorable as it is the ORIGIN tag
                    // the second line onwards conveniently have the number as
                    // the [0] tuple, and sequence string as [1] so all we have
                    // to do is concat the [1] parts and then strip out spaces,
                    // and replace '.' and '~' with '-' for our parser.
                	StringBuffer seq = new StringBuffer();
                    for (int i = 1 ; i < section.size(); i++) seq.append(((String[])section.get(i))[1]);
                	String seqData = seq.toString().replaceAll("\\s+","").replaceAll("[\\.|~]","-").toUpperCase();

                    S sequence = (S)sequenceCreator.getSequence(seqData, sequenceIndex);
                    headerParser.parseHeader(header, sequence);
                    sequence.setOriginalHeader(seqName);
                    sequence.setAccession(new AccessionID(accession));
                    
                    sequences.put(sequence.getAccession().getID(),sequence);
                }
            } while (!sectionKey.equals(END_SEQUENCE_TAG));
        } catch(RuntimeException e) {
            throw new ParserException("Bad sequence section");
        }
        return sequences;
        
        /*
        
        StringBuilder sb = new StringBuilder();
        int processedSequences=0;
        boolean keepGoing = true;

        do {
            line = line.trim(); // nice to have but probably not needed
            if (line.length() != 0) {
                if (line.startsWith(">")) {
                    if (sb.length() > 0) {//i.e. if there is already a sequence before
                    //    System.out.println("Sequence index=" + sequenceIndex);
                    	@SuppressWarnings("unchecked")
                        S sequence = (S)sequenceCreator.getSequence(sb.toString(), sequenceIndex);
                        headerParser.parseHeader(header, sequence);
                        sequences.put(sequence.getAccession().getID(),sequence);
                        processedSequences++;
//                        if (maxSequenceLength < sb.length()) {
//                            maxSequenceLength = sb.length();
//                        }
//                        sb = new StringBuilder(maxSequenceLength);
                        sb.setLength(0); //this is faster, better memory utilization (same buffer)
                    }
                    header = line.substring(1);
                } else if (line.startsWith(";")) {
                } else {
                    //mark the start of the sequence with the fileIndex before the line was read
                    if(sb.length() == 0){
                        sequenceIndex = fileIndex;
                    }
                    sb.append(line);
                }
            }
            fileIndex = br.getBytesRead();
            line = br.readLine();
			if (line == null) {//i.e. EOF
				@SuppressWarnings("unchecked")
				//    System.out.println("Sequence index=" + sequenceIndex + " " + fileIndex );
                S sequence = (S)sequenceCreator.getSequence(sb.toString(), sequenceIndex);
                headerParser.parseHeader(header, sequence);
                sequences.put(sequence.getAccession().getID(),sequence);
                processedSequences++;
                keepGoing = false;
            }
			if (max > -1 && processedSequences>=max) {
				keepGoing=false;
			}
        } while (keepGoing);
        this.line  = line;
        this.header= header;
        return sequences;
*/        
    }

	// reads an indented section, combining split lines and creating a list of
	// key->value tuples
	private List<String[]> readSection() {
		List<String[]> section = new ArrayList<String[]>();
		String line = "";
		String currKey = null;
		StringBuffer currVal = new StringBuffer();
		boolean done = false;
		int linecount = 0;

		try {
			while (!done) {
				br.mark(320);
				line = br.readLine();
				String firstSecKey = section.isEmpty() ? ""
						: ((String[]) section.get(0))[0];
				if (line != null && line.matches("\\p{Space}*")) {
					// regular expression \p{Space}* will match line
					// having only white space characters
					continue;
				}
				if (line == null
						|| (!line.startsWith(" ") && linecount++ > 0 && (!firstSecKey
								.equals(START_SEQUENCE_TAG) || line
								.startsWith(END_SEQUENCE_TAG)))) {
					// dump out last part of section
					section.add(new String[] { currKey, currVal.toString() });
					br.reset();
					done = true;
				} else {
					Matcher m = sectp.matcher(line);
					if (m.matches()) {
						// new key
						if (currKey != null)
							section.add(new String[] { currKey,
									currVal.toString() });
						// key = group(2) or group(4) or group(6) - whichever is
						// not null
						currKey = m.group(2) == null ? (m.group(4) == null ? m
								.group(6) : m.group(4)) : m.group(2);
						currVal = new StringBuffer();
						// val = group(3) if group(2) not null, group(5) if
						// group(4) not null, "" otherwise, trimmed
						currVal.append((m.group(2) == null ? (m.group(4) == null ? ""
								: m.group(5))
								: m.group(3)).trim());
					} else {
						// concatted line or SEQ START/END line?
						if (line.startsWith(START_SEQUENCE_TAG)
								|| line.startsWith(END_SEQUENCE_TAG))
							currKey = line;
						else {
							currVal.append("\n"); // newline in between lines -
													// can be removed later
							currVal.append(currKey.charAt(0) == '/' ? line
									.substring(21) : line.substring(12));
						}
					}
				}
			}
		} catch (IOException e) {
			throw new ParserException(e.getMessage());
		} catch (RuntimeException e) {
			throw new ParserException(e.getMessage());
		}
		return section;
	}

	public void close() throws IOException {
		br.close();
        isr.close();
        //If stream was created from File object then we need to close it
        if (fi != null) {
            fi.close();
        }
        this.line=this.header = null;
	}

    public static void main(String[] args) {
        try {
            String inputFile = "src/test/resources/PF00104_small.Genbank";
            FileInputStream is = new FileInputStream(inputFile);

            GenbankReader<ProteinSequence, AminoAcidCompound> GenbankReader = new GenbankReader<ProteinSequence, AminoAcidCompound>(is, new GenericGenbankHeaderParser<ProteinSequence,AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String,ProteinSequence> proteinSequences = GenbankReader.process();
            is.close();


            System.out.println(proteinSequences);

            File file = new File(inputFile);
            GenbankReader<ProteinSequence,AminoAcidCompound> GenbankProxyReader = new GenbankReader<ProteinSequence,AminoAcidCompound>(file, new GenericGenbankHeaderParser<ProteinSequence,AminoAcidCompound>(), new FileProxyProteinSequenceCreator(file, AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String,ProteinSequence> proteinProxySequences = GenbankProxyReader.process();

            for(String key : proteinProxySequences.keySet()){
                ProteinSequence proteinSequence = proteinProxySequences.get(key);
                System.out.println(key);
//                if(key.equals("Q98SJ1_CHICK/15-61")){
//                    int dummy = 1;
//                }
                System.out.println(proteinSequence.toString());

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

