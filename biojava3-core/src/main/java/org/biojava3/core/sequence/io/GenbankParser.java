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
 * @author Karl Nicholas <github:karlnicholas>
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

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava3.core.exceptions.ParserException;
import org.biojava3.core.sequence.io.template.GenbankHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

public class GenbankParser<S extends Sequence<?>, C extends Compound> {
	private String seqData = null;
    private GenbankHeaderParserInterface<S,C> headerParser;	
	private String accession = null;
	private String identifier = null;
    private InputStreamReader isr;
    private BufferedReaderBytesRead br;
    private String header;


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

//  private NCBITaxon tax = null;

	public GenbankParser(InputStream is, GenbankHeaderParserInterface<S,C> headerParser) {
		this.headerParser = headerParser;
		isr = new InputStreamReader(is);
		this.br = new BufferedReaderBytesRead(isr);
	}

	public String parse() {
	    String sectionKey = null;
	    List<String[]> section;
	    // Get an ordered list of key->value pairs in array-tuples
        do {
            section = this.readSection();
            sectionKey = ((String[])section.get(0))[0];
            if(sectionKey == null) {
                throw new ParserException("Section key was null");
            }
            // process section-by-section
            if (sectionKey.equals(LOCUS_TAG)) {
                String loc = ((String[])section.get(0))[1];
                header = loc;
                Matcher m = lp.matcher(loc);
                if (m.matches()) {
                	headerParser.setName(m.group(1));
                    accession = m.group(1); // default if no accession found
                    headerParser.setAccession(accession);
                } else {
                    throw new ParserException("Bad locus line");
                }
            } else if (sectionKey.equals(DEFINITION_TAG)) {
            	headerParser.setDescription(((String[])section.get(0))[1]);
            } else if (sectionKey.equals(ACCESSION_TAG)) {
                // if multiple accessions, store only first as accession,
                // and store rest in annotation
                String[] accs = ((String[])section.get(0))[1].split("\\s+");
                accession = accs[0].trim();
                headerParser.setAccession(accession);
            } else if (sectionKey.equals(VERSION_TAG)) {
                String ver = ((String[])section.get(0))[1];
                Matcher m = vp.matcher(ver);
                if (m.matches()) {
                    String verAcc = m.group(1);
                    if (!accession.equals(verAcc)) {
                        // the version refers to a different accession!
                        // believe the version line, and store the original
                        // accession away in the additional accession set
                        accession = verAcc;
                        headerParser.setAccession(accession);
                    }
                    if (m.group(3)!=null) headerParser.setVersion(Integer.parseInt(m.group(3)));
                    if (m.group(5)!=null) {
                        identifier = m.group(5);
                        headerParser.setIdentifier(identifier);
                    }
                } else {
                    throw new ParserException("Bad version line");
                }
            } else if (sectionKey.equals(KEYWORDS_TAG)) {
            } else if (sectionKey.equals(SOURCE_TAG)) {
                // ignore - can get all this from the first feature
            } else if (sectionKey.equals(REFERENCE_TAG) ) {
            } else if (sectionKey.equals(COMMENT_TAG) ) {
                // Set up some comments
                headerParser.setComment(((String[])section.get(0))[1]);
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
            	seqData = seq.toString().replaceAll("\\s+","").replaceAll("[\\.|~]","-").toUpperCase();
            }
        } while (!sectionKey.equals(END_SEQUENCE_TAG));
        return seqData;
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
    }

	public String getSeqData() {
		return seqData;
	}

	public S getSequence(SequenceCreatorInterface<C> sequenceCreator) {
	      @SuppressWarnings("unchecked")
	      S sequence = (S)sequenceCreator.getSequence(seqData, 0);
	      headerParser.parseHeader(header, sequence);
	      return sequence;
	}
}
