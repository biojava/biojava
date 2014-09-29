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
 * @author Richard Holland
 * @author Mark Schreiber
 * @author David Scott
 * @author Bubba Puryear
 * @author George Waldon
 * @author Deepak Sheoran
 * @author Karl Nicholas <github:karlnicholas>
 * @author Paolo Pavan 
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
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava3.core.exceptions.ParserException;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.features.AbstractFeature;
import org.biojava3.core.sequence.features.DBReferenceInfo;
import org.biojava3.core.sequence.features.Qualifier;
import org.biojava3.core.sequence.features.TextFeature;
import org.biojava3.core.sequence.io.template.SequenceParserInterface;
import org.biojava3.core.sequence.location.SequenceLocation;
import org.biojava3.core.sequence.location.SimpleLocation;
import org.biojava3.core.sequence.location.template.AbstractLocation;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

public class GenbankSequenceParser<S extends AbstractSequence<C>, C extends Compound> implements SequenceParserInterface {
	private String seqData = null;
	private GenericGenbankHeaderParser<S,C> headerParser;	
    private String header;
    private String accession;
    public LinkedHashMap<String, ArrayList<DBReferenceInfo>> mapDB;
    /**
     * this data structure collects list of features extracted from the FEATURE_TAG section
     * They are organized by list of the same type (i.e. same genbank Feature)
     * and are provided with location
     */
    private HashMap<String, ArrayList<AbstractFeature>> featureCollection;


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
    /**
     * parse a location. 
     * if group(1) is null than the feature is on the positive strand,
     * group(2) start position,
     * group(3) end position.
     */
    protected static final Pattern singleLocationPattern = Pattern.compile("^\\s?(<?)(\\d+)(?:\\.\\.(>?)(\\d+)(>?))?");
    /**
     * Decodes a split pattern.
     * Split patterns are a composition of multiple locationsString qualified by actions:
  join(location,location, ... location): The indicated elements should be joined (placed end-to-end) to form one contiguous sequence.
  order(location,location, ... location): The elements can be found in the specified order (5' to 3' direction), nothing is implied about their reasonableness
  bond(location,location...location): Found in protein files. These generally are used to describe disulfide bonds.
  complement(location,location...location): consider locations in their complement versus
 
  takes in input a comma splitted location string. The split must be done for outer level commas
  group(1) is the qualifier
  group(2) is the location string to getFeatures. In case of complex splits it will contain the nested expression
 
 Not really sure that they are not declared obsolete but they are still in several files.
     */
    //protected static final Pattern genbankSplitPattern = Pattern.compile("^\\s?(join|order|bond|complement|)\\(?([^\\)]+)\\)?");
    protected static final Pattern genbankSplitPattern = Pattern.compile("^\\s?(join|order|bond|complement|)\\(?(.+)\\)?");
    /**
     * designed to recursively split a location string in tokens.
     * Valid tokens are those divided by coma that are not inside a bracket.
     * I. e. split on the comma only if that comma has zero, or an even number of quotes ahead of it.
     */
    protected static final String locationSplitPattern = ",(?=([^\\(|\\)]+\\([^\\(|\\)]+\\))[^\\(|\\)]+)";
    
    enum complexFeaturesAppendEnum {
        FLATTEN,HIERARCHICAL;
    }
    /**
     * define the mode in which complex features should be appended
     * in FLATTEN mode their will be appended to a single feature
     * in HIERARCHICAL mode, the single mother feature will have a tree of features that will
     * reflect the construction in genbank file
     */
    private complexFeaturesAppendEnum complexFeaturesAppendMode = complexFeaturesAppendEnum.FLATTEN;

    public void setComplexFeaturesAppendMode(complexFeaturesAppendEnum complexFeaturesAppendMode) {
        this.complexFeaturesAppendMode = complexFeaturesAppendMode;
    }
    private final Qualifier partialOn5end = new Qualifier("partial on 5' end", "true");
    private final Qualifier partialOn3end = new Qualifier("partial on 3' end", "true");
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
    /**
     * these variables are used to compute the global start and end of complex features
     */
    protected Integer featureGlobalStart, featureGlobalEnd;
    
    @SuppressWarnings("unchecked")
	private String parse(BufferedReader bufferedReader) {
	    String sectionKey = null;
	    List<String[]> section;
	    // Get an ordered list of key->value pairs in array-tuples
        do {
            section = this.readSection(bufferedReader);
            sectionKey = ((String[])section.get(0))[0];
            if(sectionKey == null) {
                throw new ParserException("Section key was null");
            }
            // getFeatures section-by-section
            if (sectionKey.equals(LOCUS_TAG)) {
                String loc = ((String[])section.get(0))[1];
                header = loc;
                Matcher m = lp.matcher(loc);
                if (m.matches()) {
                	headerParser.setName(m.group(1));
                    headerParser.setAccession(m.group(1)); // default if no accession found
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
                    }
                    if (m.group(3)!=null) headerParser.setVersion(Integer.parseInt(m.group(3)));
                    if (m.group(5)!=null) {
                        headerParser.setIdentifier(m.group(5));
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
                AbstractFeature gbFeature=null;
                // starting from second line of input, start a new feature whenever we come across
                // a key that does not start with /
                //boolean skippingBond = false;
                for (int i = 1 ; i < section.size(); i++) {
                    String key = ((String[])section.get(i))[0];
                    String val = ((String[])section.get(i))[1];
                    if (key.startsWith("/")) {
                          if (gbFeature==null) throw new ParserException("Malformed GenBank file: found a qualifier without feature.");
                    	  //if(!skippingBond)
                    	  //{
                            key = key.substring(1); // strip leading slash
                            val = val.replaceAll("\\s*[\\n\\r]+\\s*"," ").trim();
                            if (val.endsWith("\"")) val = val.substring(1,val.length()-1); // strip quotes
                            // parameter on old feature
                            if (key.equals("db_xref")) {
                                Matcher m = dbxp.matcher(val);
                                if (m.matches()) {
                                    String dbname = m.group(1);
                                    String raccession = m.group(2);
                                    Qualifier xref = new DBReferenceInfo(dbname, raccession);
                                    gbFeature.addQualifier(key,xref);
                                    
                                    ArrayList<DBReferenceInfo> listDBEntry = new ArrayList<DBReferenceInfo>();
                                    listDBEntry.add((DBReferenceInfo)xref);
                                    mapDB.put(key, listDBEntry);
                                } else {
                                    throw new ParserException("Bad dbxref");
                                }
                            } else if (key.equalsIgnoreCase("organism")) {
                                Qualifier q = new Qualifier(key, val.replace('\n', ' '));
                                gbFeature.addQualifier(key, q);
                            } else {
                                if (key.equalsIgnoreCase("translation")) {
                                    // strip spaces from sequence
                                    val = val.replaceAll("\\s+","");
                                    Qualifier q = new Qualifier(key, val);
                                    gbFeature.addQualifier(key, q);
                                }
                                else {
                                    Qualifier q = new Qualifier(key, val);
                                    gbFeature.addQualifier(key, q);
                                }
                            }
                      	//}
                    } else {
                        // new feature!
                        featureGlobalStart=Integer.MAX_VALUE;
                        featureGlobalEnd=1;
                        ArrayList<AbstractFeature> ff = getFeatures(val,1);
                        
                        if (ff.size()==1){
                            gbFeature = ff.get(0);
                        } else {
                            gbFeature = new TextFeature(key,val,key,key);
                            gbFeature.setChildrenFeatures(ff);
                        }
                        gbFeature.setType(key);
                        gbFeature.setSource(val);
                        gbFeature.setShortDescription(key);
                        gbFeature.setDescription(key);
                        gbFeature.setLocation(new SequenceLocation(featureGlobalStart, featureGlobalEnd, null));
                        
                        if (!featureCollection.containsKey(key)) 
                            featureCollection.put(key, new ArrayList());
                        featureCollection.get(key).add(gbFeature);
                    }
                }
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
        
        private ArrayList<AbstractFeature> getFeatures(String string, int versus){
            Matcher m;
            ArrayList<AbstractFeature> boundedFeaturesCollection = new ArrayList();
            
            //String[] tokens = string.split(locationSplitPattern);
            List<String> tokens = splitString(string);
            for(String t: tokens){
                m = genbankSplitPattern.matcher(t);
                if (!m.find()) throw new ParserException("Cannot interpret split pattern "+t
                        +"\nin location string:"+string);
                String splitQualifier = m.group(1);
                String splitString = m.group(2);
                
                if (!splitQualifier.isEmpty()){
                    //recursive case
                    int localVersus = splitQualifier.equalsIgnoreCase("complement")? -1 : 1;
                    ArrayList subFeatures = getFeatures(splitString, versus * localVersus);
                    
                    if (complexFeaturesAppendMode==complexFeaturesAppendEnum.FLATTEN){
                        boundedFeaturesCollection.addAll(subFeatures);
                    } else {
                        if (subFeatures.size()==1){
                            boundedFeaturesCollection.addAll(subFeatures);
                        } else {
                            TextFeature motherFeature = new TextFeature(splitQualifier,"","","");
                            motherFeature.setChildrenFeatures(subFeatures);
                            boundedFeaturesCollection.add(motherFeature);
                        }
                    }
                } else {
                    //base case
                    m = singleLocationPattern.matcher(splitString);
                    if (!m.find()) throw 
                            new ParserException("Cannot interpret location pattern "+splitString
                                    +"\nin location string:"+string);
                     
                    Strand s = versus == 1? Strand.POSITIVE : Strand.NEGATIVE;
                    int start = new Integer(m.group(2));
                    int end = m.group(4) == null? start : new Integer(m.group(4));
                    
                    if (featureGlobalStart > start) featureGlobalStart = start ;
                    if (featureGlobalEnd < end) featureGlobalEnd = end;
                    
                    SequenceLocation l = new SequenceLocation(
                                    start, 
                                    end,
                                    new org.biojava3.core.sequence.DNASequence(),
                                    s
                            );
                    
                    TextFeature subFeature = new TextFeature("", "", "", "");
                    subFeature.setLocation(l);
                    if(m.group(1).equals("<")) subFeature.addQualifier("partial_end", partialOn5end);
                    if(m.group(3)!=null && (m.group(3).equals(">")||m.group(5).equals(">"))) 
                        subFeature.addQualifier("partial_end", partialOn3end);
                    
                    boundedFeaturesCollection.add(subFeature);
                    
                }
            }
            
            return boundedFeaturesCollection;
        }
        
        private List<String> splitString(String input){
            List<String> result = new ArrayList<String>();
            int start = 0;
            int openedParenthesis = 0;
            for (int current = 0; current < input.length(); current++) {
                if (input.charAt(current) == '(') openedParenthesis++;
                if (input.charAt(current) == ')') openedParenthesis--;
                boolean atLastChar = (current == input.length() - 1);
                if(atLastChar) result.add(input.substring(start));
                else if (input.charAt(current) == ',' && openedParenthesis==0) {
                    result.add(input.substring(start, current));
                    start = current + 1;
                }
            }
            return result;
        }

	// reads an indented section, combining split lines and creating a list of
	// key->value tuples
	// reads an indented section, combining split lines and creating a list of
	// key->value tuples
	// reads an indented section, combining split lines and creating a list of
	// key->value tuples
	private List<String[]> readSection(BufferedReader bufferedReader) {
		List<String[]> section = new ArrayList<String[]>();
		String line = "";
		String currKey = null;
		StringBuffer currVal = new StringBuffer();
		boolean done = false;
		int linecount = 0;
		
		try {
			while (!done) {
				bufferedReader.mark(320);
				line = bufferedReader.readLine();
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
					bufferedReader.reset();
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

	@Override
    public String getSequence(BufferedReader bufferedReader, int sequenceLength) throws IOException {
		featureCollection = new HashMap<String, ArrayList<AbstractFeature>>();
		mapDB = new LinkedHashMap<String, ArrayList<DBReferenceInfo>>();
		headerParser = new GenericGenbankHeaderParser<S,C>();

		parse(bufferedReader);
		
		return seqData;
	}

	public String getHeader() {
		return header;
	}

	public GenericGenbankHeaderParser<S,C> getSequenceHeaderParser() {
		return headerParser;
	}
	public LinkedHashMap<String, ArrayList<DBReferenceInfo>> getDatabaseReferences() {
		return mapDB;
	}

        public HashMap<String, ArrayList<AbstractFeature>> getFeatureCollection() {
            return featureCollection;
        }

        public ArrayList<String> getKeyWords() {
            return new ArrayList<String>(featureCollection.keySet());
        }

}

