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
 */

package	org.biojavax.bio.seq.io;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.Comment;
import org.biojavax.CrossRef;
import org.biojavax.DocRef;
import org.biojavax.DocRefAuthor;
import org.biojavax.Namespace;
import org.biojavax.Note;
import org.biojavax.RankedCrossRef;
import org.biojavax.RankedDocRef;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleComment;
import org.biojavax.SimpleCrossRef;
import org.biojavax.SimpleDocRef;
import org.biojavax.SimpleRankedCrossRef;
import org.biojavax.SimpleRankedDocRef;
import org.biojavax.SimpleRichAnnotation;
import org.biojavax.bio.seq.CompoundRichLocation;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.SimplePosition;
import org.biojavax.bio.seq.SimpleRichLocation;
import org.biojavax.bio.taxa.NCBITaxon;
import org.biojavax.bio.taxa.SimpleNCBITaxon;
import org.biojavax.ontology.ComparableTerm;
import org.biojavax.utils.StringTools;

/**
 * Format reader for GenBank files. This version of Genbank format will generate
 * and write RichSequence objects. Loosely Based on code from the old, deprecated,
 * org.biojava.bio.seq.io.GenbankFormat object.
 *
 * @author Richard Holland
 * @author Mark Schreiber
 * @author David Scott
 * @author Bubba Puryear
 * @author George Waldon
 * @since 1.5
 */
public class GenbankFormat extends RichSequenceFormat.HeaderlessFormat {
    
    // Register this format with the format auto-guesser.
    static {
        RichSequence.IOTools.registerFormat(GenbankFormat.class);
    }
    
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
    //  \s{0,8} word \s{1,7} value
    //  \s{21} /word = value
    //  \s{21} /word
    protected static final Pattern sectp = Pattern.compile("^(\\s{0,8}(\\S+)\\s{1,7}(.*)|\\s{21}(/\\S+?)=(.*)|\\s{21}(/\\S+))$");
    
    protected static final Pattern readableFiles = Pattern.compile(".*(g[bp]k*$|\\u002eg[bp].*)");
    protected static final Pattern headerLine = Pattern.compile("^LOCUS.*");
    
    private final static HashSet isNotQuoted = new HashSet();
    static {
        isNotQuoted.add("anticodon");
        isNotQuoted.add("citation");
        isNotQuoted.add("codon");
        isNotQuoted.add("codon_start");
        isNotQuoted.add("compare");
        isNotQuoted.add("cons_splice");
        isNotQuoted.add("direction");
        isNotQuoted.add("estimated_length");
        isNotQuoted.add("label");
        isNotQuoted.add("mod_base");
        isNotQuoted.add("number");
        isNotQuoted.add("rpt_type");
        isNotQuoted.add("rpt_unit_range");
        isNotQuoted.add("transl_except");
        isNotQuoted.add("transl_table");
    }
    
    /**
     * Implements some GenBank-specific terms.
     */
    public static class Terms extends RichSequence.Terms {        
        /**
         * Getter for the Genbank term
         * @return The genbank Term
         */
        public static ComparableTerm getGenBankTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("GenBank");
        }
    }
    
    /**
     * {@inheritDoc}
     * A file is in GenBank format if the name ends with gbk, contains the letters egb, or the first line of
     * the file starts with the word LOCUS
     */
    public boolean canRead(File file) throws IOException {
        if (readableFiles.matcher(file.getName()).matches()) return true;
        BufferedReader br = new BufferedReader(new FileReader(file));
        final String firstLine = br.readLine();
        boolean readable = firstLine!=null && headerLine.matcher(firstLine).matches();
        br.close();
        return readable;
    }
    
    /**
     * {@inheritDoc}
     * Returns an dna parser if the letters DNA or RNA appear in the first line of the file.
     * Otherwise returns a DNA tokenizer.
     */
    public SymbolTokenization guessSymbolTokenization(File file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String firstLine = br.readLine();
        boolean dna = (firstLine.indexOf("DNA") >0 || firstLine.indexOf("RNA") > 0);
        br.close();
        if (dna) return RichSequence.IOTools.getDNAParser();
        else return RichSequence.IOTools.getProteinParser();
    }
    
    /**
     * {@inheritDoc}
     * A stream is in GenBank format if the first line of the stream starts with the word LOCUS
     */
    public boolean canRead(BufferedInputStream stream) throws IOException {
        stream.mark(2000); // some streams may not support this
        BufferedReader br = new BufferedReader(new InputStreamReader(stream));
        final String firstLine = br.readLine();
        boolean readable = firstLine!=null && headerLine.matcher(firstLine).matches();
        // don't close the reader as it'll close the stream too.
        // br.close();
        stream.reset();
        return readable;
    }
    
    /**
     * {@inheritDoc}
     * Returns an dna parser if the letters DNA or RNA appear in the first line of the stream.
     * Otherwise returns a DNA tokenizer.
     */
    public SymbolTokenization guessSymbolTokenization(BufferedInputStream stream) throws IOException {
        stream.mark(2000); // some streams may not support this
        BufferedReader br = new BufferedReader(new InputStreamReader(stream));
        String firstLine = br.readLine();
        boolean dna = (firstLine.indexOf("DNA") >0 || firstLine.indexOf("RNA") > 0);
        // don't close the reader as it'll close the stream too.
        // br.close();
        stream.reset();
        if (dna) return RichSequence.IOTools.getDNAParser();
        else return RichSequence.IOTools.getProteinParser();
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean readSequence(BufferedReader reader,
            SymbolTokenization symParser,
            SeqIOListener listener)
            throws IllegalSymbolException, IOException, ParseException {
        if (!(listener instanceof RichSeqIOListener)) throw new IllegalArgumentException("Only accepting RichSeqIOListeners today");
        return this.readRichSequence(reader,symParser,(RichSeqIOListener)listener,null);
    }
    
    private String sectionKey = null;
    private NCBITaxon tax = null;
    private String organism = null;
    private String accession = null;
    private String identifier = null;
    /**
     * {@inheritDoc}
     */
    public boolean readRichSequence(BufferedReader reader,
            SymbolTokenization symParser,
            RichSeqIOListener rlistener,
            Namespace ns)
            throws IllegalSymbolException, IOException, ParseException {
        
        sectionKey = null;
        tax = null;
        organism = null;
        accession = null;
        identifier = null;
        boolean hasAnotherSequence = true;
        //boolean hasInternalWhitespace = false;
        
        rlistener.startSequence();
        
        if (ns==null) ns=RichObjectFactory.getDefaultNamespace();
        rlistener.setNamespace(ns);
        
        // Get an ordered list of key->value pairs in array-tuples
        List section = null;
        try{
            do {
                section = this.readSection(reader);
                sectionKey = ((String[])section.get(0))[0];
                if(sectionKey == null){
                    String message = ParseException.newMessage(this.getClass(), accession, identifier, "Section key was null", sectionToString(section));
                    throw new ParseException(message);
                }
                // process section-by-section
                if (sectionKey.equals(LOCUS_TAG)) {
                    String loc = ((String[])section.get(0))[1];
                    Matcher m = lp.matcher(loc);
                    if (m.matches()) {
                        rlistener.setName(m.group(1));
                        accession = m.group(1); // default if no accession found
                        rlistener.setAccession(accession);
                        if (m.group(4)!=null)
                            rlistener.addSequenceProperty(Terms.getMolTypeTerm(),m.group(4));
                        // Optional extras
                        String stranded = m.group(3);
                        if(stranded!=null && stranded.equals("ss-"))
                            stranded = "single";
                        else if(stranded!=null && stranded.equals("ms-"))
                            stranded = "mixed";
                        else if(stranded!=null && stranded.equals("ds-"))
                            stranded = "double";
                        String circular = m.group(5);
                        String fifth = m.group(6);
                        String sixth = m.group(7);
                        if (stranded!=null) rlistener.addSequenceProperty(Terms.getStrandedTerm(),stranded);
                        if (circular!=null && circular.equalsIgnoreCase("circular")) rlistener.setCircular(true);
                        if (sixth != null) {
                            rlistener.setDivision(fifth);
                            rlistener.addSequenceProperty(Terms.getDateUpdatedTerm(),sixth);
                        } else if (fifth!=null) {
                            rlistener.addSequenceProperty(Terms.getDateUpdatedTerm(),fifth);
                        }
                    } else {
                        String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad locus line", sectionToString(section));
                        throw new ParseException(message);
                    }
                } else if (sectionKey.equals(DEFINITION_TAG)) {
                    rlistener.setDescription(((String[])section.get(0))[1]);
                } else if (sectionKey.equals(ACCESSION_TAG)) {
                    // if multiple accessions, store only first as accession,
                    // and store rest in annotation
                    String[] accs = ((String[])section.get(0))[1].split("\\s+");
                    accession = accs[0].trim();
                    rlistener.setAccession(accession);
                    for (int i = 1; i < accs.length; i++) {
                        rlistener.addSequenceProperty(Terms.getAdditionalAccessionTerm(),accs[i].trim());
                    }
                } else if (sectionKey.equals(VERSION_TAG)) {
                    String ver = ((String[])section.get(0))[1];
                    Matcher m = vp.matcher(ver);
                    if (m.matches()) {
                        String verAcc = m.group(1);
                        if (!accession.equals(verAcc)) {
                            // the version refers to a different accession!
                            // believe the version line, and store the original
                            // accession away in the additional accession set
                            rlistener.addSequenceProperty(Terms.getAdditionalAccessionTerm(),accession);
                            accession = verAcc;
                            rlistener.setAccession(accession);
                        }
                        if (m.group(3)!=null) rlistener.setVersion(Integer.parseInt(m.group(3)));
                        if (m.group(5)!=null) {
                            identifier = m.group(5);
                            rlistener.setIdentifier(identifier);
                        }
                    } else {
                        String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad version line", sectionToString(section));
                        throw new ParseException(message);
                    }
                } else if (sectionKey.equals(KEYWORDS_TAG)) {
                    String val = ((String[])section.get(0))[1];
                    if (val.endsWith(".")) val = val.substring(0, val.length()-1); // chomp dot
                    val = val.replace('\n',' '); //remove newline
                    String[] kws = val.split(";");
                    
                    for (int i = 0; i < kws.length; i++) {
                        String kw = kws[i].trim();
                        if (kw.length()==0) continue;
                        rlistener.addSequenceProperty(Terms.getKeywordTerm(), kw);
                    }
                } else if (sectionKey.equals(SOURCE_TAG)) {
                    // ignore - can get all this from the first feature
                } else if (sectionKey.equals(REFERENCE_TAG) && !this.getElideReferences()) {
                    // first line of section has rank and location
                    int ref_rank;
                    List baseRangeList=null;
                    String ref = ((String[])section.get(0))[1];
                    Matcher m = refp.matcher(ref);
                    if (m.matches()) {
                        ref_rank = Integer.parseInt(m.group(1));
                        if (m.group(3) != null) baseRangeList=buildBaseRanges(m.group(3));
                    } else {
                        String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad reference line", sectionToString(section));
                        throw new ParseException(message);
                    }
                    // rest can be in any order
                    String authors = null;
                    String consortium = null;
                    String title = null;
                    String journal = null;
                    String medline = null;
                    String pubmed = null;
                    String remark = null;
                    for (int i = 1; i < section.size(); i++) {
                        String key = ((String[])section.get(i))[0];
                        String val = ((String[])section.get(i))[1];
                        if (key.equals(AUTHORS_TAG)) authors = val.replace('\n',' '); //see #2276
                        else if (key.equals(CONSORTIUM_TAG)) consortium = val.replace('\n',' '); //see #2276
                        else if (key.equals(TITLE_TAG)) title = val.replace('\n',' '); //see #2276
                        else if (key.equals(JOURNAL_TAG)) journal = val.replace('\n',' '); //see #2276
                        else if (key.equals(MEDLINE_TAG)) medline = val;
                        else if (key.equals(PUBMED_TAG)) pubmed = val;
                        else if (key.equals(REMARK_TAG)) remark = val.replace('\n',' '); //see #2276
                    }
                    
                    // create the docref object
                    try {
                        // Use consortium as well if present.
                        if (authors==null) authors = consortium + " (consortium)";
                        else if (consortium!=null) authors = authors + ", " + consortium + " (consortium)";
                        // Create docref.
                        DocRef dr = (DocRef)RichObjectFactory.getObject(SimpleDocRef.class,new Object[]{DocRefAuthor.Tools.parseAuthorString(authors),journal,title});
                        // assign either the pubmed or medline to the docref - medline gets priority
                        if (medline!=null) dr.setCrossref((CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{Terms.MEDLINE_KEY, medline, new Integer(0)}));
                        else if (pubmed!=null) dr.setCrossref((CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{Terms.PUBMED_KEY, pubmed, new Integer(0)}));
                        // assign the remarks
                        if (!this.getElideComments()) dr.setRemark(remark);
                        // assign the docref to the bioentry: null if no base ranges, Integers if 1 base range - the normal case, joined RichLocation if more than 1
                        RankedDocRef rdr = baseRangeList == null?new SimpleRankedDocRef(dr, null, null, ref_rank):(baseRangeList.size()==1?new SimpleRankedDocRef(dr, new Integer(((RichLocation)baseRangeList.get(0)).getMin()), new Integer(((RichLocation)baseRangeList.get(0)).getMax()), ref_rank):new SimpleRankedDocRef(dr, new CompoundRichLocation(baseRangeList), ref_rank));
                        rlistener.setRankedDocRef(rdr);
                    } catch (ChangeVetoException e) {
                        throw new ParseException(e+", accession:"+accession);
                    }
                } else if (sectionKey.equals(COMMENT_TAG) && !this.getElideComments()) {
                    // Set up some comments
                    rlistener.setComment(((String[])section.get(0))[1]);
                } else if (sectionKey.equals(FEATURE_TAG) && !this.getElideFeatures()) {
                    // starting from second line of input, start a new feature whenever we come across
                    // a key that does not start with /
                    boolean seenAFeature = false;
                    int rcrossrefCount = 0;
                    boolean skippingBond = false;
                    for (int i = 1 ; i < section.size(); i++) {
                        String key = ((String[])section.get(i))[0];
                        String val = ((String[])section.get(i))[1];
                        if (key.startsWith("/")) {
                        	  if(!skippingBond)
                        	  {
	                            key = key.substring(1); // strip leading slash
	                            val = val.replaceAll("\\s*[\\n\\r]+\\s*"," ").trim();
	                            if (val.endsWith("\"")) val = val.substring(1,val.length()-1); // strip quotes
	                            // parameter on old feature
	                            if (key.equals("db_xref")) {
	                                Matcher m = dbxp.matcher(val);
	                                if (m.matches()) {
	                                    String dbname = m.group(1);
	                                    String raccession = m.group(2);
	                                    if (dbname.equalsIgnoreCase("taxon")) {
	                                        // Set the Taxon instead of a dbxref
	                                        tax = (NCBITaxon)RichObjectFactory.getObject(SimpleNCBITaxon.class, new Object[]{Integer.valueOf(raccession)});
	                                        rlistener.setTaxon(tax);
	                                        try {
	                                            if (organism!=null) tax.addName(NCBITaxon.SCIENTIFIC,organism.replace('\n', ' '));// readSection can embed new lines
	                                        } catch (ChangeVetoException e) {
	                                            throw new ParseException(e+", accession:"+accession);
	                                        }
	                                    } else {
	                                        try {
	                                            CrossRef cr = (CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{dbname, raccession, new Integer(0)});
	                                            RankedCrossRef rcr = new SimpleRankedCrossRef(cr, ++rcrossrefCount);
	                                            rlistener.getCurrentFeature().addRankedCrossRef(rcr);
	                                        } catch (ChangeVetoException e) {
	                                            throw new ParseException(e+", accession:"+accession);
	                                        }
	                                    }
	                                } else {
	                                    String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad dbxref", sectionToString(section));
	                                    throw new ParseException(message);
	                                }
	                            } else if (key.equalsIgnoreCase("organism")) {
	                                try {
	                                    organism = val;
	                                    if (tax!=null) tax.addName(NCBITaxon.SCIENTIFIC,organism.replace('\n', ' '));// readSection can embed new lines
	                                } catch (ChangeVetoException e) {
	                                    throw new ParseException(e+", accession:"+accession);
	                                }
	                            } else {
	                                if (key.equalsIgnoreCase("translation")) {
	                                    // strip spaces from sequence
	                                    val = val.replaceAll("\\s+","");
	                                }
	                                rlistener.addFeatureProperty(RichObjectFactory.getDefaultOntology().getOrCreateTerm(key),val);
	                            }
                          	}
                        } else {
                            // new feature!
                            // end previous feature
                            if(key.equalsIgnoreCase("bond"))
                            {
                            	skippingBond = true;
                            }
                            else
                            {
                            	skippingBond = false;
                            	if (seenAFeature) {
                            		rlistener.endFeature();
                            	}
	                            // start next one, with lots of lovely info in it
	                            RichFeature.Template templ = new RichFeature.Template();
	                            templ.annotation = new SimpleRichAnnotation();
	                            templ.sourceTerm = Terms.getGenBankTerm();
	                            templ.typeTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm(key);
	                            templ.featureRelationshipSet = new TreeSet();
	                            templ.rankedCrossRefs = new TreeSet();
	                            String tidyLocStr = val.replaceAll("\\s+","");
	                            templ.location = GenbankLocationParser.parseLocation(ns, accession, tidyLocStr);
	                            rlistener.startFeature(templ);
	                            seenAFeature = true;
	                            rcrossrefCount = 0;
                            }
                            
                        }
                    }
                    
                    if (seenAFeature) {
                    	rlistener.endFeature();
                    }
                } else if (sectionKey.equals(BASE_COUNT_TAG)) {
                    // ignore - can calculate from sequence content later if needed
                } else if (sectionKey.equals(START_SEQUENCE_TAG) && !this.getElideSymbols()) {
                    // our first line is ignorable as it is the ORIGIN tag
                    // the second line onwards conveniently have the number as
                    // the [0] tuple, and sequence string as [1] so all we have
                    // to do is concat the [1] parts and then strip out spaces,
                    // and replace '.' and '~' with '-' for our parser.
                    StringBuffer seq = new StringBuffer();
                    for (int i = 1 ; i < section.size(); i++) seq.append(((String[])section.get(i))[1]);
                    try {
                        SymbolList sl = new SimpleSymbolList(symParser,
                                seq.toString().replaceAll("\\s+","").replaceAll("[\\.|~]","-"));
                        rlistener.addSymbols(symParser.getAlphabet(),
                                (Symbol[])(sl.toList().toArray(new Symbol[0])),
                                0, sl.length());
                    } catch (IllegalAlphabetException e) {
                        String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad sequence section", sectionToString(section));
                        throw new ParseException(e, message);
                    }
                }
            } while (!sectionKey.equals(END_SEQUENCE_TAG));
        }catch(RuntimeException e){
            String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad sequence section", sectionToString(section));
            throw new ParseException(e, message);
        }
        
        // Allows us to tolerate trailing whitespace without
        // thinking that there is another Sequence to follow
        while (true) {
            reader.mark(1);
            int c = reader.read();
            if (c == -1) {
                hasAnotherSequence = false;
                break;
            }
            if (Character.isWhitespace((char) c)) {
                //hasInternalWhitespace = true;
                continue;
            }
            //if (hasInternalWhitespace)
            //    System.err.println("Warning: whitespace found between sequence entries");
            reader.reset();
            break;
        }
        
        // Finish up.
        rlistener.endSequence();
        return hasAnotherSequence;
    }
    
    // reads an indented section, combining split lines and creating a list of key->value tuples
    private List readSection(BufferedReader br) throws ParseException {
        List section = new ArrayList();
        String line;
        String currKey = null;
        StringBuffer currVal = new StringBuffer();
        boolean done = false;
        int linecount = 0;
        
        try {
            while (!done) {
                br.mark(320);
                line = br.readLine();
                String firstSecKey = section.size() == 0 ? "" : ((String[])section.get(0))[0];
                if (line==null || (!line.startsWith(" ") && linecount++>0 && ( !firstSecKey.equals(START_SEQUENCE_TAG)  || line.startsWith(END_SEQUENCE_TAG)))) {
                    // dump out last part of section
                    section.add(new String[]{currKey,currVal.toString()});
                    br.reset();
                    done = true;
                } else {
                    Matcher m = sectp.matcher(line);
                    if (m.matches()) {
                        // new key
                        if (currKey!=null) section.add(new String[]{currKey,currVal.toString()});
                        // key = group(2) or group(4) or group(6) - whichever is not null
                        currKey = m.group(2)==null?(m.group(4)==null?m.group(6):m.group(4)):m.group(2);
                        currVal = new StringBuffer();
                        // val = group(3) if group(2) not null, group(5) if group(4) not null, "" otherwise, trimmed
                        currVal.append((m.group(2)==null?(m.group(4)==null?"":m.group(5)):m.group(3)).trim());
                    } else {
                        // concatted line or SEQ START/END line?
                        if (line.startsWith(START_SEQUENCE_TAG) || line.startsWith(END_SEQUENCE_TAG)) currKey = line;
                        else {
                            currVal.append("\n"); // newline in between lines - can be removed later
                            currVal.append(currKey.charAt(0)=='/'?line.substring(21):line.substring(12));
                        }
                    }
                }
            }
        } catch (IOException e) {
            String message = ParseException.newMessage(this.getClass(), accession, identifier, "", sectionToString(section));
            throw new ParseException(e, message);
        } catch (RuntimeException e){
            String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad section", sectionToString(section));
            throw new ParseException(e, message);
        }
        return section;
    }
    
    private final List buildBaseRanges(final String theBaseRangeList) throws ParseException {
        if (theBaseRangeList == null) return null;
        final List baseRangeList = new ArrayList();
        final String[] baseRange = theBaseRangeList.split(";");
        try{
        for (int r=0; r<baseRange.length; r++) {
            final Matcher rangeMatch = refRange.matcher(baseRange[r]);
            if (rangeMatch.matches()) {
                final int rangeStart = Integer.parseInt(rangeMatch.group(1));
                final int rangeEnd = Integer.parseInt(rangeMatch.group(2));
                baseRangeList.add(new SimpleRichLocation(new SimplePosition(rangeStart), new SimplePosition(rangeEnd), r));
            } else {
                String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad reference range found", theBaseRangeList);
                throw new ParseException(message);
            }
        }
        return baseRangeList;
        }catch(RuntimeException e){
            String message = ParseException.newMessage(this.getClass(), accession, identifier, "Bad base range", theBaseRangeList);
            throw new ParseException(e, message);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void	writeSequence(Sequence seq, PrintStream os) throws IOException {
        if (this.getPrintStream()==null) this.setPrintStream(os);
        this.writeSequence(seq, RichObjectFactory.getDefaultNamespace());
    }
    
    /**
     * {@inheritDoc}
     */
    public void writeSequence(Sequence seq, String format, PrintStream os) throws IOException {
        if (this.getPrintStream()==null) this.setPrintStream(os);
        if (!format.equals(this.getDefaultFormat())) throw new IllegalArgumentException("Unknown format: "+format);
        this.writeSequence(seq, RichObjectFactory.getDefaultNamespace());
    }
    
    /**
     * {@inheritDoc}
     * Namespace is ignored as Genbank has no concept of it.
     */
    public void writeSequence(Sequence seq, Namespace ns) throws IOException {
        RichSequence rs;
        try {
            if (seq instanceof RichSequence) rs = (RichSequence)seq;
            else rs = RichSequence.Tools.enrich(seq);
        } catch (ChangeVetoException e) {
            IOException e2 = new IOException("Unable to enrich sequence");
            e2.initCause(e);
            throw e2;
        }
        
        SymbolTokenization tok;
        try {
            tok = rs.getAlphabet().getTokenization("token");
        } catch (Exception e) {
            throw new RuntimeException("Unable to get alphabet tokenizer",e);
        }
        Set notes = rs.getNoteSet();
        String accession = rs.getAccession();
        StringBuffer accessions = new StringBuffer();
        accessions.append(accession);
        String stranded = "";
        String udat = "";
        String moltype = rs.getAlphabet().getName();
        if ("PROTEIN-TERM".equals(moltype) || "PROTEIN".equals(moltype)) moltype = null; //a genpept curiosity
        StringBuffer keywords = new StringBuffer();
        for (Iterator i = notes.iterator(); i.hasNext(); ) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getStrandedTerm())) {
                String value = n.getValue();
                if(value != null && value.equals("single"))
                    stranded= "ss-";
                else if(value != null && value.equals("mixed"))
                    stranded= "ms-";
            }
            else if (n.getTerm().equals(Terms.getDateUpdatedTerm())) udat=n.getValue();
            else if (n.getTerm().equals(Terms.getMolTypeTerm())) moltype=n.getValue();
            else if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) {
                accessions.append(" ");
                accessions.append(n.getValue());
            } else if (n.getTerm().equals(Terms.getKeywordTerm())) {
                if (n.getValue() != null) {
                    if (keywords.length()>0) keywords.append("; ");
                    keywords.append(n.getValue());
                }
            }
        }
        
        //adjust molecule type during format conversion
        if(moltype!=null && moltype.length()>6) {
            if(moltype.indexOf("DNA")!=-1) moltype = "DNA";
            else if(moltype.indexOf("RNA")!=-1) moltype = "RNA";
            else moltype = "NA";
        }
        
        // locus(name) + length + alpha + div + date line
        StringBuffer locusLine = new StringBuffer();
        locusLine.append(StringTools.rightPad(rs.getName(),16));//13->28=15+1=16
        locusLine.append(" ");//29
        locusLine.append(StringTools.leftPad(""+rs.length(),11));//30->40=10+1=11
        locusLine.append(" "+ (moltype==null? "aa":"bp") +" ");//41->44
        locusLine.append(StringTools.leftPad(stranded,3));//45->47=2+1=3
        locusLine.append(StringTools.rightPad(moltype==null?"":moltype,6));//48->53=5+1=6
        locusLine.append("  ");//54->55
        locusLine.append(StringTools.rightPad(rs.getCircular()?"circular":"linear",8));//56->63=7+1=8
        locusLine.append(" ");//64->64
        locusLine.append(StringTools.rightPad(rs.getDivision()==null?"":rs.getDivision(),3));//65->67=2+1=3
        locusLine.append(" ");//68->68
        locusLine.append(StringTools.rightPad(udat,11));//69->79=10+1=11
        StringTools.writeKeyValueLine(LOCUS_TAG, locusLine.toString(), 12, this.getLineWidth(), this.getPrintStream());
        
        // definition line
        StringTools.writeKeyValueLine(DEFINITION_TAG, rs.getDescription(), 12, this.getLineWidth(), this.getPrintStream());
        
        // accession line
        StringTools.writeKeyValueLine(ACCESSION_TAG, accessions.toString(), 12, this.getLineWidth(), this.getPrintStream());
        
        // version + gi line
        String version = accession+"."+rs.getVersion();
        if (rs.getIdentifier()!=null) version = version + "  GI:"+rs.getIdentifier();
        StringTools.writeKeyValueLine(VERSION_TAG, version, 12, this.getLineWidth(), this.getPrintStream());
        
        // keywords line
        keywords.append(".");
        StringTools.writeKeyValueLine(KEYWORDS_TAG, keywords.toString(), 12, this.getLineWidth()-1, this.getPrintStream());
        
        // source line (from taxon)
        //   organism line
        NCBITaxon tax = rs.getTaxon();
        if (tax!=null) {
            StringTools.writeKeyValueLine(SOURCE_TAG, (isMitochondrial(rs)?"mitochondrion ":"")+tax.getDisplayName(), 12, this.getLineWidth(), this.getPrintStream());
            StringTools.writeKeyValueLine("  "+ORGANISM_TAG, tax.getDisplayName().split("\\s+\\(")[0]+"\n"+tax.getNameHierarchy(), 12, this.getLineWidth()-1, this.getPrintStream());
        }
        
        // references - rank (bases x to y)
        for (Iterator r = rs.getRankedDocRefs().iterator(); r.hasNext(); ) {
            RankedDocRef rdr = (RankedDocRef)r.next();
            DocRef d = rdr.getDocumentReference();
            StringTools.writeKeyValueLine(REFERENCE_TAG, rdr.getRank()+((rdr.getLocation()==null || rdr.getLocation() ==RichLocation.EMPTY_LOCATION)?"": (moltype==null? "  (residues ":"  (bases ")+makeBaseRange(rdr)+")"), 12, this.getLineWidth(), this.getPrintStream());
            // Any authors that were in the input as CONSRTM tags will
            // be merged into the AUTHORS tag on output.
            StringTools.writeKeyValueLine("  "+AUTHORS_TAG, d.getAuthors(), 12, this.getLineWidth()-1, this.getPrintStream());
            StringTools.writeKeyValueLine("  "+TITLE_TAG, d.getTitle(), 12, this.getLineWidth(), this.getPrintStream());
            StringTools.writeKeyValueLine("  "+JOURNAL_TAG, d.getLocation(), 12, this.getLineWidth(), this.getPrintStream());
            CrossRef c = d.getCrossref();
            if (c!=null) StringTools.writeKeyValueLine(StringTools.leftPad(c.getDbname(),9), c.getAccession(), 12, this.getLineWidth(), this.getPrintStream());
            StringTools.writeKeyValueLine("  "+REMARK_TAG, d.getRemark(), 12, this.getLineWidth(), this.getPrintStream());
        }
        
        // comments - if any
        Set comments = rs.getComments();
        if (!comments.isEmpty()) {
            StringBuffer sb = new StringBuffer();
            for (Iterator i = comments.iterator(); i.hasNext(); ) {
                Comment c = (SimpleComment)i.next();
                sb.append(c.getComment());
                if (i.hasNext()) sb.append("\n");
            }
            StringTools.writeKeyValueLine(COMMENT_TAG, sb.toString(), 12, this.getLineWidth(), this.getPrintStream());
        }
        
        this.getPrintStream().println(FEATURE_TAG+"             Location/Qualifiers");
        // feature_type     location
        for (Iterator i = rs.getFeatureSet().iterator(); i.hasNext(); ) {
            RichFeature f = (RichFeature)i.next();
            StringTools.writeKeyValueLine("     "+f.getTypeTerm().getName(), GenbankLocationParser.writeLocation((RichLocation)f.getLocation()), 21, this.getLineWidth()-1, ",", this.getPrintStream());
            for (Iterator j = f.getNoteSet().iterator(); j.hasNext(); ) {
                Note n = (Note)j.next();
                // /key="val" or just /key if val==""
                if (n.getValue()==null || n.getValue().length()==0) StringTools.writeKeyValueLine("", "/"+n.getTerm().getName(), 21, this.getLineWidth(), this.getPrintStream());
                else if (isNotQuoted(n)) {// doesn't have the value enclosed in quotes
                    StringTools.writeKeyValueLine("", "/"+n.getTerm().getName()+"="+n.getValue(), 21, this.getLineWidth(), this.getPrintStream());
                } else if (n.getTerm().getName().equals("translation")) {
                    StringTools.writeKeyValueLine("", "/"+n.getTerm().getName()+"=\""+n.getValue()+"\"", 21, this.getLineWidth()-1, this.getPrintStream());
                } else {
                    StringTools.writeKeyValueLine("", "/"+n.getTerm().getName()+"=\""+n.getValue()+"\"", 21, this.getLineWidth(), this.getPrintStream());
                }
            }
            // add-in to source feature only organism and db_xref="taxon:xyz" where present
            if (f.getType().equals("source") && tax!=null) {
                String displayName = tax.getDisplayName();
                if (displayName.indexOf('(')>-1) displayName = displayName.substring(0, displayName.indexOf('(')).trim();
                StringTools.writeKeyValueLine("", "/organism=\""+displayName+"\"", 21, this.getLineWidth()-1, this.getPrintStream());// AF252370 fits in exactly 80 - but is wrapped
                for (Iterator j = f.getRankedCrossRefs().iterator(); j.hasNext(); ) {
                    RankedCrossRef rcr = (RankedCrossRef)j.next();
                    CrossRef cr = rcr.getCrossRef();
                    StringTools.writeKeyValueLine("", "/db_xref=\""+cr.getDbname()+":"+cr.getAccession()+"\"", 21, this.getLineWidth(), this.getPrintStream());
                }
                StringTools.writeKeyValueLine("", "/db_xref=\"taxon:"+tax.getNCBITaxID()+"\"", 21, this.getLineWidth(), this.getPrintStream());
            } else {
                // add-in other dbxrefs where present
                for (Iterator j = f.getRankedCrossRefs().iterator(); j.hasNext(); ) {
                    RankedCrossRef rcr = (RankedCrossRef)j.next();
                    CrossRef cr = rcr.getCrossRef();
                    StringTools.writeKeyValueLine("", "/db_xref=\""+cr.getDbname()+":"+cr.getAccession()+"\"", 21, this.getLineWidth(), this.getPrintStream());
                }
            }
        }
        
        //BASE COUNT obsolete in Genbank flatfile format since October 2003
        //if (rs.getAlphabet()==AlphabetManager.alphabetForName("DNA")) {
        //    // BASE COUNT     1510 a   1074 c    835 g   1609 t
        //    int aCount = 0;
        //    int cCount = 0;
        //    int gCount = 0;
        //    int tCount = 0;
        //    int oCount = 0;
        //    for (int i = 1; i <= rs.length(); i++) {
        //        char c;
        //        try {
        //            c = tok.tokenizeSymbol(rs.symbolAt(i)).charAt(0);
        //        } catch (Exception e) {
        //            throw new RuntimeException("Unable to get symbol at position "+i,e);
        //        }
        //        switch (c) {
        //            case 'a': case 'A':
        //                aCount++;
        //                break;
        //            case 'c': case 'C':
        //                cCount++;
        //                break;
        //            case 'g': case 'G':
        //                gCount++;
        //                break;
        //            case 't': case 'T':
        //                tCount++;
        //                break;
        //            default:
        //                oCount++;
        //        }
        //    }
        //
        //    this.getPrintStream().print(BASE_COUNT_TAG_FULL+"    ");
        //    this.getPrintStream().print(aCount + " a   ");
        //    this.getPrintStream().print(cCount + " c   ");
        //    this.getPrintStream().print(gCount + " g   ");
        //    this.getPrintStream().print(tCount + " t    ");
        //    this.getPrintStream().println(oCount + " others");
        //}
        
        this.getPrintStream().println(START_SEQUENCE_TAG);
        // sequence stuff
        Symbol[] syms = (Symbol[])rs.toList().toArray(new Symbol[0]);
        int lines = 0;
        int symCount = 0;
        for (int i = 0; i < syms.length; i++) {
            if (symCount % 60 == 0) {
                if (lines > 0) this.getPrintStream().print("\n"); // newline from previous line
                int lineNum = (lines*60) + 1;
                this.getPrintStream().print(StringTools.leftPad(""+lineNum,9));
                lines++;
            }
            if (symCount % 10 == 0) this.getPrintStream().print(" ");
            try {
                this.getPrintStream().print(tok.tokenizeSymbol(syms[i]));
            } catch (IllegalSymbolException e) {
                throw new RuntimeException("Found illegal symbol: "+syms[i]);
            }
            symCount++;
        }
        if(syms.length>0) //do not create an empty line
            this.getPrintStream().print("\n");
        this.getPrintStream().println(END_SEQUENCE_TAG);
    }
    
    /**
     * {@inheritDoc}
     */
    public String getDefaultFormat() {
        return GENBANK_FORMAT;
    }
    
    private final static boolean isMitochondrial(final RichSequence theSequence) {
        final Set featureSet = theSequence.getFeatureSet();
        final Iterator i = featureSet.iterator();
        while (i.hasNext()) {
            final RichFeature feature = (RichFeature) i.next();
            if (feature.getType().equals("source")) {
                final Set noteSet = feature.getNoteSet();
                final Iterator n = noteSet.iterator();
                while(n.hasNext()) {
                    final Note note = (Note) n.next();
                    if (note.getTerm().getName().equals("organelle")) return note.getValue().equals("mitochondrion");
                }
            }
        }
        return false;
    }
    
    private final static boolean isNotQuoted(final Note theNote) {
        return isNotQuoted(theNote.getTerm().getName(), theNote.getValue());
    }
    
    private final static boolean isNotQuoted(final String theName, final String theValue) {
        return isNotQuoted.contains(theName);
    }
    
    private final static String makeBaseRange(final RankedDocRef theReference) {
        return theReference.getLocation()==null?theReference.getStart()+" to "+theReference.getEnd():toString(theReference.getLocation());
    }
    
    private final static String toString(final RichLocation theLocation) {
        final StringBuffer list = new StringBuffer();
        final Iterator b = theLocation.blockIterator();
        while (b.hasNext()) {
            final RichLocation location = (RichLocation) b.next();
            list.append(location.getMin()+" to "+location.getMax());
            if (b.hasNext()) list.append("; ");
        }
        return list.toString();
    }
    
    /**
     * Converts the current parse section to a String. Useful for debugging.
     */
    String sectionToString(List section){
        StringBuffer parseBlock = new StringBuffer();
        for(Iterator i = section.listIterator(); i.hasNext();){
            String[] part = (String[])i.next();
            for(int x = 0; x < part.length; x++){
                parseBlock.append(part[x]);
                if(x == 0){
                    parseBlock.append("   "); //the gap will have been trimmed
                }
            }
        }
        return parseBlock.toString();
    }
}
