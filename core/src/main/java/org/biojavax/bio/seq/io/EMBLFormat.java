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
import org.biojavax.RichAnnotation;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleComment;
import org.biojavax.SimpleCrossRef;
import org.biojavax.SimpleDocRef;
import org.biojavax.SimpleDocRefAuthor;
import org.biojavax.SimpleNote;
import org.biojavax.SimpleRankedCrossRef;
import org.biojavax.SimpleRankedDocRef;
import org.biojavax.SimpleRichAnnotation;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.taxa.NCBITaxon;
import org.biojavax.bio.taxa.SimpleNCBITaxon;
import org.biojavax.ontology.ComparableTerm;
import org.biojavax.utils.StringTools;

/**
 * Format reader for EMBL files. This version of EMBL format will generate
 * and write RichSequence objects. Loosely Based on code from the old, deprecated,
 * org.biojava.bio.seq.io.EmblLikeFormat object.
 * <p>
 * This format will read both Pre-87 and 87+ versions of EMBL. It will also write
 * them both. By default, it will write the most recent version. If you want
 * an earlier one, you must specify the format by passing one of the constants
 * defined in this class to {@link #writeSequence(Sequence, String, Namespace)}.
 *
 * @author Richard Holland
 * @author Jolyon Holdstock
 * @author Mark Schreiber
 * @since 1.5
 */
public class EMBLFormat extends RichSequenceFormat.HeaderlessFormat {
    
    // Register this format with the format auto-guesser.
    static {
        RichSequence.IOTools.registerFormat(EMBLFormat.class);
    }
    
    /**
     * The name of the Pre-87 format
     */
    public static final String EMBL_PRE87_FORMAT = "EMBL_PRE87";
    
    /**
     * The name of the current format
     */
    public static final String EMBL_FORMAT = "EMBL";
    
    protected static final String LOCUS_TAG = "ID";
    protected static final String ACCESSION_TAG = "AC";
    protected static final String VERSION_TAG = "SV";
    protected static final String DEFINITION_TAG = "DE";
    protected static final String DATE_TAG = "DT";
    protected static final String DATABASE_XREF_TAG = "DR";
    protected static final String SOURCE_TAG = "OS";
    protected static final String ORGANISM_TAG = "OC";
    protected static final String ORGANELLE_TAG = "OG";
    protected static final String REFERENCE_TAG = "RN";
    protected static final String REFERENCE_POSITION_TAG = "RP";
    protected static final String REFERENCE_XREF_TAG = "RX";
    protected static final String AUTHORS_TAG = "RA";
    protected static final String CONSORTIUM_TAG = "RG";
    protected static final String TITLE_TAG = "RT";
    protected static final String LOCATOR_TAG = "RL";
    protected static final String REMARK_TAG = "RC";
    protected static final String KEYWORDS_TAG = "KW";
    protected static final String COMMENT_TAG = "CC";
    protected static final String FEATURE_HEADER_TAG = "FH";
    protected static final String FEATURE_TAG = "FT";
    protected static final String CONTIG_TAG = "CO";
    protected static final String TPA_TAG = "AH";
    protected static final String START_SEQUENCE_TAG = "SQ";
    protected static final String DELIMITER_TAG = "XX";
    protected static final String END_SEQUENCE_TAG = "//";
    
    // the date pattern
    // date (Rel. N, Created)
    // date (Rel. N, Last updated, Version M)
    protected static final Pattern dp = Pattern.compile("([^\\s]+)\\s*(\\(Rel\\.\\s+(\\d+), ([^\\)\\d]+)(\\d*)\\))?$");
    // locus line
    protected static final Pattern lp = Pattern.compile("^(\\S+);\\s+SV\\s+(\\d+);\\s+(linear|circular);\\s+(\\S+\\s?\\S+?);\\s+(\\S+);\\s+(\\S+);\\s+(\\d+)\\s+(BP|AA)\\.$");
    protected static final Pattern lpPre87 = Pattern.compile("^(\\S+)\\s+standard;\\s+(circular)?\\s*(genomic)?\\s*(\\S+);\\s+(\\S+);\\s+\\d+\\s+BP\\.$");
    // version line
    protected static final Pattern vp = Pattern.compile("^(\\S+?)\\.(\\d+)$");
    // reference position line
    protected static final Pattern rpp = Pattern.compile("^(\\d+)(-(\\d+))?,?(\\s\\d+-\\d+,?)*$");
    // dbxref line
    protected static final Pattern dbxp = Pattern.compile("^([^:]+):(\\S+)$");
    
    protected static final Pattern readableFileNames = Pattern.compile(".*\\u002e(em|dat).*");
    protected static final Pattern headerLine = Pattern.compile("^ID.*");
    
    private NCBITaxon tax = null;
    private String organism = null;
    private String accession = null;
    
    /**
     * Implements some EMBL-specific terms.
     */
    public static class Terms extends RichSequence.Terms {
        
        /**
         * Getter for the RelUpdatedRecordVersion term
         * @return The RelUpdatedRecordVersion Term
         */
        public static ComparableTerm getRelUpdatedRecordVersionTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("RelUpdatedRecordVersion");
        }
        
        /**
         * Getter for the EMBL term
         * @return The EMBL Term
         */
        public static ComparableTerm getEMBLTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("EMBL");
        }
        
        /**
         * Getter for the Ensembl-specific 'genomic' term
         * @return The genomic Term
         */
        public static ComparableTerm getGenomicTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("genomic");
        }
        
        /**
         * Getter for the Ensembl-specific 'versionLine' term
         * @return The version line Term
         */
        public static ComparableTerm getVersionLineTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("versionLine");
        }
        
        /**
         * Getter for the Ensembl-specific 'dataClass' term
         * @return The data class Term
         */
        public static ComparableTerm getDataClassTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("dataClass");
        }
    }
    
    /**
     * {@inheritDoc}
     * A file is in EMBL format if its name contains the word eem or edat, or the first line matches
     * the EMBL format for the ID line.
     */
    public boolean canRead(File file) throws IOException {
        if (readableFileNames.matcher(file.getName()).matches()) return true;
        BufferedReader br = new BufferedReader(new FileReader(file));
        String firstLine = br.readLine();
        boolean readable = firstLine!=null && headerLine.matcher(firstLine).matches() &&
                (lp.matcher(firstLine.substring(3).trim()).matches() ||
                lpPre87.matcher(firstLine.substring(3).trim()).matches()
                );
        br.close();
        return readable;
    }
    
    /**
     * {@inheritDoc}
     * Always returns a DNA tokenizer.
     */
    public SymbolTokenization guessSymbolTokenization(File file) throws IOException {
        return RichSequence.IOTools.getDNAParser();
    }
    
    /**
     * {@inheritDoc}
     * A stream is in EMBL format if its first line matches the EMBL format for the ID line.
     */
    public boolean canRead(BufferedInputStream stream) throws IOException {
        stream.mark(2000); // some streams may not support this
        BufferedReader br = new BufferedReader(new InputStreamReader(stream));
        String firstLine = br.readLine();
        boolean readable = firstLine!=null && headerLine.matcher(firstLine).matches() &&
                (lp.matcher(firstLine.substring(3).trim()).matches() ||
                lpPre87.matcher(firstLine.substring(3).trim()).matches()
                );
        // don't close the reader as it'll close the stream too.
        // br.close();
        stream.reset();
        return readable;
    }
    
    /**
     * {@inheritDoc}
     * Always returns a DNA tokenizer.
     */
    public SymbolTokenization guessSymbolTokenization(BufferedInputStream stream) throws IOException {
        return RichSequence.IOTools.getDNAParser();
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
    
    /**
     * {@inheritDoc}
     */
    public boolean readRichSequence(BufferedReader reader,
            SymbolTokenization symParser,
            RichSeqIOListener rlistener,
            Namespace ns)
            throws IllegalSymbolException, IOException, ParseException {
        tax = null;
        organism = null;
        accession = null;
        boolean hasAnotherSequence = true;
        //boolean hasInternalWhitespace = false;
        
        rlistener.startSequence();
        
        if (ns==null) ns=RichObjectFactory.getDefaultNamespace();
        rlistener.setNamespace(ns);
        
        // Get an ordered list of key->value pairs in array-tuples
        String sectionKey = null;
        do {
            List section = this.readSection(reader);
            sectionKey = ((String[])section.get(0))[0];
            if(sectionKey == null){
                
                String message = ParseException.newMessage(this.getClass(), accession, "No section key", "Not set", sectionToString(section));
                throw new ParseException(message);
            }
            // process section-by-section
            if (sectionKey.equals(LOCUS_TAG)) {
                // entryname  dataclass; [circular] molecule; division; sequencelength BP.
                String loc = ((String[])section.get(0))[1];
                Matcher m = lp.matcher(loc);
                Matcher mPre87 = lpPre87.matcher(loc);
                if (m.matches()) {
                    // first token is both name and primary accession
                    rlistener.setName(m.group(1));
                    rlistener.setAccession(m.group(1));
                    // second token is version
                    rlistener.setVersion(Integer.parseInt(m.group(2)));
                    // third token is circular/linear
                    rlistener.setCircular(m.group(3).equals("circular"));
                    // fourth token is moltype
                    rlistener.addSequenceProperty(Terms.getMolTypeTerm(),m.group(4));
                    // fifth token is data class
                    rlistener.addSequenceProperty(Terms.getDataClassTerm(),m.group(5));
                    // sixth token is taxonomic division
                    rlistener.setDivision(m.group(6));
                    // seventh token is sequence length, which is ignored
                    // as it is calculated from the sequence data later.
                } else if (mPre87.matches()) {
                    rlistener.setName(mPre87.group(1));
                    if (mPre87.group(3)!=null) {
                        // add annotation for 'genomic' (Ensembl-specific term)
                        rlistener.addSequenceProperty(Terms.getGenomicTerm(),null);
                    }
                    rlistener.addSequenceProperty(Terms.getMolTypeTerm(),mPre87.group(4));
                    rlistener.setDivision(mPre87.group(5));
                    // Optional extras
                    String circular = mPre87.group(2);
                    if (circular!=null) rlistener.setCircular(true);
                } else {
                    String message = ParseException.newMessage(this.getClass(),accession,"Not Set","Bad ID line found", sectionToString(section));
                    throw new ParseException(message);
                }
            } else if (sectionKey.equals(DEFINITION_TAG)) {
                rlistener.setDescription(((String[])section.get(0))[1]);
            } else if (sectionKey.equals(SOURCE_TAG)) {
                // only interested in organelle sub-tag
                for (int i = 1; i < section.size(); i++) {
                    sectionKey = ((String[])section.get(i))[0];
                    if (sectionKey.equals(ORGANELLE_TAG)) {
                        rlistener.addSequenceProperty(Terms.getOrganelleTerm(), ((String[])section.get(i))[1].trim());
                        break; // skip out of for loop once found
                    }
                }
            } else if (sectionKey.equals(DATE_TAG)) {
                String chunk = ((String[])section.get(0))[1].trim();
                Matcher dm = dp.matcher(chunk);
                if (dm.matches()) {
                    String date = dm.group(1);
                    String rel = dm.group(3);
                    String type = dm.group(4);
                    if (type.equals("Created")) {
                        rlistener.addSequenceProperty(Terms.getDateCreatedTerm(), date);
                        rlistener.addSequenceProperty(Terms.getRelCreatedTerm(), rel);
                    } else if (type.equals("Last updated, Version ")) {
                        rlistener.addSequenceProperty(Terms.getDateUpdatedTerm(), date);
                        rlistener.addSequenceProperty(Terms.getRelUpdatedTerm(), rel);
                        rlistener.addSequenceProperty(Terms.getRelUpdatedRecordVersionTerm(), dm.group(5));
                    } else {
                        String message = ParseException.newMessage(this.getClass(),accession,"not set", "Bad date type found",sectionToString(section));
                        throw new ParseException(message);
                    }
                } else {
                    String message = ParseException.newMessage(this.getClass(),accession,"not set", "Bad date line found",sectionToString(section));
                    throw new ParseException(message);
                    
                }
            } else if (sectionKey.equals(ACCESSION_TAG)) {
                // if multiple accessions, store only first as accession,
                // and store rest in annotation
                String[] accs = ((String[])section.get(0))[1].split(";");
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
                    rlistener.setVersion(Integer.parseInt(m.group(2)));
                } else {
                    rlistener.addSequenceProperty(Terms.getVersionLineTerm(),ver);
                }
            } else if (sectionKey.equals(KEYWORDS_TAG)) {
                String val = ((String[])section.get(0))[1];
                val = val.substring(0,val.length()-1); // chomp dot
                val = val.replace('\n',' '); //remove newline
                String[] kws = val.split(";");
                for (int i = 0; i < kws.length; i++) {
                    String kw = kws[i].trim();
                    if (kw.length()==0) continue;
                    rlistener.addSequenceProperty(Terms.getKeywordTerm(), kw);
                }
            } else if (sectionKey.equals(DATABASE_XREF_TAG)) {
                String val = ((String[])section.get(0))[1];
                val = val.substring(0,val.length()-1); // chomp dot
                // database_identifier; primary_identifier; secondary_identifier....
                String[] parts = val.split(";");
                // construct a DBXREF out of the dbname part[0] and accession part[1]
                CrossRef crossRef = (CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{parts[0].trim(),parts[1].trim(), new Integer(0)});
                // assign remaining bits of info as annotations
                for (int j = 2; j < parts.length; j++) {
                    Note note = new SimpleNote(Terms.getAdditionalAccessionTerm(),parts[j].trim(),j-1);
                    try {
                        crossRef.getRichAnnotation().addNote(note);
                    } catch (ChangeVetoException ce) {
                        String message = ParseException.newMessage(this.getClass(),accession,"not set", "Could not annotate identifier terms",sectionToString(section));
                        ParseException pe = new ParseException(message);
                        pe.initCause(ce);
                        throw pe;
                    }
                }
                RankedCrossRef rcrossRef = new SimpleRankedCrossRef(crossRef, 0);
                rlistener.setRankedCrossRef(rcrossRef);
            } else if (sectionKey.equals(REFERENCE_TAG) && !this.getElideReferences()) {
                // first line of section has rank and location
                String refrank = ((String[])section.get(0))[1];
                int ref_rank = Integer.parseInt(refrank.substring(1,refrank.length()-1));
                int ref_start = -999;
                int ref_end = -999;
                // rest can be in any order
                String consortium = null;
                String authors = "";
                String title = null;
                String locator = null;
                String pubmed = null;
                String medline = null;
                String doi = null;
                String remark = null;
                for (int i = 1; i < section.size(); i++) {
                    String key = ((String[])section.get(i))[0];
                    String val = ((String[])section.get(i))[1];
                    if (key.equals(AUTHORS_TAG)) {
                        if (val.endsWith(";")) val = val.substring(0,val.length()-1); // chomp semicolon
                        authors = val.replace('\n',' '); //see #2276
                    }
                    if (key.equals(CONSORTIUM_TAG)) {
                        if (val.endsWith(";")) val = val.substring(0,val.length()-1); // chomp semicolon
                        consortium = val.replace('\n',' '); //see #2276
                    }
                    if (key.equals(TITLE_TAG)) {
                        if (val.length()>1) {
                            if (val.endsWith(";")) val = val.substring(0,val.length()-1); // chomp semicolon
                            if (val.endsWith("\"")) val = val.substring(1,val.length()-1); // chomp quotes
                            title = val.replace('\n',' '); //see #2276
                        } else title=null; // single semi-colon indicates no title
                    }
                    if (key.equals(LOCATOR_TAG)) {
                        if (val.endsWith(".")) val = val.substring(0,val.length()-1); // chomp dot
                        locator = val.replace('\n',' '); //see #2276
                    }
                    if (key.equals(REFERENCE_XREF_TAG)) {
                        // database_identifier; primary_identifier.
                        String[] refs = val.split("\\.(\\s+|$)");
                        for (int j = 0 ; j < refs.length; j++) {
                            if (refs[j].trim().length()==0) continue;
                            String[] parts = refs[j].split(";");
                            String db = parts[0];
                            String ref = parts[1].trim();
                            if (db.equalsIgnoreCase(Terms.PUBMED_KEY)) pubmed = ref;
                            else if (db.equalsIgnoreCase(Terms.MEDLINE_KEY)) medline = ref;
                            else if (db.equalsIgnoreCase(Terms.DOI_KEY)) doi = ref;
                        }
                    }
                    if (key.equals(REMARK_TAG)) remark = val.replace('\n',' '); //see #2276
                    if (key.equals(REFERENCE_POSITION_TAG)) {
                        // only the first group is taken
                        // if we have multiple lines, only the last line is taken
                        Matcher m = rpp.matcher(val);
                        if (m.matches()) {
                            ref_start = Integer.parseInt(m.group(1));
                            if(m.group(2) != null)
                                ref_end = Integer.parseInt(m.group(3));
                        } else {
                            String message = ParseException.newMessage(this.getClass(),accession,"not set", "Bad reference line found",sectionToString(section));
                            throw new ParseException(message);
                        }
                    }
                }
                // create the docref object
                try {
                    List authSet = DocRefAuthor.Tools.parseAuthorString(authors);
                    if (consortium!=null) authSet.add(new SimpleDocRefAuthor(consortium, true, false));
                    DocRef dr = (DocRef)RichObjectFactory.getObject(SimpleDocRef.class,new Object[]{authSet,locator,title});
                    // assign either the pubmed or medline to the docref - medline gets priority, then pubmed, then doi
                    if (medline!=null) dr.setCrossref((CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{Terms.MEDLINE_KEY, medline, new Integer(0)}));
                    else if (pubmed!=null) dr.setCrossref((CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{Terms.PUBMED_KEY, pubmed, new Integer(0)}));
                    else if (doi!=null) dr.setCrossref((CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{Terms.DOI_KEY, doi, new Integer(0)}));
                    // assign the remarks
                    if (!this.getElideComments()) dr.setRemark(remark);
                    // assign the docref to the bioentry
                    RankedDocRef rdr = new SimpleRankedDocRef(dr,
                            (ref_start != -999 ? new Integer(ref_start) : null),
                            (ref_end != -999 ? new Integer(ref_end) : null),
                            ref_rank);
                    rlistener.setRankedDocRef(rdr);
                } catch (ChangeVetoException e) {
                    String message = ParseException.newMessage(this.getClass(),accession,"not set", "",sectionToString(section));
                    throw new ParseException(e, message);
                }
            } else if (sectionKey.equals(COMMENT_TAG) && !this.getElideComments()) {
                // Set up some comments
                rlistener.setComment(((String[])section.get(0))[1]);
            } else if (sectionKey.equals(FEATURE_TAG) && !this.getElideFeatures()) {
                // starting from second line of input, start a new feature whenever we come across
                // a key that does not start with /
                boolean seenAFeature = false;
                int rcrossrefCount = 0;
                for (int i = 1 ; i < section.size(); i++) {
                    String key = ((String[])section.get(i))[0];
                    String val = ((String[])section.get(i))[1];
                    if (key.startsWith("/")) {
                        key = key.substring(1); // strip leading slash
                        val = val.replaceAll("\\s*[\\n\\r]+\\s*"," ").trim();
                        if (val.startsWith("\"")) val = val.substring(1,val.length()-1); // strip quotes
                        // parameter on old feature
                        if (key.equalsIgnoreCase("db_xref")) {
                            Matcher m = dbxp.matcher(val);
                            if (m.matches()) {
                                String dbname = m.group(1);
                                String raccession = m.group(2);
                                if (dbname.equalsIgnoreCase("taxon")) {
                                    // Set the Taxon instead of a dbxref
                                    tax = (NCBITaxon)RichObjectFactory.getObject(SimpleNCBITaxon.class, new Object[]{Integer.valueOf(raccession)});
                                    rlistener.setTaxon(tax);
                                    try {
                                        if (organism!=null) tax.addName(NCBITaxon.SCIENTIFIC,organism);
                                    } catch (ChangeVetoException e) {
                                        String message = ParseException.newMessage(this.getClass(),accession,"not set", "",sectionToString(section));
                                        throw new ParseException(e, message);
                                    }
                                } else {
                                    try {
                                        CrossRef cr = (CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{dbname, raccession, new Integer(0)});
                                        RankedCrossRef rcr = new SimpleRankedCrossRef(cr, ++rcrossrefCount);
                                        rlistener.getCurrentFeature().addRankedCrossRef(rcr);
                                    } catch (ChangeVetoException e) {
                                        String message = ParseException.newMessage(this.getClass(),accession,"not set", "",sectionToString(section));
                                        throw new ParseException(e, message);
                                    }
                                }
                            } else {
                                String message = ParseException.newMessage(this.getClass(),accession,"not set", "Bad dbxref found",sectionToString(section));
                                throw new ParseException(message);
                            }
                        } else if (key.equalsIgnoreCase("organism")) {
                            try {
                                organism = val;
                                if (tax!=null) tax.addName(NCBITaxon.SCIENTIFIC,organism);
                            } catch (ChangeVetoException e) {
                                String message = ParseException.newMessage(this.getClass(),accession,"not set", "",sectionToString(section));
                                throw new ParseException(message);
                            }
                        } else {
                            if (key.equalsIgnoreCase("translation")) {
                                // strip spaces from sequence
                                val = val.replaceAll("\\s+","");
                            }
                            rlistener.addFeatureProperty(RichObjectFactory.getDefaultOntology().getOrCreateTerm(key),val);
                        }
                    } else {
                        // new feature!
                        // end previous feature
                        if (seenAFeature) rlistener.endFeature();
                        // start next one, with lots of lovely info in it
                        RichFeature.Template templ = new RichFeature.Template();
                        templ.annotation = new SimpleRichAnnotation();
                        templ.sourceTerm = Terms.getEMBLTerm();
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
                if (seenAFeature) rlistener.endFeature();
            } else if (sectionKey.equals(START_SEQUENCE_TAG) && !this.getElideSymbols()) {
                StringBuffer seq = new StringBuffer();
                for (int i = 0 ; i < section.size(); i++) seq.append(((String[])section.get(i))[1]);
                try {
                    SymbolList sl = new SimpleSymbolList(symParser,
                            seq.toString().replaceAll("\\s+","").replaceAll("[\\.|~]","-"));
                    rlistener.addSymbols(symParser.getAlphabet(),
                            (Symbol[])(sl.toList().toArray(new Symbol[0])),
                            0, sl.length());
                } catch (Exception e) {
                    String message = ParseException.newMessage(this.getClass(),accession,"not set", "Bad sequence",sectionToString(section));
                    throw new ParseException(e, message);
                }
            }
        } while (!sectionKey.equals(END_SEQUENCE_TAG));
        
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
        boolean done = false;
        
        // while not done
        try {
            while (!done) {
                // mark buffer
                br.mark(160);
                // read token
                line = br.readLine();
                if (line.length()<2) {
                    String message = ParseException.newMessage(this.getClass(),accession,"not set", "Bad line found",line);
                    throw new ParseException(message);
                }
                String token = line.substring(0,2);
                // READ SEQUENCE SECTION
                if (token.equals(START_SEQUENCE_TAG)) {
                    //      from next line, read sequence until // - leave // on stack
                    StringBuffer sb = new StringBuffer();
                    while (!done) {
                        br.mark(160);
                        line = br.readLine();
                        if (line.startsWith(END_SEQUENCE_TAG)) {
                            br.reset();
                            done = true;
                        } else {
                            //      create sequence tag->value pair to return, sans numbers
                            sb.append(line.replaceAll("\\d",""));
                        }
                    }
                    section.add(new String[]{START_SEQUENCE_TAG,sb.toString()});
                }
                // READ FEATURE TABLE SECTION
                else if (token.equals(FEATURE_HEADER_TAG)) {
                    //      create dummy feature tag->value pair and add to return set
                    section.add(new String[]{FEATURE_TAG,null});
                    //      drop next FH line
                    line = br.readLine(); // skip next line too - it is also FH
                    //      read all FT lines until XX
                    String currentTag = null;
                    StringBuffer currentVal = null;
                    while (!done) {
                        line = br.readLine();
                        if (line.startsWith(DELIMITER_TAG)) {
                            done = true;
                            // dump current tag if exists
                            if (currentTag!=null) section.add(new String[]{currentTag,currentVal.toString()});
                        } else {
                            //         FT lines:   FT   word            value
                            //         or          FT                   /word
                            //         or          FT                   /db_xref="taxon:3899....
                            //                                          ......"
                            line = line.substring(5); // chomp off "FT   "
                            if (!line.startsWith(" ")) {
                                // dump current tag if exists
                                if (currentTag!=null) section.add(new String[]{currentTag,currentVal.toString()});
                                // case 1 : word value - splits into key-value on its own
                                String[] parts = line.trim().split("\\s+");
                                currentTag = parts[0];
                                currentVal = new StringBuffer();
                                currentVal.append(parts[1]);
                            } else {
                                line = line.trim();
                                if (line.startsWith("/")) {
                                    // dump current tag if exists
                                    if (currentTag!=null) section.add(new String[]{currentTag,currentVal.toString()});
                                    // case 2 : /word[=.....]
                                    currentVal = new StringBuffer();
                                    int equalIndex = line.indexOf('=');
                                    if (equalIndex>=0) {
                                        currentTag = line.substring(0, equalIndex);
                                        currentVal.append(line.substring(equalIndex+1));
                                    } else {
                                        currentTag = line;
                                    }
                                } else {
                                    // case 3 : ...."
                                    currentVal.append("\n");
                                    currentVal.append(line);
                                }
                            }
                        }
                    }
                }
                // READ END OF SEQUENCE
                else if (token.equals(END_SEQUENCE_TAG)) {
                    section.add(new String[]{END_SEQUENCE_TAG,null});
                    done = true;
                }
                // READ DELIMITER TAG
                else if (token.equals(DELIMITER_TAG)) {
                    section.add(new String[]{DELIMITER_TAG,null});
                    done = true;
                }
                // READ THIRD PARTY ANNOTATION SECTION
                else if (token.equals(TPA_TAG)) {
                    //      exception = don't know how to do TPA yet
                    String message = ParseException.newMessage(this.getClass(),accession,"not set", "Unable to handle TPAs just yet",sectionToString(section));
                    throw new ParseException(message);
                }
                // READ CONTIG SECTION
                else if (token.equals(CONTIG_TAG)) {
                    //      exception = don't know how to do contigs yet
                    String message = ParseException.newMessage(this.getClass(),accession,"not set", "Unable to handle contig assemblies just yet",sectionToString(section));
                    throw new ParseException(message);
                }
                // READ DOCREF
                else if (token.equals(DATABASE_XREF_TAG)) {
                    section.add(new String[]{DATABASE_XREF_TAG,line.substring(5).trim()});
                    done = true;
                }
                // READ DATE
                else if (token.equals(DATE_TAG)) {
                    section.add(new String[]{DATE_TAG,line.substring(5).trim()});
                    done = true;
                }
                // READ NORMAL TAG/VALUE SECTION
                else {
                    //      rewind buffer to mark
                    br.reset();
                    //      read token/values until XX
                    String currentTag = null;
                    StringBuffer currentVal = null;
                    while (!done) {
                        line = br.readLine();
                        if (line.startsWith(DELIMITER_TAG)) {
                            done = true;
                            // dump current tag if exists
                            if (currentTag!=null) section.add(new String[]{currentTag,currentVal.toString()});
                        } else {
                            try {
                                //      merge neighbouring repeated tokens by concatting values
                                //      return tag->value pairs
                                String tag = line.substring(0,2);
                                String value = line.substring(5);
                                if (currentTag==null || !tag.equals(currentTag)) {
                                    // dump current tag if exists
                                    if (currentTag!=null) section.add(new String[]{currentTag,currentVal.toString()});
                                    // start new tag
                                    currentTag = tag;
                                    currentVal = new StringBuffer();
                                    currentVal.append(value);
                                } else {
                                    currentVal.append("\n");
                                    currentVal.append(value);
                                }
                            } catch (Exception e) {
                                String message = ParseException.newMessage(this.getClass(), accession, "not set","",sectionToString(section));
                                throw new ParseException(e, message);
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            String message = ParseException.newMessage(this.getClass(),accession,"not set", "Unable to handle TPAs just yet",sectionToString(section));
            throw new ParseException(message);
        }
        return section;
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
        this.writeSequence(seq, format, RichObjectFactory.getDefaultNamespace());
    }
    
    /**
     * {@inheritDoc}
     * Namespace is ignored as EMBL has no concept of it.
     */
    public void writeSequence(Sequence seq, Namespace ns) throws IOException {
        this.writeSequence(seq, this.getDefaultFormat(), ns);
    }
    
    /**
     * As per {@link #writeSequence(Sequence, Namespace)}, except
     * that it also takes a format parameter. This can be any of the formats
     * defined as constants in this class.
     * @param seq see {@link #writeSequence(Sequence, Namespace)}
     * @param format the format to use.
     * @param ns see {@link #writeSequence(Sequence, Namespace)}
     * @throws IOException see {@link #writeSequence(Sequence, Namespace)}
     */
    public void writeSequence(Sequence seq, String format, Namespace ns) throws IOException {
        if (!format.equals(EMBL_FORMAT) && !format.equals(EMBL_PRE87_FORMAT))
            throw new IllegalArgumentException("Format "+format+" not recognised.");
        
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
        accessions.append(";");
        String cdat = null;
        String udat = null;
        String crel = null;
        String urel = null;
        String urecv = null;
        String organelle = null;
        String versionLine = null;
        String dataClass = "STD";
        boolean genomic = false;
        String moltype = rs.getAlphabet().getName();
        for (Iterator i = notes.iterator(); i.hasNext(); ) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getDateCreatedTerm())) cdat=n.getValue();
            else if (n.getTerm().equals(Terms.getDateUpdatedTerm())) udat=n.getValue();
            else if (n.getTerm().equals(Terms.getRelCreatedTerm())) crel=n.getValue();
            else if (n.getTerm().equals(Terms.getRelUpdatedTerm())) urel=n.getValue();
            else if (n.getTerm().equals(Terms.getRelUpdatedRecordVersionTerm())) urecv=n.getValue();
            else if (n.getTerm().equals(Terms.getMolTypeTerm())) moltype=n.getValue();
            else if (n.getTerm().equals(Terms.getVersionLineTerm())) versionLine=n.getValue();
            else if (n.getTerm().equals(Terms.getGenomicTerm())) genomic = true;
            else if (n.getTerm().equals(Terms.getDataClassTerm())) dataClass = n.getValue();
            else if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) {
                accessions.append(" ");
                accessions.append(n.getValue());
                accessions.append(";");
            } else if (n.getTerm().equals(Terms.getOrganelleTerm())) organelle=n.getValue();
        }
        
        StringBuffer locusLine = new StringBuffer();
        if (format.equals(EMBL_FORMAT)) {
            // accession; SV version; circular/linear; moltype; dataclass; division; length BP.
            locusLine.append(rs.getAccession());
            locusLine.append("; SV ");
            locusLine.append(rs.getVersion());
            locusLine.append("; ");
            locusLine.append(rs.getCircular()?"circular":"linear");
            locusLine.append("; ");
            locusLine.append(moltype);
            locusLine.append("; ");
            locusLine.append(dataClass);
            locusLine.append("; ");
            locusLine.append(rs.getDivision());
            locusLine.append("; ");
            locusLine.append(rs.length());
            locusLine.append(" BP.");
        } else if (format.equals(EMBL_PRE87_FORMAT)) {
            // entryname  dataclass; [circular] molecule; division; sequencelength BP.
            locusLine.append(StringTools.rightPad(rs.getName(),9));
            locusLine.append(" standard; ");
            locusLine.append(rs.getCircular()?"circular ":"");
            // if it is Ensembl genomic, add that in too
            if (genomic==true) locusLine.append("genomic ");
            locusLine.append(moltype);
            locusLine.append("; ");
            locusLine.append(rs.getDivision()==null?"":rs.getDivision());
            locusLine.append("; ");
            locusLine.append(rs.length());
            locusLine.append(" BP.");
        }
        StringTools.writeKeyValueLine(LOCUS_TAG, locusLine.toString(), 5, this.getLineWidth(), null, LOCUS_TAG, this.getPrintStream());
        this.getPrintStream().println(DELIMITER_TAG+"   ");
        
        // accession line
        StringTools.writeKeyValueLine(ACCESSION_TAG, accessions.toString(), 5, this.getLineWidth(), null, ACCESSION_TAG, this.getPrintStream());
        this.getPrintStream().println(DELIMITER_TAG+"   ");
        
        // version line
        if (format.equals(EMBL_PRE87_FORMAT)) {
            if (versionLine!=null) StringTools.writeKeyValueLine(VERSION_TAG, versionLine, 5, this.getLineWidth(), null, VERSION_TAG, this.getPrintStream());
            else StringTools.writeKeyValueLine(VERSION_TAG, accession+"."+rs.getVersion(), 5, this.getLineWidth(), null, VERSION_TAG, this.getPrintStream());
            this.getPrintStream().println(DELIMITER_TAG+"   ");
        }
        
        // date line
        StringTools.writeKeyValueLine(DATE_TAG, (cdat==null?udat:cdat)+" (Rel. "+(crel==null?"0":crel)+", Created)", 5, this.getLineWidth(), null, DATE_TAG, this.getPrintStream());
        StringTools.writeKeyValueLine(DATE_TAG, udat+" (Rel. "+(urel==null?"0":urel)+", Last updated, Version "+(urecv==null?"0":urecv)+")", 5, this.getLineWidth(), null, DATE_TAG, this.getPrintStream());
        this.getPrintStream().println(DELIMITER_TAG+"   ");
        
        // definition line
        StringTools.writeKeyValueLine(DEFINITION_TAG, rs.getDescription(), 5, this.getLineWidth(), null, DEFINITION_TAG, this.getPrintStream());
        this.getPrintStream().println(DELIMITER_TAG+"   ");
        
        // keywords line
        StringBuffer keywords = new StringBuffer();
        for (Iterator n = notes.iterator(); n.hasNext(); ) {
            Note nt = (Note)n.next();
            if (nt.getTerm().equals(Terms.getKeywordTerm())) {
                if (keywords.length()>0) keywords.append("; ");
                keywords.append(nt.getValue());
            }
        }
        if (keywords.length()>0) {
            keywords.append(".");
            StringTools.writeKeyValueLine(KEYWORDS_TAG, keywords.toString(), 5, this.getLineWidth(), null, KEYWORDS_TAG, this.getPrintStream());
            this.getPrintStream().println(DELIMITER_TAG+"   ");
        } else {
            this.getPrintStream().println(KEYWORDS_TAG+"   .");
            this.getPrintStream().println(DELIMITER_TAG+"   ");
        }
        
        // source line (from taxon)
        //   organism line
        NCBITaxon tax = rs.getTaxon();
        if (tax!=null) {
            StringTools.writeKeyValueLine(SOURCE_TAG, tax.getDisplayName(), 5, this.getLineWidth(), null, SOURCE_TAG, this.getPrintStream());
            StringTools.writeKeyValueLine(ORGANISM_TAG, tax.getNameHierarchy(), 5, this.getLineWidth(), null, SOURCE_TAG, this.getPrintStream());
            if (organelle!=null) StringTools.writeKeyValueLine(ORGANELLE_TAG, organelle, 5, this.getLineWidth(), null, ORGANELLE_TAG, this.getPrintStream());
            this.getPrintStream().println(DELIMITER_TAG+"   ");
        }
        
        // references - rank (bases x to y)
        for (Iterator r = rs.getRankedDocRefs().iterator(); r.hasNext(); ) {
            RankedDocRef rdr = (RankedDocRef)r.next();
            DocRef d = rdr.getDocumentReference();
            // RN, RC, RP, RX, RG, RA, RT, RL
            StringTools.writeKeyValueLine(REFERENCE_TAG, "["+rdr.getRank()+"]", 5, this.getLineWidth(), null, REFERENCE_TAG, this.getPrintStream());
            StringTools.writeKeyValueLine(REMARK_TAG, d.getRemark(), 5, this.getLineWidth(), null, REMARK_TAG, this.getPrintStream());
            Integer rstart = rdr.getStart();
            if (rstart==null) rstart = new Integer(1);
            Integer rend = rdr.getEnd();
            if (rend==null) rend = new Integer(rs.length());
            StringTools.writeKeyValueLine(REFERENCE_POSITION_TAG, rstart+"-"+rend, 5, this.getLineWidth(), null, REFERENCE_POSITION_TAG, this.getPrintStream());
            CrossRef c = d.getCrossref();
            if (c!=null) StringTools.writeKeyValueLine(REFERENCE_XREF_TAG, c.getDbname()+"; "+c.getAccession()+".", 5, this.getLineWidth(), null, REFERENCE_XREF_TAG, this.getPrintStream());
            List auths = d.getAuthorList();
            for (Iterator j = auths.iterator(); j.hasNext(); ) {
                DocRefAuthor a = (DocRefAuthor)j.next();
                if (a.isConsortium()) {
                    StringTools.writeKeyValueLine(CONSORTIUM_TAG, a+";", 5, this.getLineWidth(), null, CONSORTIUM_TAG, this.getPrintStream());
                    j.remove();
                }
            }
            if (!auths.isEmpty()) StringTools.writeKeyValueLine(AUTHORS_TAG, DocRefAuthor.Tools.generateAuthorString(auths, true)+";", 5, this.getLineWidth(), null, AUTHORS_TAG, this.getPrintStream());
            else StringTools.writeKeyValueLine(AUTHORS_TAG, ";", 5, this.getLineWidth(), null, AUTHORS_TAG, this.getPrintStream());
            if (d.getTitle()!=null && d.getTitle().length()!=0) StringTools.writeKeyValueLine(TITLE_TAG, "\""+d.getTitle()+"\";", 5, this.getLineWidth(), null, TITLE_TAG, this.getPrintStream());
            else StringTools.writeKeyValueLine(TITLE_TAG, ";", 5, this.getLineWidth(), null, TITLE_TAG, this.getPrintStream());
            StringTools.writeKeyValueLine(LOCATOR_TAG, d.getLocation()+".", 5, this.getLineWidth(), null, LOCATOR_TAG, this.getPrintStream());
            this.getPrintStream().println(DELIMITER_TAG+"   ");
        }
        
        // db references - ranked
        for (Iterator r = rs.getRankedCrossRefs().iterator(); r.hasNext(); ) {
            RankedCrossRef rcr = (RankedCrossRef)r.next();
            CrossRef c = rcr.getCrossRef();
            Set noteset = c.getNoteSet();
            StringBuffer sb = new StringBuffer();
            sb.append(c.getDbname());
            sb.append("; ");
            sb.append(c.getAccession());
            boolean hasSecondary = false;
            for (Iterator i = noteset.iterator(); i.hasNext(); ) {
                Note n = (Note)i.next();
                if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) {
                    sb.append("; ");
                    sb.append(n.getValue());
                    hasSecondary = true;
                }
            }
            //if (!hasSecondary) sb.append("; -"); 
            //sb.append(".");
            if (!hasSecondary) sb.append(";");
            else sb.append(".");
            StringTools.writeKeyValueLine(DATABASE_XREF_TAG, sb.toString(), 5, this.getLineWidth(), null, DATABASE_XREF_TAG, this.getPrintStream());
        }
        if (!rs.getRankedCrossRefs().isEmpty())
            this.getPrintStream().println(DELIMITER_TAG+"   ");
        
        // comments - if any
        if (!rs.getComments().isEmpty()) {
            StringBuffer sb = new StringBuffer();
            for (Iterator i = rs.getComments().iterator(); i.hasNext(); ) {
                Comment c = (SimpleComment)i.next();
                sb.append(c.getComment());
                if (i.hasNext()) sb.append("\n");
            }
            StringTools.writeKeyValueLine(COMMENT_TAG, sb.toString(), 5, this.getLineWidth(), null, COMMENT_TAG, this.getPrintStream());
            this.getPrintStream().println(DELIMITER_TAG+"   ");
        }
        
        this.getPrintStream().println(FEATURE_HEADER_TAG+"   Key             Location/Qualifiers");
        this.getPrintStream().println(FEATURE_HEADER_TAG+"   ");
        // feature_type     location
        for (Iterator i = rs.getFeatureSet().iterator(); i.hasNext(); ) {
            RichFeature f = (RichFeature)i.next();
            StringTools.writeKeyValueLine(FEATURE_TAG+"   "+f.getTypeTerm().getName(), GenbankLocationParser.writeLocation((RichLocation)f.getLocation()), 21, this.getLineWidth(), ",", FEATURE_TAG, this.getPrintStream());
            for (Iterator j = f.getNoteSet().iterator(); j.hasNext(); ) {
                Note n = (Note)j.next();
                // /key="val" or just /key if val==""
                if (n.getValue()==null || n.getValue().length()==0) StringTools.writeKeyValueLine(FEATURE_TAG, "/"+n.getTerm().getName(), 21, this.getLineWidth(), null, FEATURE_TAG, this.getPrintStream());
                else StringTools.writeKeyValueLine(FEATURE_TAG, "/"+n.getTerm().getName()+"=\""+n.getValue()+"\"", 21, this.getLineWidth(), null, FEATURE_TAG, this.getPrintStream());
            }
            // add-in to source feature only organism and db_xref="taxon:xyz" where present
            if (f.getType().equals("source") && tax!=null) {
                String displayName = tax.getDisplayName();
                if (displayName.indexOf('(')>-1) displayName = displayName.substring(0, displayName.indexOf('(')).trim();
                StringTools.writeKeyValueLine(FEATURE_TAG, "/organism=\""+displayName+"\"", 21, this.getLineWidth(), null, FEATURE_TAG, this.getPrintStream());
                StringTools.writeKeyValueLine(FEATURE_TAG, "/db_xref=\"taxon:"+tax.getNCBITaxID()+"\"", 21, this.getLineWidth(), null, FEATURE_TAG, this.getPrintStream());
            }
            // add-in other dbxrefs where present
            for (Iterator j = f.getRankedCrossRefs().iterator(); j.hasNext(); ) {
                RankedCrossRef rcr = (RankedCrossRef)j.next();
                CrossRef cr = rcr.getCrossRef();
                StringTools.writeKeyValueLine(FEATURE_TAG, "/db_xref=\""+cr.getDbname()+":"+cr.getAccession()+"\"", 21, this.getLineWidth(), null, FEATURE_TAG, this.getPrintStream());
            }
        }
        this.getPrintStream().println(DELIMITER_TAG+"   ");
        
        // SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
        int aCount = 0;
        int cCount = 0;
        int gCount = 0;
        int tCount = 0;
        int oCount = 0;
        for (int i = 1; i <= rs.length(); i++) {
            char c;
            try {
                c = tok.tokenizeSymbol(rs.symbolAt(i)).charAt(0);
            } catch (Exception e) {
                throw new RuntimeException("Unable to get symbol at position "+i,e);
            }
            switch (c) {
                case 'a': case 'A':
                    aCount++;
                    break;
                case 'c': case 'C':
                    cCount++;
                    break;
                case 'g': case 'G':
                    gCount++;
                    break;
                case 't': case 'T':
                    tCount++;
                    break;
                default:
                    oCount++;
            }
        }
        this.getPrintStream().print(START_SEQUENCE_TAG+"   Sequence "+rs.length()+" BP; ");
        this.getPrintStream().print(aCount + " A; ");
        this.getPrintStream().print(cCount + " C; ");
        this.getPrintStream().print(gCount + " G; ");
        this.getPrintStream().print(tCount + " T; ");
        this.getPrintStream().println(oCount + " other;");
        
        // sequence stuff
        Symbol[] syms = (Symbol[])rs.toList().toArray(new Symbol[0]);
        int lineLen = 0;
        int symCount = 0;
        this.getPrintStream().print("    ");
        for (int i = 0; i < syms.length; i++) {
            if (symCount % 60 == 0 && symCount>0) {
                this.getPrintStream().print(StringTools.leftPad(""+symCount,10));
                this.getPrintStream().print("\n    ");
                lineLen = 0;
            }
            if (symCount % 10 == 0) {
                this.getPrintStream().print(" ");
                lineLen++;
            }
            try {
                this.getPrintStream().print(tok.tokenizeSymbol(syms[i]));
            } catch (IllegalSymbolException e) {
                throw new RuntimeException("Found illegal symbol: "+syms[i]);
            }
            symCount++;
            lineLen++;
        }
        this.getPrintStream().print(StringTools.leftPad(""+symCount,(66-lineLen)+10));
        this.getPrintStream().print("\n");
        this.getPrintStream().println(END_SEQUENCE_TAG);
    }
    
    /**
     * {@inheritDoc}
     */
    public String getDefaultFormat() {
        return EMBL_FORMAT;
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

