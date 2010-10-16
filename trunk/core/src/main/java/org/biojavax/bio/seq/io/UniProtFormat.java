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
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.proteomics.MassCalc;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.ontology.Term;
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
import org.biojavax.utils.CRC64Checksum;
import org.biojavax.utils.StringTools;

/**
 * Format reader for UniProt files. This version of UniProt format will generate
 * and write RichSequence objects. Loosely Based on code from the old, deprecated,
 * org.biojava.bio.seq.io.EMBLLikeFormat object. Since 1.7, the parser reads the
 * International Protein Index (IPI) pseudo-Uniprot format.
 *
 * @author Richard Holland
 * @author Mark Schreiber
 * @author George Waldon
 * @since 1.5
 */
public class UniProtFormat extends RichSequenceFormat.HeaderlessFormat {
    
    // Register this format with the format auto-guesser.
    static {
        RichSequence.IOTools.registerFormat(UniProtFormat.class);
    }
    
    /**
     * The name of this format
     */
    public static final String UNIPROT_FORMAT = "UniProt";

    private static final String SUBFORMAT_UNIPROT = "UniProt";
    private static final String SUBFORMAT_IPI = "IPI";
    
    protected static final String LOCUS_TAG = "ID";
    protected static final String ACCESSION_TAG = "AC";
    protected static final String DEFINITION_TAG = "DE";
    protected static final String DATE_TAG = "DT";
    protected static final String SOURCE_TAG = "OS";
    protected static final String ORGANELLE_TAG = "OG";
    protected static final String ORGANISM_TAG = "OC";
    protected static final String TAXON_TAG = "OX";
    protected static final String GENE_TAG = "GN";
    protected static final String DATABASE_XREF_TAG = "DR";
    protected static final String PROTEIN_EXIST_TAG = "PE";
    protected static final String REFERENCE_TAG = "RN";
    protected static final String RP_LINE_TAG = "RP";
    protected static final String REFERENCE_XREF_TAG = "RX";
    protected static final String AUTHORS_TAG = "RA";
    protected static final String CONSORTIUM_TAG = "RG";
    protected static final String TITLE_TAG = "RT";
    protected static final String LOCATION_TAG = "RL";
    protected static final String RC_LINE_TAG = "RC";
    protected static final String KEYWORDS_TAG = "KW";
    protected static final String COMMENT_TAG = "CC";
    protected static final String FEATURE_TAG = "FT";
    protected static final String START_SEQUENCE_TAG = "SQ";
    protected static final String END_SEQUENCE_TAG = "//";
    
    // locus line for uniprot format
    protected static final Pattern lp_uniprot = Pattern.compile("^((\\S+)_(\\S+))\\s+(\\S+);\\s+(PRT)?;?\\s*\\d+\\s+AA\\.$");
    // locus line for IPI format
    protected static final Pattern lp_ipi = Pattern.compile("^((\\S+)\\.(\\d+))\\s+(IPI);\\s+(PRT)?;?\\s*\\d+\\s+AA\\.$");
    // RP line parser
    protected static final Pattern rppat = Pattern.compile("SEQUENCE OF (\\d+)-(\\d+)");
    // date lineDT for uniprot
    // date, integrated into UniProtKB/database_name.
    // date, sequence version x.
    // date, entry version x.
    protected static final Pattern dp_uniprot = Pattern.compile("([^,]+),([^\\d\\.]+)(\\d+)?\\.$");
    // date lineDT for IPI
    // date (xxx, Created)
    // date (xxx, Last sequence update)
    protected static final Pattern dp_ipi = Pattern.compile("([^\\(]+)\\(([^,]+),([^\\)]+)\\)$");
    // feature line
    protected static final Pattern fp = Pattern.compile("^\\s*([\\d?<]+\\s+[\\d?>]+)(\\s+(.*))?$");
    
    protected static final Pattern headerLine = Pattern.compile("^ID.*");
    
    /**
     * Implements some UniProt-specific terms.
     */
    public static class Terms extends RichSequence.Terms {        
        private static String GENENAME_KEY = "Name";
        private static String GENESYNONYM_KEY = "Synonyms";
        private static String ORDLOCNAME_KEY = "OrderedLocusNames";
        private static String ORFNAME_KEY = "ORFNames";
        
        /**
         * Getter for the UniProt term
         * @return The UniProt Term
         */
        public static ComparableTerm getUniProtTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("UniProt");
        }
        
        /**
         * Getter for the UniProt combined database term
         * @return The combined database for UniProt Term
         */
        public static ComparableTerm getUniProtDBNameTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("UniProt database name");
        }
        
        /**
         * Getter for the protein exists term
         * @return The protein exists Term
         */
        public static ComparableTerm getProteinExistsTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("UniProt protein exists");
        }
    }
    
    /**
     * {@inheritDoc}
     * A file is in UniProt format if the first line matches the UniProt format for the ID line.
     */
    public boolean canRead(File file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String firstLine = br.readLine();
        boolean readable = firstLine!=null && headerLine.matcher(firstLine).matches() && 
                (lp_uniprot.matcher(firstLine.substring(3).trim()).matches() ||
                lp_ipi.matcher(firstLine.substring(3).trim()).matches());
        br.close();
        return readable;
    }
    
    /**
     * {@inheritDoc}
     * Always returns a protein tokenizer.
     */
    public SymbolTokenization guessSymbolTokenization(File file) throws IOException {
        return RichSequence.IOTools.getProteinParser();
    }
    
    /**
     * {@inheritDoc}
     * A stream is in UniProt format if the first line matches the UniProt format for the ID line.
     */
    public boolean canRead(BufferedInputStream stream) throws IOException {
        stream.mark(2000); // some streams may not support this
        BufferedReader br = new BufferedReader(new InputStreamReader(stream));
        String firstLine = br.readLine();
        boolean readable = firstLine!=null && headerLine.matcher(firstLine).matches() && 
                (lp_uniprot.matcher(firstLine.substring(3).trim()).matches()
                || lp_ipi.matcher(firstLine.substring(3).trim()).matches());
        // don't close the reader as it'll close the stream too.
        // br.close();
        stream.reset();
        return readable;
    }
    
    /**
     * {@inheritDoc}
     * Always returns a protein tokenizer.
     */
    public SymbolTokenization guessSymbolTokenization(BufferedInputStream stream) throws IOException {
        return RichSequence.IOTools.getProteinParser();
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
    
    private String accession = null;
    
    /**
     * {@inheritDoc}
     */
    public boolean readRichSequence(BufferedReader reader,
            SymbolTokenization symParser,
            RichSeqIOListener rlistener,
            Namespace ns)
            throws IllegalSymbolException, IOException, ParseException {
        
        boolean hasAnotherSequence = true;
        //boolean hasInternalWhitespace = false;

        String subformat = SUBFORMAT_UNIPROT;
        
        rlistener.startSequence();
        
        if (ns==null) ns=RichObjectFactory.getDefaultNamespace();
        rlistener.setNamespace(ns);
        
        // Get an ordered list of key->value pairs in array-tuples
        String sectionKey = null;
        NCBITaxon tax = null;
        accession = null;
        List section = null;
        try{
            do {
                
                section = this.readSection(reader);
                sectionKey = ((String[])section.get(0))[0];
                if(sectionKey == null){
                    String message = ParseException.newMessage(this.getClass(),accession, "", "Section key was null", sectionToString(section));
                    throw new ParseException(message);
                }
                // process section-by-section
                if (sectionKey.equals(LOCUS_TAG)) {
                    // entryname  dataclass; moltype; sequencelength AA.
                    String loc = ((String[])section.get(0))[1];
                    Matcher m = lp_uniprot.matcher(loc);
                    if (m.matches()) {
                        rlistener.setName(m.group(2));
                        rlistener.setDivision(m.group(3));
                        if (m.groupCount() > 4){
                            rlistener.addSequenceProperty(Terms.getDataClassTerm(),m.group(4));
                            rlistener.addSequenceProperty(Terms.getMolTypeTerm(),m.group(5));
                        }else{
                            rlistener.addSequenceProperty(Terms.getDataClassTerm(), m.group(4));
                            rlistener.addSequenceProperty(Terms.getMolTypeTerm(), "");
                        }
                    } else {
                        m = lp_ipi.matcher(loc);
                        if (m.matches()) {
                            subformat = SUBFORMAT_IPI;
                            rlistener.setName(m.group(2));
                            rlistener.setVersion(Integer.parseInt(m.group(3)));
                            rlistener.addSequenceProperty(Terms.getDataClassTerm(), m.group(4));
                            rlistener.addSequenceProperty(Terms.getMolTypeTerm(),m.group(5));
                        } else {
                            String message = ParseException.newMessage(this.getClass(),accession, "", "Bad ID line", sectionToString(section));
                            throw new ParseException(message);
                        }
                    }
                } else if (sectionKey.equals(DEFINITION_TAG)) {
                    String val = ((String[])section.get(0))[1];
                    if (val.endsWith(".")) val = val.substring(0, val.length()-1); // chomp dot
                    rlistener.setDescription(val);
                } else if (sectionKey.equals(SOURCE_TAG)) {
                    // use SOURCE_TAG and TAXON_TAG values
                    String sciname = null;
                    String comname = null;
                    List synonym = new ArrayList();
                    int taxid = 0;
                    for (int i = 0; i < section.size(); i++) {
                        String tag = ((String[])section.get(i))[0];
                        String value = ((String[])section.get(i))[1].trim();
                        if (tag.equals(SOURCE_TAG)) {
                            if (value.endsWith(".")) value = value.substring(0,value.length()-1); // chomp trailing dot
                            String[] parts = value.split("\\(");
                            sciname = parts[0].trim();
                            if (parts.length>1) {
                                comname = parts[1].trim();
                                if (comname.endsWith(")")) comname = comname.substring(0,comname.length()-1); // chomp trailing bracket
                                if (parts.length>2) {
                                    // synonyms
                                    for (int j = 2 ; j < parts.length; j++) {
                                        String syn = parts[j].trim();
                                        if (syn.endsWith(")")) syn = syn.substring(0,syn.length()-1); // chomp trailing bracket
                                        synonym.add(syn);
                                    }
                                }
                            }
                        } else if (tag.equals(TAXON_TAG)) {
                            String[] parts = value.split(";");
                            for (int j = 0; j < parts.length; j++) {
                                String[] bits = parts[j].split("=");
                                if (bits[0].equals("NCBI_TaxID")) {
                                    String[] morebits = bits[1].split(",");
                                    taxid = Integer.parseInt(morebits[0].trim());
                                }
                            }
                        } else if (tag.equals(ORGANELLE_TAG)) {
                            if (value.endsWith(".")) value = value.substring(0,value.length()-1); // chomp trailing dot
                            String[] parts = value.split(";");
                            for (int j = 0; j < parts.length; j++) {
                                parts[j]=parts[j].trim();
                                rlistener.addSequenceProperty(Terms.getOrganelleTerm(),parts[j]);
                            }
                        }
                    }
                    // Set the Taxon
                    tax = (NCBITaxon)RichObjectFactory.getObject(SimpleNCBITaxon.class, new Object[]{new Integer(taxid)});
                    rlistener.setTaxon(tax);
                    try {
                        if (sciname!=null) tax.addName(NCBITaxon.SCIENTIFIC,sciname);
                        if (comname!=null) tax.addName(NCBITaxon.COMMON,comname);
                        for (Iterator j = synonym.iterator(); j.hasNext(); ) tax.addName(NCBITaxon.SYNONYM, (String)j.next());
                    } catch (ChangeVetoException e) {
                        throw new ParseException(e);
                    }
                } else if (sectionKey.equals(DATE_TAG)) {
                    String chunk = ((String[])section.get(0))[1];
                    if(subformat.equals(SUBFORMAT_UNIPROT)) {
                        Matcher dm = dp_uniprot.matcher(chunk);
                        if (dm.matches()) {
                            String date = dm.group(1).trim();
                            String type = dm.group(2).trim();
                            String rel = dm.group(3);
                            if (rel!=null) rel = rel.trim();
                            if (type.startsWith("integrated into UniProtKB")) {
                                String dbname = type.split("/")[1];
                                rlistener.addSequenceProperty(Terms.getDateCreatedTerm(), date);
                                rlistener.addSequenceProperty(Terms.getUniProtDBNameTerm(), dbname);
                            } else if (type.equalsIgnoreCase("sequence version")) {
                                if (rel==null){
                                    String message = ParseException.newMessage(this.getClass(),accession, "", "Version missing for "+type, sectionToString(section));
                                    throw new ParseException(message);
                                }
                                rlistener.addSequenceProperty(Terms.getDateUpdatedTerm(), date);
                                rlistener.setVersion(Integer.parseInt(rel));
                            } else if (type.equalsIgnoreCase("entry version")) {
                                if (rel==null) {
                                    String message = ParseException.newMessage(this.getClass(),accession, "", "Version missing for "+type, sectionToString(section));
                                    throw new ParseException(message);
                                }
                                rlistener.addSequenceProperty(Terms.getDateAnnotatedTerm(), date);
                                rlistener.addSequenceProperty(Terms.getRelAnnotatedTerm(), rel);
                            } else {
                                String message = ParseException.newMessage(this.getClass(),accession, "", "Bad date type "+type, sectionToString(section));
                                throw new ParseException(message);
                            }
                        } else {
                            String message = ParseException.newMessage(this.getClass(),accession, "", "Bad date line", sectionToString(section));
                            throw new ParseException(message);
                        }
                    } else if(subformat.equals(SUBFORMAT_IPI)) {
                        Matcher dm = dp_ipi.matcher(chunk);
                        if (dm.matches()) {
                            String date = dm.group(1).trim();
                            String type = dm.group(3).trim();
                            if(type.equals("Created")) {
                                rlistener.addSequenceProperty(Terms.getDateCreatedTerm(), date);
                            } else if(type.equals("Last sequence update")) {
                                rlistener.addSequenceProperty(Terms.getDateUpdatedTerm(), date);
                            } else {
                                String message = ParseException.newMessage(this.getClass(),accession, "", "Bad date type "+type, sectionToString(section));
                                throw new ParseException(message);
                            }
                        } else {
                            String message = ParseException.newMessage(this.getClass(),accession, "", "Bad date line", sectionToString(section));
                            throw new ParseException(message);
                        }
                    } else {
                            String message = ParseException.newMessage(this.getClass(),accession, "", "Unknown date line format", sectionToString(section));
                            throw new ParseException(message);
                    }
                } else if (sectionKey.equals(ACCESSION_TAG)) {
                    // if multiple accessions, store only first as accession,
                    // and store rest in annotation
                    String[] accs = ((String[])section.get(0))[1].split(";");
                    if(accs.length>0) accession = accs[0].trim(); else accession = "";
                    rlistener.setAccession(accession);
                    for (int i = 1; i < accs.length; i++) {
                        rlistener.addSequenceProperty(Terms.getAdditionalAccessionTerm(),accs[i].trim());
                    }
                } else if (sectionKey.equals(PROTEIN_EXIST_TAG)) {
                    String val = ((String[])section.get(0))[1];
                    if (val.endsWith(";")) val = val.substring(0, val.length()-1); // chomp semicolon
                    rlistener.addSequenceProperty(Terms.getProteinExistsTerm(),val.trim());
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
                } else if (sectionKey.equals(GENE_TAG)) {
                    String[] genes = ((String[])section.get(0))[1].split("\\s+(or|and)\\s+");
                    for (int geneID = 0; geneID < genes.length; geneID++) {
                        String[] parts = genes[geneID].split(";");
                        for (int j = 0; j < parts.length; j++) {
                            String[] moreparts = parts[j].split("=");
                            String[] values = moreparts[1].split(",");
                            // nasty hack - we really should have notes on the gene object itself... if such a thing existed...
                            if (moreparts[0].trim().equals(Terms.GENENAME_KEY)) rlistener.addSequenceProperty(Terms.getGeneNameTerm(),geneID+":"+values[0].trim());
                            else if (moreparts[0].trim().equals(Terms.GENESYNONYM_KEY)) {
                                for (int k = 0; k < values.length; k++) rlistener.addSequenceProperty(Terms.getGeneSynonymTerm(),geneID+":"+values[k].trim());
                            } else if (moreparts[0].trim().equals(Terms.ORDLOCNAME_KEY)) {
                                for (int k = 0; k < values.length; k++) rlistener.addSequenceProperty(Terms.getOrderedLocusNameTerm(),geneID+":"+values[k].trim());
                            } else if (moreparts[0].trim().equals(Terms.ORFNAME_KEY)) {
                                for (int k = 0; k < values.length; k++) rlistener.addSequenceProperty(Terms.getORFNameTerm(),geneID+":"+values[k].trim());
                            }
                        }
                    }
                } else if (sectionKey.equals(DATABASE_XREF_TAG)) {
                    // database_identifier; primary_identifier; secondary_identifier....
                    String val = ((String[])section.get(0))[1];
                    if (val.endsWith(".")) val = val.substring(0, val.length()-1); // chomp dot
                    String[] parts = val.split(";");
                    // construct a DBXREF out of the dbname part[0] and accession part[1]
                    String dbname = parts[0].trim();
                    String acc = parts[1].trim();
                    CrossRef crossRef = (CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{dbname,acc,new Integer(0)});
                    // assign remaining bits of info as additional accession annotations
                    for (int j = 2; j < parts.length; j++) {
                        ComparableTerm t = (ComparableTerm)Terms.getAdditionalAccessionTerm();
                        Note note = new SimpleNote(t,parts[j].trim(),j-1);
                        try {
                            crossRef.getRichAnnotation().addNote(note);
                        } catch (ChangeVetoException ce) {
                            ParseException pe = new ParseException("Could not annotate additional accession terms");
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
                    // rest can be in any order
                    String authors = null;
                    String consortium = null;
                    String title = null;
                    String locator = null;
                    String pubmed = null;
                    String medline = null;
                    String doi = null;
                    String remark = null;
                    Integer rstart = null;
                    Integer rend = null;
                    for (int i = 1; i < section.size(); i++) {
                        String key = ((String[])section.get(i))[0];
                        String val = ((String[])section.get(i))[1];
                        //System.err.println(key+": "+val);
                        if (key.equals(AUTHORS_TAG)) {
                            if (val.endsWith(";")) val = val.substring(0, val.length()-1); // chomp semicolon
                            authors = val.replace('\n',' '); //see #2276
                        }
                        if (key.equals(CONSORTIUM_TAG)) {
                            if (val.endsWith(";")) val = val.substring(0, val.length()-1); // chomp semicolon
                            consortium = val.replace('\n',' '); //see #2276
                        }
                        if (key.equals(TITLE_TAG)) {
                            if (val.endsWith(";")) val = val.substring(0, val.length()-1); // chomp semicolon
                            if (val.endsWith("\"")) val = val.substring(1, val.length()-1); // chomp quotes
                            title = val.replace('\n',' '); //see #2276
                        }
                        if (key.equals(LOCATION_TAG)) {
                            if (val.endsWith(".")) val = val.substring(0, val.length()-1); // chomp dot
                            locator = val.replace('\n',' '); //see #2276
                        }
                        if (key.equals(REFERENCE_XREF_TAG)) {
                            // database_identifier=primary_identifier;
                            String[] refs = val.split(";");
                            for (int j = 0 ; j < refs.length; j++) {
                                if (refs[j].trim().length()==0) continue;
                                String[] parts = refs[j].split("=");
                                if ( parts.length <2) {
                                	// some DOI lines look like this and are causing problems:
                                	//DOI=10.1002/(SICI)1097-0215(19990702)82:1<137::AID-IJC23>3.0.CO;2-F;ignoring
                                	System.err.println("warning: problems while parsing: " + val);
                                	continue;
                                }
                                String db = parts[0].trim();
                                String ref = parts[1].trim();
                                if (db.equalsIgnoreCase(Terms.PUBMED_KEY)) pubmed = ref;
                                else if (db.equalsIgnoreCase(Terms.MEDLINE_KEY)) medline = ref;
                                else if (db.equalsIgnoreCase(Terms.DOI_KEY)) doi = ref;
                            }
                        }
                        if (key.equals(RP_LINE_TAG)) {
                            if (val.endsWith(".")) val = val.substring(0, val.length()-1); // chomp dot
                            remark = val.replace('\n',' '); //see #2276
                            // Try to use it to find the location of the reference, if we have one.
                            Matcher m = rppat.matcher(val);
                            if (m.matches()) {
                                rstart = Integer.valueOf(m.group(1));
                                rend = Integer.valueOf(m.group(2));
                            }
                        }
                        if (key.equals(RC_LINE_TAG)) {
                            // Split into key=value pairs separated by semicolons and terminated with semicolon.
                            String[] parts = val.split(";");
                            for (int j = 0; j < parts.length; j++) {
                                String[] subparts = parts[j].split("=");
                                // get term for first section
                                String termName = subparts[0].trim();
                                Term t;
                                if (termName.equalsIgnoreCase(Terms.SPECIES_KEY)) t = Terms.getSpeciesTerm();
                                else if (termName.equalsIgnoreCase(Terms.STRAIN_KEY)) t = Terms.getStrainTerm();
                                else if (termName.equalsIgnoreCase(Terms.TISSUE_KEY)) t = Terms.getTissueTerm();
                                else if (termName.equalsIgnoreCase(Terms.TRANSPOSON_KEY)) t = Terms.getTransposonTerm();
                                else if (termName.equalsIgnoreCase(Terms.PLASMID_KEY)) t = Terms.getPlasmidTerm();
                                else {
                                    String message = ParseException.newMessage(this.getClass(),accession, "", "Invalid RC term found: "+termName, sectionToString(section));
                                    throw new ParseException(message);
                                }
                                // assign notes using term and rank:second section as value
                                // nasty hack - we really should have notes on the reference itself.
                                rlistener.addSequenceProperty(t, ref_rank+":"+subparts[1].trim());
                            }
                        }
                    }
                    
                    // create the docref object
                    try {
                        List auths = null;
                        if(authors != null) auths = DocRefAuthor.Tools.parseAuthorString(authors);
                        if (consortium!=null){
                            if(auths == null) auths = new ArrayList();
                            auths.add(new SimpleDocRefAuthor(consortium,true,false));
                        }
                        DocRef dr = (DocRef)RichObjectFactory.getObject(SimpleDocRef.class,new Object[]{auths,locator,title});
                        // assign either the pubmed or medline to the docref - medline gets priority, then pubmed, then doi
                        if (medline!=null) dr.setCrossref((CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{Terms.MEDLINE_KEY, medline, new Integer(0)}));
                        else if (pubmed!=null) dr.setCrossref((CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{Terms.PUBMED_KEY, pubmed, new Integer(0)}));
                        else if (doi!=null) dr.setCrossref((CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{Terms.DOI_KEY, doi, new Integer(0)}));
                        // assign the remarks
                        if (!this.getElideComments()) dr.setRemark(remark);
                        // assign the docref to the bioentry
                        RankedDocRef rdr = new SimpleRankedDocRef(dr,rstart,rend,ref_rank);
                        rlistener.setRankedDocRef(rdr);
                    } catch (ChangeVetoException e) {
                        throw new ParseException(e);
                    }
                } else if (sectionKey.equals(COMMENT_TAG) && !this.getElideComments()) {
                    // Set up some comments
                    String val = ((String[])section.get(0))[1];
                    if (UniProtCommentParser.isParseable(val)) rlistener.setComment(val);
                    else {
                        // copyright message
                        rlistener.addSequenceProperty(Terms.getCopyrightTerm(), val);
                    }
                } else if (sectionKey.equals(FEATURE_TAG) && !this.getElideFeatures()) {
                    // starting from second line of input, start a new feature whenever we come across
                    // a key that does not start with /
                    boolean seenAFeature = false;
                    for (int i = 1 ; i < section.size(); i++) {
                        String key = ((String[])section.get(i))[0];
                        String val = ((String[])section.get(i))[1];
                        val = val.replaceAll("\\s*[\\n\\r]+\\s*", " ").trim();
                        if (val.endsWith(".")) val = val.substring(0,val.length()-1); // chomp dot
                        if (key.startsWith("/")) {
                            key = key.substring(1); // strip leading slash
                            if (key.equals("FTId")) rlistener.addFeatureProperty(Terms.getFTIdTerm(),val);
                            else {
                                // should never happen - but here just in case
                                rlistener.addFeatureProperty(RichObjectFactory.getDefaultOntology().getOrCreateTerm(key),val);
                            }
                        } else {
                            // new feature!
                            // end previous feature
                            if (seenAFeature) rlistener.endFeature();
                            // start next one, with lots of lovely info in it
                            RichFeature.Template templ = new RichFeature.Template();
                            templ.annotation = new SimpleRichAnnotation();
                            templ.sourceTerm = Terms.getUniProtTerm();
                            templ.typeTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm(key);
                            templ.featureRelationshipSet = new TreeSet();
                            templ.rankedCrossRefs = new TreeSet();
                            String desc = null;
                            Matcher m = fp.matcher(val);
                            if (m.matches()) {
                                String loc = m.group(1);
                                desc = m.group(3);
                                templ.location = UniProtLocationParser.parseLocation(loc);
                            } else {
                                String message = ParseException.newMessage(this.getClass(),accession, "", "Bad feature value: "+val, sectionToString(section));
                                throw new ParseException(message);
                            }
                            rlistener.startFeature(templ);
                            if (desc!=null && desc.length()>0) rlistener.addFeatureProperty(Terms.getFeatureDescTerm(),desc);
                            seenAFeature = true;
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
                    } catch (IllegalAlphabetException e) {
                        String message = ParseException.newMessage(this.getClass(),accession, "", "", sectionToString(section));
                        throw new ParseException(e, message);
                    }
                }
            } while (!sectionKey.equals(END_SEQUENCE_TAG));
        }catch (RuntimeException e){
            String message = ParseException.newMessage(this.getClass(),accession, "", "", sectionToString(section));
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
            //System.err.println("Warning: whitespace found between sequence entries");
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
                    String message = ParseException.newMessage(this.getClass(),accession, "", "Bad line found: "+line, sectionToString(section));
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
                            sb.append(line);
                        }
                    }
                    section.add(new String[]{START_SEQUENCE_TAG,sb.toString()});
                }
                // READ COMMENT SECTION
                else if (token.equals(COMMENT_TAG)) {
                    // read from first line till next that begins with "CC   -!-"
                    StringBuffer currentVal = new StringBuffer();
                    boolean wasMisc = false;
                    if (!line.startsWith(COMMENT_TAG+"   -!-")) wasMisc = true;
                    currentVal.append(line.substring(5));
                    while (!done) {
                        br.mark(160);
                        line = br.readLine();
                        if (((!wasMisc) && line.charAt(5)!=' ') || !line.startsWith("C") || line.startsWith(COMMENT_TAG+"   -!-")) {
                            br.reset();
                            done = true;
                            // dump current tag if exists
                            section.add(new String[]{COMMENT_TAG,currentVal.toString()});
                        } else {
                            currentVal.append("\n");
                            currentVal.append(line.substring(5));
                        }
                    }
                }
                // READ FEATURE TABLE SECTION
                else if (token.equals(FEATURE_TAG)) {
                    br.reset();
                    //      read all FT lines until first non-FT starting line
                    String currentTag = null;
                    StringBuffer currentVal = new StringBuffer();
                    section.add(new String[]{FEATURE_TAG,null});
                    while (!done) {
                        br.mark(160);
                        line = br.readLine();
                        if (!line.startsWith(FEATURE_TAG)) {
                            br.reset();
                            done = true;
                            // dump current tag if exists
                            if (currentTag!=null) section.add(new String[]{currentTag,currentVal.toString()});
                        } else {
                            //         FT lines:   FT   KEY_NAME     x      x        description
                            //         or:         FT                                ....
                            //         or          FT                                /FTId=899.
                            line = line.substring(5); // chomp off "FT   "
                            if (!line.startsWith(" ")) {
                                // dump current tag if exists
                                if (currentTag!=null) section.add(new String[]{currentTag,currentVal.toString()});
                                // case 1 : word value - splits into key-value based on first 8 chars
                                currentTag = line.substring(0,8).trim();
                                currentVal = new StringBuffer();
                                currentVal.append(line.substring(8).trim());
                            } else {
                                line = line.trim();
                                if (line.startsWith("/")) {
                                    // dump current tag if exists
                                    if (currentTag!=null) section.add(new String[]{currentTag,currentVal.toString()});
                                    // case 3 : /word=.....
                                    currentVal = new StringBuffer();
                                    int equalIndex = line.indexOf('=');
                                    if (equalIndex>=0) {
                                        currentTag = line.substring(0, equalIndex);
                                        currentVal.append(line.substring(equalIndex+1));
                                    } else {
                                        currentTag = line;
                                    }
                                } else {
                                    // case 2 : ...."
                                    currentVal.append("\n");
                                    currentVal.append(line);
                                }
                            }
                        }
                    }
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
                // READ END OF SEQUENCE
                else if (token.equals(END_SEQUENCE_TAG)) {
                    section.add(new String[]{END_SEQUENCE_TAG,null});
                    done = true;
                }
                // READ NORMAL TAG/VALUE SECTION
                else {
                    //      rewind buffer to mark
                    br.reset();
                    //      read token/values until first with non-same first character
                    //      exceptions: DE/DT, and RN...RN
                    String currentTag = null;
                    char currentTagStart = '\0';
                    StringBuffer currentVal = null;
                    while (!done) {
                        br.mark(160);
                        line = br.readLine();
                        if (currentTagStart=='\0') currentTagStart = line.charAt(0);
                        if (!line.startsWith(""+currentTagStart) ||
                                (currentTagStart=='D' && currentTag!=null && !line.startsWith(""+currentTag)) ||
                                (currentTagStart=='R' && currentTag!=null && line.startsWith("RN"))) {
                            br.reset();
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
                                String message = ParseException.newMessage(this.getClass(),accession, "", "", sectionToString(section));
                                throw new ParseException(e, message);
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            String message = ParseException.newMessage(this.getClass(),accession, "", "", sectionToString(section));
            throw new ParseException(e, message);
        } catch (RuntimeException e){
            String message = ParseException.newMessage(this.getClass(),accession, "", "", sectionToString(section));
            throw new ParseException(e, message);
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
        if (!format.equals(this.getDefaultFormat())) throw new IllegalArgumentException("Unknown format: "+format);
        this.writeSequence(seq, RichObjectFactory.getDefaultNamespace());
    }
    
    /**
     * {@inheritDoc}
     * Namespace is ignored as UniProt has no concept of it.
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
        accessions.append(";");
        String cdat = null;
        String udat = null;
        String adat = null;
        String dbname = "?";
        String arel = null;
        String organelle = null;
        String protExists = null;
        String dataclass = "STANDARD";
        String copyright = null;
        Map speciesRecs = new TreeMap();
        Map strainRecs = new TreeMap();
        Map tissueRecs = new TreeMap();
        Map transpRecs = new TreeMap();
        Map plasmidRecs = new TreeMap();
        Map genenames = new TreeMap();
        Map genesynonyms = new TreeMap();
        Map orfnames = new TreeMap();
        Map ordlocnames = new TreeMap();
        for (Iterator i = notes.iterator(); i.hasNext(); ) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getDateCreatedTerm())) cdat=n.getValue();
            else if (n.getTerm().equals(Terms.getDateUpdatedTerm())) udat=n.getValue();
            else if (n.getTerm().equals(Terms.getDateAnnotatedTerm())) adat=n.getValue();
            else if (n.getTerm().equals(Terms.getUniProtDBNameTerm())) dbname=n.getValue();
            else if (n.getTerm().equals(Terms.getProteinExistsTerm())) protExists=n.getValue();
            else if (n.getTerm().equals(Terms.getRelAnnotatedTerm())) arel=n.getValue();
            else if (n.getTerm().equals(Terms.getDataClassTerm())) dataclass = n.getValue();
            else if (n.getTerm().equals(Terms.getCopyrightTerm())) copyright = n.getValue();
            else if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) {
                accessions.append(" ");
                accessions.append(n.getValue());
                accessions.append(";");
            } else if (n.getTerm().equals(Terms.getOrganelleTerm())) organelle = (organelle==null?"":organelle+"; ")+n.getValue();
            // use the nasty hack to split the reference rank away from the actual value in this field
            else if (n.getTerm().equals(Terms.getGeneNameTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                genenames.put(refID, ref.substring(colon+1)); // map of id -> string as only one name per gene
            } else if (n.getTerm().equals(Terms.getGeneSynonymTerm())) {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                if (genesynonyms.get(refID)==null) genesynonyms.put(refID, new ArrayList());
                ((List)genesynonyms.get(refID)).add(ref.substring(colon+1));
            } else if (n.getTerm().equals(Terms.getOrderedLocusNameTerm())) {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                if (ordlocnames.get(refID)==null) ordlocnames.put(refID, new ArrayList());
                ((List)ordlocnames.get(refID)).add(ref.substring(colon+1));
            } else if (n.getTerm().equals(Terms.getORFNameTerm())) {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                if (orfnames.get(refID)==null) orfnames.put(refID, new ArrayList());
                ((List)orfnames.get(refID)).add(ref.substring(colon+1));
            }
            // use the nasty hack to split the reference rank away from the actual value in this field
            // we'll end up with a bunch in key 0 for those which did not come from us. We ignore these for now.
            else if (n.getTerm().equals(Terms.getSpeciesTerm())) {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                if (speciesRecs.get(refID)==null) speciesRecs.put(refID, new ArrayList());
                ((List)speciesRecs.get(refID)).add(ref.substring(colon+1));
            } else if (n.getTerm().equals(Terms.getStrainTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                if (strainRecs.get(refID)==null) strainRecs.put(refID, new ArrayList());
                ((List)strainRecs.get(refID)).add(ref.substring(colon+1));
            } else if (n.getTerm().equals(Terms.getTissueTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                if (tissueRecs.get(refID)==null) tissueRecs.put(refID, new ArrayList());
                ((List)tissueRecs.get(refID)).add(ref.substring(colon+1));
            } else if (n.getTerm().equals(Terms.getTransposonTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                if (transpRecs.get(refID)==null) transpRecs.put(refID, new ArrayList());
                ((List)transpRecs.get(refID)).add(ref.substring(colon+1));
            } else if (n.getTerm().equals(Terms.getPlasmidTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                if (plasmidRecs.get(refID)==null) plasmidRecs.put(refID, new ArrayList());
                ((List)plasmidRecs.get(refID)).add(ref.substring(colon+1));
            }
        }
        
        // entryname  dataclass; [circular] molecule; division; sequencelength BP.
        StringBuffer locusLine = new StringBuffer();
        locusLine.append(StringTools.rightPad(rs.getName()+"_"+rs.getDivision(),12));
        locusLine.append(" ");
        locusLine.append(StringTools.leftPad(dataclass,19));
        //locusLine.append(";      PRT; "); //Uniprot no longer uses the PRT;
        locusLine.append("; ");
        locusLine.append(StringTools.leftPad(""+rs.length(),11));
        locusLine.append(" AA.");
        StringTools.writeKeyValueLine(LOCUS_TAG, locusLine.toString(), 5, this.getLineWidth(), null, LOCUS_TAG, this.getPrintStream());
        
        // accession line
        StringTools.writeKeyValueLine(ACCESSION_TAG, accessions.toString(), 5, this.getLineWidth(), null, ACCESSION_TAG, this.getPrintStream());
        
        // date line
        StringTools.writeKeyValueLine(DATE_TAG, (cdat==null?udat:cdat)+", integrated into UniProtKB/"+dbname+".", 5, this.getLineWidth(), null, DATE_TAG, this.getPrintStream());
        StringTools.writeKeyValueLine(DATE_TAG, udat+", sequence version "+rs.getVersion()+".", 5, this.getLineWidth(), null, DATE_TAG, this.getPrintStream());
        StringTools.writeKeyValueLine(DATE_TAG, (adat==null?udat:adat)+", entry version "+(arel==null?"0":arel)+".", 5, this.getLineWidth(), null, DATE_TAG, this.getPrintStream());
        
        // definition line
        StringTools.writeKeyValueLine(DEFINITION_TAG, rs.getDescription()+".", 5, this.getLineWidth(), null, DEFINITION_TAG, this.getPrintStream());
        
        // gene line
        for (Iterator i = genenames.keySet().iterator(); i.hasNext(); ) {
            Integer geneid = (Integer)i.next();
            String genename = (String)genenames.get(geneid);
            List synonyms = (List)genesynonyms.get(geneid);
            List orfs = (List)orfnames.get(geneid);
            List ordlocs = (List)ordlocnames.get(geneid);
            
            StringBuffer gnline = new StringBuffer();
            gnline.append(Terms.GENENAME_KEY);
            gnline.append("=");
            gnline.append(genename);
            gnline.append("; ");
            
            if (synonyms!=null) {
                gnline.append(Terms.GENESYNONYM_KEY);
                gnline.append("=");
                for (Iterator j = synonyms.iterator(); j.hasNext(); ) {
                    gnline.append((String)j.next());
                    if (j.hasNext()) gnline.append(", ");
                }
                gnline.append("; ");
            }
            if (ordlocs!=null) {
                gnline.append(Terms.ORDLOCNAME_KEY);
                gnline.append("=");
                for (Iterator j = ordlocs.iterator(); j.hasNext(); ) {
                    gnline.append((String)j.next());
                    if (j.hasNext()) gnline.append(", ");
                }
                gnline.append("; ");
            }
            if (orfs!=null) {
                gnline.append(Terms.ORFNAME_KEY);
                gnline.append("=");
                for (Iterator j = orfs.iterator(); j.hasNext(); ) {
                    gnline.append((String)j.next());
                    if (j.hasNext()) gnline.append(", ");
                }
                gnline.append("; ");
            }
            
            StringTools.writeKeyValueLine(GENE_TAG, gnline.toString(), 5, this.getLineWidth(), null, GENE_TAG, this.getPrintStream());
            
            if (i.hasNext()) StringTools.writeKeyValueLine(GENE_TAG, "and", 5, this.getLineWidth(), null, GENE_TAG, this.getPrintStream());
        }
        
        // source line (from taxon)
        //   organism line
        NCBITaxon tax = rs.getTaxon();
        if (tax!=null) {
            StringBuffer source = new StringBuffer();
            source.append(tax.getDisplayName());
            for (Iterator j = tax.getNames(NCBITaxon.SYNONYM).iterator(); j.hasNext(); ) {
                source.append(" (");
                source.append((String)j.next());
                source.append(")");
            }
            source.append(".");
            StringTools.writeKeyValueLine(SOURCE_TAG, source.toString(), 5, this.getLineWidth(), null, SOURCE_TAG, this.getPrintStream());
            if (organelle!=null) StringTools.writeKeyValueLine(ORGANELLE_TAG, organelle+".", 5, this.getLineWidth(), null, ORGANELLE_TAG, this.getPrintStream());
            StringTools.writeKeyValueLine(ORGANISM_TAG, tax.getNameHierarchy(), 5, this.getLineWidth(), null, ORGANISM_TAG, this.getPrintStream());
            StringTools.writeKeyValueLine(TAXON_TAG, "NCBI_TaxID="+tax.getNCBITaxID()+";", 5, this.getLineWidth(), this.getPrintStream());
        }
        
        // references - rank (bases x to y)
        for (Iterator r = rs.getRankedDocRefs().iterator(); r.hasNext(); ) {
            RankedDocRef rdr = (RankedDocRef)r.next();
            DocRef d = rdr.getDocumentReference();
            // RN, RP, RC, RX, RG, RA, RT, RL
            StringTools.writeKeyValueLine(REFERENCE_TAG, "["+rdr.getRank()+"]", 5, this.getLineWidth(), null, REFERENCE_TAG, this.getPrintStream());
            if (d.getRemark()!=null)
                StringTools.writeKeyValueLine(RP_LINE_TAG, d.getRemark()+".", 5, this.getLineWidth(), null, RP_LINE_TAG, this.getPrintStream());
            // Print out ref position if present
            if (rdr.getStart()!=null && rdr.getEnd()!=null && d.getRemark()!=null && !rppat.matcher(d.getRemark()).matches()) StringTools.writeKeyValueLine(RP_LINE_TAG, "SEQUENCE OF "+rdr.getStart()+"-"+rdr.getEnd()+".", 5, this.getLineWidth(), null, RP_LINE_TAG, this.getPrintStream());
            // RC lines
            StringBuffer rcline = new StringBuffer();
            Integer rank = new Integer(rdr.getRank());
            if (speciesRecs.get(rank)!=null) {
                rcline.append(Terms.SPECIES_KEY);
                rcline.append("=");
                for (Iterator i = ((List)speciesRecs.get(rank)).iterator(); i.hasNext(); ) {
                    rcline.append((String)i.next());
                    if (i.hasNext()) rcline.append(", ");
                }
                rcline.append("; ");
            }
            if (strainRecs.get(rank)!=null) {
                rcline.append(Terms.STRAIN_KEY);
                rcline.append("=");
                for (Iterator i = ((List)strainRecs.get(rank)).iterator(); i.hasNext(); ) {
                    rcline.append((String)i.next());
                    if (i.hasNext()) rcline.append(", ");
                }
                rcline.append("; ");
            }
            if (tissueRecs.get(rank)!=null) {
                rcline.append(Terms.TISSUE_KEY);
                rcline.append("=");
                for (Iterator i = ((List)tissueRecs.get(rank)).iterator(); i.hasNext(); ) {
                    rcline.append((String)i.next());
                    if (i.hasNext()) rcline.append(", ");
                }
                rcline.append("; ");
            }
            if (transpRecs.get(rank)!=null) {
                rcline.append(Terms.TRANSPOSON_KEY);
                rcline.append("=");
                for (Iterator i = ((List)transpRecs.get(rank)).iterator(); i.hasNext(); ) {
                    rcline.append((String)i.next());
                    if (i.hasNext()) rcline.append(", ");
                }
                rcline.append("; ");
            }
            if (plasmidRecs.get(rank)!=null) {
                rcline.append(Terms.PLASMID_KEY);
                rcline.append("=");
                for (Iterator i = ((List)plasmidRecs.get(rank)).iterator(); i.hasNext(); ) {
                    rcline.append((String)i.next());
                    if (i.hasNext()) rcline.append(", ");
                }
                rcline.append("; ");
            }
            // print the rcline
            if (rcline.length()>0) StringTools.writeKeyValueLine(RC_LINE_TAG, rcline.toString(), 5, this.getLineWidth(), null, RC_LINE_TAG, this.getPrintStream());
            // Deal with RX and rest
            CrossRef c = d.getCrossref();
            if (c!=null) StringTools.writeKeyValueLine(REFERENCE_XREF_TAG, c.getDbname()+"="+c.getAccession()+";", 5, this.getLineWidth(), null, REFERENCE_XREF_TAG, this.getPrintStream());
            List auths = d.getAuthorList();
            for (Iterator j = auths.iterator(); j.hasNext(); ) {
                DocRefAuthor a = (DocRefAuthor)j.next();
                if (a.isConsortium()) {
                    StringTools.writeKeyValueLine(CONSORTIUM_TAG, a.getName()+";", 5, this.getLineWidth(), null, CONSORTIUM_TAG, this.getPrintStream());
                    j.remove();
                }
            }
            if (!auths.isEmpty()) StringTools.writeKeyValueLine(AUTHORS_TAG, DocRefAuthor.Tools.generateAuthorString(auths, false)+";", 5, this.getLineWidth(), null, AUTHORS_TAG, this.getPrintStream());
            if (d.getTitle()!=null && d.getTitle().length()!=0) StringTools.writeKeyValueLine(TITLE_TAG, "\""+d.getTitle()+"\";", 5, this.getLineWidth(), null, TITLE_TAG, this.getPrintStream());
            StringTools.writeKeyValueLine(LOCATION_TAG, d.getLocation()+".", 5, this.getLineWidth(), null, LOCATION_TAG, this.getPrintStream());
        }
        
        // comments - if any
        if (!rs.getComments().isEmpty()) {
            for (Iterator i = rs.getComments().iterator(); i.hasNext(); ) {
                Comment c = (SimpleComment)i.next();
                String text = c.getComment().trim();
                if (text.length()>3 && text.substring(0,3).equals("-!-")) StringTools.writeKeyValueLine(COMMENT_TAG, text, 5, this.getLineWidth(), null, COMMENT_TAG, this.getPrintStream());
                else StringTools.writeKeyValueLine(COMMENT_TAG, text, 5, this.getLineWidth(), null, COMMENT_TAG, this.getPrintStream());
            }
        }
        
        // copyright - if any
        if (copyright!=null)
            StringTools.writeKeyValueLine(COMMENT_TAG, copyright, 5, this.getLineWidth(), null, COMMENT_TAG, this.getPrintStream());
        
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
            if (!hasSecondary) sb.append("; -");
            sb.append(".");
            StringTools.writeKeyValueLine(DATABASE_XREF_TAG, sb.toString(), 5, this.getLineWidth(), null, DATABASE_XREF_TAG, this.getPrintStream());
        }
        
        // protein exists line
        if (protExists!=null) {
            StringTools.writeKeyValueLine(PROTEIN_EXIST_TAG, protExists+";", 5, this.getLineWidth(), null, PROTEIN_EXIST_TAG, this.getPrintStream());
        }
        
        // keywords line
        String keywords = null;
        for (Iterator n = notes.iterator(); n.hasNext(); ) {
            Note nt = (Note)n.next();
            if (nt.getTerm().equals(Terms.getKeywordTerm())) {
                if (keywords==null) keywords = nt.getValue();
                else keywords = keywords+"; "+nt.getValue();
            }
        }
        if (keywords!=null) {
            StringTools.writeKeyValueLine(KEYWORDS_TAG, keywords+".", 5, this.getLineWidth(), null, KEYWORDS_TAG, this.getPrintStream());
        }
        
        // feature_type     location
        for (Iterator i = rs.getFeatureSet().iterator(); i.hasNext(); ) {
            RichFeature f = (RichFeature)i.next();
            String desc = "";
            String ftid = null;
            for (Iterator j = f.getNoteSet().iterator(); j.hasNext(); ) {
                Note n = (Note)j.next();
                if (n.getTerm().equals(Terms.getFTIdTerm())) ftid = n.getValue();
                else if (n.getTerm().equals(Terms.getFeatureDescTerm())) desc = n.getValue();
            }
            String kw = f.getTypeTerm().getName();
            String leader = StringTools.rightPad(kw,8)+" "+UniProtLocationParser.writeLocation((RichLocation)f.getLocation());
            if(desc.length()==0)
                this.getPrintStream().println(FEATURE_TAG+"   "+leader); //see #2277
            else
                StringTools.writeKeyValueLine(FEATURE_TAG+"   "+leader, desc+".", 34, this.getLineWidth(), null, FEATURE_TAG, this.getPrintStream());
            if (ftid!=null) StringTools.writeKeyValueLine(FEATURE_TAG, "/FTId="+ftid+".", 34, this.getLineWidth(), null, FEATURE_TAG, this.getPrintStream());
        }
        
        // sequence header
        int mw = 0;
        try {
            mw = (int)MassCalc.getMolecularWeight(rs);
        } catch (IllegalSymbolException e) {
            throw new RuntimeException("Found illegal symbol", e);
        }
        CRC64Checksum crc = new CRC64Checksum();
        String seqstr = rs.seqString();
        crc.update(seqstr.getBytes(),0,seqstr.length());
        this.getPrintStream().print(START_SEQUENCE_TAG+"   SEQUENCE  "+StringTools.leftPad(""+rs.length(),4)+" AA;  ");
        this.getPrintStream().print(StringTools.leftPad(""+mw,5)+" MW;  ");
        this.getPrintStream().println(crc+" CRC64;");
        
        // sequence stuff
        Symbol[] syms = (Symbol[])rs.toList().toArray(new Symbol[0]);
        int symCount = 0;
        this.getPrintStream().print("    ");
        for (int i = 0; i < syms.length; i++) {
            if (symCount % 60 == 0 && symCount>0) {
                this.getPrintStream().print("\n    ");
            }
            if (symCount % 10 == 0) {
                this.getPrintStream().print(" ");
            }
            try {
                this.getPrintStream().print(tok.tokenizeSymbol(syms[i]));
            } catch (IllegalSymbolException e) {
                throw new RuntimeException("Found illegal symbol: "+syms[i]);
            }
            symCount++;
        }
        this.getPrintStream().print("\n");
        this.getPrintStream().println(END_SEQUENCE_TAG);
    }
    
    /**
     * {@inheritDoc}
     */
    public String getDefaultFormat() {
        return UNIPROT_FORMAT;
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