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
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.parsers.ParserConfigurationException;

import org.biojava.bio.proteomics.MassCalc;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;
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
import org.biojavax.SimpleCrossRef;
import org.biojavax.SimpleDocRef;
import org.biojavax.SimpleDocRefAuthor;
import org.biojavax.SimpleNamespace;
import org.biojavax.SimpleNote;
import org.biojavax.SimpleRankedCrossRef;
import org.biojavax.SimpleRankedDocRef;
import org.biojavax.SimpleRichAnnotation;
import org.biojavax.bio.seq.Position;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.io.UniProtCommentParser.Event;
import org.biojavax.bio.seq.io.UniProtCommentParser.Interaction;
import org.biojavax.bio.seq.io.UniProtCommentParser.Isoform;
import org.biojavax.bio.taxa.NCBITaxon;
import org.biojavax.bio.taxa.SimpleNCBITaxon;
import org.biojavax.ontology.ComparableOntology;
import org.biojavax.ontology.ComparableTerm;
import org.biojavax.ontology.SimpleComparableOntology;
import org.biojavax.utils.CRC64Checksum;
import org.biojavax.utils.StringTools;
import org.biojavax.utils.XMLTools;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Format reader for UniProtXML files. This version of UniProtXML format will generate
 * and write RichSequence objects. Loosely Based on code from the old, deprecated,
 * org.biojava.bio.seq.io.GenbankXmlFormat object.
 *
 * Understands http://www.ebi.uniprot.org/support/docs/uniprot.xsd
 *
 * @author Alan Li (code based on his work)
 * @author Richard Holland
 * @since 1.5
 */
public class UniProtXMLFormat extends RichSequenceFormat.BasicFormat {
                        
    // Register this format with the format auto-guesser.
    static {
        RichSequence.IOTools.registerFormat(UniProtXMLFormat.class);
    }
    
    /**
     * The name of this format
     */
    public static final String UNIPROTXML_FORMAT = "UniProtXML";
    
    protected static final String ENTRY_GROUP_TAG = "uniprot";
    protected static final String ENTRY_TAG = "entry";
    protected static final String ENTRY_VERSION_ATTR = "version";
    protected static final String ENTRY_NAMESPACE_ATTR = "dataset";
    protected static final String ENTRY_CREATED_ATTR = "created";
    protected static final String ENTRY_UPDATED_ATTR = "modified";
    protected static final String COPYRIGHT_TAG = "copyright";
    
    protected static final String ACCESSION_TAG = "accession";
    protected static final String NAME_TAG = "name";
    protected static final String TEXT_TAG = "text";
    
    protected static final String REF_ATTR = "ref";
    protected static final String TYPE_ATTR = "type";
    protected static final String KEY_ATTR = "key";
    protected static final String ID_ATTR = "id";
    protected static final String EVIDENCE_ATTR = "evidence";
    protected static final String VALUE_ATTR = "value";
    protected static final String STATUS_ATTR = "value";
    protected static final String NAME_ATTR = "name";
    
    protected static final String PROTEIN_TAG = "protein";
    protected static final String PROTEIN_TYPE_ATTR = "type";
    
    protected static final String DOMAIN_TAG = "domain";
    protected static final String COMPONENT_TAG = "component";
    protected static final String GENE_TAG = "gene";
    protected static final String ORGANISM_TAG = "organism";
    protected static final String DBXREF_TAG = "dbReference";
    protected static final String PROPERTY_TAG = "property";
    protected static final String LINEAGE_TAG = "lineage";
    protected static final String TAXON_TAG = "taxon";
    protected static final String GENELOCATION_TAG = "geneLocation";
    protected static final String GENELOCATION_NAME_TAG = "name";
    
    protected static final String REFERENCE_TAG = "reference";
    protected static final String CITATION_TAG = "citation";
    protected static final String TITLE_TAG = "title";
    protected static final String EDITOR_LIST_TAG = "editorList";
    protected static final String AUTHOR_LIST_TAG = "authorList";
    protected static final String PERSON_TAG = "person";
    protected static final String CONSORTIUM_TAG = "consortium";
    protected static final String LOCATOR_TAG = "locator";
    protected static final String RP_LINE_TAG = "scope";
    protected static final String RC_LINE_TAG = "source";
    protected static final String RC_SPECIES_TAG = "species";
    protected static final String RC_TISSUE_TAG = "tissue";
    protected static final String RC_TRANSP_TAG = "transposon";
    protected static final String RC_STRAIN_TAG = "strain";
    protected static final String RC_PLASMID_TAG = "plasmid";
    
    protected static final String COMMENT_TAG = "comment";
    protected static final String COMMENT_MASS_ATTR = "mass";
    protected static final String COMMENT_ERROR_ATTR = "error";
    protected static final String COMMENT_METHOD_ATTR = "method";
    protected static final String COMMENT_LOCTYPE_ATTR = "locationType";
    
    protected static final String COMMENT_ABSORPTION_TAG = "absorption";
    protected static final String COMMENT_ABS_MAX_TAG = "max";
    protected static final String COMMENT_KINETICS_TAG = "kinetics";
    protected static final String COMMENT_KIN_KM_TAG = "KM";
    protected static final String COMMENT_KIN_VMAX_TAG = "VMax";
    protected static final String COMMENT_PH_TAG = "phDependence";
    protected static final String COMMENT_REDOX_TAG = "redoxPotential";
    protected static final String COMMENT_TEMPERATURE_TAG = "temperatureDependence";
    protected static final String COMMENT_LINK_TAG = "link";
    protected static final String COMMENT_LINK_URI_ATTR = "uri";
    protected static final String COMMENT_EVENT_TAG = "event";
    protected static final String COMMENT_ISOFORM_TAG = "isoform";
    protected static final String COMMENT_INTERACTANT_TAG = "interactant";
    protected static final String COMMENT_INTERACT_INTACT_ATTR = "intactId";
    protected static final String COMMENT_INTERACT_LABEL_TAG = "label";
    protected static final String COMMENT_ORGANISMS_TAG = "organismsDiffer";
    protected static final String COMMENT_EXPERIMENTS_TAG = "experiments";
    
    protected static final String NOTE_TAG = "note";
    protected static final String KEYWORD_TAG = "keyword";
    protected static final String PROTEIN_EXISTS_TAG = "proteinExistence";
    protected static final String ID_TAG = "id";
    
    protected static final String FEATURE_TAG = "feature";
    protected static final String FEATURE_DESC_ATTR = "description";
    protected static final String FEATURE_ORIGINAL_TAG = "original";
    protected static final String FEATURE_VARIATION_TAG = "variation";
    
    protected static final String EVIDENCE_TAG = "evidence";
    protected static final String EVIDENCE_CATEGORY_ATTR = "category";
    protected static final String EVIDENCE_ATTRIBUTE_ATTR = "attribute";
    protected static final String EVIDENCE_DATE_ATTR = "date";
    
    protected static final String LOCATION_TAG = "location";
    protected static final String LOCATION_SEQ_ATTR = "sequence";
    protected static final String LOCATION_BEGIN_TAG = "begin";
    protected static final String LOCATION_END_TAG = "end";
    protected static final String LOCATION_POSITION_ATTR = "position";
    protected static final String LOCATION_POSITION_TAG = "position";
    
    protected static final String SEQUENCE_TAG = "sequence";
    protected static final String SEQUENCE_VERSION_ATTR = "version";
    protected static final String SEQUENCE_LENGTH_ATTR = "length";
    protected static final String SEQUENCE_MASS_ATTR = "mass";
    protected static final String SEQUENCE_CHECKSUM_ATTR = "checksum";
    protected static final String SEQUENCE_MODIFIED_ATTR = "modified";
    
    // RP line parser
    protected static final Pattern rppat = Pattern.compile("SEQUENCE OF (\\d+)-(\\d+)");
    
    protected static final Pattern xmlSchema = Pattern.compile(".*http://www\\.uniprot\\.org/support/docs/uniprot\\.xsd.*");
    
    /**
     * Implements some UniProtXML-specific terms.
     */
    public static class Terms extends RichSequence.Terms {        
        public static final String CONTAINS_PREFIX = "Contains:";
        public static final String INCLUDES_PREFIX = "Includes:";
        
        public static final String GENENAME_KEY = "primary";
        public static final String GENESYNONYM_KEY = "synonym";
        public static final String ORDLOCNAME_KEY = "ordered locus";
        public static final String ORFNAME_KEY = "ORF";
        
        public static final String NCBI_TAXON_KEY = "NCBI Taxonomy";
        public static final String COMMON_NAME_KEY = "common";
        public static final String FULL_NAME_KEY = "full";
        public static final String SCIENTIFIC_NAME_KEY = "scientific";
        public static final String SYNONYM_NAME_KEY = "synonym";
        public static final String ABBREV_NAME_KEY = "abbreviation";
        
        public static final String LOC_FUZZY_START_KEY = "less than";
        public static final String LOC_FUZZY_END_KEY = "greater than";
        
        // Ontology for uniprot keywords (because they have identifiers, aaargh...)
        private static ComparableOntology uniprotKWOnto = null;
        
        /**
         * Getter for the protein exists term
         * @return The protein exists Term
         */
        public static ComparableTerm getProteinExistsTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("UniProt protein exists");
        }
        
        /**
         * Getter for the private uniprot ontology.
         * @return the ontology.
         */
        public static ComparableOntology getUniprotKWOnto() {
            return (ComparableOntology)RichObjectFactory.getObject(SimpleComparableOntology.class, new Object[]{"uniprot_kw"});
        }
        
        /**
         * Getter for the UniProtXML term
         * @return The UniProtXML Term
         */
        public static ComparableTerm getUniProtXMLTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("UniProtXML");
        }
        
        /**
         * Getter for the protein type term
         * @return The protein type Term
         */
        public static ComparableTerm getProteinTypeTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("protein_type");
        }
        
        /**
         * Getter for the evidence category term
         * @return The evidence category Term
         */
        public static ComparableTerm getEvidenceCategoryTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("evidence_category");
        }
        
        /**
         * Getter for the evidence type term
         * @return The evidence type Term
         */
        public static ComparableTerm getEvidenceTypeTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("evidence_type");
        }
        
        /**
         * Getter for the evidence date term
         * @return The evidence date Term
         */
        public static ComparableTerm getEvidenceDateTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("evidence_date");
        }
        
        /**
         * Getter for the evidence attr term
         * @return The evidence attr Term
         */
        public static ComparableTerm getEvidenceAttrTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("evidence_attr");
        }
        
        /**
         * Getter for the feature ref term
         * @return The feature ref Term
         */
        public static ComparableTerm getFeatureRefTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("feature_ref");
        }
        
        /**
         * Getter for the feature status term
         * @return The feature status Term
         */
        public static ComparableTerm getFeatureStatusTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("feature_status");
        }
        
        /**
         * Getter for the feature original term
         * @return The feature original Term
         */
        public static ComparableTerm getFeatureOriginalTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("feature_original");
        }
        
        /**
         * Getter for the feature variation term
         * @return The feature variation Term
         */
        public static ComparableTerm getFeatureVariationTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("feature_variation");
        }
        
        /**
         * Getter for the location seq term
         * @return The location seq Term
         */
        public static ComparableTerm getLocationSequenceTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("locseq");
        }
    }    
    
    /**
     * {@inheritDoc}
     * A file is in UniProtXML format if the second XML line contains the phrase "http://www.uniprot.org/support/docs/uniprot.xsd".
     */
    @Override
    public boolean canRead(File file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        br.readLine(); // skip first line
        String secondLine = br.readLine();
        boolean readable = secondLine!=null && xmlSchema.matcher(secondLine).matches(); // check on second line
        br.close();
        return readable;
    }
    
    /**
     * {@inheritDoc}
     * Always returns a protein tokenizer.
     */
    @Override
    public SymbolTokenization guessSymbolTokenization(File file) throws IOException {
        return RichSequence.IOTools.getProteinParser();
    }
    
    /**
     * {@inheritDoc}
     * A stream is in UniProtXML format if the second XML line contains the phrase "http://www.uniprot.org/support/docs/uniprot.xsd".
     */
    public boolean canRead(BufferedInputStream stream) throws IOException {
        stream.mark(2000); // some streams may not support this
        BufferedReader br = new BufferedReader(new InputStreamReader(stream));
        br.readLine(); // skip first line
        String secondLine = br.readLine();
        boolean readable = secondLine!=null && xmlSchema.matcher(secondLine).matches(); // check on second line
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
    
    /**
     * {@inheritDoc}
     * If namespace is null, then the namespace of the sequence in the fasta is used.
     * If the namespace is null and so is the namespace of the sequence in the fasta,
     * then the default namespace is used.
     */
    public boolean readRichSequence(BufferedReader reader,
            SymbolTokenization symParser,
            RichSeqIOListener rlistener,
            Namespace ns)
            throws IllegalSymbolException, IOException, ParseException {
        
        Pattern copyright = Pattern.compile(".*<"+COPYRIGHT_TAG+".*");
        
        try {
            rlistener.startSequence();                   
            DefaultHandler m_handler = new UniProtXMLHandler(this,symParser,rlistener,ns);
            boolean hasMore=XMLTools.readXMLChunk(reader, m_handler, ENTRY_TAG);
            // deal with copyright chunk
            reader.mark(10000);
            String line = reader.readLine();
            reader.reset();
            if (copyright.matcher(line).matches()) XMLTools.readXMLChunk(reader, m_handler, COPYRIGHT_TAG);
            // all done!
            rlistener.endSequence();
            return hasMore;
        } catch (ParserConfigurationException e) {
            throw new ParseException(e);
        } catch (SAXException e) {
            throw new ParseException(e);
        }
    }
    
    private PrintWriter pw;
    private XMLWriter xml;
    
    /**
     * {@inheritDoc}
     */
    public void beginWriting() throws IOException {
        // make an XML writer
        pw = new PrintWriter(this.getPrintStream());
        xml = new PrettyXMLWriter(pw);
        xml.printRaw("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>");
        xml.openTag(ENTRY_GROUP_TAG);
        xml.attribute("xmlns","http://uniprot.org/uniprot");
        xml.attribute("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance");
        xml.attribute("xsi:schemaLocation","http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd");
    }
    
    /**
     * {@inheritDoc}
     */
    public void finishWriting() throws IOException {
        xml.closeTag(ENTRY_GROUP_TAG);
        pw.flush();
    }
    
    /**
     * {@inheritDoc}
     */
    public void	writeSequence(Sequence seq, PrintStream os) throws IOException {
        if (this.getPrintStream()==null) this.setPrintStream(this.getPrintStream());
        this.writeSequence(seq, RichObjectFactory.getDefaultNamespace());
    }
    
    /**
     * {@inheritDoc}
     */
    public void writeSequence(Sequence seq, String format, PrintStream os) throws IOException {
        if (this.getPrintStream()==null) this.setPrintStream(this.getPrintStream());
        if (!format.equals(this.getDefaultFormat())) throw new IllegalArgumentException("Unknown format: "+format);
        this.writeSequence(seq, RichObjectFactory.getDefaultNamespace());
    }
    
    /**
     * {@inheritDoc}
     * If namespace is null, then the sequence's own namespace is used.
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
        
        int key = 1;
        
        Set notes = rs.getNoteSet();
        List accessions = new ArrayList();
        List kws = new ArrayList();
        String cdat = null;
        String udat = null;
        String arel = null;
        String adat = null;
        String copyright = null;
        String proteinType = null;
        String proteinExists = null;
        Map genenames = new TreeMap();
        Map genesynonyms = new TreeMap();
        Map orfnames = new TreeMap();
        Map ordlocnames = new TreeMap();
        Set evidenceIDs = new TreeSet();
        Set organelles = new TreeSet();
        Map evcats = new TreeMap();
        Map evtypes = new TreeMap();
        Map evdates = new TreeMap();
        Map evattrs = new TreeMap();
        Map speciesRecs = new TreeMap();
        Map strainRecs = new TreeMap();
        Map tissueRecs = new TreeMap();
        Map transpRecs = new TreeMap();
        Map plasmidRecs = new TreeMap();
        for (Iterator i = notes.iterator(); i.hasNext();) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getDateCreatedTerm())) cdat=n.getValue();
            else if (n.getTerm().equals(Terms.getDateUpdatedTerm())) udat=n.getValue();
            else if (n.getTerm().equals(Terms.getRelAnnotatedTerm())) arel=n.getValue();
            else if (n.getTerm().equals(Terms.getDateAnnotatedTerm())) adat=n.getValue();
            else if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) accessions.add(n.getValue());
            else if (n.getTerm().equals(Terms.getOrganelleTerm())) organelles.add(n.getValue());
            else if (n.getTerm().equals(Terms.getKeywordTerm())) {
                ComparableTerm t = Terms.getUniprotKWOnto().getOrCreateTerm(n.getValue());
                try {
                    if (t.getIdentifier()==null || t.getIdentifier().length()==0) t.setIdentifier("UNKNOWN");
                } catch (ChangeVetoException ce) {
                    IOException e = new IOException("Failed to assign keyword identifier");
                    e.initCause(ce);
                    throw e;
                }
                kws.add(t);
            } else if (n.getTerm().equals(Terms.getCopyrightTerm())) copyright=n.getValue();
            else if (n.getTerm().equals(Terms.getProteinTypeTerm())) proteinType=n.getValue();
            else if (n.getTerm().equals(Terms.getProteinExistsTerm())) proteinExists=n.getValue();
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
            else if (n.getTerm().equals(Terms.getEvidenceCategoryTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                evcats.put(refID, ref.substring(colon+1)); // map of id -> string as only one name per gene
                evidenceIDs.add(refID);
            } else if (n.getTerm().equals(Terms.getEvidenceTypeTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                evtypes.put(refID, ref.substring(colon+1)); // map of id -> string as only one name per gene
                evidenceIDs.add(refID);
            } else if (n.getTerm().equals(Terms.getEvidenceDateTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                evdates.put(refID, ref.substring(colon+1)); // map of id -> string as only one name per gene
                evidenceIDs.add(refID);
            } else if (n.getTerm().equals(Terms.getEvidenceAttrTerm()))  {
                String ref = n.getValue();
                int colon = ref.indexOf(':');
                Integer refID = new Integer(0);
                if (colon>=1) refID = new Integer(ref.substring(0,colon));
                evattrs.put(refID, ref.substring(colon+1)); // map of id -> string as only one name per gene
                evidenceIDs.add(refID);
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
        
        xml.openTag(ENTRY_TAG);
        xml.attribute(ENTRY_VERSION_ATTR,""+(arel==null?""+rs.getVersion():arel));
        xml.attribute(ENTRY_NAMESPACE_ATTR,(ns==null?rs.getNamespace().getName():ns.getName()));
        xml.attribute(ENTRY_CREATED_ATTR,cdat);
        xml.attribute(ENTRY_UPDATED_ATTR,(adat==null?cdat:adat)); // annotation update
        
        xml.openTag(ACCESSION_TAG);
        xml.print(rs.getAccession());
        xml.closeTag(ACCESSION_TAG);
        
        xml.openTag(NAME_TAG);
        xml.print(rs.getName());
        xml.closeTag(NAME_TAG);
        
        xml.openTag(PROTEIN_TAG);
        if (proteinType!=null) xml.attribute(TYPE_ATTR,proteinType);
        String desc = rs.getDescription().trim(); // this is only going to make sense if it was a UniProt seq to start with
        if (desc.endsWith(".")) desc = desc.substring(0, desc.length()-1); // chomp trailing dot
        String[] parts = desc.split("\\[");
        for (int j = 0 ; j < parts.length; j++) {
            if (parts[j].startsWith(Terms.CONTAINS_PREFIX)) {
                // contains section
                String chunk = parts[j].substring(Terms.CONTAINS_PREFIX.length()+1).trim();
                if (chunk.endsWith("]")) chunk = chunk.substring(0, chunk.length()-1); // chomp trailing ]
                String[] moreparts = chunk.split(";");
                for (int k = 0; k < moreparts.length; k++) {
                    xml.openTag(DOMAIN_TAG);
                    String[] names = moreparts[k].split("\\(");
                    for (int l = 0; l < names.length; l++) {
                        String name = names[l].trim();
                        if (name.endsWith(")")) name = name.substring(0,name.length()-1); // chomp trailing )
                        xml.openTag(NAME_TAG);
                        xml.print(name);
                        xml.closeTag(NAME_TAG);
                    }
                    xml.closeTag(DOMAIN_TAG);
                }
            } else if (parts[j].startsWith(Terms.INCLUDES_PREFIX)) {
                // includes section
                String chunk = parts[j].substring(Terms.INCLUDES_PREFIX.length()+1).trim();
                if (chunk.endsWith("]")) chunk = chunk.substring(0, chunk.length()-1); // chomp trailing ]
                String[] moreparts = chunk.split(";");
                for (int k = 0; k < moreparts.length; k++) {
                    xml.openTag(COMPONENT_TAG);
                    String[] names = moreparts[k].split("\\(");
                    for (int l = 0; l < names.length; l++) {
                        String name = names[l].trim();
                        if (name.endsWith(")")) name = name.substring(0,name.length()-1); // chomp trailing )
                        xml.openTag(NAME_TAG);
                        xml.print(name);
                        xml.closeTag(NAME_TAG);
                    }
                    xml.closeTag(COMPONENT_TAG);
                }
            } else {
                // plain names
                String[] names = parts[j].split("\\(");
                for (int l = 0; l < names.length; l++) {
                    String name = names[l].trim();
                    if (name.endsWith(")")) name = name.substring(0,name.length()-1); // chomp trailing )
                    xml.openTag(NAME_TAG);
                    xml.print(name);
                    xml.closeTag(NAME_TAG);
                }
            }
        }
        xml.closeTag(PROTEIN_TAG);
        
        // gene line
        for (Iterator i = genenames.keySet().iterator(); i.hasNext(); ) {
            Integer geneid = (Integer)i.next();
            String genename = (String)genenames.get(geneid);
            List synonyms = (List)genesynonyms.get(geneid);
            List orfs = (List)orfnames.get(geneid);
            List ordlocs = (List)ordlocnames.get(geneid);
            
            xml.openTag(GENE_TAG);
            
            xml.openTag(NAME_TAG);
            xml.attribute(TYPE_ATTR,Terms.GENENAME_KEY);
            xml.print(genename);
            xml.closeTag(NAME_TAG);
            
            if (synonyms!=null) {
                for (Iterator j = synonyms.iterator(); j.hasNext(); ) {
                    xml.openTag(NAME_TAG);
                    xml.attribute(TYPE_ATTR,Terms.GENESYNONYM_KEY);
                    xml.print((String)j.next());
                    xml.closeTag(NAME_TAG);
                }
            }
            if (ordlocs!=null) {
                for (Iterator j = synonyms.iterator(); j.hasNext(); ) {
                    xml.openTag(NAME_TAG);
                    xml.attribute(TYPE_ATTR,Terms.ORDLOCNAME_KEY);
                    xml.print((String)j.next());
                    xml.closeTag(NAME_TAG);
                }
            }
            if (orfs!=null) {
                for (Iterator j = synonyms.iterator(); j.hasNext(); ) {
                    xml.openTag(NAME_TAG);
                    xml.attribute(TYPE_ATTR,Terms.ORFNAME_KEY);
                    xml.print((String)j.next());
                    xml.closeTag(NAME_TAG);
                }
            }
            
            xml.closeTag(GENE_TAG);
        }
        
        // source line (from taxon)
        //   organism line
        NCBITaxon tax = rs.getTaxon();
        if (tax!=null) {
            xml.openTag(ORGANISM_TAG);
            xml.attribute(KEY_ATTR,""+(key++));
            
            for (Iterator i = tax.getNameClasses().iterator(); i.hasNext(); ) {
                String nameclass = (String)i.next();
                String ournameclass = Terms.COMMON_NAME_KEY;
                if (nameclass.equalsIgnoreCase(Terms.FULL_NAME_KEY)) ournameclass = NCBITaxon.EQUIVALENT;
                else if (nameclass.equalsIgnoreCase(Terms.SCIENTIFIC_NAME_KEY)) ournameclass = NCBITaxon.SCIENTIFIC;
                else if (nameclass.equalsIgnoreCase(Terms.SYNONYM_NAME_KEY)) ournameclass = NCBITaxon.SYNONYM;
                else if (nameclass.equalsIgnoreCase(Terms.ABBREV_NAME_KEY)) ournameclass = NCBITaxon.ACRONYM;
                for (Iterator j = tax.getNames(nameclass).iterator(); j.hasNext(); ) {
                    xml.openTag(NAME_TAG);
                    xml.attribute(TYPE_ATTR,ournameclass);
                    xml.print((String)j.next());
                    xml.closeTag(NAME_TAG);
                }
            }
            
            xml.openTag(DBXREF_TAG);
            xml.attribute(KEY_ATTR,""+(key++));
            xml.attribute(TYPE_ATTR,Terms.NCBI_TAXON_KEY);
            xml.attribute(ID_ATTR,""+tax.getNCBITaxID());
            xml.closeTag(DBXREF_TAG);
            
            String h = tax.getNameHierarchy();
            h = h.substring(0, h.length()-1); // chomp dot
            String[] hierarch = h.split(";");
            xml.openTag(LINEAGE_TAG);
            for (int j = 0; j < hierarch.length; j++) {
                xml.openTag(TAXON_TAG);
                xml.print(hierarch[j].trim());
                xml.closeTag(TAXON_TAG);
            }
            xml.closeTag(LINEAGE_TAG);
            
            xml.closeTag(ORGANISM_TAG);
        }
        
        // gene location line (organelle)
        for (Iterator i = organelles.iterator(); i.hasNext(); ) {
            String org = (String)i.next();
            xml.openTag(GENELOCATION_TAG);
            if (org.startsWith("Plasmid")) {
                xml.attribute(TYPE_ATTR,"plasmid");
                String[] subparts = org.split(",");
                for (int j = 0; j < parts.length; j++) {
                    org = subparts[j].trim();
                    if (org.startsWith("and")) org = org.substring(3).trim();
                    org = org.substring("Plasmid".length()).trim();
                    xml.openTag(GENELOCATION_NAME_TAG);
                    xml.attribute(STATUS_ATTR,"known");
                    xml.print(org);
                    xml.closeTag(GENELOCATION_NAME_TAG);
                }
            } else {
                xml.attribute(TYPE_ATTR,org.toLowerCase()); // uniprotxml must have lower case
            }
            xml.closeTag(GENELOCATION_TAG);
        }
        
        // docrefs
        for (Iterator i = rs.getRankedDocRefs().iterator(); i.hasNext(); ) {
            RankedDocRef rdr = (RankedDocRef)i.next();
            DocRef dr = rdr.getDocumentReference();
            
            xml.openTag(REFERENCE_TAG);
            xml.attribute(KEY_ATTR,""+(key++));
            
            xml.openTag(CITATION_TAG);
            xml.attribute(TYPE_ATTR,"journal article"); // faking it i know
            
            if (dr.getTitle()!=null) {
                xml.openTag(TITLE_TAG);
                xml.print(dr.getTitle());
                xml.closeTag(TITLE_TAG);
            }
            
            List auths = new ArrayList(dr.getAuthorList());
            List editors = new ArrayList(auths);
            for (final Iterator j = editors.iterator(); j.hasNext(); ) {
                DocRefAuthor a = (DocRefAuthor)j.next();
                if (!a.isEditor())
                    j.remove();
                else
                    auths.remove(a);
            }
            if (!editors.isEmpty()) {
                xml.openTag(EDITOR_LIST_TAG);
                for (Iterator j = editors.iterator(); j.hasNext(); ) {
                    DocRefAuthor a = (DocRefAuthor)j.next();
                    if (a.isEditor()) {
                        if (a.isConsortium()) {
                            xml.openTag(CONSORTIUM_TAG);
                            xml.attribute(NAME_ATTR,a.getName());
                            xml.closeTag(CONSORTIUM_TAG);
                        } else {
                            xml.openTag(PERSON_TAG);
                            xml.attribute(NAME_ATTR,a.getName());
                            xml.closeTag(PERSON_TAG);
                        }
                    }
                }
                xml.closeTag(EDITOR_LIST_TAG);
            }
            if (!auths.isEmpty()) {
                xml.openTag(AUTHOR_LIST_TAG);
                for (Iterator j = auths.iterator(); j.hasNext(); ) {
                    DocRefAuthor a = (DocRefAuthor)j.next();
                    if (a.isConsortium()) {
                        xml.openTag(CONSORTIUM_TAG);
                        xml.attribute(NAME_ATTR,a.getName());
                        xml.closeTag(CONSORTIUM_TAG);
                    } else {
                        xml.openTag(PERSON_TAG);
                        xml.attribute(NAME_ATTR,a.getName());
                        xml.closeTag(PERSON_TAG);
                    }
                }
                xml.closeTag(AUTHOR_LIST_TAG);
            }
            
            xml.openTag(LOCATOR_TAG);
            xml.print(dr.getLocation());
            xml.closeTag(LOCATOR_TAG);
            
            CrossRef cr = dr.getCrossref();
            if (cr!=null) {
                xml.openTag(DBXREF_TAG);
                xml.attribute(TYPE_ATTR,cr.getDbname());
                xml.attribute(ID_ATTR,cr.getAccession());
                xml.attribute(KEY_ATTR,""+(key++));
                if (!cr.getNoteSet().isEmpty()) {
                    for (Iterator j = cr.getNoteSet().iterator(); j.hasNext(); ) {
                        Note n = (Note)j.next();
                        xml.openTag(PROPERTY_TAG);
                        xml.attribute(TYPE_ATTR,n.getTerm().getName());
                        xml.attribute(VALUE_ATTR,n.getValue());
                        xml.closeTag(PROPERTY_TAG);
                    }
                }
                xml.closeTag(DBXREF_TAG);
            }
            
            xml.closeTag(CITATION_TAG);
            
            // RP
            xml.openTag(RP_LINE_TAG);
            xml.print(dr.getRemark());
            xml.closeTag(RP_LINE_TAG);
            // Print out ref position if present
            if (rdr.getStart()!=null && rdr.getEnd()!=null && !rppat.matcher(dr.getRemark()).matches()) {
                xml.openTag(RP_LINE_TAG);
                xml.print("SEQUENCE OF "+rdr.getStart()+"-"+rdr.getEnd()+".");
                xml.closeTag(RP_LINE_TAG);
            }
            
            // RC
            boolean rcOpened = false;
            Integer rank = new Integer(rdr.getRank());
            if (speciesRecs.get(rank)!=null) {
                if (!rcOpened) {
                    xml.openTag(RC_LINE_TAG);
                    rcOpened = true;
                }
                for (Iterator j = ((List)speciesRecs.get(rank)).iterator(); j.hasNext(); ) {
                    xml.openTag(RC_SPECIES_TAG);
                    xml.print((String)j.next());
                    xml.closeTag(RC_SPECIES_TAG);
                }
            }
            if (strainRecs.get(rank)!=null) {
                if (!rcOpened) {
                    xml.openTag(RC_LINE_TAG);
                    rcOpened = true;
                }
                for (Iterator j = ((List)strainRecs.get(rank)).iterator(); j.hasNext(); ) {
                    xml.openTag(RC_STRAIN_TAG);
                    xml.print((String)j.next());
                    xml.closeTag(RC_STRAIN_TAG);
                }
            }
            if (tissueRecs.get(rank)!=null) {
                if (!rcOpened) {
                    xml.openTag(RC_LINE_TAG);
                    rcOpened = true;
                }
                for (Iterator j = ((List)tissueRecs.get(rank)).iterator(); j.hasNext(); ) {
                    xml.openTag(RC_TISSUE_TAG);
                    xml.print((String)j.next());
                    xml.closeTag(RC_TISSUE_TAG);
                }
            }
            if (transpRecs.get(rank)!=null) {
                if (!rcOpened) {
                    xml.openTag(RC_LINE_TAG);
                    rcOpened = true;
                }
                for (Iterator j = ((List)transpRecs.get(rank)).iterator(); j.hasNext(); ) {
                    xml.openTag(RC_TRANSP_TAG);
                    xml.print((String)j.next());
                    xml.closeTag(RC_TRANSP_TAG);
                }
            }
            if (plasmidRecs.get(rank)!=null) {
                if (!rcOpened) {
                    xml.openTag(RC_LINE_TAG);
                    rcOpened = true;
                }
                for (Iterator j = ((List)plasmidRecs.get(rank)).iterator(); j.hasNext(); ) {
                    xml.openTag(RC_PLASMID_TAG);
                    xml.print((String)j.next());
                    xml.closeTag(RC_PLASMID_TAG);
                }
            }
            if (rcOpened)
                xml.closeTag(RC_LINE_TAG);
            
            xml.closeTag(REFERENCE_TAG);
        }
        
        // comments
        for (Iterator i = rs.getComments().iterator(); i.hasNext(); ) {
            // use UniProtCommentParser to convert each text comment from string to object
            // do not print unconvertible ones (eg. no -!- on text)
            Comment c = (Comment)i.next();
            if (UniProtCommentParser.isParseable(c)) {
                // otherwise parse and display appropriately
                UniProtCommentParser ucp = new UniProtCommentParser();
                try {
                    ucp.parseComment(c);
                } catch (ParseException ce) {
                    IOException e = new IOException("Failed to parse comment when outputting");
                    e.initCause(ce);
                    throw e;
                }
                String type = ucp.getCommentType();
                String xtype = type.toLowerCase(); // uniprotxml requires lower case
                if (type.equalsIgnoreCase(UniProtCommentParser.PTM)) xtype = "posttranslational modification";
                else if (type.equalsIgnoreCase(UniProtCommentParser.DATABASE)) xtype = "online information";
                
                xml.openTag(COMMENT_TAG);
                xml.attribute(TYPE_ATTR,xtype);
                
                // database comment
                if (type.equalsIgnoreCase(UniProtCommentParser.DATABASE)) {
                    xml.attribute(NAME_ATTR,ucp.getDatabaseName());
                    
                    xml.openTag(COMMENT_LINK_TAG);
                    xml.attribute(COMMENT_LINK_URI_ATTR,ucp.getUri());
                    xml.closeTag(COMMENT_LINK_TAG);
                }
                // mass spec
                else if (type.equalsIgnoreCase(UniProtCommentParser.MASS_SPECTROMETRY)) {
                    xml.attribute(COMMENT_MASS_ATTR,""+ucp.getMolecularWeight());
                    if (ucp.getMolWeightError()!=null) xml.attribute(COMMENT_ERROR_ATTR,""+ucp.getMolWeightError());
                    xml.attribute(COMMENT_METHOD_ATTR,""+ucp.getMolWeightMethod());
                    
                    xml.openTag(LOCATION_TAG);
                    xml.openTag(LOCATION_BEGIN_TAG);
                    xml.attribute(LOCATION_POSITION_ATTR,""+ucp.getMolWeightRangeStart());
                    xml.closeTag(LOCATION_BEGIN_TAG);
                    xml.openTag(LOCATION_END_TAG);
                    xml.attribute(LOCATION_POSITION_ATTR,""+ucp.getMolWeightRangeEnd());
                    xml.closeTag(LOCATION_END_TAG);
                    xml.closeTag(LOCATION_TAG);
                }
                // interaction
                else if (type.equalsIgnoreCase(UniProtCommentParser.INTERACTION)) {
                    // UniProt flat allows for multiple interactions per comment, but
                    // UniProtXML only allows for a single one. So, we have to open/close
                    // and write additional comments as necessary.
                    for (Iterator j = ucp.getInteractions().iterator(); j.hasNext(); ) {
                        // process comment
                        Interaction interact = (Interaction)j.next();
                        
                        xml.openTag(COMMENT_INTERACTANT_TAG);
                        xml.attribute(COMMENT_INTERACT_INTACT_ATTR,interact.getFirstIntActID());
                        xml.closeTag(COMMENT_INTERACTANT_TAG);
                        
                        xml.openTag(COMMENT_INTERACTANT_TAG);
                        xml.attribute(COMMENT_INTERACT_INTACT_ATTR,interact.getSecondIntActID());
                        xml.openTag(ID_TAG);
                        xml.print(interact.getID());
                        xml.closeTag(ID_TAG);
                        if (interact.getLabel()!=null) {
                            xml.openTag(COMMENT_INTERACT_LABEL_TAG);
                            xml.print(interact.getLabel());
                            xml.closeTag(COMMENT_INTERACT_LABEL_TAG);
                        }
                        xml.closeTag(COMMENT_INTERACTANT_TAG);
                        
                        xml.openTag(COMMENT_ORGANISMS_TAG);
                        xml.print(interact.isOrganismsDiffer()?"true":"false");
                        xml.closeTag(COMMENT_ORGANISMS_TAG);
                        
                        xml.openTag(COMMENT_EXPERIMENTS_TAG);
                        xml.print(""+interact.getNumberExperiments());
                        xml.closeTag(COMMENT_EXPERIMENTS_TAG);
                        
                        // if has next, close and open next comment tag
                        if (j.hasNext()) {
                            xml.closeTag(COMMENT_TAG);
                            xml.openTag(COMMENT_TAG);
                            xml.attribute(TYPE_ATTR,xtype);
                        }
                    }
                }
                // alternative products
                else if (type.equalsIgnoreCase(UniProtCommentParser.ALTERNATIVE_PRODUCTS)) {
                    for (Iterator j = ucp.getEvents().iterator(); j.hasNext(); ) {
                        Event event = (Event)j.next();
                        xml.openTag(COMMENT_EVENT_TAG);
                        xml.attribute(TYPE_ATTR,event.getType().toLowerCase()); // uniprotxml requires lowercase
                        xml.closeTag(COMMENT_EVENT_TAG);
                    }
                    for (Iterator j = ucp.getIsoforms().iterator(); j.hasNext(); ) {
                        Isoform isoform = (Isoform)j.next();
                        xml.openTag(COMMENT_ISOFORM_TAG);
                        for (Iterator k = isoform.getIsoIDs().iterator(); k.hasNext(); ) {
                            xml.openTag(ID_TAG);
                            xml.print((String)k.next());
                            xml.closeTag(ID_TAG);
                        }
                        for (Iterator k = isoform.getNames().iterator(); k.hasNext(); ) {
                            xml.openTag(NAME_TAG);
                            xml.print((String)k.next());
                            xml.closeTag(NAME_TAG);
                        }
                        xml.openTag(SEQUENCE_TAG);
                        xml.attribute(TYPE_ATTR,isoform.getSequenceType().toLowerCase());
                        if (isoform.getSequenceType().equalsIgnoreCase("Described")) {
                            xml.attribute(REF_ATTR,isoform.getSequenceRef());
                        }
                        xml.closeTag(SEQUENCE_TAG);
                        xml.openTag(NOTE_TAG);
                        xml.print(isoform.getNote());
                        xml.closeTag(NOTE_TAG);
                        xml.closeTag(COMMENT_ISOFORM_TAG);
                    }
                }
                // biophysicoblahblah stuff
                else if (type.equalsIgnoreCase(UniProtCommentParser.BIOPHYSICOCHEMICAL_PROPERTIES)) {
                    if (ucp.getAbsorptionNote()!=null) {
                        xml.openTag(COMMENT_ABSORPTION_TAG);
                        xml.openTag(COMMENT_ABS_MAX_TAG);
                        xml.print(ucp.getAbsorptionMax());
                        xml.closeTag(COMMENT_ABS_MAX_TAG);
                        xml.openTag(TEXT_TAG);
                        xml.print(ucp.getAbsorptionNote());
                        xml.closeTag(TEXT_TAG);
                        xml.closeTag(COMMENT_ABSORPTION_TAG);
                    }
                    if (ucp.getKineticsNote()!=null) {
                        xml.openTag(COMMENT_KINETICS_TAG);
                        for (Iterator j = ucp.getKMs().iterator(); j.hasNext(); ) {
                            xml.openTag(COMMENT_KIN_KM_TAG);
                            xml.print((String)j.next());
                            xml.closeTag(COMMENT_KIN_KM_TAG);
                        }
                        for (Iterator j = ucp.getVMaxes().iterator(); j.hasNext(); ) {
                            xml.openTag(COMMENT_KIN_VMAX_TAG);
                            xml.print((String)j.next());
                            xml.closeTag(COMMENT_KIN_VMAX_TAG);
                        }
                        xml.openTag(TEXT_TAG);
                        xml.print(ucp.getKineticsNote());
                        xml.closeTag(TEXT_TAG);
                        xml.closeTag(COMMENT_KINETICS_TAG);
                    }
                    if (ucp.getPHDependence()!=null) {
                        xml.openTag(COMMENT_PH_TAG);
                        xml.print(ucp.getPHDependence());
                        xml.closeTag(COMMENT_PH_TAG);
                    }
                    if (ucp.getRedoxPotential()!=null) {
                        xml.openTag(COMMENT_REDOX_TAG);
                        xml.print(ucp.getRedoxPotential());
                        xml.closeTag(COMMENT_REDOX_TAG);
                    }
                    if (ucp.getTemperatureDependence()!=null) {
                        xml.openTag(COMMENT_TEMPERATURE_TAG);
                        xml.print(ucp.getTemperatureDependence());
                        xml.closeTag(COMMENT_TEMPERATURE_TAG);
                    }
                }
                // all other comments
                else {
                    xml.openTag(TEXT_TAG);
                    xml.print(ucp.getText());
                    xml.closeTag(TEXT_TAG);
                }
                
                // finish comment up
                if (ucp.getNote()!=null) {
                    xml.openTag(NOTE_TAG);
                    xml.print(ucp.getNote());
                    xml.closeTag(NOTE_TAG);
                }
                
                xml.closeTag(COMMENT_TAG);
            }
        }
        
        // xrefs
        for (Iterator i = rs.getRankedCrossRefs().iterator(); i.hasNext(); ) {
            RankedCrossRef rcr = (RankedCrossRef)i.next();
            CrossRef cr = rcr.getCrossRef();
            
            xml.openTag(DBXREF_TAG);
            String dbname = cr.getDbname();
            xml.attribute(TYPE_ATTR,dbname);
            xml.attribute(ID_ATTR,cr.getAccession());
            xml.attribute(KEY_ATTR,""+(key++));
            if (!cr.getNoteSet().isEmpty()) {
                int acccount = 2;
                for (Iterator j = cr.getNoteSet().iterator(); j.hasNext(); ) {
                    Note n = (Note)j.next();
                    if (n.getTerm().equals(Terms.getAdditionalAccessionTerm()) && !n.getValue().equals("-")) {
                        xml.openTag(PROPERTY_TAG);
                        String name = n.getTerm().getName();
                        if (acccount==2) {
                            // SECONDARY IDENTIFIER
                            if (dbname.equalsIgnoreCase("HIV") ||
                                    dbname.equalsIgnoreCase("INTERPRO") ||
                                    dbname.equalsIgnoreCase("PANTHER") ||
                                    dbname.equalsIgnoreCase("PFAM") ||
                                    dbname.equalsIgnoreCase("PIR") ||
                                    dbname.equalsIgnoreCase("PRINTS") ||
                                    dbname.equalsIgnoreCase("PRODOM") ||
                                    dbname.equalsIgnoreCase("REBASE") ||
                                    dbname.equalsIgnoreCase("SMART") ||
                                    dbname.equalsIgnoreCase("TIGRFAMS")) {
                                // the secondary identifier is the entry name.
                                name = "entry name";
                            } else if (dbname.equalsIgnoreCase("PDB")) {
                                // the secondary identifier is the structure determination method, which is controlled vocabulary that currently includes: X-ray(for X-ray crystallography), NMR(for NMR spectroscopy), EM(for electron microscopy and cryo-electron diffraction), Fiber(for fiber diffraction), IR(for infrared spectroscopy), Model(for predicted models) and Neutron(for neutron diffraction).
                                name = "structure determination method";
                            } else if (dbname.equalsIgnoreCase("DICTYBASE") ||
                                    dbname.equalsIgnoreCase("ECOGENE") ||
                                    dbname.equalsIgnoreCase("FLYBASE") ||
                                    dbname.equalsIgnoreCase("HGNC") ||
                                    dbname.equalsIgnoreCase("MGI") ||
                                    dbname.equalsIgnoreCase("RGD") ||
                                    dbname.equalsIgnoreCase("SGD") ||
                                    dbname.equalsIgnoreCase("STYGENE") ||
                                    dbname.equalsIgnoreCase("SUBTILIST") ||
                                    dbname.equalsIgnoreCase("WORMBASE") ||
                                    dbname.equalsIgnoreCase("ZFIN")) {
                                // the secondary identifier is the gene designation. If the gene designation is not available, a dash('-') is used.
                                name = "gene designation";
                            } else if (dbname.equalsIgnoreCase("GO")) {
                                // the second identifier is a 1-letter abbreviation for one of the 3 ontology aspects, separated from the GO term by a column. If the term is longer than 46 characters, the first 43 characters are indicated followed by 3 dots('...'). The abbreviations for the 3 distinct aspects of the ontology are P(biological Process), F(molecular Function), and C(cellular Component).
                                name = "term";
                            } else if (dbname.equalsIgnoreCase("HAMAP")) {
                                // the secondary identifier indicates if a domain is 'atypical' and/or 'fused', otherwise the field is empty('-').
                                name = "domain";
                            } else if (dbname.equalsIgnoreCase("ECO2DBASE")) {
                                // the secondary identifier is the latest release number or edition of the database that has been used to derive the cross-reference.
                                name = "release number";
                            } else if (dbname.equalsIgnoreCase("SWISS-2DPAGE") ||
                                    dbname.equalsIgnoreCase("HSC-2DPAGE")) {
                                // the secondary identifier is the species or tissue of origin.
                                name = "organism name";
                            } else if (dbname.equalsIgnoreCase("ENSEMBL")) {
                                // the secondary identifier is the species of origin.
                                name = "organism name";
                            } else if (dbname.equalsIgnoreCase("PIRSF")) {
                                // the secondary identifier is the protein family name.
                                name = "protein family name";
                            } else if (dbname.equalsIgnoreCase("AARHUS") ||
                                    dbname.equalsIgnoreCase("GHENT-2DPAGE")) {
                                // the secondary identifier is either 'IEF' (for isoelectric focusing) or 'NEPHGE' (for non-equilibrium pH gradient electrophoresis).
                                name = "secondary identifier";
                            } else if (dbname.equalsIgnoreCase("WORMPEP")) {
                                // the secondary identifier is a number attributed by the C.elegans genome-sequencing project to that protein.
                                name = "C.elegans number";
                            } else if (dbname.equalsIgnoreCase("AGD") ||
                                    dbname.equalsIgnoreCase("ANU-2DPAGE") ||
                                    dbname.equalsIgnoreCase("COMPLUYEAST-2DPAGE") ||
                                    dbname.equalsIgnoreCase("ECHOBASE") ||
                                    dbname.equalsIgnoreCase("GENEDB_SPOMBE") ||
                                    dbname.equalsIgnoreCase("GERMONLINE") ||
                                    dbname.equalsIgnoreCase("GLYCOSUITEDB") ||
                                    dbname.equalsIgnoreCase("GRAMENE") ||
                                    dbname.equalsIgnoreCase("H-INVDB") ||
                                    dbname.equalsIgnoreCase("INTACT") ||
                                    dbname.equalsIgnoreCase("LEGIOLIST") ||
                                    dbname.equalsIgnoreCase("LEPROMA") ||
                                    dbname.equalsIgnoreCase("LISTILIST") ||
                                    dbname.equalsIgnoreCase("MAIZEDB") ||
                                    dbname.equalsIgnoreCase("MEROPS") ||
                                    dbname.equalsIgnoreCase("MIM") ||
                                    dbname.equalsIgnoreCase("MYPULIST") ||
                                    dbname.equalsIgnoreCase("OGP") ||
                                    dbname.equalsIgnoreCase("PHCI-2DPAGE") ||
                                    dbname.equalsIgnoreCase("PHOSSITE") ||
                                    dbname.equalsIgnoreCase("PHOTOLIST") ||
                                    dbname.equalsIgnoreCase("PMMA-2DPAGE") ||
                                    dbname.equalsIgnoreCase("RAT-HEART-2DPAGE") ||
                                    dbname.equalsIgnoreCase("REACTOME") ||
                                    dbname.equalsIgnoreCase("SAGALIST") ||
                                    dbname.equalsIgnoreCase("SIENA-2DPAGE") ||
                                    dbname.equalsIgnoreCase("TAIR") ||
                                    dbname.equalsIgnoreCase("TIGR") ||
                                    dbname.equalsIgnoreCase("TRANSFAC") ||
                                    dbname.equalsIgnoreCase("TUBERCULIST")) {
                                // the secondary identifier is not used and a dash('-') is stored in that field.
                                // should never get here - I hope!
                            } else if (dbname.equalsIgnoreCase("HSSP")) {
                                // the secondary identifier is the entry name of the PDB structure related to that of the entry in which the HSSP cross-reference is present.
                                name = "entry name";
                            } else if (dbname.equalsIgnoreCase("GENEFARM")) {
                                // the secondary identifier is the gene family identifier. If the gene family identifier is not available, a dash('-') is used.
                                name = "gene family";
                            } else if (dbname.equalsIgnoreCase("SMR")) {
                                // the secondary identifier indicates the range(s) relevant to the structure model(s).
                                name = "range";
                            } else if (dbname.equalsIgnoreCase("EMBL") ||
                                    dbname.equalsIgnoreCase("DDBJ") ||
                                    dbname.equalsIgnoreCase("GENBANK")) {
                                // PROTEIN_ID; STATUS_IDENTIFIER; MOLECULE_TYPE
                                name = "protein id";
                            } else if (dbname.equalsIgnoreCase("PROSITE")) {
                                // ENTRY_NAME; STATUS.
                                name = "entry name";
                            }
                        } else if (acccount==3) {
                            // TERTIARY IDENTIFIER
                            if (dbname.equalsIgnoreCase("HAMAP") ||
                                    dbname.equalsIgnoreCase("PANTHER") ||
                                    dbname.equalsIgnoreCase("PFAM") ||
                                    dbname.equalsIgnoreCase("PIRSF") ||
                                    dbname.equalsIgnoreCase("PRODOM") ||
                                    dbname.equalsIgnoreCase("SMART") ||
                                    dbname.equalsIgnoreCase("TIGRFAMS")) {
                                // the tertiary identifier is the number of hits found in the sequence.
                                name = "number of hits";
                            } else if (dbname.equalsIgnoreCase("GO")) {
                                // the tertiary identifier is a 3-character GO evidence code. The meaning of the evidence codes is: IDA=inferred from direct assay, IMP=inferred from mutant phenotype, IGI=inferred from genetic interaction, IPI=inferred from physical interaction, IEP=inferred from expression pattern, TAS=traceable author statement, NAS=non-traceable author statement, IC=inferred by curator, ISS=inferred from sequence or structural similarity.
                                name = "evidence";
                            } else if (dbname.equalsIgnoreCase("PDB")) {
                                // the tertiary identifier indicates the chain(s) and the corresponding range, of which the structure has been determined. If the range is unknown, a dash is given rather than the range positions(e.g. 'A/B=-.'), if the chains and the range is unknown, a dash is used.
                                name = "chains";
                            } else if (dbname.equalsIgnoreCase("EMBL") ||
                                    dbname.equalsIgnoreCase("DDBJ") ||
                                    dbname.equalsIgnoreCase("GENBANK")) {
                                // PROTEIN_ID; STATUS_IDENTIFIER; MOLECULE_TYPE
                                name = "status identifier";
                            } else if (dbname.equalsIgnoreCase("PROSITE")) {
                                // ENTRY_NAME; STATUS.
                                name = "status";
                            }
                        } else {
                            // QUATERNARY AND ADDITIONAL
                            if (dbname.equalsIgnoreCase("EMBL") ||
                                    dbname.equalsIgnoreCase("DDBJ") ||
                                    dbname.equalsIgnoreCase("GENBANK")) {
                                // PROTEIN_ID; STATUS_IDENTIFIER; MOLECULE_TYPE
                                name = "molecule type";
                            }
                        }
                        xml.attribute(TYPE_ATTR,name);
                        xml.attribute(VALUE_ATTR,n.getValue());
                        xml.closeTag(PROPERTY_TAG);
                        acccount++;
                    }
                }
            }
            xml.closeTag(DBXREF_TAG);
        }
        
        // protein exists
        xml.openTag(PROTEIN_EXISTS_TAG);
        xml.attribute(TYPE_ATTR,proteinExists);
        xml.closeTag(PROTEIN_EXISTS_TAG);
        
        // keywords
        for (Iterator j = kws.iterator(); j.hasNext(); ) {
            ComparableTerm t = (ComparableTerm)j.next();
            xml.openTag(KEYWORD_TAG);
            xml.attribute(ID_ATTR,t.getIdentifier());
            xml.print(t.getName());
            xml.closeTag(KEYWORD_TAG);
        }
        
        // features
        for (Iterator i = rs.getFeatureSet().iterator(); i.hasNext(); ) {
            RichFeature f = (RichFeature)i.next();
            String descr = null;
            String ftid = null;
            String ref = null;
            String status = null;
            String original = null;
            String locseq = null;
            List variation = new ArrayList();
            for (Iterator j = f.getNoteSet().iterator(); j.hasNext(); ) {
                Note n = (Note)j.next();
                if (n.getTerm().equals(Terms.getFTIdTerm())) ftid = n.getValue();
                else if (n.getTerm().equals(Terms.getFeatureDescTerm())) descr = n.getValue();
                else if (n.getTerm().equals(Terms.getFeatureStatusTerm())) status = n.getValue();
                else if (n.getTerm().equals(Terms.getFeatureRefTerm())) ref = n.getValue();
                else if (n.getTerm().equals(Terms.getFeatureOriginalTerm())) original = n.getValue();
                else if (n.getTerm().equals(Terms.getFeatureVariationTerm())) variation.add(n.getValue());
                else if (n.getTerm().equals(Terms.getLocationSequenceTerm())) locseq = n.getValue();
            }
            
            xml.openTag(FEATURE_TAG);
            
            xml.attribute(TYPE_ATTR,f.getTypeTerm().getName()); // TODO : need to translate from UniProt flatfile format names?
            if (ftid!=null) xml.attribute(ID_ATTR,ftid);
            if (descr!=null) xml.attribute(FEATURE_DESC_ATTR,descr);
            if (ref!=null) xml.attribute(REF_ATTR,ref);
            if (status!=null) xml.attribute(STATUS_ATTR,status);
            if (original!=null) {
                xml.openTag(FEATURE_ORIGINAL_TAG);
                xml.print(original.trim());
                xml.closeTag(FEATURE_ORIGINAL_TAG);
            }
            for (Iterator j = variation.iterator(); j.hasNext(); ) {
                xml.openTag(FEATURE_VARIATION_TAG);
                xml.print(((String)j.next()).trim());
                xml.closeTag(FEATURE_VARIATION_TAG);
            }
            
            xml.openTag(LOCATION_TAG);
            if (locseq!=null) xml.attribute(LOCATION_SEQ_ATTR,locseq.trim());
            RichLocation rl = (RichLocation)f.getLocation();
            if (rl.getMinPosition().equals(rl.getMaxPosition())) {
                // point position
                xml.openTag(LOCATION_POSITION_TAG);
                if (rl.getMinPosition().getFuzzyStart() || rl.getMaxPosition().getFuzzyStart()) xml.attribute(STATUS_ATTR,"less than");
                else if (rl.getMinPosition().getFuzzyEnd() || rl.getMaxPosition().getFuzzyEnd()) xml.attribute(STATUS_ATTR,"greater than");
                xml.attribute(LOCATION_POSITION_ATTR,""+rl.getMin());
                xml.closeTag(LOCATION_POSITION_TAG);
            } else {
                // range position
                // begin
                xml.openTag(LOCATION_BEGIN_TAG);
                Position begin = rl.getMinPosition();
                if (begin.getFuzzyStart()) xml.attribute(STATUS_ATTR,"less than");
                else if (begin.getFuzzyEnd()) xml.attribute(STATUS_ATTR,"greater than");
                xml.attribute(LOCATION_POSITION_ATTR,""+begin.getStart());
                xml.closeTag(LOCATION_BEGIN_TAG);
                // end
                xml.openTag(LOCATION_END_TAG);
                Position end = rl.getMaxPosition();
                if (end.getFuzzyStart()) xml.attribute(STATUS_ATTR,"less than");
                else if (end.getFuzzyEnd()) xml.attribute(STATUS_ATTR,"greater than");
                xml.attribute(LOCATION_POSITION_ATTR,""+end.getEnd());
                xml.closeTag(LOCATION_END_TAG);
            }
            xml.closeTag(LOCATION_TAG);
            
            xml.closeTag(FEATURE_TAG);
        }
        
        // evidence
        for (Iterator i = evidenceIDs.iterator(); i.hasNext(); ) {
            Integer evidenceID = (Integer)i.next();
            String cat = (String)evcats.get(evidenceID);
            String type = (String)evtypes.get(evidenceID);
            String date = (String)evdates.get(evidenceID);
            String attr = (String)evattrs.get(evidenceID);
            
            xml.openTag(EVIDENCE_TAG);
            xml.attribute(KEY_ATTR,""+(key++));
            xml.attribute(EVIDENCE_CATEGORY_ATTR,cat);
            xml.attribute(EVIDENCE_DATE_ATTR,date);
            xml.attribute(TYPE_ATTR,type);
            if (attr!=null) xml.attribute(EVIDENCE_ATTRIBUTE_ATTR,attr);
            xml.closeTag(EVIDENCE_TAG);
        }
        
        // sequence
        int mw = 0;
        try {
            mw = (int)MassCalc.getMolecularWeight(rs);
        } catch (IllegalSymbolException e) {
            throw new RuntimeException("Found illegal symbol", e);
        }
        CRC64Checksum crc = new CRC64Checksum();
        String seqstr = rs.seqString();
        crc.update(seqstr.getBytes(),0,seqstr.length());
        xml.openTag(SEQUENCE_TAG);
        xml.attribute(SEQUENCE_VERSION_ATTR,""+rs.getVersion());
        xml.attribute(SEQUENCE_LENGTH_ATTR,""+rs.length());
        xml.attribute(SEQUENCE_MASS_ATTR,""+mw);
        xml.attribute(SEQUENCE_CHECKSUM_ATTR,""+crc);
        xml.attribute(SEQUENCE_MODIFIED_ATTR,(udat==null?cdat:udat)); // sequence update
        String[] lines = StringTools.wordWrap(rs.seqString(), "\\s+", this.getLineWidth());
        for (int i = 0; i < lines.length; i ++) xml.println(lines[i]);
        xml.closeTag(SEQUENCE_TAG);
        
        // close entry
        xml.closeTag(ENTRY_TAG);
        
        // copyright (if present)
        if (copyright!=null) {
            xml.openTag(COPYRIGHT_TAG);
            xml.println(copyright);
            xml.closeTag(COPYRIGHT_TAG);
        }
        
        pw.flush();
    }
    
    /**
     * {@inheritDoc}
     */
    public String getDefaultFormat() {
        return UNIPROTXML_FORMAT;
    }
    
    // SAX event handler for parsing http://www.ebi.uniprot.org/support/docs/uniprot.xsd
    private class UniProtXMLHandler extends DefaultHandler {
        
        private RichSequenceFormat parent;
        private SymbolTokenization symParser;
        private RichSeqIOListener rlistener;
        private Namespace ns;
        private StringBuffer m_currentString;
        
        private NCBITaxon tax;
        private RichFeature.Template templ;
        private StringBuffer proteinDesc;
        private boolean firstNameInProteinGroup;
        private boolean firstDomainInProteinGroup;
        private boolean firstComponentInProteinGroup;
        private int currGene;
        private String geneNameClass;
        private String organismNameClass;
        private Map currNames = new TreeMap();
        private StringBuffer organelleDesc;
        private List currDBXrefs = new ArrayList();
        private List currComments = new ArrayList();
        private String currRefLocation;
        private List currRefAuthors;
        private String currRefTitle;
        private int currRefStart;
        private int currRefEnd;
        private int currRefRank;
        private String currPersonIs;
        private int currRCID;
        private int currEvID;
        private String currKWID;
        private UniProtCommentParser currUCParser;
        private Interaction currUCParserInteract;
        private Event currUCParserEvent;
        private Isoform currUCParserIsoform;
        private String currLocIsFor;
        private String currTextIsFor;
        private String currNoteIsFor;
        private String currSeqIsFor;
        private String currIDIsFor;
        private String currNameIsFor;
        private int interactantCount;
        private StringBuffer currLocStr;
        private int featNoteRank;
        
        // construct a new handler that will populate the given list of sequences
        private UniProtXMLHandler(RichSequenceFormat parent,
                SymbolTokenization symParser,
                RichSeqIOListener rlistener,
                Namespace ns) {
            this.parent = parent;
            this.symParser = symParser;
            this.rlistener = rlistener;
            this.ns = ns;
            this.m_currentString = new StringBuffer();
        }
        
        
        // process an opening tag
        @Override
        public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
            
            if (qName.equals(ENTRY_TAG)) {
                try {
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(ENTRY_NAMESPACE_ATTR) && this.ns==null) ns=(Namespace)RichObjectFactory.getObject(SimpleNamespace.class,new Object[]{val});
                        else if (name.equals(ENTRY_VERSION_ATTR)) rlistener.addSequenceProperty(Terms.getRelAnnotatedTerm(), val);
                        else if (name.equals(ENTRY_CREATED_ATTR)) rlistener.addSequenceProperty(Terms.getDateCreatedTerm(), val);
                        else if (name.equals(ENTRY_UPDATED_ATTR)) rlistener.addSequenceProperty(Terms.getDateAnnotatedTerm(), val);
                    }
                    if (this.ns==null) ns=RichObjectFactory.getDefaultNamespace();
                    rlistener.setNamespace(ns);
                } catch (ParseException e) {
                    throw new SAXException(e);
                }
                this.currNameIsFor = "ENTRY";
                this.currSeqIsFor = "ENTRY";
                this.currGene = 0;
                this.currNames.clear();
                this.currRefRank = 0;
                this.currRCID = 0;
                this.currEvID = 0;
            }
            
            else if (qName.equals(PROTEIN_TAG)) {
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i).trim();
                    String val = attributes.getValue(i).trim();
                    try {
                        if (name.equals(TYPE_ATTR)) rlistener.addSequenceProperty(Terms.getProteinTypeTerm(),val);
                    } catch (ParseException e) {
                        throw new SAXException(e);
                    }
                }
                this.proteinDesc = new StringBuffer();
                this.currNameIsFor = "PROTEIN";
                this.firstNameInProteinGroup = true;
                this.firstDomainInProteinGroup = true;
                this.firstComponentInProteinGroup = true;
            } else if (qName.equals(DOMAIN_TAG)) {
                if (!this.firstComponentInProteinGroup) proteinDesc.append("]");
                if (this.firstDomainInProteinGroup) proteinDesc.append(" ["+Terms.CONTAINS_PREFIX);
                else proteinDesc.append(";");
                this.firstDomainInProteinGroup = false;
                this.firstNameInProteinGroup = true;
            } else if (qName.equals(COMPONENT_TAG)) {
                if (!this.firstDomainInProteinGroup) proteinDesc.append("]");
                if (this.firstComponentInProteinGroup) proteinDesc.append(" ["+Terms.INCLUDES_PREFIX);
                else proteinDesc.append(";");
                this.firstComponentInProteinGroup = false;
                this.firstNameInProteinGroup = true;
            }
            
            else if (qName.equals(GENE_TAG)) {
                this.currGene++;
                this.currNameIsFor="GENE";
            }
            
            else if (qName.equals(NAME_TAG)) {
                if (this.currNameIsFor.equals("GENE")) {
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(TYPE_ATTR)) this.geneNameClass=val;
                    }
                }
                
                else if (this.currNameIsFor.equals("ORGANISM")) {
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(TYPE_ATTR)) this.organismNameClass=val;
                    }
                }
            }
            
            else if (qName.equals(ORGANISM_TAG)) {
                this.currNameIsFor="ORGANISM";
            }
            
            else if (qName.equals(DBXREF_TAG)) {
                if (this.currNameIsFor.equals("ORGANISM")) {
                    Integer taxID = null;
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(ID_ATTR)) taxID = Integer.valueOf(val);
                    }
                    try {
                        tax = (NCBITaxon)RichObjectFactory.getObject(SimpleNCBITaxon.class, new Object[]{taxID});
                        rlistener.setTaxon(tax);
                        for (Iterator j = currNames.keySet().iterator(); j.hasNext(); ) {
                            String nameClass = (String)j.next();
                            Set nameSet = (Set)this.currNames.get(nameClass);
                            try {
                                for (Iterator k = nameSet.iterator(); k.hasNext(); ) {
                                    String name = (String)k.next();
                                    tax.addName(nameClass,name);
                                }
                            } catch (ChangeVetoException ce) {
                                throw new ParseException(ce);
                            }
                        }
                    } catch (ParseException e) {
                        throw new SAXException(e);
                    }
                    this.currNames.clear();
                }
                
                else {
                    String type = null;
                    String id = null;
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(ID_ATTR)) id = val;
                        else if (name.equals(TYPE_ATTR)) type = val;
                    }
                    CrossRef dbx = (CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{type, id, new Integer(0)});
                    this.currDBXrefs.add(dbx);
                }
            } else if (qName.equals(PROPERTY_TAG)) {
                String id = null;
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(VALUE_ATTR)) id = val;
                }
                Note note = new SimpleNote(Terms.getAdditionalAccessionTerm(),id,1);
                try {
                    int last = this.currDBXrefs.size();
                    ((CrossRef)this.currDBXrefs.get(last-1)).getRichAnnotation().addNote(note);
                } catch (ChangeVetoException ce) {
                    SAXException pe = new SAXException("Could not annotate identifier terms");
                    pe.initCause(ce);
                    throw pe;
                }
            }
            
            else if (qName.equals(GENELOCATION_TAG)) {
                this.currNameIsFor = "ORGANELLE";
                this.organelleDesc = new StringBuffer();
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(TYPE_ATTR)) {
                        val = val.toUpperCase().charAt(0)+val.substring(1); // init caps for flat format
                        if (!val.equals("Plasmid")) this.organelleDesc.append(val);
                    }
                }
            }
            
            else if (qName.equals(REFERENCE_TAG) && !this.parent.getElideReferences()) {
                this.currRefLocation = null;
                this.currRefAuthors = new ArrayList();
                this.currRefTitle = null;
                this.currDBXrefs.clear();
                this.currComments.clear();
                this.currRefRank++;
                this.currRefStart = -999;
                this.currRefEnd = -999;
            } else if (qName.equals(CITATION_TAG) && !this.parent.getElideReferences()) {
                StringBuffer currRef = new StringBuffer();
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    // combine everything except type into a fake reference to use if locator is a no-show
                    if (!name.equals(TYPE_ATTR)) {
                        if (currRef.length()>0) currRef.append(" ");
                        currRef.append(val);
                    }
                }
                this.currRefLocation = currRef.toString();
            } else if (qName.equals(EDITOR_LIST_TAG)) {
                this.currPersonIs = "EDITOR";
            } else if (qName.equals(AUTHOR_LIST_TAG)) {
                this.currPersonIs = "AUTHOR";
            }  else if (qName.equals(PERSON_TAG)) {
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(NAME_ATTR)) {
                        if (this.currPersonIs.equals("AUTHOR")) currRefAuthors.add(new SimpleDocRefAuthor(val, false, false));
                        else if (this.currPersonIs.equals("EDITOR")) currRefAuthors.add(new SimpleDocRefAuthor(val, false, true));
                    }
                }
            } else if (qName.equals(CONSORTIUM_TAG)) {
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(NAME_ATTR)) {
                        if (this.currPersonIs.equals("AUTHOR")) currRefAuthors.add(new SimpleDocRefAuthor(val, true, false));
                        else if (this.currPersonIs.equals("EDITOR")) currRefAuthors.add(new SimpleDocRefAuthor(val, true, true));
                    }
                }
            }  else if (qName.equals(RC_LINE_TAG)) {
                this.currRCID++;
            }
            
            else if (qName.equals(PROTEIN_EXISTS_TAG)) {
                try {
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(TYPE_ATTR)) rlistener.addSequenceProperty(Terms.getProteinExistsTerm(),val);
                    }
                } catch (ParseException e) {
                    SAXException pe = new SAXException("Could not annotate protein exists terms");
                    pe.initCause(e);
                    throw pe;
                }
            }
            
            else if (qName.equals(KEYWORD_TAG)) {
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(ID_ATTR)) this.currKWID = val;
                }
            }
            
            else if (qName.equals(EVIDENCE_TAG)) {
                this.currEvID++;
                try {
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(EVIDENCE_CATEGORY_ATTR)) rlistener.addSequenceProperty(Terms.getEvidenceCategoryTerm(),val);
                        else if (name.equals(EVIDENCE_DATE_ATTR)) rlistener.addSequenceProperty(Terms.getEvidenceDateTerm(),val);
                        else if (name.equals(TYPE_ATTR)) rlistener.addSequenceProperty(Terms.getEvidenceTypeTerm(),val);
                        else if (name.equals(EVIDENCE_ATTRIBUTE_ATTR)) rlistener.addSequenceProperty(Terms.getEvidenceAttrTerm(),val);
                    }
                } catch (ParseException e) {
                    SAXException pe = new SAXException("Could not annotate evidence terms");
                    pe.initCause(e);
                    throw pe;
                }
            }
            
            else if (qName.equals(LOCATION_TAG)) {
                this.currLocStr = new StringBuffer();
                if (this.currLocIsFor.equals("FEATURE")) {
                    try {
                        for (int i = 0; i < attributes.getLength(); i++) {
                            String name = attributes.getQName(i);
                            String val = attributes.getValue(i);
                            if (name.equals(LOCATION_SEQ_ATTR)) {
                                Note note = new SimpleNote(Terms.getLocationSequenceTerm(), val, this.featNoteRank++);
                                ((RichAnnotation)templ.annotation).addNote(note);
                            }
                        }
                    } catch (ChangeVetoException e) {
                        SAXException pe = new SAXException("Could not create location terms");
                        pe.initCause(e);
                        throw pe;
                    }
                }
            } else if (qName.equals(LOCATION_BEGIN_TAG) || qName.equals(LOCATION_END_TAG) || qName.equals(LOCATION_POSITION_TAG)) {
                StringBuffer pos = new StringBuffer();
                pos.append(" "); // space between start and end
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(STATUS_ATTR)) {
                        if (val.equals("less than")) pos.append("<");
                        else if (val.equals("greater than")) pos.append(">");
                    } else if (name.equals(LOCATION_POSITION_ATTR)) {
                        pos.append(val);
                    }
                }
                this.currLocStr.append(pos.toString());
                if (qName.equals(LOCATION_POSITION_TAG)) currLocStr.append(pos.toString()); // fake it as begin=end
            }
            
            else if (qName.equals(FEATURE_TAG) && !this.parent.getElideFeatures()) {
                this.featNoteRank = 1;
                templ = new RichFeature.Template();
                templ.annotation = new SimpleRichAnnotation();
                templ.sourceTerm = Terms.getUniProtXMLTerm();
                templ.featureRelationshipSet = new TreeSet();
                templ.rankedCrossRefs = new TreeSet();
                try {
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(TYPE_ATTR)) {
                            templ.typeTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm(val);
                        } else if (name.equals(ID_ATTR)) {
                            Note note = new SimpleNote(Terms.getFTIdTerm(), val, this.featNoteRank++);
                            ((RichAnnotation)templ.annotation).addNote(note);
                        } else if (name.equals(FEATURE_DESC_ATTR)) {
                            Note note = new SimpleNote(Terms.getFeatureDescTerm(), val, this.featNoteRank++);
                            ((RichAnnotation)templ.annotation).addNote(note);
                        } else if (name.equals(STATUS_ATTR)) {
                            Note note = new SimpleNote(Terms.getFeatureStatusTerm(), val, this.featNoteRank++);
                            ((RichAnnotation)templ.annotation).addNote(note);
                        } else if (name.equals(REF_ATTR)) {
                            Note note = new SimpleNote(Terms.getFeatureRefTerm(), val, this.featNoteRank++);
                            ((RichAnnotation)templ.annotation).addNote(note);
                        }
                    }
                } catch (ChangeVetoException e) {
                    SAXException pe = new SAXException("Could not create location terms");
                    pe.initCause(e);
                    throw pe;
                }
                this.currLocStr = new StringBuffer();
                this.currLocIsFor = "FEATURE";
            }
            
            else if (qName.equals(COMMENT_TAG)) {
                this.currUCParser = new UniProtCommentParser();
                this.currUCParser.setInteractions(new ArrayList());
                this.currUCParser.setEvents(new ArrayList());
                this.currUCParser.setIsoforms(new ArrayList());
                this.currUCParser.setKMs(new ArrayList());
                this.currUCParser.setVMaxes(new ArrayList());
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i).trim();
                    String val = attributes.getValue(i).trim();
                    if (name.equals(TYPE_ATTR)) {
                        String type = val.toUpperCase(); // easier to check this way, plus flat uniprot requires it
                        if (type.equals("POSTTRANSLATIONAL MODIFICATION")) type="PTM";
                        else if (type.equals("ONLINE INFORMATION")) type="DATABASE";
                        currUCParser.setCommentType(type);
                    } else if (name.equals(COMMENT_MASS_ATTR)) this.currUCParser.setMolecularWeight(Integer.parseInt(val));
                    else if (name.equals(COMMENT_ERROR_ATTR)) this.currUCParser.setMolWeightError(Integer.valueOf(val));
                    else if (name.equals(COMMENT_METHOD_ATTR)) this.currUCParser.setMolWeightMethod(val);
                    else if (name.equals(NAME_ATTR)) this.currUCParser.setDatabaseName(val);
                }
                this.currLocIsFor="COMMENT";
                this.currTextIsFor="COMMENT";
                this.currNoteIsFor="COMMENT";
                this.interactantCount = 0;
            } else if (qName.equals(COMMENT_ABSORPTION_TAG)) {
                this.currTextIsFor="ABSORPTION";
            } else if (qName.equals(COMMENT_KINETICS_TAG)) {
                this.currTextIsFor="KINETICS";
            } else if (qName.equals(COMMENT_LINK_TAG)) {
                this.currTextIsFor="KINETICS";
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(COMMENT_LINK_URI_ATTR)) this.currUCParser.setUri(val);
                }
            } else if (qName.equals(COMMENT_EVENT_TAG)) {
                this.currUCParserEvent = new Event();
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(TYPE_ATTR)) {
                        val = val.toUpperCase().charAt(0)+val.substring(1); // make first letter upper case for flat uniprot
                        this.currUCParserEvent.setType(val);
                    }
                }
                currUCParser.getEvents().add(currUCParserEvent);
            } else if (qName.equals(COMMENT_ISOFORM_TAG)) {
                this.currUCParserIsoform = new Isoform();
                this.currUCParser.getIsoforms().add(currUCParserIsoform);
                this.currUCParserEvent.setNamedIsoforms(this.currUCParser.getIsoforms().size());
                this.currNameIsFor="ISOFORM";
                this.currNoteIsFor="ISOFORM";
                this.currSeqIsFor="ISOFORM";
                this.currIDIsFor="ISOFORM";
            } else if (qName.equals(COMMENT_INTERACTANT_TAG)) {
                this.currIDIsFor="INTERACTION";
                this.interactantCount++;
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(COMMENT_INTERACT_INTACT_ATTR)) {
                        if (this.interactantCount%2==1) {
                            this.currUCParserInteract = new Interaction();
                            this.currUCParserInteract.setFirstIntActID(val);
                            this.currUCParser.getInteractions().add(this.currUCParserInteract);
                        }
                        else this.currUCParserInteract.setSecondIntActID(val);
                    }
                }
            }
            
            else if (qName.equals(SEQUENCE_TAG)) {
                if (this.currSeqIsFor.equals("ENTRY")) {
                    try {
                        for (int i = 0; i < attributes.getLength(); i++) {
                            String name = attributes.getQName(i);
                            String val = attributes.getValue(i);
                            if (name.equals(SEQUENCE_MODIFIED_ATTR)) {
                                rlistener.addSequenceProperty(Terms.getDateUpdatedTerm(),val);
                            }
                            else if (name.equals(SEQUENCE_VERSION_ATTR))
                                rlistener.setVersion(Integer.parseInt(val));
                        }
                    } catch (ParseException e) {
                        SAXException pe = new SAXException("Could not set sequence properties");
                        pe.initCause(e);
                        throw pe;
                    }
                }
                
                else if (this.currSeqIsFor.equals("ISOFORM")) {
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(TYPE_ATTR)) {
                            val = val.toUpperCase().charAt(0)+val.substring(1); // init caps for flat uniprot
                            this.currUCParserIsoform.setSequenceType(val);
                        } else if (name.equals(REF_ATTR)) {
                            this.currUCParserIsoform.setSequenceRef(val);
                        }
                    }
                }
            }
        }
        
        // process a closing tag - we will have read the text already
        @Override
        public void endElement(String uri, String localName, String qName) throws SAXException {
            String val = this.m_currentString.toString().trim();
            
            try {
                if (qName.equals(COPYRIGHT_TAG)) {
                    rlistener.addSequenceProperty(Terms.getCopyrightTerm(),val);
                }
                
                else if (qName.equals(ACCESSION_TAG)) {
                    rlistener.setAccession(val);
                } else if (qName.equals(NAME_TAG)) {
                    if (this.currNameIsFor.equals("ENTRY")) rlistener.setName(val);
                    
                    else if (this.currNameIsFor.equals("PROTEIN")) {
                        if (this.firstNameInProteinGroup) {
                            proteinDesc.append(" ");
                            proteinDesc.append(val);
                        } else {
                            proteinDesc.append(" (");
                            proteinDesc.append(val);
                            proteinDesc.append(")");
                        }
                        this.firstNameInProteinGroup = false;
                    }
                    
                    else if (this.currNameIsFor.equals("GENE")) {
                        if (this.geneNameClass.equals(Terms.GENENAME_KEY)) rlistener.addSequenceProperty(Terms.getGeneNameTerm(), this.currGene+":"+val);
                        else if (this.geneNameClass.equals(Terms.GENESYNONYM_KEY)) rlistener.addSequenceProperty(Terms.getGeneSynonymTerm(), this.currGene+":"+val);
                        else if (this.geneNameClass.equals(Terms.ORDLOCNAME_KEY)) rlistener.addSequenceProperty(Terms.getOrderedLocusNameTerm(), this.currGene+":"+val);
                        else if (this.geneNameClass.equals(Terms.ORFNAME_KEY)) rlistener.addSequenceProperty(Terms.getORFNameTerm(), this.currGene+":"+val);
                    }
                    
                    else if (this.currNameIsFor.equals("ORGANISM")) {
                        String ournameclass = NCBITaxon.COMMON;
                        if (this.organismNameClass.equals(Terms.ABBREV_NAME_KEY)) ournameclass = NCBITaxon.ACRONYM;
                        else if (this.organismNameClass.equals(Terms.FULL_NAME_KEY)) ournameclass = NCBITaxon.EQUIVALENT;
                        else if (this.organismNameClass.equals(Terms.SCIENTIFIC_NAME_KEY)) ournameclass = NCBITaxon.SCIENTIFIC;
                        else if (this.organismNameClass.equals(Terms.SYNONYM_NAME_KEY)) ournameclass = NCBITaxon.SYNONYM;
                        if (!this.currNames.containsKey(ournameclass)) this.currNames.put(ournameclass,new TreeSet());
                        ((Set)this.currNames.get(ournameclass)).add(val);
                    }
                    
                    else if (this.currNameIsFor.equals("ORGANELLE")) {
                        this.organelleDesc.append(", Plasmid ");
                        this.organelleDesc.append(val);
                    }
                    
                    else if (this.currNameIsFor.equals("ISOFORM")) {
                        this.currUCParserIsoform.getNames().add(val);
                    }
                }
                
                else if (qName.equals(PROTEIN_TAG)) {
                    if (!this.firstDomainInProteinGroup || !this.firstComponentInProteinGroup) this.proteinDesc.append("]");
                    this.proteinDesc.append(".");
                    rlistener.setDescription(this.proteinDesc.toString());
                }
                
                else if (qName.equals(ORGANISM_TAG)) {
                    this.currNameIsFor="";
                }
                
                else if (qName.equals(GENELOCATION_TAG)) {
                    String total = this.organelleDesc.toString().substring(3); // chomp leading ", "
                    int lastComma = total.lastIndexOf(',');
                    if (lastComma>-1) {
                        this.organelleDesc.insert(lastComma+1," and");
                        total = this.organelleDesc.toString();
                    }
                    rlistener.addSequenceProperty(Terms.getOrganelleTerm(), total);
                }
                
                else if (qName.equals(RC_SPECIES_TAG)) {
                    rlistener.addSequenceProperty(Terms.getSpeciesTerm(), this.currRCID+":"+val);
                } else if (qName.equals(RC_TISSUE_TAG)) {
                    rlistener.addSequenceProperty(Terms.getTissueTerm(), this.currRCID+":"+val);
                } else if (qName.equals(RC_TRANSP_TAG)) {
                    rlistener.addSequenceProperty(Terms.getTransposonTerm(), this.currRCID+":"+val);
                } else if (qName.equals(RC_PLASMID_TAG)) {
                    rlistener.addSequenceProperty(Terms.getPlasmidTerm(), this.currRCID+":"+val);
                }
                
                else if (qName.equals(TITLE_TAG)) {
                    this.currRefTitle = val;
                } else if (qName.equals(LOCATOR_TAG)) {
                    this.currRefLocation = val;
                } else if (qName.equals(RP_LINE_TAG)) {
                    this.currComments.add(val);
                    // Try to use it to find the location of the reference, if we have one.
                    Matcher m = rppat.matcher(val);
                    if (m.matches()) {
                        this.currRefStart = Integer.parseInt(m.group(1));
                        this.currRefEnd = Integer.parseInt(m.group(2));
                    }
                } else if (qName.equals(REFERENCE_TAG) && !this.parent.getElideReferences()) {
                    // do the crossrefs
                    CrossRef useForDocRef = null;
                    for (Iterator j = this.currDBXrefs.iterator(); j.hasNext();) {
                        CrossRef dbx = (CrossRef)j.next();
                        RankedCrossRef rdbx = new SimpleRankedCrossRef(dbx,0);
                        rlistener.setRankedCrossRef(rdbx);
                        if (useForDocRef==null) useForDocRef = dbx;
                        else {
                            // medline gets priority, then pubmed - if multiple, use last
                            if (dbx.getDbname().equalsIgnoreCase(Terms.MEDLINE_KEY) || 
                                    (dbx.getDbname().equalsIgnoreCase(Terms.PUBMED_KEY) && 
                                    !useForDocRef.getDbname().equalsIgnoreCase(Terms.MEDLINE_KEY))) {
                                useForDocRef = dbx;
                            }
                        }
                    }
                    // do the comment - can only be one in this object model
                    String currRefRemark = null;
                    if (currComments.size()>0) currRefRemark = (String)currComments.iterator().next();
                    // create the docref object
                    try {
                        DocRef dr = (DocRef)RichObjectFactory.getObject(SimpleDocRef.class,new Object[]{currRefAuthors,currRefLocation,currRefTitle});
                        // assign the pubmed or medline to the docref - medline gets priority
                        if (useForDocRef!=null) dr.setCrossref(useForDocRef);
                        // assign the remarks
                        dr.setRemark(currRefRemark);
                        // assign the docref to the bioentry
                        RankedDocRef rdr = new SimpleRankedDocRef(dr,
                                (currRefStart != -999 ? new Integer(currRefStart) : null),
                                (currRefEnd != -999 ? new Integer(currRefEnd) : null),
                                currRefRank);
                        rlistener.setRankedDocRef(rdr);
                    } catch (ChangeVetoException e) {
                        throw new ParseException(e);
                    }
                    currDBXrefs.clear();
                    currComments.clear();
                }
                
                // keywords
                else if (qName.equals(KEYWORD_TAG)) {
                    // create and persist term
                    ComparableTerm t = Terms.getUniprotKWOnto().getOrCreateTerm(val);
                    try {
                        t.setIdentifier(currKWID);
                    } catch (ChangeVetoException e) {
                        throw new ParseException(e);
                    }
                    rlistener.addSequenceProperty(Terms.getKeywordTerm(), val);
                }
                
                else if (qName.equals(LOCATION_TAG)) {
                    if (currLocIsFor.equals("FEATURE")) {
                        templ.location = UniProtLocationParser.parseLocation(currLocStr.toString());
                    } else if (currLocIsFor.equals("COMMENT")) {
                        Location l = UniProtLocationParser.parseLocation(currLocStr.toString());
                        this.currUCParser.setMolWeightRangeStart(l.getMin());
                        this.currUCParser.setMolWeightRangeEnd(l.getMax());
                    }
                }
                
                else if (qName.equals(FEATURE_TAG)) {
                    // start the feature from the template we built
                    rlistener.startFeature(templ);
                    // end the feature
                    rlistener.endFeature();
                } else if (qName.equals(FEATURE_ORIGINAL_TAG)) {
                    try {
                        Note note = new SimpleNote(Terms.getFeatureOriginalTerm(), val, featNoteRank++);
                        ((RichAnnotation)templ.annotation).addNote(note);
                    } catch (ChangeVetoException e) {
                        SAXException pe = new SAXException("Could not create location terms");
                        pe.initCause(e);
                        throw pe;
                    }
                } else if (qName.equals(FEATURE_VARIATION_TAG)) {
                    try {
                        Note note = new SimpleNote(Terms.getFeatureVariationTerm(), val, featNoteRank++);
                        ((RichAnnotation)templ.annotation).addNote(note);
                    } catch (ChangeVetoException e) {
                        SAXException pe = new SAXException("Could not create location terms");
                        pe.initCause(e);
                        throw pe;
                    }
                }
                
                else if (qName.equals(COMMENT_TAG)) {
                    rlistener.setComment(currUCParser.generate());
                } else if (qName.equals(TEXT_TAG)) {
                    if (this.currTextIsFor.equals("COMMENT")) currUCParser.setText(val);
                    else if (this.currTextIsFor.equals("ABSORPTION")) currUCParser.setAbsorptionNote(val);
                    else if (this.currTextIsFor.equals("KINETICS")) currUCParser.setKineticsNote(val);
                } else if (qName.equals(COMMENT_ABS_MAX_TAG)) {
                    currUCParser.setAbsorptionMax(val);
                } else if (qName.equals(COMMENT_KIN_KM_TAG)) {
                    currUCParser.getKMs().add(val);
                } else if (qName.equals(COMMENT_KIN_VMAX_TAG)) {
                    currUCParser.getVMaxes().add(val);
                } else if (qName.equals(COMMENT_PH_TAG)) {
                    currUCParser.setPHDependence(val);
                } else if (qName.equals(COMMENT_REDOX_TAG)) {
                    currUCParser.setRedoxPotential(val);
                } else if (qName.equals(COMMENT_TEMPERATURE_TAG)) {
                    currUCParser.setTemperatureDependence(val);
                } else if (qName.equals(COMMENT_ORGANISMS_TAG)) {
                    if (val.equalsIgnoreCase("true")) currUCParserInteract.setOrganismsDiffer(true);
                    else currUCParserInteract.setOrganismsDiffer(false);
                } else if (qName.equals(COMMENT_EXPERIMENTS_TAG)) {
                    currUCParserInteract.setNumberExperiments(Integer.parseInt(val));
                } else if (qName.equals(NOTE_TAG)) {
                    if (currNoteIsFor.equals("COMMENT")) currUCParser.setNote(val);
                    else if (currNoteIsFor.equals("ISOFORM")) currUCParser.setNote(val);
                } else if (qName.equals(COMMENT_EVENT_TAG)) {
                    currUCParserEvent.setComment(val);
                } else if (qName.equals(COMMENT_ISOFORM_TAG)) {
                    this.currSeqIsFor = "ENTRY";
                    this.currNoteIsFor = "COMMENT";
                } else if (qName.equals(ID_TAG)) {
                    if (currIDIsFor.equals("ISOFORM")) currUCParserIsoform.getIsoIDs().add(val);
                    else if (currIDIsFor.equals("INTERACTION")) currUCParserInteract.setID(val);
                } else if (qName.equals(COMMENT_INTERACT_LABEL_TAG)) {
                    currUCParserInteract.setLabel(val);
                }
                
                else if (qName.equals(SEQUENCE_TAG)) {
                    if (this.currSeqIsFor.equals("ENTRY") && !this.parent.getElideSymbols()) {
                        try {
                            SymbolList sl = new SimpleSymbolList(symParser,
                                    val.replaceAll("\\s+","").replaceAll("[\\.|~]","-"));
                            rlistener.addSymbols(symParser.getAlphabet(),
                                    (Symbol[])(sl.toList().toArray(new Symbol[0])),
                                    0, sl.length());
                        } catch (Exception e) {
                            throw new ParseException(e);
                        }
                    }
                }
                
                else if (qName.equals(ENTRY_TAG)) {
                    // do the comments
                    for (Iterator j = currComments.iterator(); j.hasNext();) {
                        rlistener.setComment((String)j.next());
                    }
                    // do the crossrefs
                    for (Iterator j = currDBXrefs.iterator(); j.hasNext();) {
                        CrossRef dbx = (CrossRef)j.next();
                        RankedCrossRef rdbx = new SimpleRankedCrossRef(dbx, 0);
                        rlistener.setRankedCrossRef(rdbx);
                    }
                    // end the sequence
                    currComments.clear();
                    currDBXrefs.clear();
                }
                
            } catch (ParseException e) {
                throw new SAXException(e);
            }
            
            // drop old string
            this.m_currentString.setLength(0);
        }
        
        // process text inside tags
        @Override
        public void characters(char[] ch, int start, int length) {
            this.m_currentString.append(ch, start, length);
        }
    }
}

