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
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;

import javax.xml.parsers.ParserConfigurationException;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.IllegalSymbolException;
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
import org.biojavax.SimpleNote;
import org.biojavax.SimpleRankedCrossRef;
import org.biojavax.SimpleRankedDocRef;
import org.biojavax.SimpleRichAnnotation;
import org.biojavax.bio.seq.Position;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichLocation.Strand;
import org.biojavax.bio.taxa.NCBITaxon;
import org.biojavax.bio.taxa.SimpleNCBITaxon;
import org.biojavax.ontology.ComparableTerm;
import org.biojavax.utils.StringTools;
import org.biojavax.utils.XMLTools;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Format reader for EMBLxml files. This version of EMBLxml format will generate
 * and write RichSequence objects. Loosely Based on code from the old, deprecated,
 * org.biojava.bio.seq.io.GenbankXmlFormat object.
 *
 * Understands http://www.ebi.ac.uk/embl/dtd/EMBL_Services_V1.1.dtd
 *
 * @author Alan Li (code based on his work)
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 */
public class EMBLxmlFormat extends RichSequenceFormat.BasicFormat {
    
    // Register this format with the format auto-guesser.
    static {
        RichSequence.IOTools.registerFormat(EMBLxmlFormat.class);
    }
    
    /**
     * The name of this format
     */
    public static final String EMBLXML_FORMAT = "EMBLxml";
    
    protected static final String ENTRY_GROUP_TAG = "EMBL_Services";
    protected static final String ENTRY_TAG = "entry";
    protected static final String ENTRY_ACCESSION_ATTR = "accession";
    protected static final String ENTRY_TAX_DIVISION_ATTR = "taxonomicDivision";
    protected static final String ENTRY_DATACLASS_ATTR = "dataClass";
    protected static final String ENTRY_CREATED_ATTR = "created";
    protected static final String ENTRY_RELCREATED_ATTR = "releaseCreated";
    protected static final String ENTRY_UPDATED_ATTR = "lastUpdated";
    protected static final String ENTRY_RELUPDATED_ATTR = "releaseLastUpdated";
    protected static final String ENTRY_VER_ATTR = "version";
    protected static final String ENTRY_SUBACC_ATTR = "submitterAccession";
    protected static final String ENTRY_SUBVER_ATTR = "submitterVersion";
    protected static final String ENTRY_SUBWGSVER_ATTR = "submitterWgsVersion";
    protected static final String ENTRY_STATUS_ATTR = "status";
    protected static final String ENTRY_STATUS_DATE_ATTR = "statusDate";
    
    protected static final String SEC_ACC_TAG = "secondaryAccession";
    protected static final String PROJ_ACC_TAG = "projectAccession";
    protected static final String DESC_TAG = "description";
    protected static final String KEYWORD_TAG = "keyword";
    protected static final String REFERENCE_TAG = "reference";
    
    protected static final String CITATION_TAG = "citation";
    protected static final String CITATION_ID_ATTR = "id";
    protected static final String CITATION_TYPE_ATTR = "type";
    protected static final String CITATION_DATE_ATTR = "date";
    protected static final String CITATION_NAME_ATTR = "name";
    protected static final String CITATION_VOL_ATTR = "volume";
    protected static final String CITATION_ISSUE_ATTR = "issue";
    protected static final String CITATION_FIRST_ATTR = "first";
    protected static final String CITATION_LAST_ATTR = "last";
    protected static final String CITATION_PUB_ATTR = "publisher";
    protected static final String CITATION_PATENT_ATTR = "patentNumber";
    protected static final String CITATION_INSTITUTE_ATTR = "institute";
    protected static final String CITATION_YEAR_ATTR = "year";
    
    protected static final String DBREFERENCE_TAG = "dbreference";
    protected static final String DBREF_DB_ATTR = "db";
    protected static final String DBREF_PRIMARY_ATTR = "primary";
    protected static final String DBREF_SEC_ATTR = "secondary";
    
    protected static final String CONSORTIUM_TAG = "consortium";
    protected static final String TITLE_TAG = "title";
    protected static final String EDITOR_TAG = "editor";
    protected static final String AUTHOR_TAG = "author";
    protected static final String PATENT_TAG = "patentApplicant";
    protected static final String LOCATOR_TAG = "locator";
    
    protected static final String CITATION_LOCATION_TAG = "citationLocation";
    protected static final String REF_POS_BEGIN_ATTR = "begin";
    protected static final String REF_POS_END_ATTR = "end";
    
    protected static final String COMMENT_TAG = "comment";
    
    protected static final String FEATURE_TAG = "feature";
    protected static final String FEATURE_NAME_ATTR = "name";
    
    protected static final String ORGANISM_TAG = "organism";
    protected static final String SCINAME_TAG = "scientificName";
    protected static final String COMNAME_TAG = "preferredCommonName";
    protected static final String TAXID_TAG = "taxId";
    protected static final String LINEAGE_TAG = "lineage";
    protected static final String TAXON_TAG = "taxon";
    protected static final String ORGANELLE_TAG = "organelle";
    
    protected static final String QUALIFIER_TAG = "qualifier";
    protected static final String QUALIFIER_NAME_ATTR = "name";
    
    protected static final String LOCATION_TAG = "location";
    protected static final String LOCATION_TYPE_ATTR = "type";
    protected static final String LOCATION_COMPL_ATTR = "complement";
    
    protected static final String LOCATION_ELEMENT_TAG = "locationElement";
    protected static final String LOC_ELEMENT_TYPE_ATTR = "type";
    protected static final String LOC_ELEMENT_ACC_ATTR = "accession";
    protected static final String LOC_ELEMENT_VER_ATTR = "version";
    protected static final String LOC_ELEMENT_COMPL_ATTR = "complement";
    
    protected static final String BASEPOSITION_TAG = "basePosition";
    protected static final String BASEPOSITION_TYPE_ATTR = "type";
    
    protected static final String CONTIG_TAG = "contig";
    protected static final String SEQUENCE_TAG = "sequence";
    protected static final String SEQUENCE_TYPE_ATTR = "type";
    protected static final String SEQUENCE_LENGTH_ATTR = "length";
    protected static final String SEQUENCE_TOPOLOGY_ATTR = "topology";
    protected static final String SEQUENCE_VER_ATTR = "version";
    
    protected static final Pattern xmlSchema = Pattern.compile(".*http://www\\.ebi\\.ac\\.uk/schema/EMBL_schema\\.xsd.*");
    
    /**
     * Implements some EMBLxml-specific terms.
     */
    public static class Terms extends RichSequence.Terms {        
        /**
         * Getter for the SubmitterAccession term
         * @return The SubmitterAccession Term
         */
        public static ComparableTerm getSubmitterAccessionTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("SubmitterAccession");
        }
        
        /**
         * Getter for the SubmitterVersion term
         * @return The SubmitterVersion Term
         */
        public static ComparableTerm getSubmitterVersionTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("SubmitterVersion");
        }
        
        /**
         * Getter for the SubmitterWgsVersion term
         * @return The SubmitterWgsVersion Term
         */
        public static ComparableTerm getSubmitterWgsVersionTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("SubmitterWgsVersion");
        }
        
        /**
         * Getter for the Status term
         * @return The Status Term
         */
        public static ComparableTerm getStatusTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("Status");
        }
        
        /**
         * Getter for the StatusDate term
         * @return The StatusDate Term
         */
        public static ComparableTerm getStatusDateTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("StatusDate");
        }
        
        /**
         * Getter for the ProjectAccession term
         * @return The ProjectAccession Term
         */
        public static ComparableTerm getProjectAccessionTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("ProjectAccession");
        }
        
        /**
         * Getter for the EMBLxml term
         * @return The EMBLxml Term
         */
        public static ComparableTerm getEMBLxmlTerm() {
            return RichObjectFactory.getDefaultOntology().getOrCreateTerm("EMBLxml");
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
     * A file is in EMBLxml format if the second XML line contains the phrase "http://www.ebi.ac.uk/schema/EMBL_schema.xsd".
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
     * Always returns a DNA tokenizer.
     */
    @Override
    public SymbolTokenization guessSymbolTokenization(File file) throws IOException {
        return RichSequence.IOTools.getDNAParser();
    }
    
    /**
     * {@inheritDoc}
     * A stream is in EMBLxml format if the second XML line contains the phrase "http://www.ebi.ac.uk/schema/EMBL_schema.xsd".
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
        
        try {
            DefaultHandler m_handler = new EMBLxmlHandler(this,symParser,rlistener,ns);
            return XMLTools.readXMLChunk(reader, m_handler, ENTRY_TAG);
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
        xml.attribute("xmlns:ebi", "http://www.ebi.ac.uk/embl/schema");
        xml.attribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
        xml.attribute("xsi:noNamespaceSchemaLocation","http://www.ebi.ac.uk/embl/schema/EMBL_Services_V1.1.xsd");
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
     * Namespace is ignored as EMBLxml has no concept of it.
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
                
        Set notes = rs.getNoteSet();
        List accessions = new ArrayList();
        List projAccessions = new ArrayList();
        List kws = new ArrayList();
        List organelles = new ArrayList();
        String cdat = null;
        String udat = null;
        String crel = null;
        String urel = null;
        String dataClass = null;
        String moltype = rs.getAlphabet().getName();
        String subWgsVer = null;
        String subAcc = null;
        String subVer = null;
        String status = null;
        String statusDate = null;
        for (Iterator i = notes.iterator(); i.hasNext();) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getDateCreatedTerm())) cdat=n.getValue();
            else if (n.getTerm().equals(Terms.getDateUpdatedTerm())) udat=n.getValue();
            else if (n.getTerm().equals(Terms.getRelCreatedTerm())) crel=n.getValue();
            else if (n.getTerm().equals(Terms.getRelUpdatedTerm())) urel=n.getValue();
            else if (n.getTerm().equals(Terms.getMolTypeTerm())) moltype=n.getValue();
            else if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) accessions.add(n.getValue());
            else if (n.getTerm().equals(Terms.getProjectAccessionTerm())) projAccessions.add(n.getValue());
            else if (n.getTerm().equals(Terms.getOrganelleTerm())) organelles.add(n.getValue());
            else if (n.getTerm().equals(Terms.getKeywordTerm())) kws.add(n.getValue());
            else if (n.getTerm().equals(Terms.getDataClassTerm())) dataClass = n.getValue();
            else if (n.getTerm().equals(Terms.getSubmitterAccessionTerm())) subAcc = n.getValue();
            else if (n.getTerm().equals(Terms.getSubmitterVersionTerm())) subVer = n.getValue();
            else if (n.getTerm().equals(Terms.getSubmitterWgsVersionTerm())) subWgsVer = n.getValue();
            else if (n.getTerm().equals(Terms.getStatusTerm())) status = n.getValue();
            else if (n.getTerm().equals(Terms.getStatusDateTerm())) statusDate = n.getValue();
        }
        
        xml.openTag(ENTRY_TAG);
        xml.attribute(ENTRY_ACCESSION_ATTR,rs.getAccession());
        xml.attribute(ENTRY_TAX_DIVISION_ATTR,rs.getDivision());
        xml.attribute(ENTRY_DATACLASS_ATTR,dataClass);
        xml.attribute(ENTRY_CREATED_ATTR,cdat==null?udat:cdat);
        xml.attribute(ENTRY_RELCREATED_ATTR,crel==null?"0":crel);
        xml.attribute(ENTRY_UPDATED_ATTR,udat);
        xml.attribute(ENTRY_RELUPDATED_ATTR,urel==null?"0":urel);
        xml.attribute(ENTRY_VER_ATTR,""+rs.getVersion());
        if (subAcc!=null)
            xml.attribute(ENTRY_SUBACC_ATTR,subAcc);
        if (subVer!=null)
            xml.attribute(ENTRY_SUBVER_ATTR,subVer);
        if (subWgsVer!=null)
            xml.attribute(ENTRY_SUBWGSVER_ATTR,subWgsVer);
        if (status!=null)
            xml.attribute(ENTRY_STATUS_ATTR,status);
        if (statusDate!=null)
            xml.attribute(ENTRY_STATUS_DATE_ATTR,statusDate);
        
        for (Iterator i = accessions.iterator(); i.hasNext(); ) {
            xml.openTag(SEC_ACC_TAG);
            xml.print((String)i.next());
            xml.closeTag(SEC_ACC_TAG);
        }
        
        for (Iterator i = projAccessions.iterator(); i.hasNext(); ) {
            xml.openTag(PROJ_ACC_TAG);
            xml.print((String)i.next());
            xml.closeTag(PROJ_ACC_TAG);
        }
        
        xml.openTag(DESC_TAG);
        xml.print(rs.getDescription());
        xml.closeTag(DESC_TAG);
        
        for (Iterator i = kws.iterator(); i.hasNext(); ) {
            xml.openTag(KEYWORD_TAG);
            xml.print((String)i.next());
            xml.closeTag(KEYWORD_TAG);
        }
        
        for (Iterator i = rs.getRankedDocRefs().iterator(); i.hasNext(); ) {
            RankedDocRef rdr = (RankedDocRef)i.next();
            DocRef dr = rdr.getDocumentReference();
            
            xml.openTag(REFERENCE_TAG);
            
            xml.openTag(CITATION_TAG);
            xml.attribute(CITATION_ID_ATTR,""+rdr.getRank());
            xml.attribute(CITATION_TYPE_ATTR,"journal article");
            
            CrossRef cr = dr.getCrossref();
            if (cr!=null) {
                xml.openTag(DBREFERENCE_TAG);
                xml.attribute(DBREF_DB_ATTR,cr.getDbname());
                xml.attribute(DBREF_PRIMARY_ATTR,cr.getAccession());
                if (!cr.getNoteSet().isEmpty()) {
                    for (Iterator j = cr.getNoteSet().iterator(); j.hasNext(); ) {
                        Note n = (Note)j.next();
                        if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) {
                            xml.attribute(DBREF_SEC_ATTR,n.getValue());
                            break;
                        }
                    }
                }
                xml.closeTag(DBREFERENCE_TAG);
            }
            
            List auths = dr.getAuthorList();
            
            for (Iterator j = auths.iterator(); j.hasNext(); ) {
                DocRefAuthor a = (DocRefAuthor)j.next();
                if (a.isConsortium()) {
                    xml.openTag(CONSORTIUM_TAG);
                    xml.print(a.getName());
                    xml.closeTag(CONSORTIUM_TAG);
                    j.remove();
                }
            }
            
            if (dr.getTitle()!=null) {
                xml.openTag(TITLE_TAG);
                xml.print(dr.getTitle());
                xml.closeTag(TITLE_TAG);
            }
            
            for (Iterator j = auths.iterator(); j.hasNext(); ) {
                DocRefAuthor a = (DocRefAuthor)j.next();
                if (a.isEditor()) {
                    xml.openTag(EDITOR_TAG);
                    xml.print(a.getName());
                    xml.closeTag(EDITOR_TAG);
                } else {
                    xml.openTag(AUTHOR_TAG);
                    xml.print(a.getName());
                    xml.closeTag(AUTHOR_TAG);
                }
            }
            
            xml.openTag(LOCATOR_TAG);
            xml.print(dr.getLocation());
            xml.closeTag(LOCATOR_TAG);
            xml.closeTag(CITATION_TAG);
            
            xml.openTag(CITATION_LOCATION_TAG);     
            Integer rstart = rdr.getStart();
            if (rstart==null) rstart = new Integer(1);
            Integer rend = rdr.getEnd();
            if (rend==null) rend = new Integer(rs.length());
            xml.attribute(REF_POS_BEGIN_ATTR,""+rstart);
            xml.attribute(REF_POS_END_ATTR,""+rend);
            if (dr.getRemark()!=null) {
                xml.openTag(COMMENT_TAG);
                xml.print(dr.getRemark());
                xml.closeTag(COMMENT_TAG);
            }
            xml.closeTag(CITATION_LOCATION_TAG);
            
            xml.closeTag(REFERENCE_TAG);
        }
        
        for (Iterator i = rs.getRankedCrossRefs().iterator(); i.hasNext(); ) {
            RankedCrossRef rcr = (RankedCrossRef)i.next();
            CrossRef cr = rcr.getCrossRef();
            
            xml.openTag(DBREFERENCE_TAG);
            xml.attribute(DBREF_DB_ATTR,cr.getDbname());
            xml.attribute(DBREF_PRIMARY_ATTR,cr.getAccession());
            
            if (!cr.getNoteSet().isEmpty()) {
                for (Iterator j = cr.getNoteSet().iterator(); j.hasNext(); ) {
                    Note n = (Note)j.next();
                    if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) {
                        xml.attribute(DBREF_SEC_ATTR,n.getValue());
                        break;
                    }
                }
            }
            
            xml.closeTag(DBREFERENCE_TAG);
        }
        
        for (Iterator i = rs.getComments().iterator(); i.hasNext(); ) {
            xml.openTag(COMMENT_TAG);
            xml.println(((Comment)i.next()).getComment());
            xml.closeTag(COMMENT_TAG);
        }
        
        NCBITaxon tax = rs.getTaxon();
        for (Iterator i = rs.getFeatureSet().iterator(); i.hasNext(); ) {
            RichFeature f = (RichFeature)i.next();
            xml.openTag(FEATURE_TAG);
            xml.attribute(FEATURE_NAME_ATTR,f.getTypeTerm().getName());
            
            // display organism on source feature only
            if (f.getTypeTerm().getName().equals("source") && tax!=null) {
                xml.openTag(ORGANISM_TAG);
                
                String[] parts = tax.getDisplayName().split("(\\(|\\))");
                xml.openTag(SCINAME_TAG);
                xml.print(parts[0].trim());
                xml.closeTag(SCINAME_TAG);
                if (parts.length>1) {
                    xml.openTag(COMNAME_TAG);
                    xml.print(parts[1].trim());
                    xml.closeTag(COMNAME_TAG);
                }
                
                xml.openTag(TAXID_TAG);
                xml.print(""+tax.getNCBITaxID());
                xml.closeTag(TAXID_TAG);
                
                String hierarchy = tax.getNameHierarchy();
                hierarchy = hierarchy.substring(0,hierarchy.length()-1); // chomp "."
                if (hierarchy.length()>0) {
                    parts = hierarchy.split(";");
                    xml.openTag(LINEAGE_TAG);
                    for (int j = 0; j < parts.length; j++) {
                        xml.openTag(TAXON_TAG);
                        xml.print(parts[j].trim());
                        xml.closeTag(TAXON_TAG);
                    }
                    xml.closeTag(LINEAGE_TAG);
                }
                
                for (final Iterator j = organelles.iterator(); j.hasNext(); ) {
                    final String organelle = (String)j.next();
                    xml.openTag(ORGANELLE_TAG);
                    xml.print(organelle);
                    xml.closeTag(ORGANELLE_TAG);
                }
                
                xml.closeTag(ORGANISM_TAG);
            }
            
            for (Iterator j = f.getRankedCrossRefs().iterator(); j.hasNext(); ) {
                RankedCrossRef rcr = (RankedCrossRef)j.next();
                CrossRef cr = rcr.getCrossRef();
                
                xml.openTag(DBREFERENCE_TAG);
                xml.attribute(DBREF_DB_ATTR,cr.getDbname());
                xml.attribute(DBREF_PRIMARY_ATTR,cr.getAccession());
                
                if (!cr.getNoteSet().isEmpty()) {
                    for (Iterator k = cr.getNoteSet().iterator(); k.hasNext(); ) {
                        Note n = (Note)k.next();
                        if (n.getTerm().equals(Terms.getAdditionalAccessionTerm())) {
                            xml.attribute(DBREF_SEC_ATTR,n.getValue());
                            break;
                        }
                    }
                }
                
                xml.closeTag(DBREFERENCE_TAG);
            }
            
            for (Iterator j = f.getNoteSet().iterator(); j.hasNext();) {
                Note n = (Note)j.next();
                xml.openTag(QUALIFIER_TAG);
                xml.attribute(QUALIFIER_NAME_ATTR,n.getTerm().getName());
                if (n.getValue()!=null && !n.getValue().equals("")) {
                	if (n.getTerm().getName().equalsIgnoreCase("translation")) {
                		String[] lines = StringTools.wordWrap(n.getValue(), "\\s+", this.getLineWidth());
                		for (int k = 0; k < lines.length; k++) xml.println(lines[k]);
                	} else {	
                		xml.print(n.getValue());
                	}
                }	
                xml.closeTag(QUALIFIER_TAG);
            }
            
            // make it easy for ourselves by flattening into a single compound location
            RichLocation rle = (RichLocation)f.getLocation();
            Collection locElements = RichLocation.Tools.flatten(rle);
            xml.openTag(LOCATION_TAG);
            xml.attribute(LOCATION_TYPE_ATTR,(locElements.size()>1?rle.getTerm().getName():"single"));
            xml.attribute(LOCATION_COMPL_ATTR,"false");
            for (Iterator j = locElements.iterator(); j.hasNext(); ) {
                RichLocation rl = (RichLocation)j.next();
                xml.openTag(LOCATION_ELEMENT_TAG);
                
                if (rl.getStrand().equals(Strand.NEGATIVE_STRAND)) {
                    xml.attribute(LOC_ELEMENT_COMPL_ATTR,"true");
                } else {
                    xml.attribute(LOC_ELEMENT_COMPL_ATTR,"false");
                }
                
                if (rl.getCrossRef()!=null) {
                    xml.attribute(LOC_ELEMENT_ACC_ATTR,rl.getCrossRef().getAccession());
                    xml.attribute(LOC_ELEMENT_VER_ATTR,""+rl.getCrossRef().getVersion());
                }
                
                Position start = rl.getMinPosition();
                // EMBLxml does not support fuzzy locations so we only ever
                // use the start coordinate.
                
                // output first base only
                xml.attribute(LOC_ELEMENT_TYPE_ATTR,"site");
                    
                xml.openTag(BASEPOSITION_TAG);
                if (start.getFuzzyStart()) xml.attribute(BASEPOSITION_TYPE_ATTR,"<");
                else if (start.getFuzzyEnd()) xml.attribute(BASEPOSITION_TYPE_ATTR,"<");
                else xml.attribute(BASEPOSITION_TYPE_ATTR,"simple");
                xml.print(""+start.getStart());
                xml.closeTag(BASEPOSITION_TAG);
                
                xml.closeTag(LOCATION_ELEMENT_TAG);
            }

            xml.closeTag(LOCATION_TAG);
            
            xml.closeTag(FEATURE_TAG);
        }
        
        xml.openTag(SEQUENCE_TAG);
        xml.attribute(SEQUENCE_TYPE_ATTR,moltype);
        xml.attribute(SEQUENCE_LENGTH_ATTR,""+rs.length());
        xml.attribute(SEQUENCE_TOPOLOGY_ATTR,rs.getCircular()?"circular":"linear");
        xml.attribute(SEQUENCE_VER_ATTR,""+rs.getSeqVersion().intValue());
        String[] lines = StringTools.wordWrap(rs.seqString(), "\\s+", this.getLineWidth());
        for (int i = 0; i < lines.length; i ++) xml.println(lines[i]);
        xml.closeTag(SEQUENCE_TAG);
        
        xml.closeTag(ENTRY_TAG);
        
        pw.flush();
    }
    
    /**
     * {@inheritDoc}
     */
    public String getDefaultFormat() {
        return EMBLXML_FORMAT;
    }
    
    // SAX event handler for parsing http://www.ebi.ac.uk/embl/Documentation/DTD/EMBL_dtd.txt
    private class EMBLxmlHandler extends DefaultHandler {
        
        private RichSequenceFormat parent;
        private SymbolTokenization symParser;
        private RichSeqIOListener rlistener;
        private Namespace ns;
        private StringBuffer m_currentString;
        
        private NCBITaxon tax;
        private String accession;
        private RichFeature.Template templ;
        private String currFeatQual;
        private String currRefLocation;
        private List currRefAuthors;
        private String currRefTitle;
        private Map currNames = new TreeMap();
        private int currRefStart;
        private int currRefEnd;
        private int currRefRank;
        private int currLocBrackets;
        private int currLocElemBrackets;
        private StringBuffer currLocStr;
        private String currBaseType;
        private boolean firstBase; // oooh err!
        private boolean firstLocationElement;
        private List currDBXrefs = new ArrayList();
        private List currComments = new ArrayList();
        private Map currQuals = new LinkedHashMap();
        
        // construct a new handler that will populate the given list of sequences
        private EMBLxmlHandler(RichSequenceFormat parent,
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
                    rlistener.startSequence();
                    if (ns==null) ns=RichObjectFactory.getDefaultNamespace();
                    rlistener.setNamespace(ns);
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(ENTRY_ACCESSION_ATTR)) {
                            accession = val;
                            rlistener.setAccession(accession);
                            rlistener.setName(accession);
                        } else if (name.equals(ENTRY_TAX_DIVISION_ATTR)) rlistener.setDivision(val);
                        else if (name.equals(ENTRY_DATACLASS_ATTR)) rlistener.addSequenceProperty(Terms.getDataClassTerm(),val);
                        else if (name.equals(ENTRY_CREATED_ATTR)) rlistener.addSequenceProperty(Terms.getDateCreatedTerm(),val);
                        else if (name.equals(ENTRY_UPDATED_ATTR)) rlistener.addSequenceProperty(Terms.getDateUpdatedTerm(),val);
                        else if (name.equals(ENTRY_RELCREATED_ATTR)) rlistener.addSequenceProperty(Terms.getRelCreatedTerm(),val);
                        else if (name.equals(ENTRY_RELUPDATED_ATTR)) rlistener.addSequenceProperty(Terms.getRelUpdatedTerm(),val);
                        else if (name.equals(ENTRY_VER_ATTR)) rlistener.setVersion(Integer.parseInt(val));
                        else if (name.equals(ENTRY_SUBACC_ATTR)) rlistener.addSequenceProperty(Terms.getSubmitterAccessionTerm(),val);
                        else if (name.equals(ENTRY_SUBVER_ATTR)) rlistener.addSequenceProperty(Terms.getSubmitterVersionTerm(),val);
                        else if (name.equals(ENTRY_SUBWGSVER_ATTR)) rlistener.addSequenceProperty(Terms.getSubmitterWgsVersionTerm(),val);
                        else if (name.equals(ENTRY_STATUS_ATTR)) rlistener.addSequenceProperty(Terms.getStatusTerm(),val);
                        else if (name.equals(ENTRY_STATUS_DATE_ATTR)) rlistener.addSequenceProperty(Terms.getStatusDateTerm(),val);
                    }
                    currNames.clear();
                    currComments.clear();
                    currDBXrefs.clear();
                } catch (ParseException e) {
                    throw new SAXException(e);
                }
            }
            
            else if (qName.equals(REFERENCE_TAG) && !this.parent.getElideReferences()) {
                currRefLocation = null;
                currRefAuthors = new ArrayList();
                currRefTitle = null;
                currRefStart = -999;
                currRefEnd = -999;
                currRefRank = 0;
                currDBXrefs.clear();
                currComments.clear();
            } else if (qName.equals(CITATION_LOCATION_TAG) && !this.parent.getElideReferences()) {
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(REF_POS_BEGIN_ATTR)) currRefStart = Integer.parseInt(val);
                    else if (name.equals(REF_POS_END_ATTR)) currRefEnd = Integer.parseInt(val);
                }
            } else if (qName.equals(CITATION_TAG) && !this.parent.getElideReferences()) {
                StringBuffer currRef = new StringBuffer();
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(CITATION_ID_ATTR)) currRefRank = Integer.parseInt(val);
                    // combine everything else into a fake reference to use if locator is a no-show
                    else if (!name.equals(CITATION_TYPE_ATTR)) {
                        if (currRef.length()>0) currRef.append(" ");
                        currRef.append(val);
                    }
                }
                currRefLocation = currRef.toString();
            }
            
            else if (qName.equals(DBREFERENCE_TAG)) {
                String db = null;
                String primary = null;
                String secondary = null;
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(DBREF_DB_ATTR)) db = val;
                    else if (name.equals(DBREF_PRIMARY_ATTR)) primary = val;
                    else if (name.equals(DBREF_SEC_ATTR)) secondary = val;
                }
                CrossRef dbx = (CrossRef)RichObjectFactory.getObject(SimpleCrossRef.class,new Object[]{db, primary, new Integer(0)});
                if (secondary!=null) {
                    Note note = new SimpleNote(Terms.getAdditionalAccessionTerm(),secondary,0);
                    try {
                        dbx.getRichAnnotation().addNote(note);
                    } catch (ChangeVetoException ce) {
                        SAXException pe = new SAXException("Could not annotate identifier terms");
                        pe.initCause(ce);
                        throw pe;
                    }
                }
                currDBXrefs.add(dbx);
            }
            
            else if (qName.equals(FEATURE_TAG) && !this.parent.getElideFeatures()) {
                templ = new RichFeature.Template();
                templ.annotation = new SimpleRichAnnotation();
                templ.sourceTerm = Terms.getEMBLxmlTerm();
                templ.featureRelationshipSet = new TreeSet();
                templ.rankedCrossRefs = new TreeSet();
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(FEATURE_NAME_ATTR)) templ.typeTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm(val);
                }
                currLocStr = new StringBuffer();
                currDBXrefs.clear();
                currQuals.clear();
            } else if (qName.equals(QUALIFIER_TAG) && !this.parent.getElideFeatures()) {
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(QUALIFIER_NAME_ATTR)) currFeatQual = val;
                }
            } else if (qName.equals(LOCATION_TAG) && !this.parent.getElideFeatures()) {
                currLocBrackets = 0;
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(LOCATION_TYPE_ATTR) && !val.equalsIgnoreCase("single")) {
                        // open a bracket just in case
                        currLocStr.append(val);
                        currLocStr.append("(");
                        currLocBrackets++;
                    } else if (name.equals(LOCATION_COMPL_ATTR) && val.equalsIgnoreCase("true")) {
                        currLocStr.append("complement");
                        currLocStr.append("(");
                        currLocBrackets++;
                    }
                }
                firstLocationElement = true;
            } else if (qName.equals(LOCATION_ELEMENT_TAG) && !this.parent.getElideFeatures()) {
                String currAcc = null;
                String currVer = null;
                if (!firstLocationElement) currLocStr.append(",");
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(LOCATION_COMPL_ATTR) && val.equalsIgnoreCase("true")) {
                        currLocStr.append("complement");
                        currLocStr.append("(");
                        currLocElemBrackets++;
                    } else if (name.equals(LOC_ELEMENT_ACC_ATTR)) currAcc = val;
                    else if (name.equals(LOC_ELEMENT_VER_ATTR)) currVer = val;
                }
                if (currAcc!=null) {
                    currLocStr.append(currAcc);
                    if (currVer!=null) {
                        currLocStr.append(".");
                        currLocStr.append(currVer);
                    }
                    currLocStr.append(":");
                }
                firstBase = true;
            } else if (qName.equals(BASEPOSITION_TAG) && !this.parent.getElideFeatures()) {
                for (int i = 0; i < attributes.getLength(); i++) {
                    String name = attributes.getQName(i);
                    String val = attributes.getValue(i);
                    if (name.equals(BASEPOSITION_TYPE_ATTR)) currBaseType = val;
                }
            }
            
            else if (qName.equals(CONTIG_TAG))  {
                String message = ParseException.newMessage(this.getClass(),accession,"not set", "Unable to handle contig assemblies just yet", qName);
                ParseException e = new ParseException(message);
                SAXException pe = new SAXException("Could not set contig properties");
                pe.initCause(e);
                throw pe;
            }
            
            else if (qName.equals(SEQUENCE_TAG)) {
                try {
                    for (int i = 0; i < attributes.getLength(); i++) {
                        String name = attributes.getQName(i);
                        String val = attributes.getValue(i);
                        if (name.equals(SEQUENCE_TYPE_ATTR)) rlistener.addSequenceProperty(Terms.getMolTypeTerm(),val);
                        else if (name.equals(SEQUENCE_VER_ATTR)) rlistener.setSeqVersion(val);
                        else if (name.equals(SEQUENCE_TOPOLOGY_ATTR) && val.equalsIgnoreCase("circular")) rlistener.setCircular(true);
                    }
                } catch (ParseException e) {
                    SAXException pe = new SAXException("Could not set sequence properties");
                    pe.initCause(e);
                    throw pe;
                }
            }
        }
        
        // process a closing tag - we will have read the text already
        @Override
        public void endElement(String uri, String localName, String qName) throws SAXException {
            String val = this.m_currentString.toString().trim();
            
            try {
                if (qName.equals(SEC_ACC_TAG)) {
                    rlistener.addSequenceProperty(Terms.getAdditionalAccessionTerm(),val);
                } else if (qName.equals(PROJ_ACC_TAG)) {
                    rlistener.addSequenceProperty(Terms.getProjectAccessionTerm(),val);
                } else if (qName.equals(ORGANELLE_TAG)) {
                    rlistener.addSequenceProperty(Terms.getOrganelleTerm(),val);
                } else if (qName.equals(DESC_TAG)) {
                    rlistener.setDescription(val);
                } else if (qName.equals(KEYWORD_TAG)) {
                    rlistener.addSequenceProperty(Terms.getKeywordTerm(), val);
                } else if (qName.equals(COMMENT_TAG)) {
                    currComments.add(val);
                }
                
                else if (qName.equals(TITLE_TAG)) {
                    currRefTitle = val;
                } else if (qName.equals(AUTHOR_TAG)) {
                    currRefAuthors.add(new SimpleDocRefAuthor(val,false,false));
                } else if (qName.equals(EDITOR_TAG)) {
                    currRefAuthors.add(new SimpleDocRefAuthor(val,false,true));
                } else if (qName.equals(CONSORTIUM_TAG)) {
                    currRefAuthors.add(new SimpleDocRefAuthor(val,true,false));
                } else if (qName.equals(LOCATOR_TAG)) {
                    currRefLocation = val;
                } else if (qName.equals(REFERENCE_TAG) && !this.parent.getElideReferences()) {
                    // do the crossrefs
                    CrossRef useForDocRef = null;
                    for (Iterator j = currDBXrefs.iterator(); j.hasNext();) {
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
                    // do the comment - will only be one, if any
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
                
                else if (qName.equals(LOCATION_TAG) && !this.parent.getElideFeatures()) {
                    while (currLocBrackets-->0) currLocStr.append(")"); // close the location groups
                    String tidyLocStr = currLocStr.toString().replaceAll("\\s+","");
                    templ.location = GenbankLocationParser.parseLocation(ns, accession, tidyLocStr);
                } else if (qName.equals(LOCATION_ELEMENT_TAG) && !this.parent.getElideFeatures()) {
                    while (currLocElemBrackets-->0) currLocStr.append(")"); // close the location groups
                    firstLocationElement = false;
                } else if (qName.equals(BASEPOSITION_TAG) && !this.parent.getElideFeatures()) {
                    if (!firstBase) currLocStr.append("..");
                    // left angle bracket, right angle bracket, simple, fuzzy
                    if (currBaseType.equals("<")) {
                        currLocStr.append("<");
                        currLocStr.append(val);
                    } else if (currBaseType.equals(">")) {
                        currLocStr.append(val);
                        currLocStr.append(">");
                    } else if (currBaseType.equalsIgnoreCase("simple")) {
                        currLocStr.append(val);
                    }
                    firstBase = false;
                } else if (qName.equals(QUALIFIER_TAG) && !this.parent.getElideFeatures()) {
                    currQuals.put(currFeatQual,val);
                } else if (qName.equals(FEATURE_TAG) && !this.parent.getElideFeatures()) {
                    // start the feature
                    rlistener.startFeature(templ);
                    // assign qualifiers
                    for (Iterator j = currQuals.keySet().iterator(); j.hasNext(); ) {
                        String qualName = (String)j.next();
                        String qualVal = (String)currQuals.get(qualName);
                        if (qualName.equalsIgnoreCase("translation")) {
                            // strip spaces from sequence
                            qualVal = qualVal.replaceAll("\\s+","");
                        }
                        rlistener.addFeatureProperty(RichObjectFactory.getDefaultOntology().getOrCreateTerm(qualName),qualVal);
                    }
                    // do the crossrefs
                    int rcrossrefCount = 0;
                    for (Iterator j = currDBXrefs.iterator(); j.hasNext();) {
                        CrossRef dbx = (CrossRef)j.next();
                        RankedCrossRef rdbx = new SimpleRankedCrossRef(dbx, ++rcrossrefCount);
                        try {
                            rlistener.getCurrentFeature().addRankedCrossRef(rdbx);
                        } catch (ChangeVetoException ce) {
                            throw new ParseException(ce);
                        }
                    }
                    // end the feature
                    rlistener.endFeature();
                    currDBXrefs.clear();
                }
                
                else if (qName.equals(TAXID_TAG)) {
                    tax = (NCBITaxon)RichObjectFactory.getObject(SimpleNCBITaxon.class, new Object[]{Integer.valueOf(val)});
                    rlistener.setTaxon(tax);
                    for (Iterator j = currNames.keySet().iterator(); j.hasNext(); ) {
                        String nameClass = (String)j.next();
                        Set nameSet = (Set)currNames.get(nameClass);
                        try {
                            for (Iterator k = nameSet.iterator(); k.hasNext(); ) {
                                String name = (String)k.next();
                                tax.addName(nameClass,name);
                            }
                        } catch (ChangeVetoException ce) {
                            throw new ParseException(ce);
                        }
                    }
                    currNames.clear();
                } else if (qName.equals(SCINAME_TAG)) {
                    try {
                        if (tax==null) {
                            if (!currNames.containsKey(NCBITaxon.SCIENTIFIC)) currNames.put(NCBITaxon.SCIENTIFIC,new TreeSet());
                            ((Set)currNames.get(NCBITaxon.SCIENTIFIC)).add(val);
                        } else {
                            tax.addName(NCBITaxon.SCIENTIFIC,val);
                        }
                    } catch (ChangeVetoException ce) {
                        throw new ParseException(ce);
                    }
                } else if (qName.equals(COMNAME_TAG)) {
                    try {
                        if (tax==null) {
                            if (!currNames.containsKey(NCBITaxon.COMMON)) currNames.put(NCBITaxon.COMMON,new TreeSet());
                            ((Set)currNames.get(NCBITaxon.COMMON)).add(val);
                        } else {
                            tax.addName(NCBITaxon.COMMON,val);
                        }
                    } catch (ChangeVetoException ce) {
                        throw new ParseException(ce);
                    }
                }
                
                else if (qName.equals(SEQUENCE_TAG) && !this.parent.getElideSymbols()) {
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
                    rlistener.endSequence();
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

