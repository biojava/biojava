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

package org.biojava.bio.seq.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.ListIterator;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.symbol.IllegalSymbolException;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Format reader for GenBank XML files.
 * 
 * @author Alan Li - alanli[at]xyworks.com
 * @deprecated Use org.biojavax.bio.seq.io.INSDseqFormat
 */

public class GenbankXmlFormat extends GenbankFormat
{
    private SAXParser m_xmlParser;
    private GenbankXmlHandler m_handler;
    private boolean m_parsed;
    private int m_sequenceIndex;
    
    public GenbankXmlFormat()
    {
        m_parsed = false;
    }
    
    public boolean readSequence( BufferedReader reader,
                                 SymbolTokenization symParser,
                                 SeqIOListener listener)
    throws IllegalSymbolException, IOException, ParseException
    {
        // only parse the sequence only once
        if( ! m_parsed )
        {
            SAXParserFactory factory = SAXParserFactory.newInstance();
            factory.setValidating( true );
            try
            {
                m_xmlParser = factory.newSAXParser();
            }
            catch( ParserConfigurationException ex )
            {
                throw new ParseException( ex );
            }
            catch( SAXException ex )
            {
                throw new ParseException( ex );
            }
            
            InputSource source = new InputSource( reader );
            m_handler = new GenbankXmlHandler();
                    
            try
            {
                m_xmlParser.parse( source, m_handler );
            }
            catch( SAXException ex )
            {
                throw new ParseException( ex );
            }
            
            m_parsed = true;
        }
        
        // after parsing the input into one or more sequences, pass only 
        // one sequence per readSequence()
        ArrayList sequences = m_handler.getSequences();     
        listener.startSequence();
        GenbankXmlSequence sequence = 
            (GenbankXmlSequence) sequences.get( m_sequenceIndex++ );
        generateSequenceForListener( sequence, listener, symParser );
        listener.endSequence();
        
        return ( m_sequenceIndex < sequences.size() );
    }
   
    // helper method to handle the genbank xml objects to the SeqIOListener
    private void generateSequenceForListener( GenbankXmlSequence sequence, SeqIOListener listener,
                                              SymbolTokenization symParser )
        throws ParseException, IllegalSymbolException
    {
        listener.addSequenceProperty( LOCUS_TAG, sequence.getLocus() );
        listener.addSequenceProperty( SIZE_TAG, sequence.getLength() );
        
        String strandedness = sequence.getStrandedness();
        if( strandedness != null )
            listener.addSequenceProperty( STRAND_NUMBER_TAG, 
                                          convertToStrandednessName( strandedness ) );
        
        String moltype = sequence.getMolType();
        if( moltype != null )
            listener.addSequenceProperty( TYPE_TAG, 
                                          convertToMolTypeName( moltype ) );
        String topology = sequence.getTopology();
        if( topology != null )
            listener.addSequenceProperty( CIRCULAR_TAG,
                                          convertToTopologyName( topology ) );
        
        listener.addSequenceProperty( DIVISION_TAG, sequence.getDivision() );
        listener.addSequenceProperty( DATE_TAG, sequence.getUpdateDate() );
        listener.addSequenceProperty( DEFINITION_TAG, sequence.getDefinition() );
        
        String primaryAccession = sequence.getPrimaryAccession();
        if( primaryAccession != null )
            listener.addSequenceProperty( ACCESSION_TAG, primaryAccession );
        
        String accessionVersion = sequence.getAccessionVersion();
        if( accessionVersion != null )
            listener.addSequenceProperty( VERSION_TAG, accessionVersion );
        
        handleOtherSequenceIds( sequence, listener );  
        handleKeywords( sequence, listener );
        
        listener.addSequenceProperty( SOURCE_TAG, sequence.getSource() );
        listener.addSequenceProperty( ORGANISM_TAG, sequence.getOrganism() + sequence.getTaxonomy() );
        
        handleReferences( sequence, listener );
        
        String comment = sequence.getComment();
        if( comment != null )
            listener.addSequenceProperty( COMMENT_TAG, comment );
        
        handleFeatures( sequence, listener );
        
        String seq = sequence.getSequence();
        if( seq != null && ! getElideSymbols() )
        {
            StreamParser streamParser = symParser.parseStream( listener );
            streamParser.characters( seq.toCharArray(), 0, seq.length() );
        }
    }

    private void handleOtherSequenceIds( GenbankXmlSequence sequence, SeqIOListener listener ) 
        throws ParseException
    {
        ArrayList otherSeqIds = sequence.getOtherSequencesIds();
        ListIterator iter = otherSeqIds.listIterator();
        while( iter.hasNext() )
        {
            String seqId = (String) iter.next();
            if( seqId.startsWith( "gi|" ) )
                listener.addSequenceProperty( GI_TAG, seqId.substring( 3 ) );
        }
    }

    private void handleKeywords( GenbankXmlSequence sequence, SeqIOListener listener ) 
        throws ParseException
    {
        ArrayList keywords = sequence.getKeywords();
        ListIterator iter = keywords.listIterator();
        while( iter.hasNext() )
        {
            String keyword = (String) iter.next();
            listener.addSequenceProperty( "KEYWORDS", keyword );
        }
    }

    private void handleReferences( GenbankXmlSequence sequence, SeqIOListener listener ) throws ParseException
    {
        ArrayList references = sequence.getReferences();
        ListIterator iter = references.listIterator();
        while( iter.hasNext() )
        {
            GenbankXmlReference reference = (GenbankXmlReference) iter.next();
            
            listener.addSequenceProperty( REFERENCE_TAG, reference.getReference() );
            
            ArrayList authors = reference.getAuthors();
            ListIterator authorIter = authors.listIterator();
            StringBuffer theAuthors = new StringBuffer( "" );
            while( authorIter.hasNext() )
            {
                String author = (String) authorIter.next();
                theAuthors.append( author );
                if( authorIter.hasNext() )
                    theAuthors.append( ", " );
            }
            listener.addSequenceProperty( AUTHORS_TAG, theAuthors.toString() );
            
            String title = reference.getTitle();
            if( title != null )
                listener.addSequenceProperty( TITLE_TAG, title );
            
            listener.addSequenceProperty( JOURNAL_TAG, reference.getJournal() );
            
            // Added by RichardH to include Medline/Pubmed refs.
            String pubmed = reference.getPubmed();
            if (pubmed!=null) listener.addSequenceProperty(PUBMED_TAG, pubmed);
            String medline = reference.getMedline();
            if (medline!=null) listener.addSequenceProperty(MEDLINE_TAG, medline);
        }
    }
    
    private void handleFeatures( GenbankXmlSequence sequence, SeqIOListener listener ) throws ParseException
    {
        ArrayList features = sequence.getFeatures();
        ListIterator iter = features.listIterator();
        while( iter.hasNext() )
        {
            GenbankXmlFeature feature = (GenbankXmlFeature) iter.next();
            String key = feature.getKey();
            int keyLength = key.length();
            String featureString = feature.getKey() + createBlankString( 16 - keyLength ) + 
                                   feature.getLocation();
            listener.addSequenceProperty( FEATURE_FLAG, featureString );
            
            ArrayList qualifiers = feature.getQualifiers();
            ListIterator qualifiersIter = qualifiers.listIterator();
            while( qualifiersIter.hasNext() )
            {
                GenbankXmlQualifier qualifier = (GenbankXmlQualifier) qualifiersIter.next();
                String qualifierString = "                /" + qualifier.getName() + "=" +
                                         "\"" + qualifier.getValue() + "\"";
                listener.addSequenceProperty( FEATURE_FLAG, qualifierString );
            }
        }
    }
    
    public String getDefaultFormat()
    {
        return "GenbankXml";
    }
 
    private String createBlankString( int size )
    {
        StringBuffer sb = new StringBuffer();
        for( int i = 0; i < size; i++ )
            sb.append( ' ' );
        return sb.toString();
    }
   
    private String convertToStrandednessName( String strandednessIndex ) throws ParseException
    {
        int i = Integer.parseInt( strandednessIndex );
        switch( i )
        {
            case 0:
                return "not-set";
            case 1:
                return "single-stranded";
            case 2:
                return "double-stranded";
            case 3:
                return "mixed-stranded";
            default:
                throw new ParseException( "Unknown strandedness: " + strandednessIndex );
        }
    }
    
    private String convertToMolTypeName( String moltypeIndex ) throws ParseException
    {
        int i = Integer.parseInt( moltypeIndex );
        switch( i )
        {
            case 0:
                return "nucleic-acid";
            case 1:
                return "dna";
            case 2:
                return "rna";
            case 3:
                return "trna";
            case 4:
                return "rrna";
            case 5:
                return "mrna";
            case 6:
                return "urna";
            case 7:
                return "snrna";
            case 8:
                return "snorna";
            case 9:
                return "peptide";
            default:
                throw new ParseException( "Unknown molecule type: " + i );
        }
    }
    
    private String convertToTopologyName( String topologyIndex ) throws ParseException
    {
        int i = Integer.parseInt( topologyIndex );
        switch( i )
        {
            case 1:
                return "linear";
            case 2:
                return "circular";
            default:
                throw new ParseException( "Unknown topology: " + i );
        }
    }
     
    // SAX event handler extended to parse the Genbank XML format
    // see details about GenbankXML format at http://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.mod
    private static class GenbankXmlHandler extends DefaultHandler
    {
        private GenbankXmlSequence m_currentSequence;
        private ArrayList m_sequences;
        private StringBuffer m_currentString;
        
        private GenbankXmlHandler()
        {
            m_sequences = new ArrayList();
            m_currentString = new StringBuffer();
        }
        
        public void startElement( String uri, String localName, String qName, 
                                  Attributes attributes )
        {
            if( qName.equals( "GBSeq" ) )
                m_currentSequence = new GenbankXmlSequence();
            else if( qName.equals( "GBReference" ) )
                m_currentSequence.addNewReference();
            else if( qName.equals( "GBFeature" ) )
                m_currentSequence.addNewFeature();
            else if( qName.equals( "GBQualifier" ) )
                m_currentSequence.getCurrentFeature().addNewQualifier();
            else if( qName.equals( "GBInterval" ) )
                m_currentSequence.getCurrentFeature().addNewInterval();
        }
        
        public void endElement( String uri, String localName, String qName ) throws SAXException
        {
            if( qName.equals( "GBSet" ) )
                return;     // do nothing
            else if( qName.equals( "GBSeq" ) )
                m_sequences.add( m_currentSequence );
            else if( qName.equals( "GBSeq_locus" ) )
                m_currentSequence.setLocus( m_currentString.toString() );
            else if( qName.equals( "GBSeq_length" ) )
                m_currentSequence.setLength( m_currentString.toString() );
            else if( qName.equals( "GBSeq_strandedness" ) )
                m_currentSequence.setStrandedness( m_currentString.toString() );
            else if( qName.equals( "GBSeq_moltype" ) )
                m_currentSequence.setMolType( m_currentString.toString() );
            else if( qName.equals( "GBSeq_topology" ) )
                m_currentSequence.setTopology( m_currentString.toString() );
            else if( qName.equals( "GBSeq_division" ) )
                m_currentSequence.setDivision( m_currentString.toString() );
            else if( qName.equals( "GBSeq_update-date" ) )
                m_currentSequence.setUpdateDate( m_currentString.toString() );
            else if( qName.equals( "GBSeq_create-date") )
                m_currentSequence.setCreateDate( m_currentString.toString() );
            else if( qName.equals( "GBSeq_update-release" ) )
                m_currentSequence.setUpdateRelease( m_currentString.toString() );
            else if( qName.equals( "GBSeq_create-release") )
                m_currentSequence.setCreateRelease( m_currentString.toString() );
            else if( qName.equals( "GBSeq_definition" ) )
                m_currentSequence.setDefinition( m_currentString.toString() );
            else if( qName.equals( "GBSeq_primary-accession" ) )
                m_currentSequence.setPrimaryAccession( m_currentString.toString() );
            else if( qName.equals( "GBSeq_entry-version" ) )
                m_currentSequence.setEntryVersion( m_currentString.toString() );
            else if( qName.equals( "GBSeq_accession-version" ) )
                m_currentSequence.setAccessionVersion( m_currentString.toString() );
            else if( qName.equals( "GBSeq_other-seqids" ) )
                return;     // nothing to do
            else if( qName.equals( "GBSeqid" ) )
                m_currentSequence.addOtherSequenceId( m_currentString.toString() );
            else if( qName.equals( "GBSeq_secondary-accessions" ) )
                return;     // nothing to do
            else if( qName.equals( "GBSecondary-accn" ) )
                m_currentSequence.addSecondaryAccession( m_currentString.toString() );
            else if( qName.equals( "GBSeq_keywords" ) )
                return;     // nothing to do
            else if( qName.equals( "GBKeyword" ) )
                m_currentSequence.addKeyword( m_currentString.toString() );
            else if( qName.equals( "GBSeq_segment" ) )
                m_currentSequence.setSegment( m_currentString.toString() );
            else if( qName.equals( "GBSeq_source" ) )
                m_currentSequence.setSource( m_currentString.toString() );
            else if( qName.equals( "GBSeq_organism" ) )
                m_currentSequence.setOrganism( m_currentString.toString() );
            else if( qName.equals( "GBSeq_taxonomy" ) )
                m_currentSequence.setTaxonomy( m_currentString.toString() );
            else if( qName.equals( "GBSeq_references" ) )
                return;     // nothing to do
            else if( qName.equals( "GBReference" ) )
                return;     // nothing to do
            else if( qName.equals( "GBReference_reference" ) )
                m_currentSequence.getCurrentReference().setReference( m_currentString.toString() );
            else if( qName.equals( "GBReference_authors" ) )
                return;     //nothing to do
            else if( qName.equals( "GBAuthor" ) )
                m_currentSequence.getCurrentReference().addAuthor( m_currentString.toString() );
            else if( qName.equals( "GBReference_consortium" ) )
                m_currentSequence.getCurrentReference().setConsortium( m_currentString.toString() );
            else if( qName.equals( "GBReference_title" ) )
                m_currentSequence.getCurrentReference().setTitle( m_currentString.toString() );
            else if( qName.equals( "GBReference_journal" ) )
                m_currentSequence.getCurrentReference().setJournal( m_currentString.toString() );
            else if( qName.equals( "GBReference_medline" ) )
                m_currentSequence.getCurrentReference().setMedline( m_currentString.toString() );
            else if( qName.equals( "GBReference_pubmed" ) )
                m_currentSequence.getCurrentReference().setPubmed( m_currentString.toString() );
            else if( qName.equals( "GBReference_remark" ) )
                m_currentSequence.getCurrentReference().setRemark( m_currentString.toString() );
            else if( qName.equals( "GBSeq_comment" ) )
                m_currentSequence.setComment( m_currentString.toString() );
            else if( qName.equals( "GBSeq_primary" ) )
                m_currentSequence.setPrimary( m_currentString.toString() );
            else if( qName.equals( "GBSeq_source-db" ) )
                m_currentSequence.setSourceDb( m_currentString.toString() );
            else if( qName.equals( "GBSeq_database-reference" ) )
                m_currentSequence.setDatabaseReference( m_currentString.toString() );
            else if( qName.equals( "GBSeq_feature-table" ) )
                return;     // do nothing
            else if( qName.equals( "GBFeature" ) )
                return;     // do nothing
            else if( qName.equals( "GBFeature_key" ) )
                m_currentSequence.getCurrentFeature().setKey( m_currentString.toString() );
            else if( qName.equals( "GBFeature_location" ) )
                m_currentSequence.getCurrentFeature().setLocation( m_currentString.toString() );
            else if( qName.equals( "GBFeature_intervals" ) )
                return;     // do nothing
            else if( qName.equals( "GBInterval" ) )
                return;     // do nothing
            else if( qName.equals( "GBInterval_from" ) )
                m_currentSequence.getCurrentFeature().getCurrentInterval().setFrom( m_currentString.toString() );
            else if( qName.equals( "GBInterval_to" ) )
                m_currentSequence.getCurrentFeature().getCurrentInterval().setTo( m_currentString.toString() );
            else if( qName.equals( "GBInterval_point" ) )
                m_currentSequence.getCurrentFeature().getCurrentInterval().setPoint( m_currentString.toString() );
            else if( qName.equals( "GBInterval_accession" ) )
                m_currentSequence.getCurrentFeature().getCurrentInterval().setAccession( m_currentString.toString() );
            else if( qName.equals( "GBFeature_quals" ) )
                return;     // do nothing
            else if( qName.equals( "GBQualifier" ) )
                return;     // do nothing
            else if( qName.equals( "GBQualifier_name" ) )
                m_currentSequence.getCurrentFeature().getCurrentQualifier().setName( m_currentString.toString() );
            else if( qName.equals( "GBQualifier_value" ) )
                m_currentSequence.getCurrentFeature().getCurrentQualifier().setValue( m_currentString.toString() );
            else if( qName.equals( "GBSeq_sequence" ) )
                m_currentSequence.setSequence( m_currentString.toString() );
            else if( qName.equals( "GBSeq_contig" ) )
                m_currentSequence.setContig( m_currentString.toString() );
            else
                throw new SAXException( "Unrecognized tag: " + qName );
            
            m_currentString.delete( 0, m_currentString.length() );
        }
        
        public void characters( char[] ch, int start, int length )
        {
            m_currentString.append( ch, start, length );
        }
        
        private ArrayList getSequences()
        {
            return m_sequences;
        }
    }
    
    // see GBSeq at http://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.mod
    private static class GenbankXmlSequence
    {
        private String m_locus;
        private String m_length;
        private String m_strandedness;
        private String m_moltype;
        private String m_topology;
        private String m_division;
        private String m_updateDate;
        private String m_createDate;
        private String m_updateRelease;
        private String m_createRelease;
        private String m_definition;
        private String m_primaryAccession;
        private String m_entryVersion;
        private String m_accessionVersion;
        private ArrayList m_otherSeqIds;
        private ArrayList m_secondaryAccessions;
        private ArrayList m_keywords;
        private String m_segment;
        private String m_source;
        private String m_organism;
        private String m_taxonomy;
        private ArrayList m_references;
        private GenbankXmlReference m_currReference;
        private String m_comment;
        private String m_primary;
        private String m_sourceDb;
        private String m_databaseReference;
        private ArrayList m_features;
        private GenbankXmlFeature m_currFeature;
        private String m_sequence;
        private String m_contig;
        
        private GenbankXmlSequence()
        {
            m_otherSeqIds = new ArrayList();
            m_secondaryAccessions = new ArrayList();
            m_keywords = new ArrayList();
            m_references = new ArrayList();
            m_features = new ArrayList();
        }
        
        private void setLocus( String locus )
        {
            m_locus = locus;
        }
        
        private String getLocus()
        {
            return m_locus;
        }
        
        private void setLength( String length )
        {
            m_length = length;
        }
        
        private String getLength()
        {
            return m_length;
        }
        
        private void setStrandedness( String strandedness )
        {
            m_strandedness = strandedness;
        }
        
        private String getStrandedness()
        {
            return m_strandedness;
        }
        
        private void setMolType( String moltype )
        {
            m_moltype = moltype;
        }
        
        private String getMolType()
        {
            return m_moltype;
        }
        
        private void setTopology( String topology )
        {
            m_topology = topology;
        }
        
        private String getTopology()
        {
            return m_topology;
        }
        
        private void setDivision( String division )
        {
            m_division = division;
        }
        
        private String getDivision()
        {
            return m_division;
        }
        
        private void setUpdateDate( String updateDate )
        {
            m_updateDate = updateDate;
        }
        
        private String getUpdateDate()
        {
            return m_updateDate;
        }
        
        private void setCreateDate( String createDate )
        {
            m_createDate = createDate;
        }
        
        String getCreateDate()
        {
            return m_createDate;
        }
        
        private void setUpdateRelease( String updateRelease )
        {
            m_updateRelease = updateRelease;
        }
        
        String getUpdateRelease()
        {
            return m_updateRelease;
        }
        
        private void setCreateRelease( String createRelease )
        {
            m_createRelease = createRelease;
        }
        
        String getCreateRelease()
        {
            return m_createRelease;
        }
        
        private void setDefinition( String definition )
        {
            m_definition = definition;
        }
        
        private String getDefinition()
        {
            return m_definition;
        }
        
        private void setPrimaryAccession( String primaryAccession )
        {
            m_primaryAccession = primaryAccession;
        }
        
        private String getPrimaryAccession()
        {
            return m_primaryAccession;
        }
        
        private void setEntryVersion( String entryVersion )
        {
            m_entryVersion = entryVersion;
        }
        
        String getEntryVersion()
        {
            return m_entryVersion;
        }
        
        private void setAccessionVersion( String accessionVersion )
        {
            m_accessionVersion = accessionVersion;
        }
        
        private String getAccessionVersion()
        {
            return m_accessionVersion;
        }
        
        private void addOtherSequenceId( String seqId )
        {
            m_otherSeqIds.add( seqId );
        }
        
        private ArrayList getOtherSequencesIds()
        {
            return new ArrayList( m_otherSeqIds );
        }
        
        private void addSecondaryAccession( String secondaryAccession )
        {
            m_secondaryAccessions.add( secondaryAccession );
        }
        
        ArrayList getSecondaryAccessions()
        {
            return new ArrayList( m_secondaryAccessions );
        }
        
        private void addKeyword( String keyword )
        {
            m_keywords.add( keyword );
        }
        
        private ArrayList getKeywords()
        {
            return new ArrayList( m_keywords );
        }
        
        private void setSegment( String segment )
        {
            m_segment = segment;
        }
        
        String getSegment()
        {
            return m_segment;
        }
        
        private void setSource( String source )
        {
            m_source = source;
        }
        
        private String getSource()
        {
            return m_source;
        }
        
        private void setOrganism( String organism )
        {
            m_organism = organism;
        }
        
        private String getOrganism()
        {
            return m_organism;
        }
        
        private void setTaxonomy( String taxonomy )
        {
            m_taxonomy = taxonomy;
        }
        
        private String getTaxonomy()
        {
            return m_taxonomy;
        }
        
        private void addNewReference()
        {
            m_currReference = new GenbankXmlReference();
            m_references.add( m_currReference );
        }
        
        private GenbankXmlReference getCurrentReference()
        {
            return m_currReference;
        }
        
        private ArrayList getReferences()
        {
            return new ArrayList( m_references );
        }
        
        private void setComment( String comment )
        {
            m_comment = comment;
        }
        
        private String getComment()
        {
            return m_comment;
        }
        
        private void setPrimary( String primary )
        {
            m_primary = primary;
        }
        
        String getPrimary()
        {
            return m_primary;
        }
        
        private void setSourceDb( String sourceDb )
        {
            m_sourceDb = sourceDb;
        }
        
        String getSourceDb()
        {
            return m_sourceDb;
        }
        
        private void setDatabaseReference( String databaseReference )
        {
            m_databaseReference = databaseReference;
        }
        
        String getDatabaseReference()
        {
            return m_databaseReference;
        }
        
        private void addNewFeature()
        {
            m_currFeature = new GenbankXmlFeature();
            m_features.add( m_currFeature );
        }
        
        private GenbankXmlFeature getCurrentFeature()
        {
            return m_currFeature;
        }
        
        private ArrayList getFeatures()
        {
            return new ArrayList( m_features );
        }
        
        private void setSequence( String sequence )
        {
            m_sequence = sequence;
        }
        
        private String getSequence()
        {
            return m_sequence;
        }
        
        private void setContig( String contig )
        {
            m_contig = contig;
        }
        
        String getContig()
        {
            return m_contig;
        }
    }
    
    // see GBReference at http://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.mod
    private static class GenbankXmlReference
    {
        private String m_reference;
        private ArrayList m_authors;
        private String m_consortium;
        private String m_title;
        private String m_journal;
        private String m_medline;
        private String m_pubmed;
        private String m_remark;
        
        private GenbankXmlReference()
        {
            m_authors = new ArrayList();
        }
        
        private void setReference( String reference )
        {
            m_reference = reference;
        }
        
        private String getReference()
        {
            return m_reference;
        }
        
        private void addAuthor( String author )
        {
            m_authors.add( author );
        }
        
        private ArrayList getAuthors()
        {
            return new ArrayList( m_authors );
        }
        
        private void setConsortium( String consortium )
        {
            m_consortium = consortium;
        }
        
        String getConsortium()
        {
            return m_consortium;
        }
        
        private void setTitle( String title )
        {
            m_title = title;
        }
        
        private String getTitle()
        {
            return m_title;
        }
        
        private void setJournal( String journal )
        {
            m_journal = journal;
        }
        
        private String getJournal()
        {
            return m_journal;
        }
        
        private void setMedline( String medline )
        {
            m_medline = medline;
        }
        
        private String getMedline()
        {
            return m_medline;
        }
        
        private void setPubmed( String pubmed )
        {
            m_pubmed = pubmed;
        }
        
        private String getPubmed()
        {
            return m_pubmed;
        }
        
        private void setRemark( String remark )
        {
            m_remark = remark;
        }
        
        String getRemark()
        {
            return m_remark;
        }
    }
    
    // see GBFeature at http://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.mod
    private static class GenbankXmlFeature
    {
        private String m_key;
        private String m_location;
        private ArrayList m_intervals;
        private GenbankXmlInterval m_currInterval;
        private ArrayList m_qualifiers;
        private GenbankXmlQualifier m_currQualifier;
        
        private GenbankXmlFeature()
        {
            m_intervals = new ArrayList();
            m_qualifiers = new ArrayList();
        }
        
        private void setKey( String key )
        {
            m_key = key;
        }
        
        private String getKey()
        {
            return m_key;
        }
        
        private void setLocation( String location )
        {
            m_location = location;
        }
        
        private String getLocation()
        {
            return m_location;
        }
        
        private void addNewInterval()
        {
            m_currInterval = new GenbankXmlInterval();
            m_intervals.add( m_currInterval );
        }
        
        private GenbankXmlInterval getCurrentInterval()
        {
            return m_currInterval;
        }
        
        ArrayList getIntervals()
        {
            return new ArrayList( m_intervals );
        }
        
        private void addNewQualifier()
        {
            m_currQualifier = new GenbankXmlQualifier();
            m_qualifiers.add( m_currQualifier );
        }
        
        private GenbankXmlQualifier getCurrentQualifier()
        {
            return m_currQualifier;
        }
        
        private ArrayList getQualifiers()
        {
            return new ArrayList( m_qualifiers );
        }
    }
    
    // see GBQualifier at http://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.mod
    private static class GenbankXmlQualifier
    {
        private String m_name;
        private String m_value;
        
        private void setName( String name )
        {
            m_name = name;
        }
        
        private String getName()
        {
            return m_name;
        }
        
        private void setValue( String value )
        {
            m_value = value;
        }
        
        private String getValue()
        {
            return m_value;
        }
    }
    
    // see GBInterval at http://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.mod
    private static class GenbankXmlInterval
    {
        private String m_from;
        private String m_to;
        private String m_point;
        private String m_accession;
        
        private void setFrom( String from )
        {
            m_from = from;
        }
        
        String getFrom()
        {
            return m_from;
        }
        
        private void setTo( String to )
        {
            m_to = to;
        }
        
        String getTo()
        {
            return m_to;
        }
        
        private void setPoint( String point )
        {
            m_point = point;
        }
        
        String getPoint()
        {
            return m_point;
        }
        
        private void setAccession( String accession )
        {
            m_accession = accession;
        }
        
        String getAccession()
        {
            return m_accession;
        }
    }
}