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

package org.biojavax.bio.seq.io;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.io.ChunkedSymbolListFactory;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SimpleSymbolListFactory;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListFactory;
import org.biojava.ontology.InvalidTermException;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.Comment;
import org.biojavax.Namespace;
import org.biojavax.Note;
import org.biojavax.RankedCrossRef;
import org.biojavax.RankedDocRef;
import org.biojavax.RichAnnotation;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleComment;
import org.biojavax.SimpleNote;
import org.biojavax.SimpleRichAnnotation;
import org.biojavax.bio.BioEntryRelationship;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.SimpleRichFeature;
import org.biojavax.bio.seq.SimpleRichFeatureRelationship;
import org.biojavax.bio.seq.SimpleRichSequence;
import org.biojavax.bio.taxa.NCBITaxon;
import org.biojavax.ontology.ComparableTerm;


/**
 * Constructs BioEntry objects by listening to events.
 * @author Richard Holland
 * @author George Waldon
 * @since 1.5
 */
public class SimpleRichSequenceBuilder extends RichSeqIOAdapter implements RichSequenceBuilder {
    
    private RichAnnotation notes = new SimpleRichAnnotation();
    
    /**
     * Creates a new instance of SimpleRichSequenceBuilder using a SimpleSymbolListFactory
     * with a threshold of zero.
     */
    public SimpleRichSequenceBuilder() {
        this(new SimpleSymbolListFactory(),0);
    }
        
    /**
     * Creates a new instance of SimpleRichSequenceBuilder with the
     * desired symbollistfactory and threshold of zero.
     * @param factory the symbollistfactory to use from the start.
     */
    public SimpleRichSequenceBuilder(SymbolListFactory factory) {
        this(factory,0);
    }
    
    /**
     * Creates a new instance of SimpleRichSequenceBuilder with the
     * desired symbollistfactory and threshold.
     * @param factory the symbollistfactory to use.
     * @param threshold the threshold at which the specified symbollistfactory
     * should come into use. If <=0, it will be used from the start.
     */
    public SimpleRichSequenceBuilder(SymbolListFactory factory, int threshold) {
        this.reset();
        this.factory = factory;
        this.threshold = threshold;
    }
    
    /**
     * Sets the sequence info back to default values, ie. in order to start
     * constructing a new sequence from scratch.
     */
    private void reset() {
        try{
            this.version = 0;
            this.versionSeen = false;
            this.seqVersion = 0.0;
            this.seqVersionSeen = false;
            this.accession = null;
            this.description = null;
            this.division = null;
            this.identifier = null;
            this.name = null;
            this.crossRefs.clear();
            this.symbols = null;
            this.namespace = null;
            this.taxon = null;
            this.seqPropCount = 1;   //annotation rank
            this.referenceCount = 1; //doc reference rank
            this.commentRank = 1;    //comment rank
            this.featureRank = 1;    //feature rank
            this.featPropCount = 1;  //feature annotation rank
            this.comments.clear();
            this.relations.clear();
            this.references.clear();
            this.rootFeatures.clear();
            this.featureStack.clear();
            this.allFeatures.clear();
            this.notes.clear();
        }catch(ChangeVetoException ex){
            throw new BioError("A ChangeListener should not have been applied", ex);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public void setVersion(int version) throws ParseException {
        if (this.versionSeen) throw new ParseException("Current BioEntry already has a version");
        else {
            try {
                this.version = version;
                this.versionSeen = true;
            } catch (NumberFormatException e) {
                throw new ParseException("Could not parse version as an integer");
            }
        }
    }
    private int version;
    private boolean versionSeen;
    
    /**
     * {@inheritDoc}
     * NOT IMPLEMENTED
     */
    public void setURI(String uri) throws ParseException {
        throw new ParseException("We don't understand URIs");
    }
    
    /**
     * {@inheritDoc}
     */
    public void setSeqVersion(String seqVersion) throws ParseException {
        if (this.seqVersionSeen) throw new ParseException("Current BioEntry already has a sequence version");
        if (seqVersion==null) this.seqVersion=0.0;
        else {
            try {
                this.seqVersion = Double.parseDouble(seqVersion);
                this.seqVersionSeen = true;
            } catch (NumberFormatException e) {
                throw new ParseException("Could not parse sequence version as a double");
            }
        }
    }
    private double seqVersion = 0.0;
    private boolean seqVersionSeen;
    
    /**
     * {@inheritDoc}
     * The last accession passed to this routine will always be the one used.
     */
    public void setAccession(String accession) throws ParseException {
        if (accession==null) throw new ParseException("Accession cannot be null");
        this.accession = accession;
    }
    private String accession;
    
    /**
     * {@inheritDoc}
     */
    public void setDescription(String description) throws ParseException {
        if (this.description!=null) throw new ParseException("Current BioEntry already has a description");
        this.description = description;
    }
    private String description;
    
    /**
     * {@inheritDoc}
     */
    public void setDivision(String division) throws ParseException {
        if (division==null) throw new ParseException("Division cannot be null");
        if (this.division!=null) throw new ParseException("Current BioEntry already has a division");
        this.division = division;
    }
    private String division;
    
    /**
     * {@inheritDoc}
     */
    public void setIdentifier(String identifier) throws ParseException {
        if (identifier==null) throw new ParseException("Identifier cannot be null");
        if (this.identifier!=null) throw new ParseException("Current BioEntry already has a identifier");
        this.identifier = identifier;
    }
    private String identifier;
    
    /**
     * {@inheritDoc}
     */
    public void setName(String name) throws ParseException {
        if (name==null) throw new ParseException("Name cannot be null");
        if (this.name!=null) throw new ParseException("Current BioEntry already has a name");
        this.name = name;
    }
    private String name;
    
    /**
     * {@inheritDoc}
     */
    public void setRankedCrossRef(RankedCrossRef ref) throws ParseException {
        if (ref==null) throw new ParseException("Reference cannot be null");
        ref.setRank(crossRefsRank++);
        this.crossRefs.add(ref);
    }
    private Set crossRefs = new TreeSet();
    private int crossRefsRank = 1;
    
    /**
     * {@inheritDoc}
     */
    public void addSymbols(Alphabet alpha, Symbol[] syms, int start, int length) throws IllegalAlphabetException {
        if (this.symbols==null) {
            if (threshold<=0) {
                this.symbols = new ChunkedSymbolListFactory(this.factory);
            } else {
                this.symbols = new ChunkedSymbolListFactory(this.factory,threshold);
            }
        }
        this.symbols.addSymbols(alpha, syms, start, length);
    }
    private SymbolListFactory factory;
    private int threshold;
    private ChunkedSymbolListFactory symbols;
    
    /**
     * {@inheritDoc}
     */
    public void setComment(String comment) throws ParseException {
        if (comment==null) throw new ParseException("Comment cannot be null");
        this.comments.add(new SimpleComment(comment,commentRank++));
    }
    private Set comments = new TreeSet();
    private int commentRank = 1;
    
    /**
     * {@inheritDoc}
     */
    public void setNamespace(Namespace namespace) throws ParseException {
        if (namespace==null) throw new ParseException("Namespace cannot be null");
        if (this.namespace!=null) throw new ParseException("Current BioEntry already has a namespace");
        this.namespace = namespace;
    }
    private Namespace namespace;
    
    /**
     * {@inheritDoc}
     */
    public void startFeature(Feature.Template templ) throws ParseException {
        try {
            RichFeature f = new SimpleRichFeature(featureHolder,templ);
            f.setRank(this.featureRank++);
            this.allFeatures.add(f);
            if (this.featureStack.size() == 0) this.rootFeatures.add(f);
            else {
                RichFeature parent = (RichFeature)this.featureStack.get(this.featureStack.size() - 1);
                parent.addFeatureRelationship(
                        new SimpleRichFeatureRelationship(parent, f, SimpleRichFeatureRelationship.getContainsTerm(), 0)
                        );
            }
            this.featPropCount = 1; //reset feature anotation rank
            this.featureStack.add(f);
        } catch (ChangeVetoException e) {
            throw new ParseException(e);
        } catch (InvalidTermException e) {
            throw new ParseException(e);
        }
    }
    private FeatureHolder featureHolder = new SimpleFeatureHolder();
    private Set rootFeatures = new TreeSet();
    private List allFeatures = new ArrayList();
    private List featureStack = new ArrayList();
    private int featureRank = 1;
    
    /**
     * {@inheritDoc}
     */ 
    public RichFeature getCurrentFeature() throws ParseException {
        if (this.featureStack.size()==0) throw new ParseException("Not currently within a feature");
        else return (RichFeature)this.featureStack.get(this.featureStack.size()-1);
    }
    
    /**
     * {@inheritDoc}
     */
    public void setTaxon(NCBITaxon taxon) throws ParseException {
        if (taxon==null) throw new ParseException("Taxon cannot be null");
        if (this.taxon!=null){
            if(! this.taxon.equals(taxon)){
                System.err.println(
                        "Warning: attempted to set taxon twice with different values. Keeping first value. "+
                        "old value (retained): "+this.taxon+" new value: "+taxon+", accession: <"+accession+">, version:"+version);
            }
        }
        this.taxon = taxon;
    }
    private NCBITaxon taxon;
    
    /**
     * {@inheritDoc}
     */
    public void setRelationship(BioEntryRelationship relationship) throws ParseException {
        if (relationship==null) throw new ParseException("Relationship cannot be null");
        this.relations.add(relationship);
    }
    private Set relations = new TreeSet();
    
    /**
     * {@inheritDoc}
     */
    public void setRankedDocRef(RankedDocRef ref) throws ParseException {
        if (ref==null) throw new ParseException("Reference cannot be null");
        ref.setRank(referenceCount++);
        this.references.add(ref);
    }
    private Set references = new TreeSet();
    private int referenceCount = 1;
    
    /**
     * {@inheritDoc}
     */
    public void startSequence() throws ParseException {
        this.reset();
    }
    
    /**
     * {@inheritDoc}
     */
    public void addFeatureProperty(Object key, Object value) throws ParseException {
        if (this.featureStack.size() == 0) throw new ParseException("Assertion failed: Not within a feature");
        if (!(key instanceof ComparableTerm)) key = RichObjectFactory.getDefaultOntology().getOrCreateTerm(key.toString());
        if ((value != null) && !(value instanceof String)) value = value.toString();
        RichFeature f = this.getCurrentFeature();
        try {
            Note n = new SimpleNote((ComparableTerm)key,(String)value,this.featPropCount++);
            f.getRichAnnotation().addNote(n);
        } catch (ChangeVetoException e) {
            throw new ParseException(e);
        }
    }
    int featPropCount = 1;
    
    /**
     * {@inheritDoc}
     */
    public void addSequenceProperty(Object key, Object value) throws ParseException {
        if (!(key instanceof ComparableTerm)) key = RichObjectFactory.getDefaultOntology().getOrCreateTerm(key.toString());
        if (value!=null && !(value instanceof String)) value = value.toString();
        try {
            Note n;
            if (value==null) n = new SimpleNote((ComparableTerm)key,null,this.seqPropCount++);
            else n = new SimpleNote((ComparableTerm)key,(String)value,this.seqPropCount++);
            this.notes.addNote(n);
        } catch (ChangeVetoException e) {
            throw new ParseException(e);
        }
    }
    int seqPropCount = 1;
    
    /**
     * {@inheritDoc}
     */
    public void endFeature() throws ParseException {
        if (this.featureStack.size() == 0) throw new ParseException("Assertion failed: Not within a feature");
        this.featureStack.remove(this.featureStack.size() - 1);
    }
    
    /**
     * {@inheritDoc}
     */
    public void endSequence() throws ParseException {
        if (this.name==null) throw new ParseException("Name has not been supplied");
        if (this.namespace==null) throw new ParseException("Namespace has not been supplied");
        if (this.accession==null) throw new ParseException("No accessions have been supplied");
    }
        
    /**
     * {@inheritDoc}
     */
    public void setCircular(boolean circular) throws ParseException { this.circular = circular; }
    private boolean circular = false;
    
    /**
     * {@inheritDoc}
     */
    public Sequence makeSequence() throws BioException {
        this.endSequence(); // Check our input.
        // make our basic object
        SymbolList syms = this.symbols==null?SymbolList.EMPTY_LIST:this.symbols.makeSymbolList();
        RichSequence rs = new SimpleRichSequence(this.namespace,this.name,this.accession,this.version,syms,new Double(this.seqVersion));
        // set misc stuff
        try {
            // set features
            for (Iterator i = this.allFeatures.iterator(); i.hasNext(); ){
                RichFeature f = (RichFeature)i.next();
                f.setParent(rs);
                if (f.getName()==null || f.getName().length()==0) f.setName(rs.getAccession()+"#"+f.getRank()); // dummy feature name for use in GBrowse
            }
            rs.setDescription(this.description);
            rs.setDivision(this.division);
            rs.setIdentifier(this.identifier);
            rs.setTaxon(this.taxon);
            rs.setCircular(this.circular);
            rs.setFeatureSet(this.rootFeatures);
            for (Iterator i = this.crossRefs.iterator(); i.hasNext(); ) rs.addRankedCrossRef((RankedCrossRef)i.next());
            for (Iterator i = this.relations.iterator(); i.hasNext(); ) rs.addRelationship((BioEntryRelationship)i.next());
            for (Iterator i = this.references.iterator(); i.hasNext(); ) rs.addRankedDocRef((RankedDocRef)i.next());
            for (Iterator i = this.comments.iterator(); i.hasNext(); ) rs.addComment((Comment)i.next());
            // set annotations
            rs.setNoteSet(this.notes.getNoteSet());
        } catch (Exception e) {
            throw new ParseException(e); // Convert them all to parse exceptions.
        }
        // return the object
        return rs;
    }
    
    /**
     * {@inheritDoc}
     */
    public RichSequence makeRichSequence() throws BioException { return (RichSequence)this.makeSequence(); }
}
