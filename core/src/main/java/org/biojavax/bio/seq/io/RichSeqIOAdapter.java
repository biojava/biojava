/*
 * RichSeqIOAdapter.java
 *
 * Created on October 6, 2005, 4:53 PM
 */

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

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojavax.Namespace;
import org.biojavax.RankedCrossRef;
import org.biojavax.RankedDocRef;
import org.biojavax.bio.BioEntryRelationship;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.taxa.NCBITaxon;

/**
 * This class implements all methods of RichSeqIOListener and takes no action.
 * It should be overridden to implement custom listeners that only listen for
 * a small subset of events.
 *
 * @author Mark Schreiber
 * @since 1.5
 */
public class RichSeqIOAdapter implements RichSeqIOListener {
    
    /**
     * This is a dummy feature. It is returned by the method 
     * {@link #getCurrentFeature() getCurrentFeature()}. Access is provided so
     * you can override it.
     */
    protected RichFeature emptyFeature;
    
    /** Creates a new instance of RichSeqIOAdapter */
    public RichSeqIOAdapter() {
        emptyFeature = RichFeature.Tools.makeEmptyFeature();
    }
    
    
    
    public void setAccession(String accession) throws ParseException{}
    public void setIdentifier(String identifier) throws ParseException{}
    public void setDivision(String division) throws ParseException{}
    public void setDescription(String description)throws ParseException{}
    public void setVersion(int version) throws ParseException{}
    public void setSeqVersion(String version) throws ParseException{}
    public void setComment(String comment) throws ParseException{}
    public void setRankedDocRef(RankedDocRef ref) throws ParseException{}
    public void setTaxon(NCBITaxon taxon) throws ParseException{}
    public void setNamespace(Namespace namespace) throws ParseException{}
    public void setRelationship(BioEntryRelationship relationship) throws ParseException{}
    public void setRankedCrossRef(RankedCrossRef crossRef) throws ParseException{}
    public void setURI(String uri) throws ParseException{}
    public RichFeature getCurrentFeature() throws ParseException{return this.emptyFeature;}
    public void setCircular(boolean circular) throws ParseException{}
    public void addFeatureProperty(Object key, Object value) throws ParseException{}
    public void endFeature() throws ParseException{}
    public void startFeature(Feature.Template templ) throws ParseException{
        this.emptyFeature = RichFeature.Tools.makeEmptyFeature();
    }
    public void addSequenceProperty(Object key, Object value) throws ParseException{}
    public void addSymbols(Alphabet alpha, Symbol[] syms, int start, int length)
    throws IllegalAlphabetException{}
    public void setName(String name) throws ParseException{}
    public void endSequence() throws ParseException{}
    public void startSequence() throws ParseException{}
}
