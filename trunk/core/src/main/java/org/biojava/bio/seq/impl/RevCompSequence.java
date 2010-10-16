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

package org.biojava.bio.seq.impl;

import java.util.Iterator;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceTools;
import org.biojava.bio.seq.projection.ProjectedFeatureHolder;
import org.biojava.bio.seq.projection.Projection;
import org.biojava.bio.seq.projection.TranslateFlipContext;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * A reverse complement view onto <code>Sequence</code> interface.
 * <p>
 * All features of the underlying sequence are reflected onto the RevCompSequence using a ProjectedFeatureHolder</p>
 * calling createFeature() on a RevCompSequence creates a feature on the underlying
 * sequence. Non-Stranded features will return the reverse compemented view of the sequence
 * when getSymbols() is called that is to say if you get what you expect as if your RevCompSequence
 * was a regular Sequence.
 *
 * @author David Waring
 * @author Thomas Down
 */
public class RevCompSequence
    extends SimpleSequence
{
    private ProjectedFeatureHolder pfh;
    protected Sequence origSeq;
    
    
    /**
    *  URN, Name and Annotation are copied as is from the original Sequence, unless you use the
    *  the other contructor that sets these.
    */
    
    public RevCompSequence(Sequence seq)
        throws IllegalAlphabetException
    {
        this(seq,seq.getURN(),seq.getName(),seq.getAnnotation());
    }
    
    
    public RevCompSequence(Sequence seq, String urn, String name, Annotation annotation)throws IllegalAlphabetException {
        super(DNATools.reverseComplement(seq),urn,name,annotation);
        pfh = new ProjectedFeatureHolder(new TranslateFlipContext(this,seq,seq.length()+1,true));
        origSeq = seq;
    }
    

    // SymbolList stuff
    /**
    * edit() will try to edit the underlying Sequence. So if it is editable this will be too
    * <p>Since I have not seen and editable Sequence I have not tested this </p>
    *
    */
    public void edit(Edit e)throws ChangeVetoException,IndexOutOfBoundsException{
        int pos = (this.length() - (e.pos + e.length)) + 2;
        Edit newE = null;
        try {
            newE = new Edit (pos,e.length,DNATools.reverseComplement(e.replacement));
            origSeq.edit(newE);
        }catch (IllegalAlphabetException iae){
            throw new BioError("Error while editing RevCompSequence " + iae.getMessage());
        }
        
    }
    
    // Sequence stuff
    public Iterator features(){
        return pfh.features();
    }
    
    public int countFeatures(){
        return pfh.countFeatures();
    }
    
    public FeatureHolder filter(FeatureFilter ff) {
        return pfh.filter(ff);
    }
    
    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        return pfh.filter(ff, recurse);
    }
    
    /**
    * containsFeature() will return true if this seq contains the feature in question, or
    * if if the original (non reverse complement) sequence contains the feature;
    */ 
    
    public boolean containsFeature(Feature f) {
        return pfh.containsFeature(f) || origSeq.containsFeature(f);
    }
    
    public void removeFeature(Feature f)
            throws ChangeVetoException, BioException {
        pfh.removeFeature(f);
    }

    /**
    * createFeature() will call createFeature() on the underlying Sequence.
    * returns the feature as it will be projected onto the reverse complement sequence 
    * not the actual feature that was created.
    *
    */
    public Feature createFeature(Feature.Template ft) throws ChangeVetoException,BioException{
    	 return pfh.getContext().createFeature(ft);
    }
    
    /**
    * getFeatureFromOriginal() Since you can not create a feature on a projectedFeature at this time, I am 
    * including this method so that you can get the corresponding feature from the original sequence.
    * (which is not projected) and do something with that such as createFeature().
    */
    
    public Feature getFeatureFromOriginal(Feature f){
        return ((Projection) f).getViewedFeature();
    }
    
    /**
    * clone() should make a complete copy of the Sequence with  all features (and children) and return
    * a SimpleSequence that is unconnected from the original sequence.
    */
    
    public Object clone(){
        SymbolList sl = new SimpleSymbolList(this);
        Sequence newSeq = new SimpleSequence(sl,this.getURN(),this.getName(),this.getAnnotation());
        try{
            SequenceTools.addAllFeatures(newSeq, this);
        } catch ( BioException e){
            throw new BioError( "Error while cloning RevCompSequenece: " + e.getMessage());
        } catch (ChangeVetoException cve) {
            throw new BioError("Couldn't modify newly created SimpleSequence", cve);
        }
            
        return newSeq;
        
    }     
}