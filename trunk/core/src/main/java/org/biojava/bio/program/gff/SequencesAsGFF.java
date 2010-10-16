/* -*- c-basic-offset: 2; indent-tabs-mode: nil -*- */
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

package org.biojava.bio.program.gff;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.Location;

/**
 * Turns a sequence database into a GFF event stream.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Len Trigg
 */
public class SequencesAsGFF {
  /**
   * The <span class="type">FeatureFilter</span> for selecting features to
   * report as <span class="type">GFFRecord</span>s.
   */
  private FeatureFilter filter = FeatureFilter.all;
  
  /**
   * Whether or not to recurse through the features during searching.
   */
  private boolean recurse = false;
  
  /**
   * Whether or not non-contiguous features should be broken into blocks
   * 
   * @since 1.4
   */
  
  private boolean shatter = false;
  
  private boolean generateSequenceHeader = true;
  
  /**
   * Specify whether features with non-contiguous locations should be broken
   * up such that a GFF feature line is emitted for each contiguous block.
   * 
   * @param b
   * @since 1.4
   */
  
  public void setShatter(boolean b) {
      this.shatter = b;
  }
  
  /**
   * Determine if features with non-contiguous locations will be broken into
   * multiple GFF records.
   * 
   * @since 1.4
   */
  
  public boolean getShatter() {
      return shatter;
  }
  
  /**
   * Specify whether a per-sequence header line, giving the length of the
   * sequence, should be generated.
   *
   * @since 1.4
   */
   
   public void setGenerateSequenceHeader(boolean b) {
       this.generateSequenceHeader = b;
   }
   
   /**
    * Discover if per-sequence header lines will be generated.
    *
    * @since 1.4
    */
   
   public boolean getGenerateSequenceHeader() {
       return generateSequenceHeader;
   }
  
  /**
   * Return the current <span class="type">FeatureFilter</span>.
   * <p>
   * This is the object that will accept or reject individual features.
   *
   * @return the current <span class="type">FeatureFilter</span>
   */
  public FeatureFilter getFeatureFilter() {
    return filter;
  }
  
  /**
   * Replace the current <span class="type">FeatureFilter</span> with
   * <span class="arg">filter</span>.
   *
   * @param filter  the new <span class="type">FeatureFilter</span>
   */
  public void setFeatureFilter(FeatureFilter filter) {
    this.filter = filter;
  }
  
  /**
   * Return whether features will be filtered recursively or not.
   *
   * @return whether or not to recurse
   */
  public boolean getRecurse() {
    return recurse;
  }
  
  /**
   * Set whether features will be filtered recursively to
   * <span class="arg">recurse</span>.
   *
   * @param recurse  <span class="kw">true</span> if you want to recurse,
   *                 <span class="kw">false</span> otherwise
   */
  public void setRecurse(boolean recurse) {
    this.recurse = recurse;
  }

  /**
   * Emit any per-sequence header information.
   * The default implementation emits sequence-region comment lines.
   *
   * @since 1.4
   */
  
  protected void doPreProcessSequence(
    Sequence seq,
    GFFDocumentHandler handler,
    String id
  )
    throws BioException
  {
      if (generateSequenceHeader) {
          handler.commentLine("#sequence-region " + id + " 1 " + seq.length());
      }
  }
  
  /**
   * Internal method to process an individual <span class="type">Sequence</span>.
   *
   * @param seq  the <span class="type">Sequence</span> to GFFify
   * @param handler the <span class="type">GFFDocumentHandler</span> that will
   *                receive the GFF for all suitable features within
   *                <span class="arg">seq</span>
   * @param id the value of the <span class="method">seqName</span> field in any
   *           <span class="type">GFFRecord</span>s produced
   */
  protected void doProcessSequence(Sequence seq,
                                   GFFDocumentHandler handler,
                                   String id) 
    throws BioException 
  {
    Iterator fi = seq.filter(getFeatureFilter(), getRecurse()).features();
      
    while (fi.hasNext()) {
      doProcessFeature((Feature) fi.next(), handler, id);
    }
  }


  /**
   * Internal method to process an individual <span class="type">Feature</span>.
   *
   * @param feature  the <span class="type">Feature</span> to GFFify
   * @param handler the <span class="type">GFFDocumentHandler</span> that will
   *                receive the GFF for this feature
   * @param id the value of the <span class="method">seqName</span> field in any
   *           <span class="type">GFFRecord</span>s produced
   */
  protected void doProcessFeature(Feature feature,
                                  GFFDocumentHandler handler,
                                  String id) 
    throws BioException 
  {
    SimpleGFFRecord record = createGFFRecord(feature, id);
    if (shatter && !feature.getLocation().isContiguous()) {
        for (Iterator si = feature.getLocation().blockIterator(); si.hasNext(); ) {
            Location shatterBloc = (Location) si.next();
            record.setStart(shatterBloc.getMin());
            record.setEnd(shatterBloc.getMax());
            handler.recordLine(record);
        }
    } else {
        handler.recordLine(record);
    }
  }


  /**
   * Internal method to create a <span class="type">GFFRecord</span>
   * from an individual <span class="type">Feature</span>.
   *
   * @param feature  the <span class="type">Feature</span> to GFFify
   * @param id the value of the <span class="method">seqName</span> field in any
   *           <span class="type">GFFRecord</span>s produced
   */
  protected SimpleGFFRecord createGFFRecord(Feature feature,
                                            String id) 
    throws BioException {
    
    SimpleGFFRecord record = new SimpleGFFRecord();
    record.setSeqName(id);
    record.setSource(feature.getSource());
    record.setFeature(feature.getType());
    Location loc = feature.getLocation();
    record.setStart(loc.getMin());
    record.setEnd(loc.getMax());
    record.setScore(GFFTools.NO_SCORE);
    record.setStrand(StrandedFeature.UNKNOWN);
    if (feature instanceof StrandedFeature) {
      StrandedFeature sf = (StrandedFeature) feature;
      if (sf.getStrand() == StrandedFeature.POSITIVE) {
        record.setStrand(StrandedFeature.POSITIVE);
      } else if (sf.getStrand() == StrandedFeature.NEGATIVE) {
        record.setStrand(StrandedFeature.NEGATIVE);
      }
    }
    record.setFrame(GFFTools.NO_FRAME);
    Map fMap = feature.getAnnotation().asMap();
    Map fMap2 = new HashMap();
    for (Iterator ki = fMap.keySet().iterator(); ki.hasNext(); ) {
      Object key = ki.next();
      Object value = fMap.get(key);
      String keyS = key.toString();
      List valueList;
      if (value instanceof Collection) {
        valueList = new ArrayList((Collection) value);
      } else {
        //valueList = Collections.singletonList(value); 1.3?
        valueList = new ArrayList();
        valueList.add(value);
      }
      for (int i = 0; i < valueList.size(); i++) {
        Object o = valueList.get(i);
        valueList.set(i, o.toString());
      }
      fMap2.put(keyS, valueList);
    }
    record.setGroupAttributes(fMap2);
    record.setComment(null);        

    return record;
  }


  /**
   * Process an individual <span class="type">Sequence</span>, informing
   * <span class="arg">handler</span> of any suitable features.
   *
   * @param seq  the <span class="type">Sequence</span> to GFFify
   * @param handler the <span class="type">GFFDocumentHandler</span> that will
   *                receive the GFF for all suitable features within
   *                <span class="arg">seq</span>
   */
  public void processSequence(Sequence seq, GFFDocumentHandler handler) 
  throws BioException {
    handler.startDocument(seq.getURN());
    doPreProcessSequence(seq, handler, seq.getName());
    doProcessSequence(seq, handler, seq.getName());
    handler.endDocument();
  }

  /**
   * Process all <span class="type">Sequence</span>s within a
   * <span class="type">SequenceDB</span>, informing
   * <span class="arg">handler</span> of any suitable features.
   *
   * @param seqDB  the <span class="type">SequenceDB</span> to GFFify
   * @param handler the <span class="type">GFFDocumentHandler</span> that will
   *                receive the GFF for all suitable features within
   *                <span class="arg">seqDB</span>
   */
  public void processDB(SequenceDB seqDB, GFFDocumentHandler handler)
  throws BioException {
    handler.startDocument("unknown:SequenceDB:" + seqDB.getName());
    for(Iterator i = seqDB.ids().iterator(); i.hasNext(); ) {
      String id = (String) i.next();
      Sequence seq = seqDB.getSequence(id);
      doPreProcessSequence(seq, handler, id);
    }
    for(Iterator i = seqDB.ids().iterator(); i.hasNext(); ) {
      String id = (String) i.next();
      Sequence seq = seqDB.getSequence(id);
      doProcessSequence(seq, handler, id);
    }
    handler.endDocument();
  }
}
