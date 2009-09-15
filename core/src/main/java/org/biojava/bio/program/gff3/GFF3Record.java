package org.biojava.bio.program.gff3;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.ontology.OntoTools;
import org.biojava.ontology.Term;
import org.biojava.utils.AbstractChangeable;

/**
 * A record in a GFF3 formatted file.
 *
 * @author Matthew Pocock
 */
public interface GFF3Record
extends Annotatable {
  public String getSequenceID();
  
  public Term getSource();
  
  public Term getType();
  
  public int getStart();
  
  public int getEnd();
  
  public double getScore();
  
  public StrandedFeature.Strand getStrand();
  
  public int getPhase();
  
  public static final class Impl
  extends AbstractChangeable
  implements GFF3Record {
    
    private String sequenceID;
    private Term source;
    private Term type;
    private int start;
    private int end;
    private double score;
    private StrandedFeature.Strand strand;
    private int phase;
    private Annotation annotation;
    
    public Impl() {
      // do nothing much - initialize us with uninformative data
      sequenceID = null;
      source = OntoTools.ANY;
      type = OntoTools.ANY;
      start = Integer.MAX_VALUE;
      end = Integer.MIN_VALUE;
      score = GFFTools.NO_SCORE;
      strand = StrandedFeature.UNKNOWN;
      phase = GFFTools.NO_FRAME;
    }
    
    public Impl(GFF3Record rec) {
      this.sequenceID = rec.getSequenceID();
      this.source = rec.getSource();
      this.type = rec.getType();
      this.start = rec.getStart();
      this.end = rec.getEnd();
      this.score = rec.getScore();
      this.strand = rec.getStrand();
      this.phase = rec.getPhase();
    }
    
    public String getSequenceID() {
      return this.sequenceID;
    }
    
    public void setSequenceID(String sequenceID) {
      this.sequenceID = sequenceID;
    }
    
    public Term getSource() {
      return this.source;
    }
    
    public void setSource(Term source) {
      this.source = source;
    }
    
    public Term getType() {
      return this.type;
    }
    
    public void setType(Term type) {
      this.type = type;
    }
    
    public int getStart() {
      return this.start;
    }
    
    public void setStart(int start) {
      this.start = start;
    }
    
    public int getEnd() {
      return this.end;
    }
    
    public void setEnd(int end) {
      this.end = end;
    }
    
    public double getScore() {
      return this.score;
    }
    
    public void setScore(double score) {
      this.score = score;
    }
    
    public StrandedFeature.Strand getStrand() {
      return this.strand;
    }
    
    public void setStrand(StrandedFeature.Strand strand) {
      this.strand = strand;
    }
    
    public int getPhase() {
      return this.phase;
    }
    
    public void setPhase(int phase) {
      this.phase = phase;
    }
    
    public Annotation getAnnotation() {
      if(annotation == null) {
        annotation = new SmallAnnotation();
      }
      
      return annotation;
    }
  }
}
