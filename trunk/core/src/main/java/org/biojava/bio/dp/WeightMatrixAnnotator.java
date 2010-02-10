package org.biojava.bio.dp;

import java.io.Serializable;

import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceAnnotator;
import org.biojava.bio.seq.impl.ViewSequence;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;

/**
 * Annotates a sequence with hits to a weight-matrix.
 *
 * <p>
 * This SequenceAnnotator implementation returns a new
 * ViewSequence wrapping the underlying Sequence
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Tanya Vavouri
 */
public class WeightMatrixAnnotator implements SequenceAnnotator,
    Serializable {
  private WeightMatrix matrix;
  private double threshold;
  private final ScoreType scoreType;
  private String wmID;

  public Sequence annotate(Sequence seq) throws IllegalAlphabetException,
      BioException, ChangeVetoException {
    seq = new ViewSequence(seq);

    int cols = matrix.columns();
    Feature.Template template = new Feature.Template();
    template.source = "WeightMatrixAnnotator";
    template.type = wmID;
    for (int offset = 1;
         offset <= seq.length() - cols + 1;
         offset++) {
      double score = DP.scoreWeightMatrix(matrix, seq, scoreType, offset);
      double q = Math.exp(score);
      if (q >= threshold) {
        template.location = new RangeLocation(offset, offset + cols - 1);
        SimpleAnnotation ann = new SimpleAnnotation();
        ann.setProperty("score", new Double(q));
        ann.setProperty("weightMatrix", matrix);
        template.annotation = ann;
        seq.createFeature(template);
      }
    }
    return seq;
  }

  /**
   * Create a new annotator that uses the PROBABILITY score type and an ID
   for the weight matrix.
   *
   * @param wm        the weight matrix
   * @param threshold the threshold
   * @param wmID the weight matrix ID
   */
  public WeightMatrixAnnotator(WeightMatrix wm, ScoreType scoreType,
                               double threshold, String wmID) {
    this.matrix = wm;
    this.threshold = threshold;
    this.scoreType = ScoreType.PROBABILITY;
    this.wmID = wmID;
  }

  /**
   * Create a new annotator that uses PROBABILITY score type.
   *
   * @param wm a <code>WeightMatrix</code> value
   * @param threshold a <code>double</code> value
   */
  public WeightMatrixAnnotator(WeightMatrix wm, double threshold) {
    this(wm, ScoreType.PROBABILITY, threshold, "hit");
  }

  /**
   * Create a new annotator that uses a specific score type.
   *
   * @param wm        the weigth matrix
   * @param scoreType the score type
   * @param threshold the threshold
   * @since 1.4
   */
  public WeightMatrixAnnotator(WeightMatrix wm, ScoreType scoreType,
                               double threshold) {
    this.matrix = wm;
    this.scoreType = scoreType;
    this.threshold = threshold;
    this.wmID = "hit";
  }

  /**
   * Get the value of the weight matrix id.
   * @return value of the weight matrix id.
   */
  public String getWeightMatrixID() {
    return wmID;
  }

  /**
   * Set the weight matrix id.
   * @param id  Value to assign to the weight matrix id.
   */
  public void setWeightMatrixID(String id) {
    this.wmID = id;
  }

}