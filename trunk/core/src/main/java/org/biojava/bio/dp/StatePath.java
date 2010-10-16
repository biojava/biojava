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


package org.biojava.bio.dp;

import org.biojava.bio.alignment.Alignment;

/**
 * Extends the Alignment interface so that it is explicitly used to represent
 * a state path through an HMM, and the associated emitted sequence and
 * likelihoods.
 * <p>
 * A state path should have the following structure:
 * <bq>
 * STATES -> list of all states used by the machine
 * <br>
 * SCORES -> list of step-wise scores for each state (transition + emission)
 * <br>
 * SEQUENCE -> sequence emitted by the machine 
 * </bq>
 * The sequence emitted by the machine will be some function of the sequences
 * that were aligned to the machine, and the state-path taken. Whenever the
 * state used is a non-emitting state, this emitted sequence is a gap. Whenever
 * it is an emission state, it is the symbol matched by that state. This is
 * modeled by the following nesting:
 * <bq><pre>
 * SEQUENCE
 *   -> Gapped view (gap inserted for every position aligned with a dot-state
 *     -> Sequence emitted by emission states as Alignment
 *       label_n = input_SymbolList_n
 *         -> gapped view of SymbolList_n
 * </pre></bq>
 * A multi-head HMM (2 or more) emits a single sequence that is
 * an alignment of the input sequences with gaps added. In this case, the
 * emitted sequence should be an Alignment object with labels being the input
 * sequences, and the associated SymbolList objects being gapped views. For the
 * sake of least-suprise, single-head HMMs should emit an alignment of one
 * sequence, where the label is the input sequence, and the associated
 * SymbolList is also the input sequence.
 *
 * <p>
 * I think that this scheme keeps the emitted alignment as close as possible to
 * a sensible path through the sequence coordinate space, while making this
 * gappable adapts this to the same co-ordinate system as the HMM state-path
 * space.
 * </p>
 *
 * @author Matthew Pocock
 */
public interface StatePath extends Alignment {
  /**
   * Alignment label for the emitted sequence.
   */
  public static final Object SEQUENCE = "SEQUENCE";

  /**
   * Alignment label for the state path.
   */
  public static final Object STATES   = "STATES";

  /**
   * Alignment label for the likelyhood at each step.
   */
  public static final Object SCORES   = "SCORES";
  
  /**
   * Return the overall score for this state-path and it's emissions.
   *
   * @return the score
   */
  public double getScore();
}
