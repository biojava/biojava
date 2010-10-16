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
package org.biojava.bio.seq.homol;
 
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.utils.ChangeType;

/**
 * <p>
 * Signifies that two or more features are homologous.
 * </p> 
 *
 * <p> Blast hits or local multiple-sequence alignments can be
 * represented as a set of features on sequences that have an
 * alignment. The features will probably implement
 * HomologyFeature.
 * </p>
 *
 * @author Matthew Pocock
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 * @since 1.2
 */
public interface Homology {
  /**
   * Signals that the alignment describing the homologous sequences
   * has changed. For implementations which implement
   * <code>Changeable</code>.
   */
  public static final ChangeType ALIGNMENT =
      new ChangeType("The alignment has been changed",
                     "org.biojava.bio.seq.homol.Homology",
                     "ALIGNMENT");

  /**
   * Retrieve the set of features that mark homologous regions.
   *
   * @return the FeatureHolder containing each homologous region
   */
  FeatureHolder getFeatures();
  /**
   * Retrieve the Alignment that specifies how the homologous regions are
   * aligned. The labels of the alignment are the HomologyFeature objects.
   *
   * @return the Alignment between the HomologyFeatures
   */
  Alignment getAlignment();
}
