 
/**
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

package org.biojava.bio.seq.io.game;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.RangeLocation;
/**
 * An interface that can be tested for by nested handlers
 * when trying to do a callback.
 * This one allows an element that represents a transcript
 * to gain information about locations of nested elements
 * that represent exons/introns.
 *
 * @author David Huen
 * @since 1.8
 */
public interface GAMETranscriptCallbackItf {

/**
 * Allows nesting class that manages a transcript template
 * to gain information about its extent from nested
 * elements that represent exons.
 */
  public void reportExon(RangeLocation range, StrandedFeature.Strand strand);
}

