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
 * Created on June 7, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.template;

import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

public interface AlignedSequence<C extends Compound> extends Sequence<C> {

    int getAlignmentIndexAt(int sequenceIndex);

    int getEnd();

    Location getLocationInAlignment();

    int getNumGaps();

    Sequence<C> getOriginalSequence();

    int getOverlapCount(); // if !isCircular() ? == 1 : >= 1

    int getSequenceIndexAt(int alignmentIndex);

    int getStart();

    boolean isCircular();

}
