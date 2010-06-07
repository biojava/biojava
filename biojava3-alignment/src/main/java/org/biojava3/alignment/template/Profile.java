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

import java.util.List;

import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

public interface Profile<S extends Sequence<C>, C extends Compound> extends Iterable<S> {

    AlignedSequence<C> getAlignedSequence(int index);

    AlignedSequence<C> getAlignedSequence(S sequence); // will find either aligned or original sequences

    List<AlignedSequence<C>> getAlignedSequences(); // unmodifiable unless class implements MutableProfile

    List<AlignedSequence<C>> getAlignedSequences(int... indices); // useful for views

    List<AlignedSequence<C>> getAlignedSequences(S... sequences); // useful for views

    C getCompoundAt(int index, int alignmentIndex);

    C getCompoundAt(S sequence, int alignmentIndex); // will find either aligned or original sequences

    List<C> getCompoundsAt(int alignmentIndex); // useful for views

    CompoundSet<C> getCompoundSet();

    int[] getIndicesAt(int alignmentIndex); // useful for views

    int getIndexOf(C compound);

    int getLastIndexOf(C compound);

    int getLength(); // number of columns

    int getSize(); // number of rows ... if !isCircular() ? == number of sequences : >= number of sequences

    ProfileView<S, C> getSubProfile(Location location); // only include sequences that overlap Location

    boolean isCircular(); // if so, sequences longer than length() return multiple compounds at any location

    String toString(); // simple view: each sequence on 1 line

    String toString(int width); // formatted view: show start and end indices of profile and sequences, limited line length

}
