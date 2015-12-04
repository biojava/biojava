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
 * Created on July 12, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment;

import org.biojava.nbio.core.alignment.SimpleSequencePair;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.PairInProfileScorer;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Implements an algorithm which computes a score for a sequence alignment pair picked from an alignment
 * {@link Profile}.  The reported score is the number of alignment columns which have identical {@link Compound}s.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class FractionalIdentityInProfileScorer<S extends Sequence<C>, C extends Compound>
        extends FractionalIdentityScorer<S, C> implements PairInProfileScorer<S, C> {

    private Profile<S, C> profile;

    /**
     * Creates a fractional identity scorer for an aligned pair of sequences in the given alignment profile.
     *
     * @param profile alignment profile containing pair of sequences
     * @param query index in the profile of the first sequence of the pair
     * @param target index in the profile of the second sequence of the pair
     */
    public FractionalIdentityInProfileScorer(Profile<S, C> profile, int query, int target) {
        super(new SimpleSequencePair<S, C>(profile.getAlignedSequence(query), profile.getAlignedSequence(target)));
        this.profile = profile;
    }

    @Override
    public Profile<S, C> getProfile() {
        return profile;
    }

}
