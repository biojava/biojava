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
 * Created on June 22, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.core.alignment;

import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.ProfilePair;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.List;

/**
 * Implements a data structure for the results of the alignment of a pair of {@link Profile}s.
 *
 * @author Mark Chapman
 * @author Paolo Pavan
 * @param <S> each element of an alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class SimpleProfilePair<S extends Sequence<C>, C extends Compound> extends SimpleProfile<S, C>
        implements ProfilePair<S, C> {

    private Profile<S, C> query, target;

    /**
     * Creates a pair profile for the given profiles.
     *
     * @param query the first profile of the pair
     * @param target the second profile of the pair
     * @param sx lists whether the query profile aligns a {@link Compound} or gap at each index of the alignment
     * @param sy lists whether the target profile aligns a {@link Compound} or gap at each index of the alignment
     * @throws IllegalArgumentException if alignments differ in size or given profiles do not fit in alignments
     */
    public SimpleProfilePair(Profile<S, C> query, Profile<S, C> target, List<Step> sx, List<Step> sy) {
        super(query, target, sx, sy);
        this.query = query;
        this.target = target;
    }

    @Override
    public Profile<S, C> getQuery() {
        return query;
    }

    @Override
    public Profile<S, C> getTarget() {
        return target;
    }

}
