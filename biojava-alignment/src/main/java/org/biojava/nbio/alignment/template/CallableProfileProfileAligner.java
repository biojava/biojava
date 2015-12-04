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
 * Created on July 9, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment.template;

import org.biojava.nbio.core.alignment.template.ProfilePair;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.concurrent.Callable;

/**
 * Implements a concurrency wrapper for a {@link ProfileProfileAligner}.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the {@link Profile} pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class CallableProfileProfileAligner<S extends Sequence<C>, C extends Compound>
        implements Callable<ProfilePair<S, C>> {

    private ProfileProfileAligner<S, C> ppa;

    /**
     * Creates a profile-profile alignment task for simplified parallel execution.
     *
     * @param ppa already initialized profile-profile aligner
     */
    public CallableProfileProfileAligner(ProfileProfileAligner<S, C> ppa) {
        this.ppa = ppa;
    }

    @Override
    public ProfilePair<S, C> call() {
        return ppa.getPair();
    }

}
