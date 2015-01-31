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

import org.biojava.nbio.alignment.Alignments.PairInProfileScorerType;
import org.biojava.nbio.alignment.Alignments.ProfileProfileAlignerType;
import org.biojava.nbio.alignment.template.AbstractScorer;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.alignment.template.RescoreRefiner;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

public class StandardRescoreRefiner<S extends Sequence<C>, C extends Compound> extends AbstractScorer
        implements RescoreRefiner<S, C> {

    private PairInProfileScorerType pips;
    private ProfileProfileAlignerType ppa;

    public StandardRescoreRefiner(PairInProfileScorerType pips, ProfileProfileAlignerType ppa) {
        this.pips = pips;
        this.ppa = ppa;
    }

    // methods for RescoreRefiner

    @Override
    public PairInProfileScorerType getPairInProfileScorer() {
        return pips;
    }

    @Override
    public ProfileProfileAlignerType getProfileProfileAligner() {
        return ppa;
    }

    // methods for Aligner

    @Override
    public long getComputationTime() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public Profile<S, C> getProfile() {
        // TODO Auto-generated method stub
        return null;
    }

    // methods for Scorer

    @Override
    public double getMaxScore() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public double getMinScore() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public double getScore() {
        // TODO Auto-generated method stub
        return 0;
    }

}
