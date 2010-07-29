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

package org.biojava3.alignment;

import org.biojava3.alignment.Alignments.PairInProfileScorerType;
import org.biojava3.alignment.Alignments.ProfileProfileAlignerType;
import org.biojava3.alignment.template.AbstractScorer;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.RescoreRefiner;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

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
    public int getMaxScore() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int getMinScore() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int getScore() {
        // TODO Auto-generated method stub
        return 0;
    }

}
