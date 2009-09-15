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

package org.biojava.bio.molbio;

import org.biojava.bio.seq.StrandedFeature;

/**
 * <code>RestrictionSite</code> represents the recognition site of a
 * restriction enzyme.
 *
 * @author Keith James
 * @since 1.3
 */
public interface RestrictionSite extends StrandedFeature
{
    /**
     * <code>getPosition</code> returns the common, forward strand cut
     * site. Note that some enzymes cut in more than one
     * position. Such supplementary sites may be calculated by
     * retrieving the <code>RestrictionEnzyme</code> instance and
     * using its methods to calculate the position.
     *
     * @return an <code>int</code> indicating the base immediately
     * before the cleavage site on the forward strand.
     */
    public int getPosition();

    /**
     * <code>getEnzyme</code> returns the enzyme which cuts at this
     * site. A sequence which is the target for several different
     * enzymes is expected to have a corresponding
     * <code>RestrictionSite</code> feature for each.
     *
     * @return a <code>RestrictionEnzyme</code>.
     */
    public RestrictionEnzyme getEnzyme();

    /**
     * <code>Template</code> for construction of
     * <code>RestrictionSite</code>s.
     */
    public static class Template extends StrandedFeature.Template
    {
        /**
         * <code>enzyme</code> <code>RestrictionEnzyme</code> field.
         */
        public RestrictionEnzyme enzyme;
    }
}
