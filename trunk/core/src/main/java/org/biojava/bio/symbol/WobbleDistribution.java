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

package org.biojava.bio.symbol;

import java.util.Set;

import org.biojava.bio.dist.Distribution;

/**
 * an object to return statistics about
 * the frequency of the wobble base
 * in a set of synonymous codons.
 *
 * @author David Huen
 * @since 1.3
 */
public interface WobbleDistribution
{
    /**
     * returns the residue encoded by this WobbleDistribution
     */
    public Symbol getResidue();

    /**
     * returns Set containing the nonWobbleBases that
     * occur in codons that encode this residue
     */
    public Set getNonWobbleBases();

    /**
     * returns the frequency with which
     * synonymous codons start with a
     * specified pair of bases.
     */
    public Distribution getFrequencyOfNonWobbleBases();

    /**
     * returns the frequency of a specific
     * wobble base in a set of synonymous
     * codons that start with the same two
     * bases. (sums to one over this set).
     */
    public Distribution getWobbleFrequency(Symbol nonWobbleBases);
}

