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

import java.util.Collections;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.seq.ProteinTools;

/**
 * Simple no-frills implementation of a WobbleDistribution object.
 *
 * @author David Huen
 * @since 1.3
 */
class SimpleWobbleDistribution
    implements WobbleDistribution
{
    Symbol residue;
    Set nonWobbleBases;
    Distribution nonWobbleDist;
    Map wobbleDist;

    SimpleWobbleDistribution(Symbol residue, Set nonWobbleBases, Distribution nonWobbleDist, Map wobbleDist)
        throws NullPointerException, IllegalAlphabetException
    {
        if ((residue == null)
            || (nonWobbleBases == null)
            || (nonWobbleDist == null) 
            || (wobbleDist == null)
            ) throw new NullPointerException();

        if (!ProteinTools.getTAlphabet().contains(residue)) throw new IllegalAlphabetException();

        this.residue = residue;
        this.nonWobbleBases = nonWobbleBases;
        this.nonWobbleDist = nonWobbleDist;
        this.wobbleDist = wobbleDist;
    }

    public Symbol getResidue()
    {
        return residue;
    }

    public Set getNonWobbleBases()
    {
        return Collections.unmodifiableSet(nonWobbleBases);
    }

    public Distribution getFrequencyOfNonWobbleBases()
    {
        return nonWobbleDist;
    }

    public Distribution getWobbleFrequency(Symbol nonWobbleBases)
    {
        return (Distribution) wobbleDist.get(nonWobbleBases);
    }
}
