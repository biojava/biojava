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

import org.biojava.bio.dist.Distribution;

public interface CodonPref
{
    /**
     * Distributions on codon preferences can be based on
     * i) codon frequency over all codons
     * ii) per residue: codon fraction over all codons (returns a Distribution)
     * iii) per residue: codon fraction per dinucleotide + wobble distribution (returns a WobbleDistribution)
     */

    /**
     * get name of object
     */
    public String getName();

    /**
     * get the name of the genetic code
     */
    public String getGeneticCodeName();

    /**
     * the genetic code that this codon
     * preference is based on.
     */
    public ManyToOneTranslationTable getGeneticCode();

    /**
     * returns a Distribution giving the 
     * frequency of codons (sums to one over the 
     * totality of codons).
     */
    public Distribution getFrequency();

    /**
     * returns a Distribution giving the
     * frequency of synonymous codons.
     * (sums to one over the total number
     * of codons that encode that residue).
     */
    public Distribution getFrequencyForSynonyms(Symbol residue) 
        throws IllegalSymbolException;

    /**
     * returns a WobbleDistribution for
     * a specified residue.
     */
    public WobbleDistribution getWobbleDistributionForSynonyms(Symbol residue) 
        throws IllegalSymbolException;
}

