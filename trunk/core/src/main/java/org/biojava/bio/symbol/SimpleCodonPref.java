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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.IndexedCount;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.RNATools;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeVetoException;

/**
 * a simple no-frills implementation of the
 * CodonPref object that encapsulates
 * codon preference data.
 *
 * @author David Huen
 * @author gwaldon pyrrolysine
 * @since 1.3
 */
public class SimpleCodonPref
    extends AbstractChangeable
    implements CodonPref
{
    String name;
    String geneticCodeName;
    Distribution codonPref;

    // residue-based codon preference stats
    Map codonPrefByResidue = null;

    // codon wobble-based codon preference stats
    Map wobbleDistributions;

    public SimpleCodonPref(String geneticCodeName, Distribution codonPref, String name)
        throws IllegalAlphabetException
    {
        this.name = name;
        this.geneticCodeName = geneticCodeName;
        this.codonPref = codonPref;

        // validate the Distribution
        if (codonPref.getAlphabet() != RNATools.getCodonAlphabet())
            throw new IllegalAlphabetException("codon preferences must be over codon alphabet");
    }

    public String getName()
    {
        return name;
    }

    public String getGeneticCodeName()
    {
        return geneticCodeName;
    }

    public ManyToOneTranslationTable getGeneticCode()
    {
        return RNATools.getGeneticCode(geneticCodeName);
    }

    public Distribution getFrequency()
    {
        return codonPref;
    }

    public Distribution getFrequencyForSynonyms(Symbol residue)
        throws IllegalSymbolException
    {
        if (codonPrefByResidue == null) preparePrefsByResidue();

        return (Distribution) codonPrefByResidue.get(residue);
    }

    public WobbleDistribution getWobbleDistributionForSynonyms(Symbol residue)
        throws IllegalSymbolException
    {
        if (wobbleDistributions == null) preparePrefsByWobble();

        return (WobbleDistribution) wobbleDistributions.get(residue);
    }

    private void preparePrefsByResidue()
        throws IllegalSymbolException
    {
        try {
            codonPrefByResidue = new HashMap();

            // what we want is to create residue-specific distributions

            for (Iterator residueI = ProteinTools.getTAlphabet().iterator(); residueI.hasNext(); ) {
                Symbol residue = (Symbol) residueI.next();

                
                // filter out selenocysteine!
                if (residue.getName().equals("SEC")) continue;
                // filter out pyrrolysine!
                if (residue.getName().equals("PYL")) continue;
                
                // get the synonymous codons and sum their frequencies
                double residueFreq = 0.0;
                Set synonyms = getGeneticCode().untranslate(residue);

                for (Iterator synonymI = synonyms.iterator(); synonymI.hasNext(); ) {
                     Symbol synonym = (Symbol) synonymI.next();

                    // sum frequency of synonyms for this residue
                    residueFreq += codonPref.getWeight(synonym);
                }

                // now create a new distribution over the synonyms
                Distribution residueCodonDist = DistributionFactory.DEFAULT.createDistribution(RNATools.getCodonAlphabet());

                for (Iterator synonymI = synonyms.iterator(); synonymI.hasNext(); ) {
                    Symbol synonym = (Symbol) synonymI.next();

                    // compute the probability of the current codon
                    residueCodonDist.setWeight(synonym, codonPref.getWeight(synonym)/residueFreq);            
                }

                // lock the Distribution and stash in map for later use
                residueCodonDist.addChangeListener(ChangeListener.ALWAYS_VETO);
                codonPrefByResidue.put(residue, residueCodonDist);
            }
        }
        catch (ChangeVetoException cve) {}
        catch (IllegalAlphabetException iae) {} // none of these should be thrown since the alphabet was preverified.
    }

    private void preparePrefsByWobble()
        throws IllegalSymbolException
    {
        try {
            wobbleDistributions = new HashMap();

            // what we want is to create residue-specific distributions
            FiniteAlphabet nonWobbleAlfa = CodonPrefTools.getDinucleotideAlphabet();

            for (Iterator residueI = ProteinTools.getTAlphabet().iterator(); residueI.hasNext(); ) {
                Symbol residue = (Symbol) residueI.next();

                // filter out selenocysteine!
                if (residue.getName().equals("SEC")) continue;
                // filter out pyrrolysine!
                if (residue.getName().equals("PYL")) continue;
                
                // create bins keyed on non-wobble bases
                IndexedCount nonWobbleCounts = new IndexedCount(nonWobbleAlfa);;
                IndexedCount wobbleCounts;

                Map wobbleDists = new HashMap();
                Set nonWobbleBases = new HashSet();

                // get the synonymous codons
                Set synonyms = getGeneticCode().untranslate(residue);

                for (Iterator synonymI = synonyms.iterator(); synonymI.hasNext(); ) {
                     BasisSymbol synonym = (BasisSymbol) synonymI.next();

                     // retrieve the non-wobble bases for these codons
                     List codonSymbols = synonym.getSymbols();
                     AtomicSymbol wobble = (AtomicSymbol) codonSymbols.get(2);

                     List nonWobbleSymbols = new ArrayList(2);
                     nonWobbleSymbols.add(codonSymbols.get(0));
                     nonWobbleSymbols.add(codonSymbols.get(1));
                     AtomicSymbol nonWobble = (AtomicSymbol) nonWobbleAlfa.getSymbol(nonWobbleSymbols);
                     nonWobbleBases.add(nonWobble);

                     // add counts to the appropriate Count objects
                     double codonFreq = codonPref.getWeight(synonym);
                     nonWobbleCounts.increaseCount(nonWobble, codonFreq);

                     // add Counts
                     wobbleCounts = (IndexedCount) wobbleDists.get(nonWobble);
                     if (wobbleCounts == null) {
                         wobbleCounts = new IndexedCount(RNATools.getRNA());
                         wobbleDists.put(nonWobble, wobbleCounts);
                     }

                     wobbleCounts.increaseCount(wobble, codonFreq);
                }

                // convert the accumulated Counts into Distributions
                Distribution nonWobbleDist = DistributionTools.countToDistribution(nonWobbleCounts);

                for (Iterator nonWobbleBasesI = nonWobbleBases.iterator(); nonWobbleBasesI.hasNext(); ) {
                    AtomicSymbol nonWobbleBase = (AtomicSymbol) nonWobbleBasesI.next();

                    // retrieve and replace each Count with its corresponding Distribution
                    IndexedCount count = (IndexedCount) wobbleDists.get(nonWobbleBase);

                    if (count != null) {
                        Distribution wobbleDist = DistributionTools.countToDistribution(count);
                        wobbleDists.put(nonWobbleBase, wobbleDist);
                    }
                }
                wobbleDistributions.put(residue, new SimpleWobbleDistribution(residue, nonWobbleBases, nonWobbleDist, wobbleDists));
            }
        }
        catch (IllegalAlphabetException iae) {System.err.println("unexpected IllegalAlphabetException"); }
        catch (ChangeVetoException cve) {System.err.println("unexpected ChangeVetoException"); }
    }
}

