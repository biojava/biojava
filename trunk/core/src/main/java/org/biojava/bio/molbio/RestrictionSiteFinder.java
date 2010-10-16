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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.SymbolListCharSequence;
import org.biojava.bio.symbol.RangeLocation;

/**
 * <code>RestrictionSiteFinder</code>s do the work of finding sites
 * for one <code>RestrictionEnzyme</code> in a target
 * <code>Sequence</code>. Instances are passed to a
 * <code>ThreadPool</code> in order to perform several concurrent
 * searches.
 *
 * @author Keith James
 * @since 1.3
 */
class RestrictionSiteFinder implements Runnable
{
    private Sequence target;
    private boolean findAll;
    private RestrictionEnzyme enzyme;

    /**
     * Creates a new <code>RestrictionSiteFinder</code>.
     *
     * @param enzyme a <code>RestrictionEnzyme</code> for which to
     * find sites.
     * @param findAll a <code>boolean</code> indicating whether all
     * sites should be found, including those which have recognition
     * sites within the sequence, but cut outside it.
     * @param target a <code>Sequence</code> to search.
     */
    RestrictionSiteFinder(RestrictionEnzyme enzyme,
                          boolean           findAll,
                          Sequence          target)
    {
        this.enzyme  = enzyme;
        this.findAll = findAll;
        this.target  = target;
    }

    /**
     * <code>run</code> searches for restriction sites.
     */
    public void run()
    {
        SymbolListCharSequence charSeq = new SymbolListCharSequence(target);

        try
        {
            Pattern [] patterns = RestrictionEnzymeManager.getPatterns(enzyme);

            int siteLen = enzyme.getRecognitionSite().length();
            int seqLen  = target.length();
            int usOffset = 0;
            int dsOffset = 0;

            int [] dsCut = enzyme.getDownstreamCut();
            dsOffset = Math.max(dsCut[0], dsCut[1]);

            if (enzyme.getCutType() == RestrictionEnzyme.CUT_COMPOUND)
            {
                // In coordinate space of recognition site, so
                // upstream coordinates are negative
                int [] usCut = enzyme.getUpstreamCut();
                usOffset = Math.min(usCut[0], usCut[1]);
            }

            RestrictionSite.Template t = new RestrictionSite.Template();
            t.type       = RestrictionMapper.SITE_FEATURE_TYPE;
            t.source     = RestrictionMapper.SITE_FEATURE_SOURCE;
            t.strand     = StrandedFeature.POSITIVE;
            t.annotation = RestrictionEnzymeManager.getAnnotation(enzyme);
            t.enzyme     = enzyme;

            Matcher m = patterns[0].matcher(charSeq);
            while (m.find())
            {
                int idx = m.start() + 1;

                // Cuts outside target sequence
                if (! findAll && (idx + usOffset < 0 || idx + dsOffset > seqLen))
                    continue;

                t.location = new RangeLocation(idx, idx + siteLen - 1);
                synchronized(target){
                    target.createFeature(t);
                }
            }

            // If not palindromic we have to search reverse strand too
            if (! enzyme.isPalindromic())
            {
                t.strand = StrandedFeature.NEGATIVE;
                m = patterns[1].matcher(charSeq);

                while (m.find())
                {
                    int idx = m.start() + 1;

                    // Cuts outside target sequence
                    if (! findAll && (idx + usOffset < 0 || idx + dsOffset > seqLen))
                        continue;

                    t.location = new RangeLocation(idx, idx + siteLen - 1);
                    synchronized(target){
                        target.createFeature(t);
                    }
                }
            }
        }
        catch (Exception e)
        {
            throw new BioRuntimeException("Failed to complete search for "
                                          + enzyme,e);
        }
    }
}
