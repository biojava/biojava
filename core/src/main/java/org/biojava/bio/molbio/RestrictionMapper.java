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

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceAnnotator;
import org.biojava.bio.seq.impl.ViewSequence;
import org.biojava.utils.ThreadPool;

/**
 * <p><code>RestrictionMapper</code> is a class for annotating
 * <code>Sequence</code>s with <code>Feature</code>s which represent
 * restriction sites. Calling <code>annotate(Sequence sequence)</code>
 * will annotate the <code>Sequence</code> with the sites of any
 * <code>RestrictionEnzyme</code>s which have been added to the
 * <code>RestrictionMapper</code>. The returned <code>Sequence</code>
 * is a <code>ViewSequence</code> wrapping the original.</p>
 *
 * <p>The <code>Feature</code>s created are
 * <code>RestrictionSite</code>s which have a flyweight
 * <code>Annotation</code> containing a single <code>String</code>
 * property "dbxref" whose value is "REBASE:" plus name of the enzyme
 * (e.g. EcoRI).</p>
 *
 * <p>The mapper will by default map only those sites which have both
 * their recognition sites and their cut sites within the
 * <code>Sequence</code>. This behaviour may be changed to map all
 * sites which have their recognition sites within the
 * <code>Sequence</code> using the <code>setMapAll(boolean
 * on)</code> method.</p>
 *
 * <p>The current implementation requires that
 * <code>RestrictionEnzyme</code>s to be searched must first be
 * registered with the <code>RestrictionEnzymeManager</code>.</p>
 *
 * @author Keith James
 * @since 1.3
 */
public class RestrictionMapper implements SequenceAnnotator
{
    /**
     * <code>SITE_FEATURE_SOURCE</code> the source <code>String</code>
     * used by <code>RestrictionMapper</code> when creating
     * restriction site <code>Feature</code>s. This is the
     * <code>String</code> which is returned when a
     * <code>Feature</code>'s <code>getSource()</code> method is
     * called.
     */
    public static final String SITE_FEATURE_SOURCE = "regex";

    /**
     * <code>SITE_FEATURE_TYPE</code> the type <code>String</code>
     * used by <code>RestrictionMapper</code> when creating
     * restriction site <code>Feature</code>s. This is the
     * <code>String</code> which is returned when a
     * <code>Feature</code>'s <code>getType()</code> method is called.
     */
    public static final String SITE_FEATURE_TYPE = "misc_binding";

    private List restrictionEnzymes;
    private boolean mapAll;
    private ThreadPool threadPool;

    /**
     * <p>Creates a new <code>RestrictionMapper</code> which will use
     * the specified <code>ThreadPool</code>. Do not share one pool
     * between a number of <code>RestrictionMapper</code>s because
     * <code>annotate(Sequence sequence)</code> waits for all threads
     * in the pool to finish work before returning and this will lead
     * to a race condition between mappers. One mapper could end up
     * waiting for another mapper's threads before returning.</p>
     *
     * @param threadPool a <code>ThreadPool</code>.
     */
    public RestrictionMapper(ThreadPool threadPool)
    {
        restrictionEnzymes = new ArrayList();
        mapAll = false;
        this.threadPool = threadPool;
    }

    /**
     * <code>annotate</code> adds <code>Feature</code>s which
     * represent restriction sites.
     *
     * @param sequence a <code>Sequence</code>.
     *
     * @return a <code>Sequence</code> view with restriction sites
     * marked.
     */
    public Sequence annotate(Sequence sequence)
    {
        Sequence mapped = new ViewSequence(sequence);

        for (int i = 0; i < restrictionEnzymes.size(); i++)
        {
            RestrictionEnzyme enzyme =
                (RestrictionEnzyme) restrictionEnzymes.get(i);
            threadPool.addRequest(new RestrictionSiteFinder(enzyme,
                                                            mapAll,
                                                            mapped));
        }

        // Threads will finish work and become idle
        threadPool.waitForThreads();

        return mapped;
    }

    /**
     * <code>getMapAll</code> returns whether all sites should be
     * marked, including those which have recognition sites within the
     * sequence, but cut outside it. The default is false, indicating
     * only sites which can actually be cut are mapped.
     *
     * @return a <code>boolean</code>.
     */
    public boolean getMapAll()
    {
        return mapAll;
    }

    /**
     * <code>setMapAll</code> sets whether all sites should be marked,
     * including those which have recognition sites within the
     * sequence, but cut outside it. The default is false, indicating
     * only sites which can actually be cut are mapped.
     *
     * @param on a <code>boolean</code>.
     */
    public void setMapAll(boolean on)
    {
        mapAll = on;
    }

    /**
     * <code>addEnzyme</code> adds an enzyme to be searched for in the
     * <code>Sequence</code>.
     *
     * @param enzyme a <code>RestrictionEnzyme</code>.
     */
    public void addEnzyme(RestrictionEnzyme enzyme)
    {
        if (restrictionEnzymes.contains(enzyme))
            throw new IllegalArgumentException("RestrictionMapper is already mapping '"
                                               + enzyme
                                               + "'");
        restrictionEnzymes.add(enzyme);
    }

    /**
     * <code>removeEnzyme</code> removes an enzyme from those to be
     * searched for in the <code>Sequence</code>.
     *
     * @param enzyme a <code>RestrictionEnzyme</code>.
     */
    public void removeEnzyme(RestrictionEnzyme enzyme)
    {
        if (! restrictionEnzymes.contains(enzyme))
            throw new IllegalArgumentException("RestrictionMapper is not mapping '"
                                               + enzyme
                                               + "'");

        restrictionEnzymes.remove(enzyme);
    }

    /**
     * <code>clearEnzymes</code> removes all enzymes from those to be
     * searched for in the <code>Sequence</code>.
     */
    public void clearEnzymes()
    {
        restrictionEnzymes.clear();
    }
}
