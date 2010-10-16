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

package org.biojava.bio.search;

import java.io.Serializable;

import org.biojava.utils.TriState;
import org.biojava.utils.walker.WalkerFactory;
/**
 * A SearchContentHandler class that implements filtering
 * in chains of SearchContentHandler instances.
 * <p>
 * The SearchContentHandler organise Blast-like searches
 * as a hierarchy of search/hit/subhit.  Each search
 * is conducted with a single query sequence.  Hits
 * of the query sequence are reported against different
 * target sequences.  The hit is further subdivided into
 * one of more subhits which represent the positions
 * within the target sequence that alignments of the query
 * sequence were achieved against the query sequence
 * (e.g. HSPs).
 * <p>
 * This implementation depends on the a well ordered use
 * of the SearchContentHandler interface.  In particular,
 * it requires that search/hit/subhit properties are
 * reported immediately following the associated 
 * startSearch/startHit/startSubHit call.  For example,
 * search properties should not be reported following
 * the corresponding startHit() call.
 * <p>
 * <u>Semantics of this interface</u><br>
 * BlastLikeSearchFilters test different levels of the
 * SearchContentHandler property hierarchy and each
 * filter should be seen as being applied when a full
 * set of events from that level is received.  So the
 * ByHitProperty filter is applied when endHit() is
 * called and determines whether all events received
 * between startHit() and endHit() are to be passed on
 * or discarded.
 * <p>
 * <u>Some keys used by SearchContentHandlers</u><br>
 * <u>SearchProperties</u><br>
 * <table border="1">
 * <tr>
 * <td>KEY_QUERY_ID</td>
 * <td>String. Value from setQueryID</td>
 * </tr>
 * <tr>
 * <td>queryDescription</td>
 * <td>String. FASTA description line</td>
 * </tr>
 * <tr>
 * <td>program</td>
 * <td>String. variant of BLAST used</td>
 * </tr>
 * <tr>
 * <td>version</td>
 * <td>software version</td>
 * </tr>
 * </table>
 * <p>
 * <u>HitProperties</u><br>
 * <table border="1" > 
* <tr>
 * <td>subjectId</td>
 * <td>String.  Identity of subject (target) sequence.</td>
 * </tr>
 * <tr>
 * <td>subjectSequenceLength</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>subjectDescription</td>
 * <td>String.</td>
 * </tr>
 * <tr>
 * <td></td>
 * <td></td>
 * </tr>
 * </table>
 * <p>
 * <u>SubHitProperties</u><br>
 * <table border="1" >
 * <tr>
 * <td>bitScore</td>
 * <td>String representation of real value</td>
 * </tr>
 * <tr>
 * <td>queryStrand</td>
 * <td>plus/minus</td>
 * </tr>
 * <tr>
 * <td>percentageIdentity</td>
 * <td>String representation of real value</td>
 * </tr>
 * <tr>
 * <td>querySequenceEnd</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>expectValue</td>
 * <td>String representation of real value</td>
 * </tr>
 * <tr>
 * <td>subjectStrand</td>
 * <td>plus/minus</td>
 * </tr>
 * <tr>
 * <td>subjectSequenceEnd</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>numberOfPositives</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>score</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>subjectSequence</td>
 * <td>String representation of sequence</td>
 * </tr>
 * <tr>
 * <td>alignmentSize</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>querySequenceStart</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>subjectSequenceStart</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>numberOfIdentities</td>
 * <td>String representation of integer value</td>
 * </tr>
 * <tr>
 * <td>querySequence</td>
 * <td>String representation of sequence</td>
 * </tr>
 * </table>
 *
 *
 * @author David Huen
 */
public interface BlastLikeSearchFilter
    extends Serializable
{
    public interface Node
    {
        public Object getSearchProperty(Object key);
        public Object getHitProperty(Object key);
        public Object getSubHitProperty(Object key);
    }

    public static final String KEY_QUERY_ID = "___QUERY_ID___";

    /**
     * returns a TriState indicating the current outcome
     * of evaluating this filter.  This is usually the
     * outcome saved when evaluate(FilteringContentHandler fch) was called.
     */
    public TriState accept();

    /**
     * computes the outcome of this filter on the 
     * specified node and stores it.  <b>This method
     * is only exposed to permit it to be included
     * in an interface.  Users should not use it.</b>
     */
    public void evaluate(Node fch);

    /**
     * resets the internal state of this filter including
     * any cached evaluations. <b>This method
     * is only exposed to permit it to be included
     * in an interface.  Users should not use it.</b>
     */
    public void reset();

    public abstract static class AbstractBlastLikeSearchFilter
        implements BlastLikeSearchFilter
    {
        protected TriState cachedOutcome = TriState.INDETERMINATE;
        public TriState accept() { return cachedOutcome; }
        abstract public void evaluate(Node fch);
        public void reset() { cachedOutcome = TriState.INDETERMINATE; }

        private AbstractBlastLikeSearchFilter() {}
    }

    public static final class And
    {
    static { WalkerFactory.getInstance().addTypeWithParent(And.class); }

        private AbstractBlastLikeSearchFilter filter0;
        private AbstractBlastLikeSearchFilter filter1;

        public And(
            AbstractBlastLikeSearchFilter filter0,
            AbstractBlastLikeSearchFilter filter1)
        {
            this.filter0 = filter0;
            this.filter1 = filter1;
        }

        public TriState accept()
        {
            TriState outcome0 = filter0.accept();
            TriState outcome1 = filter1.accept();

            if ((outcome0 == TriState.FALSE) || (outcome1 == TriState.FALSE))
                return TriState.FALSE;

            // neither can be false now
            if ((outcome0 == TriState.INDETERMINATE) || (outcome1 == TriState.INDETERMINATE))
                return TriState.INDETERMINATE;

            // neither is false nor indeterminate so it must be true!
            return TriState.TRUE;
        }
    }

    public static final class Or
    {
    static { WalkerFactory.getInstance().addTypeWithParent(Or.class); }

        private AbstractBlastLikeSearchFilter filter0;
        private AbstractBlastLikeSearchFilter filter1;

        public Or(
            AbstractBlastLikeSearchFilter filter0,
            AbstractBlastLikeSearchFilter filter1)
        {
            this.filter0 = filter0;
            this.filter1 = filter1;
        }

        public TriState accept()
        {
            TriState outcome0 = filter0.accept();
            TriState outcome1 = filter1.accept();

            if ((outcome0 == TriState.TRUE) || (outcome1 == TriState.TRUE))
                return TriState.TRUE;

            // neither can be false now
            if ((outcome0 == TriState.INDETERMINATE) || (outcome1 == TriState.INDETERMINATE))
                return TriState.INDETERMINATE;

            // neither is true nor indeterminate so it must be false!
            return TriState.FALSE;
        }
    }

    public static final class Not
        extends AbstractBlastLikeSearchFilter
    {
    static { WalkerFactory.getInstance().addTypeWithParent(Not.class); }

        private AbstractBlastLikeSearchFilter filter;

        public Not(AbstractBlastLikeSearchFilter filter)
        {
            this.filter = filter;
        }

        public TriState accept()
        {
            TriState outcome = filter.accept();

            if (outcome == TriState.INDETERMINATE)
                return TriState.INDETERMINATE;

            if (outcome == TriState.TRUE)
                return TriState.FALSE;
            else
                return TriState.TRUE;
        }

        public void evaluate(Node fch) {}
    }

    /**
     * Applies test to the value specified by the key in search properties.
     */
    public static final class BySearchProperty
        extends AbstractBlastLikeSearchFilter
    {
        private Object key;
        private FilterTest test;
        public BySearchProperty(String key, FilterTest test)
        {
            this.key = key;
            this.test = test;
        }

        public void evaluate(Node fch)
        {
            Object propertyValue = fch.getSearchProperty(key);

            cachedOutcome = ((propertyValue != null) && test.accept(propertyValue)) ? TriState.TRUE : TriState.FALSE;
        }
    }

    /**
     * Applies test to the value specified by the key in hit properties.
     */
    public static final class ByHitProperty
        extends AbstractBlastLikeSearchFilter
    {
        private Object key;
        private FilterTest test;
        public ByHitProperty(String key, FilterTest test)
        {
            this.key = key;
            this.test = test;
        }

        public void evaluate(Node fch)
        {
            Object propertyValue = fch.getHitProperty(key);

            cachedOutcome = ((propertyValue != null) && test.accept(propertyValue)) ? TriState.TRUE : TriState.FALSE;
        }
    }

    /**
     * Applies test to the value specified by the key in subhit properties.
     */
    public static final class BySubHitProperty
        extends AbstractBlastLikeSearchFilter
    {
        private Object key;
        private FilterTest test;
        public BySubHitProperty(String key, FilterTest test)
        {
            this.key = key;
            this.test = test;
        }

        public void evaluate(Node fch)
        {
            Object propertyValue = fch.getSubHitProperty(key);
            //cachedOutcome = ((propertyValue == null)
            //    ? TriState.INDETERMINATE 
            //    : (test.accept(propertyValue)) ? TriState.TRUE : TriState.FALSE);
            cachedOutcome = ((propertyValue != null) && test.accept(propertyValue)) ? TriState.TRUE : TriState.FALSE;
        }
    }
}

