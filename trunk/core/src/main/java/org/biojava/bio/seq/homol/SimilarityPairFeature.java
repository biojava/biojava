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

package org.biojava.bio.seq.homol;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

/**
 * <p><code>SimilarityPairFeature</code> describes a pairwise
 * similarity between two nucleotide sequences (as it extends
 * <code>StrandedFeature</code>). It is analagous to, and based on,
 * the BioPerl Bio::SeqFeature::SimilarityPair.</p>
 *
 * <p>It is different from <code>HomologyFeature</code> in that it
 * expresses a relationship between only two sequence regions (rather
 * than >= 2), with one clearly defined as the query sequence and the
 * other as the subject (database hit). These are identified by
 * constant labels in the
 * <code>Alignment</code>. <code>HomologyFeature</code> identifies the
 * related sequence regions by means of an <code>Homology</code>
 * instance which contains an <code>Alignment</code> which uses the
 * <code>HomologyFeature</code>s themselves as labels.</p>
 *
 * <p>In cases where there is no alignment available, for example when
 * MSPCrunch output or GFF have been used, the
 * <code>EmptyPairwiseAlignment</code> in the EMPTY_PAIRWISE field may
 * be used. This may also be useful if an implementation elides the
 * alignment data for some reason.</p>
 *
 * @author Keith James
 * @since 1.2
 */
public interface SimilarityPairFeature extends StrandedFeature
{
    /**
     * The sibling of this feature has altered.
     */
    public static final ChangeType SIBLING =
        new ChangeType("Sibling has altered", SimilarityPairFeature.class, "SIBLING");

    /**
     * Constant <code>QUERY_LABEL</code> is the alignment label used
     * for all query sequences.
     */
    public static final String QUERY_LABEL = "query";

    /**
     * Constant <code>SUBJECT_LABEL</code> is the alignment label used
     * for all subject sequences.
     */
    public static final String SUBJECT_LABEL = "subject";

    /**
     * Constant <code>EMPTY_PAIRWISE</code> is an empty alignment for
     * situations where there is no available alignment data or the
     * implementation does not want to create one.
     */
    public static final Alignment EMPTY_PAIRWISE = new EmptyPairwiseAlignment();

    /**
     * <code>getSibling</code> returns the sibling
     * <code>Feature</code>, query for subject and vice versa.
     *
     * @return a <code>Feature</code>.
     */
    public SimilarityPairFeature getSibling();

    /**
     * <code>setSibling</code> sets the sibling feature of the
     * pair. This is used to set the reciprocal
     * <code>SimilarityPairFeature</code> as both cannot be set using
     * the <code>Template</code>.
     */
    public void setSibling(SimilarityPairFeature sibling)
        throws ChangeVetoException;

    /**
     * <code>getAlignment</code> returns the <code>Alignment</code> of
     * two similar features.
     *
     * @return an <code>Alignment</code> value.
     */
    public Alignment getAlignment();

    /**
     * <code>getScore</code> returns the alignment score.
     *
     * @return a <code>double</code>.
     */
    public double getScore();

    /**
     * <code>Template</code> for construction of
     * <code>SimilarityPairFeature</code>s.
     */
    public static class Template extends StrandedFeature.Template
    {
        /**
         * <code>sibling</code> <code>SimilarityPairFeature</code>
         * field. May be null if the reciprocal
         * <code>SimilarityPairFeature</code> has not yet been
         * created.
         */
        public SimilarityPairFeature sibling;

        /**
         * <code>alignment</code> <code>Alignment</code> field.
         */
        public Alignment alignment;

        /**
         * <code>score</code> of the search which produced the
         * alignment.
         */
        public double score;
    }

    /**
     * <code>EmptyPairwiseAlignment</code> empty pairwise alignment
     * which has labels to empty symbol lists.
     */
    static final class EmptyPairwiseAlignment extends Unchangeable
        implements Alignment
    {
        private List labels = new ArrayList(2);

        EmptyPairwiseAlignment()
        {
            labels = new ArrayList(2);
            labels.add(QUERY_LABEL);
            labels.add(SUBJECT_LABEL);
        }

        public List getLabels()
        {
            return labels;
        }

        public Symbol symbolAt(Object label, int index)
            throws NoSuchElementException
        {
            throw new NoSuchElementException("Attempted to retrieve symbol from empty list at "
                                             + label
                                             + ":"
                                             + index);
        }

        public SymbolList symbolListForLabel(Object label)
            throws NoSuchElementException
        {
            return SymbolList.EMPTY_LIST;
        }

        public Alignment subAlignment(Set labels, Location loc)
            throws NoSuchElementException
        {
            throw new NoSuchElementException("Attempted to retrieve sub-alignment from empty list at "
                                             + labels
                                             + ":"
                                             + loc);
        }

        public int length()
        {
            return 0;
        }

        public Iterator iterator()
        {
            return Collections.EMPTY_LIST.iterator();
        }

        public SymbolList subList(int index1, int index2)
            throws IndexOutOfBoundsException
        {
            Collections.EMPTY_LIST.subList(index1 - 1, index2);
            return SymbolList.EMPTY_LIST;
        }

        public Symbol symbolAt(int index) throws IndexOutOfBoundsException
        {
            throw new IndexOutOfBoundsException("Attempted to retrieve symbol from empty list at "
                                                + index);
        }

        public Alphabet getAlphabet()
        {
            return Alphabet.EMPTY_ALPHABET;
        }

        public List toList()
        {
            return Collections.EMPTY_LIST;
        }

        public String seqString()
        {
            return "";
        }

        public String subStr(int index1, int index2)
            throws IndexOutOfBoundsException
        {
            throw new IndexOutOfBoundsException("You can not retrieve part of an empty symbol list");
        }

        public void edit(Edit edit)
            throws IndexOutOfBoundsException,IllegalAlphabetException,
                   ChangeVetoException
        {
            throw new ChangeVetoException("You can't edit the empty symbol list");
        }
      
        public Iterator symbolListIterator()
        {
            return new Alignment.SymbolListIterator(this);
        }
    }
}
