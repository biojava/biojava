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

package org.biojava.bio.seq.io;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * <p><code>GenEmblPropertyComparator</code> compares Genbank/EMBL
 * file format tags by the order in which they should appear in their
 * respective formats.</p>
 *
 * <p>EMBL tags sort before Genbank tags. This is arbitrary. Given the
 * subtle differences in the values accompanying equivalent tags in
 * these formats the two sets shouldn't be mixed anyway.</p>
 *
 * <p>Any tags which belong to neither set sort before anything
 * else.<p>
 *
 * @author Keith James
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
public final class GenEmblPropertyComparator implements Comparator
{
    public static final Comparator INSTANCE = new GenEmblPropertyComparator();

    private List tagOrder;

    private GenEmblPropertyComparator()
    {
        tagOrder = new ArrayList();
        tagOrder.add(EmblLikeFormat.ID_TAG);
        tagOrder.add(EmblLikeFormat.ACCESSION_TAG);
        tagOrder.add(EmblLikeFormat.VERSION_TAG);
        tagOrder.add(EmblLikeFormat.DATE_TAG);
        tagOrder.add(EmblLikeFormat.DEFINITION_TAG);
        tagOrder.add(EmblLikeFormat.KEYWORDS_TAG);
        tagOrder.add(EmblLikeFormat.SOURCE_TAG);
        tagOrder.add(EmblLikeFormat.ORGANISM_TAG);
        /*tagOrder.add(EmblLikeFormat.REFERENCE_TAG);
        tagOrder.add(EmblLikeFormat.COORDINATE_TAG);
        tagOrder.add(EmblLikeFormat.REF_ACCESSION_TAG);
        tagOrder.add(EmblLikeFormat.AUTHORS_TAG);
        tagOrder.add(EmblLikeFormat.TITLE_TAG);
        tagOrder.add(EmblLikeFormat.JOURNAL_TAG);*/
        tagOrder.add(ReferenceAnnotation.class);
        tagOrder.add(EmblLikeFormat.DR_TAG);//lorna:added 21.08.03
        tagOrder.add(EmblLikeFormat.COORDINATE_TAG);
        tagOrder.add(EmblLikeFormat.REF_ACCESSION_TAG);
        tagOrder.add(EmblLikeFormat.AUTHORS_TAG);
        tagOrder.add(EmblLikeFormat.TITLE_TAG);
        tagOrder.add(EmblLikeFormat.JOURNAL_TAG);
        tagOrder.add(EmblLikeFormat.REF_XREF_TAG);
        tagOrder.add(EmblLikeFormat.COMMENT_TAG);
        tagOrder.add(EmblLikeFormat.FEATURE_TAG);

        tagOrder.add(GenbankFormat.LOCUS_TAG);
        tagOrder.add(GenbankFormat.SIZE_TAG);
        tagOrder.add(GenbankFormat.STRAND_NUMBER_TAG);
        tagOrder.add(GenbankFormat.TYPE_TAG);
        tagOrder.add(GenbankFormat.CIRCULAR_TAG);
        tagOrder.add(GenbankFormat.DIVISION_TAG);
        tagOrder.add(GenbankFormat.DATE_TAG);
        tagOrder.add(GenbankFormat.DEFINITION_TAG);
        tagOrder.add(GenbankFormat.ACCESSION_TAG);
        tagOrder.add(GenbankFormat.VERSION_TAG);
        tagOrder.add(GenbankFormat.GI_TAG);
        tagOrder.add(GenbankFormat.KEYWORDS_TAG);
        tagOrder.add(GenbankFormat.SOURCE_TAG);
        tagOrder.add(GenbankFormat.ORGANISM_TAG);
        tagOrder.add(GenbankFormat.REFERENCE_TAG);
        tagOrder.add(GenbankFormat.AUTHORS_TAG);
        tagOrder.add(GenbankFormat.TITLE_TAG);
        tagOrder.add(GenbankFormat.JOURNAL_TAG);
        tagOrder.add(GenbankFormat.PUBMED_TAG);
        tagOrder.add(GenbankFormat.MEDLINE_TAG);
        tagOrder.add(GenbankFormat.COMMENT_TAG);
        tagOrder.add(GenbankFormat.FEATURE_TAG);
    }

    public int compare(Object o1, Object o2)
    {
        int index1 = tagOrder.indexOf(o1);
        int index2 = tagOrder.indexOf(o2);

        return (index1 - index2);
    }
}
