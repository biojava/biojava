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
 */

/*
 * @author Lorna Morris
 * @since 1.3.1
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
package org.biojava.bio.seq.io;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class EmblReferenceComparator implements Comparator {

    static final Comparator INSTANCE = new EmblReferenceComparator();

    private List tagOrder;

    {
        tagOrder = new ArrayList();
        tagOrder.add(EmblLikeFormat.REFERENCE_TAG);
        tagOrder.add(EmblLikeFormat.COORDINATE_TAG);
        tagOrder.add(EmblLikeFormat.REF_ACCESSION_TAG);
        tagOrder.add(EmblLikeFormat.AUTHORS_TAG);
        tagOrder.add(EmblLikeFormat.TITLE_TAG);
        tagOrder.add(EmblLikeFormat.JOURNAL_TAG);
        tagOrder.add(EmblLikeFormat.REF_XREF_TAG);
        tagOrder.add(EmblLikeFormat.SEPARATOR_TAG);
    }

    public int compare(Object o1, Object o2)
    {
        int index1 = tagOrder.indexOf(o1);
        int index2 = tagOrder.indexOf(o2);

        return (index1 - index2);
    }

}
