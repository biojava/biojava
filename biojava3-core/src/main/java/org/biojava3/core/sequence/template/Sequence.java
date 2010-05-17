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
 * Created on 01-21-2010
 *
 * @author Richard Holland
 *
 *
 */
package org.biojava3.core.sequence.template;

import java.util.List;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;

public interface Sequence<C extends Compound> extends Iterable<C> {

    public int getLength();

    public C getCompoundAt(int position);

    public int getIndexOf(C compound);

    public int getLastIndexOf(C compound);

    public String getSequenceAsString();

    public String getSequenceAsString(Integer start, Integer end, Strand strand);

    public List<C> getAsList();

    public SequenceView<C> getSubSequence(Integer start, Integer end);

    public CompoundSet<C> getCompoundSet();

    public AccessionID getAccession();

    public int countCompounds(C... compounds);
}
