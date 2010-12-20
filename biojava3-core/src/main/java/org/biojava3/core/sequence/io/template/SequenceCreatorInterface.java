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
 */

package org.biojava3.core.sequence.io.template;

import java.util.List;

import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.ProxySequenceReader;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public interface SequenceCreatorInterface<C extends Compound> {
/**
 *
 * @param sequence
 * @param index
 * @return
 */
    public AbstractSequence<C> getSequence(String sequence, long index);
 /**
  * 
  * @param proxyLoader
  * @param index
  * @return
  */
    public AbstractSequence<C> getSequence(ProxySequenceReader<C> proxyLoader, long index);
 /**
  * 
  * @param list
  * @return
  */
    public AbstractSequence<C> getSequence(List<C> list);

}
