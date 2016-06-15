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
 * Created on 3/1/2010
 * @author Andy Yates
 */
package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.template.*;

/**
 * Bare bones version of the Sequence object to be used sparingly. You should
 * really use a specialized version of {@link Sequence} which describes
 * your domain.
 */
public class BasicSequence<C extends Compound> extends AbstractSequence<C> {

	public BasicSequence(String sequence, CompoundSet<C> compoundSet) throws CompoundNotFoundException {
		super(sequence, compoundSet);
	}

	public BasicSequence(ProxySequenceReader<C> reader) {
		super(reader, reader.getCompoundSet());
	}
}
