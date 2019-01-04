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
 * Created on DATE
 *
 */
package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.sequence.template.AbstractSequence;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Used to sort sequences
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SequenceComparator implements Comparator<AbstractSequence<?>>, Serializable{
	private static final long serialVersionUID = 1;

	@Override
	public int compare(AbstractSequence<?> o1, AbstractSequence<?> o2) {
		return o1.getBioBegin() - o2.getBioBegin();
	}

}
