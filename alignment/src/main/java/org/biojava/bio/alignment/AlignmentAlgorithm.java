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

/*
 * Created on 2005-08-03
 */
package org.biojava.bio.alignment;

import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.SymbolList;

/**
 * This Interface provides methods for the alignment of bio-sequences.
 * 
 * @author Andreas Dr&auml;ger <andreas.draeger@uni-tuebingen.de>
 * @author Mark Schreiber
 */
public abstract class AlignmentAlgorithm {

	/**
	 * @param source
	 *            a SequenceIterator containing a set of sequences to be aligned
	 *            with
	 * @param subjectDB
	 *            the SequenceDB containing another set of sequences.
	 * @return a list containing the results of all single alignments performed
	 *         by this method.
	 * @throws NoSuchElementException
	 * @throws Exception
	 */
	public List<AlignmentPair> alignAll(SequenceIterator source,
			SequenceDB subjectDB) throws Exception {
		List<AlignmentPair> l = new LinkedList<AlignmentPair>();
		while (source.hasNext()) {
			Sequence query = source.nextSequence();
			// compare all the sequences of both sets.
			SequenceIterator target = subjectDB.sequenceIterator();
			while (target.hasNext())
				try {
					l.add(pairwiseAlignment(query, target.nextSequence()));
					// pairwiseAlignment(query, target.nextSequence());
				} catch (Exception exc) {
					exc.printStackTrace();
				}
		}
		return l;
	}

	/**
	 * Performs a pairwise sequence alignment of the two given sequences.
	 * 
	 * @param query
	 * @param subject
	 * @return score of the alignment or the distance.
	 * @throws Exception
	 */
	public abstract AlignmentPair pairwiseAlignment(SymbolList query, SymbolList subject)
			throws Exception;
}
