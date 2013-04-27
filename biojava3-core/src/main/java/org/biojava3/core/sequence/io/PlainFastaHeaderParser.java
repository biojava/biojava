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
package org.biojava3.core.sequence.io;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 * The plain fasta header takes everything in the header as a single entity. It
 * is useful for non-standard header formats that don't follow a single rule.<br>
 * If the user has a custom header with local data that is kept constant all
 * over the data then they can create their own implementation of a
 * FastaHeaderParserInterface
 * 
 * @author Amr AL-Hossary
 * @since 3.0.6
 */
public class PlainFastaHeaderParser<S extends AbstractSequence<C>, C extends Compound>
		implements FastaHeaderParserInterface<S, C> {

	/**
	 * Parse out the all header as one entity
	 * 
	 * @param header
	 * @return
	 */
	private String[] getHeaderValues(String header) {
		return new String[] { header };
	}

	/**
	 * Parse the header and set the values in the sequence
	 * 
	 * @param header
	 * @param sequence
	 */
	@Override
	public void parseHeader(String header, S sequence) {
		sequence.setOriginalHeader(header);
		String[] data = getHeaderValues(header);

		if (data.length == 1) {
			sequence.setAccession(new AccessionID(data[0]));
		} else {
			throw new RuntimeException(
					"No header or Some Error Occurred while reading header");
		}
	}
}
