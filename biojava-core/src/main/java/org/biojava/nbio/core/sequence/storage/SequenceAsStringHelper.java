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
package org.biojava.nbio.core.sequence.storage;

import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;

import java.util.List;

/**
 * This is a common method that can be used across multiple storage/proxy implementations to
 * handle Negative strand and other interesting elements of sequence data.
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SequenceAsStringHelper<C extends Compound> {

	/**
	 *
	 * @param parsedCompounds
	 * @param compoundSet
	 * @param bioBegin
	 * @param bioEnd
	 * @param strand
	 * @return
	 */
	public String getSequenceAsString(List<C> parsedCompounds, CompoundSet<C> compoundSet, Integer bioBegin, Integer bioEnd, Strand strand) {
		// TODO Optimise/cache.
		if(parsedCompounds.size() == 0)
			return "";
		StringBuilder builder = new StringBuilder();
		if (strand.equals(Strand.NEGATIVE)) {
			//we expect bioBegin to be bigger but could have circular case
			if (bioBegin <= bioEnd) {
				for (int index = bioEnd - 1; index >= bioBegin - 1; index--) {
					C compound = parsedCompounds.get(index);
					builder.append(compoundSet.getStringForCompound(compound));
				}
			}else{
				//go to 0 and the up
				for (int index = bioBegin - 1; index >= 0; index--) {
					C compound = parsedCompounds.get(index);
					builder.append(compoundSet.getStringForCompound(compound));
				}

				for (int index = parsedCompounds.size() - 1; index >= bioEnd - 1; index--) {
					C compound = parsedCompounds.get(index);
					builder.append(compoundSet.getStringForCompound(compound));
				}
			}
		} else {
			if (bioBegin <= bioEnd) {
				for (int index = bioBegin - 1; index <= bioEnd - 1 ; index++) {
					C compound = parsedCompounds.get(index);
					builder.append(compoundSet.getStringForCompound(compound));
				}
			}else{
				//go to 0 and the up
				for (int index = bioBegin - 1; index <=  parsedCompounds.size() - 1; index++) {
					C compound = parsedCompounds.get(index);
					builder.append(compoundSet.getStringForCompound(compound));
				}

				for (int index = 0; index <= bioEnd - 1; index++) {
					C compound = parsedCompounds.get(index);
					builder.append(compoundSet.getStringForCompound(compound));
				}
			}


		}

		return builder.toString();
	}
}
