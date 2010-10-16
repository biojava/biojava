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

package org.biojava.bio.seq;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SymbolList;

/**
 * Collects the references to translation methods in one place.  Right now this
 * is just a wrapper on RNATools
 *
 * @author Greg Cox
 */
public class GeneticCodes
{
// Static variables

// Member variables

// Constructors and initialization

    private GeneticCodes() {
    }

// Interface implementations

// Public methods
	/**
	 * Transcribe DNA into RNA.
	 *
	 * @param theList the SymbolList of DNA symbols to transcribe
	 * @return a SymbolList that is the transcribed view
	 * @throws IllegalAlphabetException if the list is not DNA
	 */
	public static SymbolList transcribe(SymbolList theList)
		throws IllegalAlphabetException
	{
		return RNATools.transcribe(theList);
	}

	/**
	 * Translate RNA into protein (with termination symbols).
	 *
	 * @param theList the SymbolList of RNA symbols to translate
	 * @return a SymbolList that is the translated view
	 * @throws IllegalAlphabetException if the list is not RNA
	 * @since 1.1
	 */
	public static SymbolList translate(SymbolList theList)
		throws IllegalAlphabetException
	{
		SymbolList returnList = null;
		try
		{
			// Assume it's a RNA list
			returnList = RNATools.translate(theList);
		}
		catch(IllegalAlphabetException iae)
		{
			// Since it isn't, lets try DNA.  If this fails, the method fails,
			// so we won't catch the exception
			SymbolList tempList = GeneticCodes.transcribe(theList);
			returnList = RNATools.translate(tempList);
		}
		return returnList;
	}

// Protected methods

// Private methods
}
