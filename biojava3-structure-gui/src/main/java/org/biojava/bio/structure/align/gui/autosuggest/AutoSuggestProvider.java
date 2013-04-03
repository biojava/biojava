/**
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
 * Created on Sep 14, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.align.gui.autosuggest;

import java.util.Vector;

/** A class that provides auto-completion suggestions for JAutoSuggest
 * 
 * @author AndreasPrlic
 *
 */
public interface AutoSuggestProvider {
	
	/** get a list of suggestions for a userInput
	 * 
	 * @param userInput
	 * @return list of suggestions
	 */
	public Vector<String> getSuggestion(String userInput);
	
	
	/** set the maximum number of suggestions to return
	 * 
	 * @param maxNrSuggestions
	 */
	public void setMaxNrSuggestions(int maxNrSuggestions);

	
	/** Get the maximun nr of suggestions 
	 * 
	 * @return maxNrSuggestions
	 */
	public int getMaxNrSuggestions();
	
	
	/** reset all suggestions
	 * 
	 */
	public void clear();
	
	
	/** Interrupt searching for suggestions
	 * 
	 */
	public void stop();
}
