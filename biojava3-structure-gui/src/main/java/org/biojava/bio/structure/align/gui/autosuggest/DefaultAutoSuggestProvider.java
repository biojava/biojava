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

public class DefaultAutoSuggestProvider implements AutoSuggestProvider {

	@Override
	public Vector<String> getSuggestion(String userInput) {
		Vector<String> data = new Vector<String>();
		
		data.add(userInput + " no AutoSuggestProvider registered yet!");
		return data;
	}

	@Override
	public void setMaxNrSuggestions(int maxNrSuggestions) {
		// TODO Auto-generated method stub

	}

	@Override
	public int getMaxNrSuggestions() {
		// TODO Auto-generated method stub
		return 0;
	}
	
	@Override
	public void clear(){
		
	}
	
	
	public void stop(){
		
	}
}
