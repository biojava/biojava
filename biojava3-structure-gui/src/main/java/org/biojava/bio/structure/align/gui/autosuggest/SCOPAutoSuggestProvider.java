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
 * Created on Sep 15, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.align.gui.autosuggest;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicBoolean;

import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;

public class SCOPAutoSuggestProvider implements AutoSuggestProvider{

	boolean DEBUG = false;
	
	int maxResults = 20;

	AtomicBoolean stop = new AtomicBoolean(false);

	@Override
	public Vector<String> getSuggestion(String userInput) {

		long timeS = System.currentTimeMillis();

		List<ScopDomain> domains = new ArrayList<ScopDomain>();

		domains = getPossibleScopDomains(userInput);

	

		// convert domains to Strings

		Vector<String> v=new Vector<String>();

		int counter = 0;
		for ( ScopDomain d : domains){
			counter ++;

			String scopId = d.getScopId();
			v.add(scopId);


			if ( counter > maxResults)
				break;
		}

		long timeE = System.currentTimeMillis();
		
		if ( DEBUG)
			System.out.println("ScopAutoSuggestProvider took " + (timeE - timeS) + " ms. to get " + v.size() + " suggestions");
		
		return v;

	}



	private List<ScopDomain> getPossibleScopDomains(String userInput) {
		
		List<ScopDomain> domains = new ArrayList<ScopDomain>();

		ScopDatabase scop = ScopFactory.getSCOP();

		if (userInput.length() ==5 && userInput.startsWith("d") && (! userInput.contains("."))) {
			userInput = userInput.substring(1);
		}
				
		if ( userInput.length() ==4){
			domains = scop.getDomainsForPDB(userInput);

		} else {
			int suni = -1;

			try {
				suni = Integer.parseInt(userInput);
			} catch (NumberFormatException e){
				//supress
			}

			if ( stop.get())
				return domains;

			if ( suni != -1)
				domains = scop.getScopDomainsBySunid(suni);

			if ( stop.get())
				return domains;

			if ( domains == null || domains.size() < 1){

				if ( userInput.length() > 5){
					// e.g. d4hhba

					domains.addAll(scop.filterByDomainName(userInput));

				}
			}

			if ( stop.get())
				return domains;

			if (DEBUG)
				System.out.println("domains: " + domains);
			
			if ( domains == null || domains.size() < 1) {
				if ( userInput.length() > 0 ){
					List<ScopDescription> descs = scop.filterByClassificationId(userInput);

					if ( descs == null || descs.size() < 1){
						descs = scop.filterByDescription(userInput); 
					}


					for (ScopDescription d : descs){
						domains.addAll(scop.getScopDomainsBySunid(d.getSunID()));
						if ( domains.size()> maxResults){
							break;
						}

						if ( stop.get())
							return domains;
					}
				}

			}

		}
		
		return domains;
	}





	@Override
	public void setMaxNrSuggestions(int maxNrSuggestions) {
		maxResults = maxNrSuggestions;

	}

	@Override
	public int getMaxNrSuggestions() {
		return maxResults;
	}

	@Override
	public void clear() {
		// TODO Auto-generated method stub

	}

	@Override
	public void stop() {
		stop.set(true);
		if (DEBUG)
			System.out.println("ScopAutoSuggestProvider got signal stop");

	}



}
