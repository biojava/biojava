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
package org.biojava.nbio.structure.domain;

import java.io.IOException;


/** A simple factory object that returns the system wide default DomainProvider
 *
 * @author andreas
 *
 */
public class DomainProviderFactory {

	private DomainProviderFactory(){

	}

	static DomainProvider domainProvider ;



	public static void setDomainProvider(DomainProvider provider){
		domainProvider = provider;

	}

	public static DomainProvider getDomainProvider() throws IOException{
		if ( domainProvider == null)
			domainProvider = new RemoteDomainProvider(true);

		return domainProvider;
	}
}
