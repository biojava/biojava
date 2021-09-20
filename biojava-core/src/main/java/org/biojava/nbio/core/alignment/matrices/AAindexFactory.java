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
package org.biojava.nbio.core.alignment.matrices;

/**
 * Factory class to get Providers for substitution matrices that are provided by
 * the AAINDEX database.
 *
 * @author Andreas Prlic
 *
 */
public class AAindexFactory {

	private static AAIndexProvider provider = null;

	/**
	 * Gets singleton instance of an {@code AAIndexProvider}, always non-null
	 * 
	 * @return
	 */
	public static AAIndexProvider getAAIndexProvider() {
		if (provider == null)
			provider = new DefaultAAIndexProvider();
		return provider;
	}

	/**
	 * Sets the singleton provider.
	 * @param provider
	 */
	public static void setAAIndexProvider(AAIndexProvider provider) {
		AAindexFactory.provider = provider;
	}

}
