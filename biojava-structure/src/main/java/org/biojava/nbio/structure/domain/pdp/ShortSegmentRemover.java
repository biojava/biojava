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
package org.biojava.nbio.structure.domain.pdp;

import java.util.List;


public class ShortSegmentRemover {

	public static void cleanup(List<Domain> domains) {

		int ndom = domains.size();

		for(int j=0;j<ndom;j++) {

			int n=0;
			boolean allshort=true;
			// count the length of segments for this domain.
			for (int i=0;i<domains.get(j).nseg;i++) {
				int seglen=(domains.get(j).getSegmentAtPos(i).getTo()-domains.get(j).getSegmentAtPos(i).getFrom()+1);
				if(seglen>=30) allshort = false;
				n+=seglen;
			}

			if(n<PDPParameters.MIN_DOMAIN_LENGTH||allshort) {
				ndom--;
				domains.remove(j);
				j--;

			}
		}
	}
}
