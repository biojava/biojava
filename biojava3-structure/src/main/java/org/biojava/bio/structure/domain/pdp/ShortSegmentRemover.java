package org.biojava.bio.structure.domain.pdp;

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
