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
package org.biojava.nbio.genome.homology;

import org.biojava.nbio.genome.query.BlastXMLQuery;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class BlastHomologyHits {

	static public LinkedHashMap<String, ArrayList<String>> getMatches(File xmlBlastHits, double ecutoff) throws Exception {
		LinkedHashMap<String, ArrayList<String>> homologyHits = new LinkedHashMap<String, ArrayList<String>>();
		BlastXMLQuery blastXMLQuery = new BlastXMLQuery(xmlBlastHits.getAbsolutePath());
		LinkedHashMap<String, ArrayList<String>> hits = blastXMLQuery.getHitsQueryDef(ecutoff);
		for (String accessionid : hits.keySet()) {
			String[] data = accessionid.split(" "); // deal with notes/comments in blast results
			String id = data[0];
			ArrayList<String> uniprotProteinHits = hits.get(accessionid);
			homologyHits.put(id, uniprotProteinHits);

		}
		return homologyHits;
	}
}
