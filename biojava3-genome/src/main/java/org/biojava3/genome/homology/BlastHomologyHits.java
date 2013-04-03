/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.genome.homology;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import org.biojava3.genome.query.BlastXMLQuery;

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
