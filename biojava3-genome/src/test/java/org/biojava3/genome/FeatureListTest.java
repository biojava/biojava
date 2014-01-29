/**
 * 
 */
package org.biojava3.genome;

import org.biojava3.genome.parsers.gff.Feature;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.Location;
import org.junit.Test;

import junit.framework.TestCase;

/**
 * @author mckeee1
 *
 */
public class FeatureListTest extends TestCase {
    @Test
    public void testAddIndex() throws Exception
    {
        FeatureList fl = new FeatureList();
        fl.add(new Feature("seqname", "source", "type", new Location(1, 2), (double)0, 0, "gene_id \"gene_id_1\"; transcript_id \"transcript_id_1\";"));
        fl.addIndex("transcript_id");
        assertEquals(1, fl.selectByAttribute("transcript_id").size());

    
        FeatureList f2 = new FeatureList();
        f2.addIndex("transcript_id");
        f2.add(new Feature("seqname", "source", "type", new Location(1, 2), (double)0, 0, "gene_id \"gene_id_1\"; transcript_id \"transcript_id_1\";"));
        assertEquals(1, f2.selectByAttribute("transcript_id").size());
    }
}
