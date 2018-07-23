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
/**
 *
 */
package org.biojava.nbio.genome;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author mckeee1
 *
 */
public class FeatureListTest {
	@Test
	public void testAddIndex() throws Exception
	{
		FeatureList fl = new FeatureList();
		fl.add(new Feature("seqname", "source", "type", new Location(1, 2), (double)0, 0, "gene_id \"gene_id_1\"; transcript_id \"transcript_id_1\";"));
		fl.addIndex("transcript_id");
		Assert.assertEquals(1, fl.selectByAttribute("transcript_id").size());


		FeatureList f2 = new FeatureList();
		f2.addIndex("transcript_id");
		f2.add(new Feature("seqname", "source", "type", new Location(1, 2), (double)0, 0, "gene_id \"gene_id_1\"; transcript_id \"transcript_id_1\";"));
		Assert.assertEquals(1, f2.selectByAttribute("transcript_id").size());
	}
}
