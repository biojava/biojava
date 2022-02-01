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
package org.biojava.nbio.genome;

import com.google.common.collect.Range;
import org.biojava.nbio.genome.util.ChromosomeMappingTools;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by andreas on 7/19/16.
 */
public class TestGenomeMapping {

	@Test
	public void testGenomeMappingToolGetCDSRanges(){

		List<Integer> lst1 = new ArrayList<>(Arrays.asList( 86346823, 86352858, 86354529));
		List<Integer> lst2 = new ArrayList<>(Arrays.asList(86348878, 86352984, 86354692));

		int cdsStart=86348749, cdsEnd=86387027;

		List<Range<Integer>> result = ChromosomeMappingTools.getCDSRegions(lst1,lst2,cdsStart,cdsEnd);

		// makes sure the first list does not get  changed;
		Assert.assertEquals(86346823, (int) lst1.get(0));


		Assert.assertEquals(86348749, (int) result.get(0).lowerEndpoint());
		Assert.assertEquals(86352858, (int) result.get(1).lowerEndpoint());
		Assert.assertEquals(86354529, (int) result.get(2).lowerEndpoint());

		Assert.assertEquals(86348878, (int) result.get(0).upperEndpoint());
		Assert.assertEquals(86352984, (int) result.get(1).upperEndpoint());
		Assert.assertEquals(86387027, (int) result.get(2).upperEndpoint());

	}

	@Test
	public void testGenomeMappingToolGetCDSRangesSERINC2(){

		List<Integer> lst1 = new ArrayList<>(Arrays.asList(31413812, 31415872, 31423692));
		List<Integer> lst2 = new ArrayList<>(Arrays.asList(31414777, 31415907, 31423854));

		int cdsStart=31423818, cdsEnd=31434199;

		List<Range<Integer>> result = ChromosomeMappingTools.getCDSRegions(lst1,lst2,cdsStart,cdsEnd);

		// makes sure the first list does not get  changed;
		Assert.assertEquals(31423818, (int) result.get(0).lowerEndpoint());

	}
}

