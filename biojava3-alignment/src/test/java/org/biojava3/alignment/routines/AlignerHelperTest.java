package org.biojava3.alignment.routines;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.biojava3.alignment.routines.AlignerHelper.Cut;
import org.biojava3.alignment.routines.AlignerHelper.Subproblem;
import org.junit.Test;

/**
 * 
 * @author Daniel Cameron
 *
 */
public class AlignerHelperTest {
	@Test
	public void getCuts_should_not_return_start_position_for_starting_anchor() {
		Cut[] cuts = AlignerHelper.getCuts(10, new Subproblem(5, 20, 10, 30), new int[] { 50, 50, 3 }, true);
		assertNotEquals(5, cuts[0].getQueryIndex());
	}
	@Test
	public void getCuts_should_return_all_positions_when_cuts_exceeds_query_size() {
		Cut[] cuts = AlignerHelper.getCuts(10, new Subproblem(5, 20, 10, 30), new int[] { 50, 50, 3 }, false);
		assertEquals(5, cuts.length);
		assertEquals(5, cuts[0].getQueryIndex());
		assertEquals(6, cuts[1].getQueryIndex());
		assertEquals(7, cuts[2].getQueryIndex());
		assertEquals(8, cuts[3].getQueryIndex());
		assertEquals(9, cuts[4].getQueryIndex());
	}
	@Test
	public void getCuts_should_return_spaced_cuts_when_query_interval_larger_than_cut_size() {
		Cut[] cuts = AlignerHelper.getCuts(3, new Subproblem(5, 20, 10, 30), new int[] { 50, 50, 3 }, false);
		assertEquals(3, cuts.length);
		assertEquals(5, cuts[0].getQueryIndex());
		assertEquals(7, cuts[1].getQueryIndex());
		assertEquals(9, cuts[2].getQueryIndex());
	}
}
