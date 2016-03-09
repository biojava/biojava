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
package org.biojava.nbio.alignment.routines;

import org.biojava.nbio.alignment.routines.AlignerHelper.Anchor;
import org.biojava.nbio.alignment.routines.AlignerHelper.Cut;
import org.biojava.nbio.alignment.routines.AlignerHelper.Subproblem;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

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
	@Test
	public void getSubproblems_should_return_score_indicies_of_alignment_subproblems() {
		List<Anchor> anchors = new ArrayList<Anchor>();
		anchors.add(new Anchor(1, 2));
		anchors.add(new Anchor(5, 5));
		List<Subproblem> problems = AlignerHelper.Subproblem.getSubproblems(anchors, 10, 15);
		assertEquals(3, problems.size());
		assertEquals(0, problems.get(0).getQueryStartIndex());
		assertEquals(0, problems.get(0).getTargetStartIndex());
		assertEquals(1, problems.get(0).getQueryEndIndex());
		assertEquals(2, problems.get(0).getTargetEndIndex());
		assertEquals(2, problems.get(1).getQueryStartIndex());
		assertEquals(3, problems.get(1).getTargetStartIndex());
		assertEquals(5, problems.get(1).getQueryEndIndex());
		assertEquals(5, problems.get(1).getTargetEndIndex());
		assertEquals(6, problems.get(2).getQueryStartIndex());
		assertEquals(6, problems.get(2).getTargetStartIndex());
		assertEquals(10, problems.get(2).getQueryEndIndex());
		assertEquals(15, problems.get(2).getTargetEndIndex());
	}
	@Test
	public void getSubproblems_should_allow_zero_anchors() {
		List<Anchor> anchors = new ArrayList<Anchor>();
		List<Subproblem> problems = AlignerHelper.Subproblem.getSubproblems(anchors, 10, 15);
		assertEquals(1, problems.size());
		assertEquals(0, problems.get(0).getQueryStartIndex());
		assertEquals(0, problems.get(0).getTargetStartIndex());
		assertEquals(10, problems.get(0).getQueryEndIndex());
		assertEquals(15, problems.get(0).getTargetEndIndex());
		assertEquals(false, problems.get(0).isStartAnchored());
	}
	@Test
	public void getSubproblems_should_allow_start_and_end_anchors() {
		List<Anchor> anchors = new ArrayList<Anchor>();
		anchors.add(new Anchor(0, 0));
		anchors.add(new Anchor(9, 14));
		List<Subproblem> problems = AlignerHelper.Subproblem.getSubproblems(anchors, 10, 15);
		assertEquals(3, problems.size());
		assertEquals(0, problems.get(0).getQueryStartIndex());
		assertEquals(0, problems.get(0).getTargetStartIndex());
		assertEquals(0, problems.get(0).getQueryEndIndex());
		assertEquals(0, problems.get(0).getTargetEndIndex());
		assertEquals(false, problems.get(0).isStartAnchored());
		assertEquals(1, problems.get(1).getQueryStartIndex());
		assertEquals(1, problems.get(1).getTargetStartIndex());
		assertEquals(9, problems.get(1).getQueryEndIndex());
		assertEquals(14, problems.get(1).getTargetEndIndex());
		assertEquals(true, problems.get(1).isStartAnchored());
		assertEquals(10, problems.get(2).getQueryStartIndex());
		assertEquals(15, problems.get(2).getTargetStartIndex());
		assertEquals(10, problems.get(2).getQueryEndIndex());
		assertEquals(15, problems.get(2).getTargetEndIndex());
		assertEquals(true, problems.get(2).isStartAnchored());
	}
	@Test
	public void getSubproblems_should_allow_adjacent_anchors() {
		List<Anchor> anchors = new ArrayList<Anchor>();
		anchors.add(new Anchor(1, 1));
		anchors.add(new Anchor(2, 3));
		List<Subproblem> problems = AlignerHelper.Subproblem.getSubproblems(anchors, 10, 15);
		assertEquals(3, problems.size());
		assertEquals(2, problems.get(1).getQueryStartIndex());
		assertEquals(2, problems.get(1).getTargetStartIndex());
		assertEquals(2, problems.get(1).getQueryEndIndex());
		assertEquals(3, problems.get(1).getTargetEndIndex());
		assertEquals(3, problems.get(2).getQueryStartIndex());
		assertEquals(4, problems.get(2).getTargetStartIndex());
		assertEquals(10, problems.get(2).getQueryEndIndex());
		assertEquals(15, problems.get(2).getTargetEndIndex());
	}
	@Test(expected=IllegalArgumentException.class)
	public void getSubproblems_should_not_allow_repeated_anchors() {
		List<Anchor> anchors = new ArrayList<Anchor>();
		anchors.add(new Anchor(1, 1));
		anchors.add(new Anchor(1, 2));
		AlignerHelper.Subproblem.getSubproblems(anchors, 10, 15);
	}
	@Test(expected=IllegalArgumentException.class)
	public void getSubproblems_should_not_allow_unalignable_anchors() {
		List<Anchor> anchors = new ArrayList<Anchor>();
		anchors.add(new Anchor(1, 2));
		anchors.add(new Anchor(2, 1));
		AlignerHelper.Subproblem.getSubproblems(anchors, 10, 15);
	}
}
