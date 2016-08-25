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

import static org.junit.Assert.*;

import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;

public class TestLocation {

	@Test
	public void testLocation() {
		// tests taken from Location.main()

		//Location p3_7= new Location( 3, 7 );
		Location p16_19= new Location( 16, 19 );
		Location p15_19= new Location( 15, 19 );
		Location p15_16= new Location( 15, 16 );
		Location p10_17= new Location( 10, 17 );
		Location p10_12= new Location( 10, 12 );
		Location p14_17= new Location( 14, 17 );
		//Location p14_14= new Location( 14, 14 );

		Location r13_17= new Location( 13, 17 );
		Location r21_25= new Location( 21, 25 );

		Location r4_7= new Location( 4, 7 );
		Location r2_5= new Location( 2, 5 );
		Location r0_3= new Location( 0, 3 );
		Location r5_8= new Location( 5, 8 );

		//distance
		assertEquals(7, L(14,14).distance( L(3,7) ));
		assertEquals(7, L(3,7).distance( L(14,14) ));
		assertEquals(3, L(1,4).distance( L(7, 10) ));

		//union
		assertEquals(p10_17, p10_12.union( p14_17 ));
		assertEquals(p10_17, p14_17.union( p10_12 ));
		assertEquals(p15_19, p15_19.union( p15_16 ));

		//intersection
		assertEquals(new Location( 21, 25 ), r13_17.union( r21_25 ).intersection( r21_25 ));


		//isBefore
		assertTrue( r2_5.isBefore( r5_8 ));
		assertTrue( !r2_5.isBefore( r4_7 ));

		//isAfter
		assertTrue(r5_8.isAfter( r2_5 ));
		assertTrue(!r5_8.isAfter( r4_7 ));

		//contains
		assertTrue(p15_19.contains( p16_19 ));

		//overlaps
		assertTrue(r2_5.overlaps( r4_7 ));
		assertTrue(r2_5.overlaps( r0_3 ));
		assertTrue(!r5_8.overlaps( r2_5 ));
		assertTrue(!r2_5.overlaps( r5_8 ));


		//prefix
		assertEquals(L(2,3), L(2,20).prefix(1));
		assertEquals(L(2,19), L(2,20).prefix(-1));
		assertEquals( L(2,10), L(2,20).prefix( L(10,12)));

		//suffix
		assertEquals(L(3,20), L(2,20).suffix(1));
		assertEquals(L(19,20), L(2,20).suffix(-1));
		assertEquals(L(12,20), L(2,20).suffix( L(10,12)));

	}
	
	@Test
	public void testLocationIntersections() {
		// One inside another
		Location r21_25 = new Location( 21, 25 );
		Location r1_100 = new Location(1, 100 );
		
		assertEquals(r21_25, r21_25.intersection( r1_100));
		assertEquals(r21_25, r1_100.intersection( r21_25));
		
		// Non overlapping
		Location r10_100 = new Location(10, 100 );
		Location r1_9 = new Location( 1, 9 );
		
		assertNull(r10_100.intersection( r1_9));
		assertNull(r1_9.intersection( new Location( 9, 10 )));
		
		// Partially overlappping
		Location r1_25 = new Location( 1, 25 );
		Location r21_100 = new Location(21, 100 );
		assertEquals(r21_25, r1_25.intersection( r21_100));
		assertEquals(r21_25, r21_100.intersection( r1_25));		
	}

	//shorthand for testing
	private static Location L( int s, int e ) {
		return new Location( s, e );
	}

}
