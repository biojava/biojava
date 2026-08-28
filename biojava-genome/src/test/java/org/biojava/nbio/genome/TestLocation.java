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

import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

class TestLocation {

	@Test
    void testLocation() {
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
		Assertions.assertEquals(7, L(14,14).distance( L(3,7) ));
		Assertions.assertEquals(7, L(3,7).distance( L(14,14) ));
		Assertions.assertEquals(3, L(1,4).distance( L(7, 10) ));

		//union
		Assertions.assertEquals(p10_17, p10_12.union( p14_17 ));
		Assertions.assertEquals(p10_17, p14_17.union( p10_12 ));
		Assertions.assertEquals(p15_19, p15_19.union( p15_16 ));

		//intersection
		Assertions.assertEquals(new Location( 21, 25 ), r13_17.union( r21_25 ).intersection( r21_25 ));


		//isBefore
		Assertions.assertTrue(r2_5.isBefore( r5_8 ));
		Assertions.assertTrue(!r2_5.isBefore( r4_7 ));

		//isAfter
		Assertions.assertTrue(r5_8.isAfter( r2_5 ));
		Assertions.assertTrue(!r5_8.isAfter( r4_7 ));

		//contains
		Assertions.assertTrue(p15_19.contains( p16_19 ));

		//overlaps
		Assertions.assertTrue(r2_5.overlaps( r4_7 ));
		Assertions.assertTrue(r2_5.overlaps( r0_3 ));
		Assertions.assertTrue(!r5_8.overlaps( r2_5 ));
		Assertions.assertTrue(!r2_5.overlaps( r5_8 ));


		//prefix
		Assertions.assertEquals(L(2,3), L(2,20).prefix(1));
		Assertions.assertEquals(L(2,19), L(2,20).prefix(-1));
		Assertions.assertEquals(L(2,10), L(2,20).prefix( L(10,12)));

		//suffix
		Assertions.assertEquals(L(3,20), L(2,20).suffix(1));
		Assertions.assertEquals(L(19,20), L(2,20).suffix(-1));
		Assertions.assertEquals(L(12,20), L(2,20).suffix( L(10,12)));

	}

	@Test
    void testLocationIntersections() {
		// One inside another
		Location r21_25 = new Location( 21, 25 );
		Location r1_100 = new Location(1, 100 );

		Assertions.assertEquals(r21_25, r21_25.intersection( r1_100));
		Assertions.assertEquals(r21_25, r1_100.intersection( r21_25));

		// Non overlapping
		Location r10_100 = new Location(10, 100 );
		Location r1_9 = new Location( 1, 9 );

		Assertions.assertNull(r10_100.intersection( r1_9));
		Assertions.assertNull(r1_9.intersection( new Location( 9, 10 )));

		// Partially overlappping
		Location r1_25 = new Location( 1, 25 );
		Location r21_100 = new Location(21, 100 );
		Assertions.assertEquals(r21_25, r1_25.intersection( r21_100));
		Assertions.assertEquals(r21_25, r21_100.intersection( r1_25));
	}

	//shorthand for testing
	private static Location L( int s, int e ) {
		return new Location( s, e );
	}

}
