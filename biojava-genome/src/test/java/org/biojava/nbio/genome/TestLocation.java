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

	//shorthand for testing
	private static Location L( int s, int e ) {
		return new Location( s, e );
	}
	
}
