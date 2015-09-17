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
 * Created on Sep 5, 2013
 * Author: andreas 
 *
 */

package org.biojava.nbio.structure.io;

import junit.framework.TestCase;
import org.biojava.nbio.structure.io.sifts.SiftsEntity;
import org.biojava.nbio.structure.io.sifts.SiftsMappingProvider;
import org.biojava.nbio.structure.io.sifts.SiftsResidue;
import org.biojava.nbio.structure.io.sifts.SiftsSegment;

import java.util.List;

public class TestSiftsParsing extends TestCase {


	public void test4DIA(){
		try {
			List<SiftsEntity> entities = SiftsMappingProvider.getSiftsMapping("4DIA");

			assertNotNull(entities);

			assertTrue(entities.size() == 1);

			for (SiftsEntity e : entities){
				//System.out.println(e.getEntityId() + " " +e.getType());


				assertTrue(e.getSegments().size() > 0 );
				for ( SiftsSegment seg: e.getSegments()) {
					assertTrue(seg.getResidues().size() > 0);
					
					for ( SiftsResidue res: seg.getResidues() ) {

						if ( res.getUniProtResName() != null ) {
							assertNotNull(res.getUniProtAccessionId() );
							assertNotNull(res.getUniProtResName());

							// test for github ticket #280
							if ( res.getUniProtPos() == 129) {

								assertTrue(res.getNotObserved());
							}

						}
					}
				}

			}

		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}


	}

	public void test4jn3(){
		try {
			List<SiftsEntity> entities = SiftsMappingProvider.getSiftsMapping("4jn3");

			assertNotNull(entities);

			assertTrue(entities.size() == 2);

			for (SiftsEntity e : entities){
				//System.out.println(e.getEntityId() + " " +e.getType());


				assertTrue(e.getSegments().size() > 0 );
				for ( SiftsSegment seg: e.getSegments()) {
					assertTrue(seg.getResidues().size() > 0);
					//System.out.println(seg.getResidues().size());
					//System.out.println(" Segment: " + seg.getSegId() + " " + seg.getStart() + " " + seg.getEnd()) ;
					//					
					for ( SiftsResidue res: seg.getResidues() ) {
						//System.out.println("  " + res);
						if ( res.getUniProtResName() != null ) {
							assertNotNull(res.getUniProtAccessionId() );
							assertNotNull(res.getUniProtResName());
							
						}
					}
				}

			}

		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}


	}

	public void test4DOU(){
		try {
			List<SiftsEntity> entities = SiftsMappingProvider.getSiftsMapping("4dou");

			assertNotNull(entities);

			assertTrue(entities.size() == 1);

			for (SiftsEntity e : entities){
				//System.out.println(e.getEntityId() + " " +e.getType());

				//	4DOU has 3 segments
				assertTrue(e.getSegments().size() == 3);
				
				// test segment 1:
				
				//SiftsSegment seg1 = e.getSegments().get(0);
				//System.out.println(" Segment: " + seg1.getSegId() + " " + seg1.getStart() + " " + seg1.getEnd() + " res. size: " + seg1.getResidues().size());
				//assertTrue(seg1.getResidues().size() == 17);
				
				for ( SiftsSegment seg: e.getSegments()) {
					assertTrue(seg.getResidues().size() > 0);
					
					//System.out.println(" Segment: " + seg.getSegId() + " " + seg.getStart() + " " + seg.getEnd() + " res. size: " + seg.getResidues().size()) ;
										
					for ( SiftsResidue res: seg.getResidues() ) {
						
						
						if ( res.getUniProtResName() != null ) {
							//System.out.println("  " + res);
							assertNotNull(res.getUniProtAccessionId() );
							assertNotNull(res.getUniProtResName());
							
						}
					}
					//break;
				}

			}

		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}


	}

	public void test4O6W(){
		try {
			List<SiftsEntity> entities = SiftsMappingProvider.getSiftsMapping("4O6W");

			assertNotNull(entities);

			assertTrue(entities.size() == 2);

			int ecount = 0;
			for (SiftsEntity e : entities){
				ecount++;

				// we only test the 2nd segment in entity #1
				if ( ecount != 1)
					continue;


				assertTrue(e.getEntityId().equals("A"));
				assertTrue(e.getType().equals("protein"));


				//	4O6W A has 2 segments
				assertTrue(e.getSegments().size() == 2);

				// test segment 2:

				SiftsSegment seg = e.getSegments().get(1);

				//SiftsSegment seg1 = e.getSegments().get(0);
				//System.out.println(" Segment: " + seg1.getSegId() + " " + seg1.getStart() + " " + seg1.getEnd() + " res. size: " + seg1.getResidues().size());
				//assertTrue(seg1.getResidues().size() == 17);

				assertTrue(seg.getResidues().size() > 0);


				for ( SiftsResidue res: seg.getResidues() ) {



					if ( res.getUniProtResName() != null ) {
						//System.out.println("  " + res);
						assertNotNull(res.getUniProtAccessionId() );
						assertNotNull(res.getUniProtResName());

					}

					if (res.getPdbResNum().equals("502")){
						
						assertTrue(res.getNotObserved());

					}
				}
				//break;
			}



		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}


}
