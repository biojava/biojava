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

package org.biojava.bio.structure.io;

import java.util.List;

import org.biojava.bio.structure.io.sifts.SiftsEntity;
import org.biojava.bio.structure.io.sifts.SiftsMappingProvider;
import org.biojava.bio.structure.io.sifts.SiftsResidue;
import org.biojava.bio.structure.io.sifts.SiftsSegment;

import junit.framework.TestCase;

public class TestSiftsParsing extends TestCase {

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

}
