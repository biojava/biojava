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

import org.biojava.nbio.structure.io.sifts.*;
import org.junit.Assert;
import org.junit.Test;

import java.io.InputStream;
import java.util.List;
import java.util.zip.GZIPInputStream;

public class TestSiftsParsing {


	@Test
	public void test4DIA() throws Exception {
		List<SiftsEntity> entities = SiftsMappingProvider.getSiftsMapping("4DIA");

		Assert.assertNotNull(entities);

		Assert.assertEquals(1, entities.size());

		for (SiftsEntity e : entities) {
			//System.out.println(e.getEntityId() + " " +e.getType());


			Assert.assertTrue(e.getSegments().size() > 0);
			for (SiftsSegment seg : e.getSegments()) {
				Assert.assertTrue(seg.getResidues().size() > 0);

				for (SiftsResidue res : seg.getResidues()) {

					if (res.getUniProtResName() != null) {
						Assert.assertNotNull(res.getUniProtAccessionId());
						Assert.assertNotNull(res.getUniProtResName());

						// test for github ticket #280
						if (res.getUniProtPos() == 129) {

							Assert.assertTrue(res.getNotObserved());
						}

					}
				}
			}

		}


	}

	@Test
	public void test4jn3() throws Exception {
		List<SiftsEntity> entities = SiftsMappingProvider.getSiftsMapping("4jn3");

		Assert.assertNotNull(entities);

		Assert.assertEquals(2, entities.size());

		for (SiftsEntity e : entities) {
			//System.out.println(e.getEntityId() + " " +e.getType());


			Assert.assertTrue(e.getSegments().size() > 0);
			for (SiftsSegment seg : e.getSegments()) {
				Assert.assertTrue(seg.getResidues().size() > 0);
				//System.out.println(seg.getResidues().size());
				//System.out.println(" Segment: " + seg.getSegId() + " " + seg.getStart() + " " + seg.getEnd()) ;
				//
				for (SiftsResidue res : seg.getResidues()) {
					//System.out.println("  " + res);
					if (res.getUniProtResName() != null) {
						Assert.assertNotNull(res.getUniProtAccessionId());
						Assert.assertNotNull(res.getUniProtResName());

					}
				}
			}

		}


	}

	@Test
	public void test4DOU() throws Exception {

		// get file from test resource folder, since the file from EBI seems to have issues (on 20170405).

		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4dou.sifts.xml.gz"));

		SiftsXMLParser parser = new SiftsXMLParser();

		parser.parseXmlFile(inStream);

		List<SiftsEntity> entities = parser.getEntities();

		Assert.assertNotNull(entities);

		Assert.assertEquals(1, entities.size());

		for (SiftsEntity e : entities) {
			//System.out.println(e.getEntityId() + " " +e.getType());

			//	4DOU has 3 segments
			Assert.assertEquals("SiftsEntity does not have 3 segments, but " + e.getSegments().size(), 3, e.getSegments().size());

			// test segment 1:

			//SiftsSegment seg1 = e.getSegments().get(0);
			//System.out.println(" Segment: " + seg1.getSegId() + " " + seg1.getStart() + " " + seg1.getEnd() + " res. size: " + seg1.getResidues().size());
			//assertTrue(seg1.getResidues().size() == 17);

			for (SiftsSegment seg : e.getSegments()) {
				Assert.assertTrue(seg.getResidues().size() > 0);

				//System.out.println(" Segment: " + seg.getSegId() + " " + seg.getStart() + " " + seg.getEnd() + " res. size: " + seg.getResidues().size()) ;

				for (SiftsResidue res : seg.getResidues()) {


					if (res.getUniProtResName() != null) {
						//System.out.println("  " + res);
						Assert.assertNotNull(res.getUniProtAccessionId());
						Assert.assertNotNull(res.getUniProtResName());

					}
				}
				//break;
			}

		}


	}

	@Test
	public void test4O6W() throws Exception {
		List<SiftsEntity> entities = SiftsMappingProvider.getSiftsMapping("4O6W");

		Assert.assertNotNull(entities);

		Assert.assertEquals(2, entities.size());

		int ecount = 0;
		for (SiftsEntity e : entities) {
			ecount++;

			// we only test the 2nd segment in entity #1
			if (ecount != 1)
				continue;


			Assert.assertEquals("A", e.getEntityId());
			Assert.assertEquals("protein", e.getType());


			//	4O6W A has 2 segments
			Assert.assertEquals(2, e.getSegments().size());

			// test segment 2:

			SiftsSegment seg = e.getSegments().get(1);

			//SiftsSegment seg1 = e.getSegments().get(0);
			//System.out.println(" Segment: " + seg1.getSegId() + " " + seg1.getStart() + " " + seg1.getEnd() + " res. size: " + seg1.getResidues().size());
			//assertTrue(seg1.getResidues().size() == 17);

			Assert.assertTrue(seg.getResidues().size() > 0);


			for (SiftsResidue res : seg.getResidues()) {


				if (res.getUniProtResName() != null) {
					//System.out.println("  " + res);
					Assert.assertNotNull(res.getUniProtAccessionId());
					Assert.assertNotNull(res.getUniProtResName());

				}

				if (res.getPdbResNum().equals("502")) {

					Assert.assertTrue(res.getNotObserved());

				}
			}
			//break;
		}


	}


}
