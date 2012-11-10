/*
 * BioJava development code
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
 * Author: Daniel Asarnow
 * Date:   2012-6-25
 */

package org.biojava.bio.structure.cath;

import junit.framework.TestCase;

import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.junit.BeforeClass;
import org.junit.Test;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Daniel Asarnow
 */
public class CathTest extends TestCase{



	@Test
	public void testCATH(){

//		CathDatabase database;
//		CathDatabase databaseCDDF;
//		CathDomain knownDomain;
//		CathNode knownNode;
//		CathSegment knownSegment;
//		List<CathSegment> knownSegments;
//		List<CathFragment> knownFragments;
//
//		System.out.println(" CathTest: Current memory: " + Runtime.getRuntime().totalMemory() / 1024/1024 + " Mb");
//		
//		database     = new CathInstallation((new UserConfiguration()).getPdbFilePath(),false,true);
//		databaseCDDF = new CathInstallation((new UserConfiguration()).getPdbFilePath(),true,false);
//		
//		System.out.println(" CathTest: After CATH init: " + Runtime.getRuntime().totalMemory() / 1024/1024 + " Mb");
//
//		knownDomain = new CathDomain();
//		knownDomain.setDomainName( "1oaiA00" );
//		knownDomain.setCATH("1.10.8.10");
//		knownDomain.setSOLID("1.1.1.1.1");
//		knownDomain.setLength(59);
//		knownDomain.setResolution(1.0);
//
//		//        knownDomain = new CathDescription();
//		knownDomain.setFormat("CDDF1.0");
//		//        knownDomain.setDomain(knownDomain);
//		knownDomain.setVersion("3.5.0");
//		DateFormat format = new SimpleDateFormat("dd-MMM-yyyy");
//		try {
//			knownDomain.setDate(format.parse("21-Sep-2011"));
//		} catch (ParseException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//			fail(e.getMessage());
//		}
//		knownDomain.setName("Nuclear RNA export factor. Chain: a. Fragment: uba domain, residues 56" +
//				"1-619. Synonym: tap, tip associating protein, mRNA export factor  tap." +
//				" Engineered: yes. Fxfg nucleoporin peptide. Chain: b. Fragment: nucleo" +
//				"porin peptide, residues 10-18. Engineered: yes");
//		knownDomain.setSource("Homo sapiens. Human. Organism_taxid: 9606.  Expressed in: escherichia " +
//				"coli. Expression_system_taxid: 562.");
//		knownDomain.setSequenceHeader(">pdb|1oaiA00");
//		knownDomain.setSequence("PTLSPEQQEMLQAFSTQSGMNLEWSQKCLQDNNWDYTRSAQAFTHLKAKGEIPEVAFMK");
//
//		knownSegment = new CathSegment();
//		knownSegment.setSegmentId(1);
//		knownSegment.setLength(59);
//		knownSegment.setStart("561");
//		knownSegment.setStop("619");
//		knownSegment.setSequenceHeader(">pdb|1oaiA00:1:1");
//		knownSegment.setSequence("PTLSPEQQEMLQAFSTQSGMNLEWSQKCLQDNNWDYTRSAQAFTHLKAKGEIPEVAFMK");
//		knownSegments = new ArrayList<CathSegment>();
//		knownSegments.add(knownSegment);
//
//		knownDomain.setSegments(knownSegments);
//
//		knownNode = new CathNode();
//		knownNode.setDescription("DNA helicase RuvA subunit, C-terminal domain");
//		knownNode.setNodeId("1.10.8.10");
//		knownNode.setRepresentative("1oaiA00");
//		knownNode.setParentId("1.10.8");
//
//		CathFragment fragment1 = new CathFragment();
//		fragment1.setFragmentId(1);
//		fragment1.setStart("0");
//		fragment1.setStop("1");
//		fragment1.setLength(2);
//		CathFragment fragment2 = new CathFragment();
//		fragment2.setFragmentId(2);
//		fragment2.setStart("209");
//		fragment2.setStop("209");
//		fragment2.setLength(1);
//
//		knownFragments = new ArrayList<CathFragment>();
//		knownFragments.add(fragment1);
//		knownFragments.add(fragment2);
//
//	
//
//	
//		
//		// test 4HHB
//		String pdbID = "4hhb";
//
//		List<CathDomain> cathDomains = database.getDomainsForPdb(pdbID);
//
//		assertTrue(cathDomains.size() == 4);
//
//		List<CathDomain> cathDomains2 =  databaseCDDF.getDomainsForPdb(pdbID);
//		assertEquals(cathDomains.size(),cathDomains2.size());
//
//		for (CathDomain domain : cathDomains){
//			StructureName n = new StructureName(domain.getDomainName());
//
//			assertTrue(n.getPdbId().equals(pdbID));
//
//			assertNotNull(domain.getSegments());
//
//			List<CathSegment> segments = domain.getSegments();
//
//			assertTrue(segments.size() > 0);
//
//			if ( n.getChainId().equals("A")) {
//				for (CathSegment s : segments){
//					assertEquals(s.getStart(),"1");
//					assertEquals(s.getStop(),"141");
//				}
//			}
//		}
//
//		// test CATH node
//		CathNode node = database.getCathNode("1.10.8.10");
//		assertEquals(knownNode.getNodeId(),node.getNodeId());
//		assertEquals(knownNode.getParentId(),node.getParentId());
//		assertEquals(knownNode.getRepresentative(),node.getRepresentative());
//		assertEquals(knownNode.getDescription(),node.getDescription());
//		assertEquals(knownNode.getCategory(),node.getCategory());
//	
//		// test CATH domain
//		
//		CathDomain domain = database.getDomainByCathId("1oaiA00");
//		System.out.println(domain);
//		assertEquals(knownDomain.getDomainName(),domain.getDomainName());
//		assertEquals(knownDomain.getLength(),domain.getLength());
//		assertEquals(knownDomain.getResolution(),domain.getResolution());
//		assertEquals(knownDomain.getCATH(),domain.getCATH());
//		assertEquals(knownDomain.getSOILD(),domain.getSOILD());
//	
//		
//		// test CATH desc
//		CathDomain descriptionCDDF = databaseCDDF.getDescriptionByCathId("1oaiA00");
//		assertEquals(knownDomain.getFormat(),descriptionCDDF.getFormat());
//		assertEquals(knownDomain.getDate(),descriptionCDDF.getDate());
//		assertEquals(knownDomain.getName(),descriptionCDDF.getName());
//		assertEquals(knownDomain.getSource(),descriptionCDDF.getSource());
//		assertEquals(knownDomain.getSequenceHeader(),descriptionCDDF.getSequenceHeader());
//		assertEquals(knownDomain.getSequence(),descriptionCDDF.getSequence());
//
//		List<CathSegment> segmentsCDDF = databaseCDDF.getDescriptionByCathId("1oaiA00").getSegments();
//		for (int i=0; i<segmentsCDDF.size(); i++) {
//			assertEquals(knownSegments.get(i).getLength(),segmentsCDDF.get(i).getLength());
//			assertEquals(knownSegments.get(i).getSegmentId(),segmentsCDDF.get(i).getSegmentId());
//			assertEquals(knownSegments.get(i).getSequenceHeader(),segmentsCDDF.get(i).getSequenceHeader());
//			assertEquals(knownSegments.get(i).getSequence(),segmentsCDDF.get(i).getSequence());
//			assertEquals(knownSegments.get(i).getStart(),segmentsCDDF.get(i).getStart());
//			assertEquals(knownSegments.get(i).getStop(),segmentsCDDF.get(i).getStop());
//		}
//	
//		
//		// test CATH segment
//		
//		List<CathSegment> segments = database.getDescriptionByCathId("1oaiA00").getSegments();
//		assertEquals(segments.size(),knownSegments.size());
//		for (int i = 0; i<segments.size(); i++) {
//			assertEquals(knownSegments.get(i).getLength(),segments.get(i).getLength());
//			assertEquals(knownSegments.get(i).getSegmentId(),segments.get(i).getSegmentId());
//			assertEquals(knownSegments.get(i).getStart(),segments.get(i).getStart());
//			assertEquals(knownSegments.get(i).getStop(),segments.get(i).getStop());
//		}
//	
//		
//		// test CATH fragments
//
//
//	
//		List<CathFragment> fragments = database.getFragmentsByPdbId("17gsA");
//		for (int i=0; i<fragments.size(); i++) {
//			assertEquals(knownFragments.get(i).getFragmentId(),fragments.get(i).getFragmentId());
//			assertEquals(knownFragments.get(i).getStart(),fragments.get(i).getStart());
//			assertEquals(knownFragments.get(i).getStop(),fragments.get(i).getStop());
//			assertEquals(knownFragments.get(i).getLength(),fragments.get(i).getLength());
//		}
//		
//		
		System.out.println(" CathTest: After CATH tests: " + Runtime.getRuntime().totalMemory() / 1024/1024 + " Mb");
	}

}
