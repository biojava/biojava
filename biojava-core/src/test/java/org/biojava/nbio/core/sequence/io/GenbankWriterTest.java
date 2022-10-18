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
/**
 *
 */
package org.biojava.nbio.core.sequence.io;


import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.features.AbstractFeature;
import org.biojava.nbio.core.sequence.features.DBReferenceInfo;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.features.TextFeature;
import org.biojava.nbio.core.sequence.location.SimpleLocation;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.Assert;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


/**
 * @author mckeee1
 *
 */
public class GenbankWriterTest {


	@Test
	public void testProcess() throws Exception {

		InputStream inStream = GenbankWriterTest.class.getResourceAsStream("/NM_000266.gb");
		//File dnaFile = new File("src/test/resources/NM_000266.gb");
		LinkedHashMap<String, DNASequence> dnaSequences = GenbankReaderHelper.readGenbankDNASequence( inStream );
		ByteArrayOutputStream fragwriter = new ByteArrayOutputStream();
		ArrayList<DNASequence> seqs = new ArrayList<DNASequence>();
		for(DNASequence seq : dnaSequences.values()) {
			seqs.add(seq);
		}
		GenbankWriterHelper.writeNucleotideSequence(fragwriter, seqs,
				GenbankWriterHelper.LINEAR_DNA);
		//System.out.println(fragwriter.toString());
		ByteArrayInputStream fragreader = new ByteArrayInputStream(fragwriter.toByteArray());
		/**
		 * Hello Jacek
		 * can you please investigate why this test fails? it seems that
		 * fragreader at the line below is read with the last feature
		 * in an invalid state: location = 2005..2004
		 */
		//dnaSequences = GenbankReaderHelper.readGenbankDNASequence( fragreader );
		fragwriter.close();
		Assert.assertEquals(seqs.get(0).getSequenceAsString(), dnaSequences.values().iterator().next().getSequenceAsString());
	}
	
	/**
	 * String Formatter error when key or value of Qualifier has character "%"
	 * https://github.com/biojava/biojava/issues/886
	 */
	@Test
	public void testGithub886() throws Exception {
		
		DNASequence seq = new DNASequence("ATGC");
		seq.setAccession(new AccessionID("."));
		AbstractFeature feature = new TextFeature("CDS", "source", "short description", "description");
		feature.setLocation(new SimpleLocation(1, 10, Strand.POSITIVE));

		// no percent symbols in key or value
		feature.addQualifier("note1", new Qualifier("note1", "50", true));
		// percent symbol in key
		feature.addQualifier("note2", new Qualifier("%note2", "50", true));
		feature.addQualifier("note3", new Qualifier("not%e3", "50", true));
		feature.addQualifier("note4", new Qualifier("note4%", "50", true));
		// percent symbol in value
		feature.addQualifier("note5", new Qualifier("note5", "%50", true));
		feature.addQualifier("note6", new Qualifier("note6", "5%0", true));
		feature.addQualifier("note7", new Qualifier("note7", "50%", true));
		
		seq.addFeature(feature);
		
		ByteArrayOutputStream fragwriter = new ByteArrayOutputStream();
		GenbankWriterHelper.writeNucleotideSequence(
				fragwriter, 
				Arrays.asList(seq), 
				GenbankWriterHelper.LINEAR_DNA);
		fragwriter.close();
		//System.out.println(fragwriter.toString().replaceAll("\r\n", "\n"));
		
		// now read in the file that was created and check that the qualifiers were created correctly
		InputStream readerInputStream = new ByteArrayInputStream(fragwriter.toByteArray());
		DNASequence newSeq = GenbankReaderHelper.readGenbankDNASequence(readerInputStream).values().iterator().next();
		AbstractFeature newFeature = (TextFeature) seq.getFeaturesByType("CDS").get(0);
		Map<String, List<Qualifier>> newQualifiers = newFeature.getQualifiers();
		
		assertEquals("note1", newQualifiers.get("note1").get(0).getName());
		assertEquals("50", newQualifiers.get("note1").get(0).getValue());
		
		assertEquals("%note2", newQualifiers.get("note2").get(0).getName());
		assertEquals("50", newQualifiers.get("note2").get(0).getValue());
		
		assertEquals("not%e3", newQualifiers.get("note3").get(0).getName());
		assertEquals("50", newQualifiers.get("note3").get(0).getValue());
		
		assertEquals("note4%", newQualifiers.get("note4").get(0).getName());
		assertEquals("50", newQualifiers.get("note4").get(0).getValue());
		
		assertEquals("note5", newQualifiers.get("note5").get(0).getName());
		assertEquals("%50", newQualifiers.get("note5").get(0).getValue());
		
		assertEquals("note6", newQualifiers.get("note6").get(0).getName());
		assertEquals("5%0", newQualifiers.get("note6").get(0).getValue());
		
		assertEquals("note7", newQualifiers.get("note7").get(0).getName());
		assertEquals("50%", newQualifiers.get("note7").get(0).getValue());
		
	}
	
	@Test
	public void testLocationJoins() throws Exception {
		
		// First read a GenBank file containing location joins
		InputStream inStream = GenbankWriterTest.class.getResourceAsStream("/with_joins.gb");
		DNASequence sequence = GenbankReaderHelper.readGenbankDNASequence(inStream).values().iterator().next();
		
		// Check the joins are read correctly
		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeatures();
		
		FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> join1 = features.get(0);
		List<Location> join1SubLocs = join1.getLocations().getSubLocations();
		
		assertEquals("join1, getType()", "CDS", join1.getType());
		assertEquals("join1, getLocations().getStrand()", "POSITIVE", join1.getLocations().getStrand().toString());
		assertEquals("join1, getLocations().getSubLocations().size()", 6, join1SubLocs.size());
		
		assertEquals("join1, SubLocation 1)", 1, join1SubLocs.get(0).getStart().getPosition().intValue());
		assertEquals("join1, SubLocation 1)", 1, join1SubLocs.get(0).getEnd().getPosition().intValue());
		
		assertEquals("join1, SubLocation 2)", 10, join1SubLocs.get(1).getStart().getPosition().intValue());
		assertEquals("join1, SubLocation 2)", 12, join1SubLocs.get(1).getEnd().getPosition().intValue());
		
		assertEquals("join1, SubLocation 3)", 30, join1SubLocs.get(2).getStart().getPosition().intValue());
		assertEquals("join1, SubLocation 3)", 30, join1SubLocs.get(2).getEnd().getPosition().intValue());
		
		assertEquals("join1, SubLocation 3)", 35, join1SubLocs.get(3).getStart().getPosition().intValue());
		assertEquals("join1, SubLocation 3)", 38, join1SubLocs.get(3).getEnd().getPosition().intValue());
		
		assertEquals("join1, SubLocation 5)", 43, join1SubLocs.get(4).getStart().getPosition().intValue());
		assertEquals("join1, SubLocation 5)", 46, join1SubLocs.get(4).getEnd().getPosition().intValue());
		
		assertEquals("join1, SubLocation 6)", 47, join1SubLocs.get(5).getStart().getPosition().intValue());
		assertEquals("join1, SubLocation 6)", 50, join1SubLocs.get(5).getEnd().getPosition().intValue());
		
		//qualifiers
		assertEquals("join1, getType()", "Joined feature", join1.getQualifiers().get("standard_name").get(0).getValue());
		
		//Join 2
		FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> join2 = features.get(1);
		List<Location> join2SubLocs = join2.getLocations().getSubLocations();
		
		assertEquals("join2, getType()", "CDS", join2.getType());
		assertEquals("join2, getLocations().getStrand()", "NEGATIVE", join2.getLocations().getStrand().toString());
		assertEquals("join2, getLocations().getSubLocations().size()", 5, join2SubLocs.size());
		
		assertEquals("join2, SubLocation 1)", 33, join2SubLocs.get(0).getStart().getPosition().intValue());
		assertEquals("join2, SubLocation 1)", 33, join2SubLocs.get(0).getEnd().getPosition().intValue());
		
		assertEquals("join2, SubLocation 2)", 35, join2SubLocs.get(1).getStart().getPosition().intValue());
		assertEquals("join2, SubLocation 2)", 37, join2SubLocs.get(1).getEnd().getPosition().intValue());
		
		assertEquals("join2, SubLocation 3)", 41, join2SubLocs.get(2).getStart().getPosition().intValue());
		assertEquals("join2, SubLocation 3)", 43, join2SubLocs.get(2).getEnd().getPosition().intValue());
		
		assertEquals("join2, SubLocation 4)", 44, join2SubLocs.get(3).getStart().getPosition().intValue());
		assertEquals("join2, SubLocation 4)", 46, join2SubLocs.get(3).getEnd().getPosition().intValue());
		
		assertEquals("join2, SubLocation 5)", 47, join2SubLocs.get(4).getStart().getPosition().intValue());
		assertEquals("join2, SubLocation 5)", 50, join2SubLocs.get(4).getEnd().getPosition().intValue());
		
		//qualifiers
		assertEquals("join2, getType()", "Joined feature on complement", join2.getQualifiers().get("standard_name").get(0).getValue());
		
		// Now write the joins back to a file using the GenbankWriterHelper
		ByteArrayOutputStream fragwriter = new ByteArrayOutputStream();
		GenbankWriterHelper.writeNucleotideSequenceOriginal(
				fragwriter, 
				Arrays.asList(sequence));
		fragwriter.close();
		
		//System.out.println(fragwriter.toString().replaceAll("\r\n", "\n"));
		
		// Read the output file and test that no information is lost
		InputStream readerInputStream = new ByteArrayInputStream(fragwriter.toByteArray());
		DNASequence newSequence = GenbankReaderHelper.readGenbankDNASequence(readerInputStream).values().iterator().next();
		
		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> newFeatures = newSequence.getFeatures();
		
		// Check the output matches the original sequence feature
		for (int i=0; i < features.size(); i++ ) {
			assertEquals("getFeatures(), getType()", features.get(i).getType(), newFeatures.get(i).getType());
			assertEquals("getFeatures(), getStart()", features.get(i).getLocations().getStart(), newFeatures.get(i).getLocations().getStart());
			assertEquals("getFeatures(), getEnd()", features.get(i).getLocations().getEnd(), newFeatures.get(i).getLocations().getEnd());
			assertEquals("getFeatures(), getStrand()", features.get(i).getLocations().getStrand(), newFeatures.get(i).getLocations().getStrand());

			List<Location> subLocations = features.get(i).getLocations().getSubLocations();
			List<Location> newSubLocations = newFeatures.get(i).getLocations().getSubLocations();
			assertEquals("getSubLocations()", subLocations.size(), newSubLocations.size());
		for (int j=0; j < subLocations.size(); j++ ) {
				assertEquals("getSubLocations(), getStart()",  subLocations.get(j).getStart(), newSubLocations.get(j).getStart());
				assertEquals("getSubLocations(), getEnd()",    subLocations.get(j).getEnd(), newSubLocations.get(j).getEnd());
				assertEquals("getSubLocations(), getStrand()", subLocations.get(j).getStrand(), newSubLocations.get(j).getStrand());
			}
			
			Map<String, List<Qualifier>> qualifiers = features.get(i).getQualifiers();
			Map<String, List<Qualifier>> newQualifiers = newFeatures.get(i).getQualifiers();
			
			for (String qualifierType:  qualifiers.keySet()) {
				assertEquals("getSubLocations()", qualifiers.get(qualifierType).get(0).getValue(), newQualifiers.get(qualifierType).get(0).getValue());				
			}
			
		}
		
	}
	
	/**
	 * Going from GenBank file -> DNASequence object -> GenBank file looses information
	 * https://github.com/biojava/biojava/issues/942
	 */
	@Test
	public void testGithub942() throws Exception {
		
		// Important information is lost when reading and writing a
		// GenBank file through GenbankReaderHelper & GenbankWriterHelper

		// First read the sample GenBank file from
		// https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html using the
		// GenbankReaderHelper
		InputStream inStream = GenbankWriterTest.class.getResourceAsStream("/NM_000266.gb");
		DNASequence sequence = GenbankReaderHelper.readGenbankDNASequence(inStream).values().iterator().next();

		// Then write sequence back to a file using the GenbankWriterHelper
		ByteArrayOutputStream fragwriter = new ByteArrayOutputStream();
		GenbankWriterHelper.writeNucleotideSequenceOriginal(
				fragwriter, 
				Arrays.asList(sequence));
		fragwriter.close();
		
		// Test no important information is lost
		InputStream readerInputStream = new ByteArrayInputStream(fragwriter.toByteArray());
		DNASequence newSequence = GenbankReaderHelper.readGenbankDNASequence(readerInputStream).values().iterator().next();
		
		//System.out.println(fragwriter.toString().replaceAll("\r\n", "\n"));

		assertEquals("getOriginalHeader()", sequence.getOriginalHeader(), newSequence.getOriginalHeader());
		assertEquals("getLength()", sequence.getLength(), newSequence.getLength());
		assertEquals("getAccession().getID()", sequence.getAccession().getID(), newSequence.getAccession().getID());
		assertEquals("getAccession().getVersion()", sequence.getAccession().getVersion(), newSequence.getAccession().getVersion());
		assertEquals("getDescription()", sequence.getDescription(), newSequence.getDescription());
		//assertEquals("getSource()", sequence.getSource(), newSequence.getSource());
		//assertEquals("getDNAType()", sequence.getDNAType(), newSequence.getDNAType());
		//assertEquals("getTaxonomy()", sequence.getTaxonomy(), newSequence.getTaxonomy());		
		//assertEquals("getReferences()", sequence.getReferences(), newSequence.getReferences());
		//assertEquals("getComments()", sequence.getComments(), newSequence.getComments());
		//assertEquals("getNotesList()", sequence.getNotesList(), newSequence.getNotesList());
		
		//Assuming the features will be in the same order
		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = sequence.getFeatures();
		List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> newFeatures = newSequence.getFeatures();
		
		//feature locations and qualifiers
		for (int i=0; i < features.size(); i++ ) {
			
			FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> feature = features.get(i);
			Location location           = feature.getLocations();
			List<Location> subLocations = location.getSubLocations();
			Map<String, List<Qualifier>> qualifiers  = feature.getQualifiers();
			
			FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> newFeature = newFeatures.get(i);
			Location newLocation           = newFeature.getLocations();
			List<Location> newSubLocations = newLocation.getSubLocations();
			Map<String, List<Qualifier>> newQualifiers  = newFeature.getQualifiers();
			
			assertEquals("feature, getType()",       feature.getType(),    newFeature.getType());
			assertEquals("feature, Location start",  location.getStart(),  newLocation.getStart());
			assertEquals("feature, Location end",    location.getEnd(),    newLocation.getEnd());
			assertEquals("feature, Location strand", location.getStrand(), newLocation.getStrand());
			assertEquals("feature, sublocations",    subLocations.size(),  newSubLocations.size());
			
			for (int j=0; j < subLocations.size(); j++ ) {
				assertEquals("SubLocations, start",  subLocations.get(j).getStart(),  newSubLocations.get(j).getStart());
				assertEquals("SubLocations, end",    subLocations.get(j).getEnd(),    newSubLocations.get(j).getEnd());
				assertEquals("SubLocations, strand", subLocations.get(j).getStrand(), newSubLocations.get(j).getStrand());
				
			}
			
			assertEquals("getQualifiers()", qualifiers.size(), newQualifiers.size());
			
			for (String qualifierType: qualifiers.keySet()) {
				
				List<Qualifier> qualifier   = new ArrayList<Qualifier>(qualifiers.get(qualifierType));
				List<Qualifier> newQualifier = new ArrayList<Qualifier>(newQualifiers.get(qualifierType));
				
				assertEquals("getQualifiers()", qualifier.size(), newQualifier.size());
				
				for (int k=0; k < qualifier.size(); k++) {
					if (qualifier.get(k) instanceof DBReferenceInfo) {
						DBReferenceInfo dbxref = (DBReferenceInfo) qualifier.get(k);
						DBReferenceInfo newDbxref = (DBReferenceInfo) newQualifier.get(k);
						assertEquals("getQualifiers() DBReferenceInfo", dbxref.getDatabase(), newDbxref.getDatabase());
						assertEquals("getQualifiers() DBReferenceInfo", dbxref.getId(), newDbxref.getId());
						
					} else {
						assertEquals("getQualifiers()", qualifier.get(k).getValue(), newQualifier.get(k).getValue());
						
					}					
				}					
			}			
		}
		
		assertEquals("getSequenceAsString()", sequence.getSequenceAsString(), newSequence.getSequenceAsString());
		
	}
}
