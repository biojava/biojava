package org.biojava3.aaproperties;

import java.util.Map;

import org.biojava3.aaproperties.PeptideProperties.SingleLetterAACode;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.Test;
import static junit.framework.Assert.*;

public class PeptidePropertiesImplTester {

	/**
	 * Test input 
	 */
	private final String peptide = "MTADGPRELLQLRAAVRHR";
	
	@Test
	public void testAAComposition() { 
		Map<AminoAcidCompound, Double> composition = PeptideProperties.getAAComposition(peptide);
		AminoAcidCompound m = AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("M");
		int plen = peptide.length();
		assertEquals(1.0/plen,composition.get(m));
		AminoAcidCompound a = AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("A");
		assertEquals(3.0/plen,composition.get(a));
		//TODO finish for every AA
		/* 
		 * TODO This is horrible:  
		 * AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("M");
		 * either AminoAcidCompoundSet should be modified or adaptor method introduced
		 */
	}
	
	/*
	 * TODO 
	 */
	@Test
	public void testEnrichment() {
		
	}
	
}
