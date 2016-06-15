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
package org.biojava.nbio.core.sequence.location;

import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;

/**
 *
 * @author Jacek Grzebyta
 */
@RunWith(Parameterized.class)
public class TargetedLocationParserTest {

	private Data request;

	public static class Data {
		private String Insdc;

		/**
		 * Parser data input. Based on that input it should be able to identity the origin and wanted target
		 * @param gi origin GI number
		 * @param originType origin compound set type
		 * @param Insdc string with INSDC notation
		 * @param compound wanted compound type {@see CompoundSet}
		 */
		public Data(String gi, CompoundSet<?> originType, String Insdc, CompoundSet<?> compound) {
			this.Insdc = Insdc;
		}
	};

	@Parameterized.Parameters
	public static Collection<Data[]> getLocations() throws Exception {


		Data[][] out = new Data[][]{
			{new Data("7525057", AminoAcidCompoundSet.getAminoAcidCompoundSet(),
					"join(complement(NC_000932.1:69611..69724),NC_000932.1:139856..140087,NC_000932.1:140625..140650)", DNACompoundSet.getDNACompoundSet())},

			{new Data("7525059", AminoAcidCompoundSet.getAminoAcidCompoundSet(),
					"NC_000932.1:72371..73897", DNACompoundSet.getDNACompoundSet())},

			{new Data("7525073", DNACompoundSet.getDNACompoundSet() ,
					"complement(NC_000932.1:84005..84283)", DNACompoundSet.getDNACompoundSet())},


			{new Data("7525012", DNACompoundSet.getDNACompoundSet(),
					"complement(9938..11461)", DNACompoundSet.getDNACompoundSet())}
		};

		return Arrays.asList(out);
	}

	public TargetedLocationParserTest(Data request) {
		this.request = request;
	}


	@Test
	public void locationTest() throws Exception {

		InsdcParser parser = new InsdcParser(DataSource.GENBANK);
		Location loc = parser.parse(request.Insdc);

		Assert.assertNotNull(loc);
		if (loc.isComplex()) {
			Assert.assertFalse(loc.getSubLocations().isEmpty());
		} else {
		}
	}
}
