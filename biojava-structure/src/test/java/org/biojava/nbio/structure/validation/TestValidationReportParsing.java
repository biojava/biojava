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
 * created at Sep 18, 2013
 * Author: ap3
 */

package org.biojava.nbio.structure.validation;

import org.junit.Test;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.fail;

public class TestValidationReportParsing {

	@Test
	public void test() {

		String[] testPDBids = new String[]{

				"3vtq",
				"3vtu",
				"3vtv",
				"3vtw",
				"3vu8",
				"3vua",
				"3vv5",
				"3vvd",
				"3vve",
				"3vvf",
				"3vw5",
				"3w1f",
				"3w5p",
				"3w5q",
				"3w5r",
				"3w5t",
				"3w9y",
				"3wcp",
				"3zjh",
				"3zji",
				"3zjj",
				"3zjm",
				"3zjn",
				"3zjo",
				"3zjp",
				"3zjq",
				"3zjr",
				"3zjs",
				"3znv",
				"3znx",
				"3znz",
				"3zoi",
				"3zoj",
				"3zpy",
		};

		for (String pdbId : testPDBids){
			testPDB(pdbId);
		}

	}

	private void testPDB(String pdbId) {
		try {
			JAXBContext ctx = JAXBContext.newInstance(new Class[] {WwPDBValidationInformation.class});

			Unmarshaller um = ctx.createUnmarshaller();

			InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/validation/"+pdbId+"-valdata.xml.gz"));

			WwPDBValidationInformation validationReport = (WwPDBValidationInformation) um.unmarshal(inStream);

			validationReport.getEntry();

//			Entry entry = validationReport.getEntry();
//			System.out.println(pdbId + " " + entry.getPDBRevisionNumber() +
//					"\t Rfree: " + entry.getDCCRfree() +
//					"\t Clashscore " + entry.getClashscore() +
//					"\t % Ramachandran outliers: "  + entry.getPercentRamaOutliers() +
//					"\t % RSRC outliers: " + entry.getPercentRSRZOutliers() );


		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

}
