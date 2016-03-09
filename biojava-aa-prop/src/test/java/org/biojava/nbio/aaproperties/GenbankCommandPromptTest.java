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
package org.biojava.nbio.aaproperties;

import org.biojava.nbio.aaproperties.CommandPrompt;
import org.junit.Test;

import java.io.File;

public class GenbankCommandPromptTest {
	@Test
	public void testAdvancedXMLExample() throws Exception{
		File output = new File(System.getProperty("java.io.tmpdir"),"modifiedTestGB.csv");
		output.deleteOnExit();
		//-i BondFeature.gb -a
		String[] args = new String[11];
		args[0] = "-i";
		args[1] = "./src/test/resources/BondFeature.gb";;
		args[2] = "-x";
		args[3] = "./src/main/resources/AdvancedAminoAcidComposition.xml";
		args[4] = "-0";
		args[5] = "0";
		args[6] = "-0";
		args[7] = "1";
		args[8] = "-a";
		args[9] = "-o";
		args[10] = output.toString();
		CommandPrompt.run(args);
	}

	@Test
	public void testExample1() throws Exception{
		File output = new File(System.getProperty("java.io.tmpdir"),"testgb.tsv");
		output.deleteOnExit();
		//-i BondFeature.gb -a
		String[] args = new String[7];
		args[0] = "-i";
		args[1] = "./src/test/resources/BondFeature.gb";;
		args[2] = "-a";
		args[3] = "-o";
		args[4] = output.toString();
		args[5] = "-f";
		args[6] = "tsv";
		CommandPrompt.run(args);
	}

	@Test
	public void testExample1WithCSV() throws Exception{
		File output = new File(System.getProperty("java.io.tmpdir"),"testgb.csv");
		output.deleteOnExit();
		//-i BondFeature.gb -a
		String[] args = new String[7];
		args[0] = "-i";
		args[1] = "./src/test/resources/BondFeature.gb";;
		args[2] = "-a";
		args[3] = "-o";
		args[4] = output.toString();
		args[5] = "-f";
		args[6] = "csv";
		CommandPrompt.run(args);
	}

	@Test
	public void testExample2() throws Exception{
		//-i BondFeature.gb -1 -3 -7
		String[] args = new String[5];
		args[0] = "-i";
		args[1] = "./src/test/resources/BondFeature.gb";
		args[2] = "-1";
		args[3] = "-3";
		args[4] = "-7";
		CommandPrompt.run(args);
	}

	@Test
	public void testExample3() throws Exception{
		//-i BondFeature.gb -0 A -0 N -1
		String[] args = new String[7];
		args[0] = "-i";
		args[1] = "./src/test/resources/BondFeature.gb";
		args[2] = "-0";
		args[3] = "A";
		args[4] = "-0";
		args[5] = "N";
		args[6] = "-1";
		CommandPrompt.run(args);
	}

	@Test
	public void testWithCases() throws Exception{
		//-i BondFeature.gb -0 A -0 N -1
		String[] args = new String[7];
		args[0] = "-i";
		args[1] = "./src/test/resources/BondFeature.gb";
		args[2] = "-0";
		args[3] = "A";
		args[4] = "-0";
		args[5] = "N";
		args[6] = "-1";
		CommandPrompt.run(args);
	}
}
