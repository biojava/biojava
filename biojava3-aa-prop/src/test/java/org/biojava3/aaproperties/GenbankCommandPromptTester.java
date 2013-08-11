package org.biojava3.aaproperties;

import org.junit.Test;

public class GenbankCommandPromptTester {
	@Test
	public void testAdvancedXMLExample() throws Exception{
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
		args[10] = "./src/test/resources/modifiedTestGB.csv";
		CommandPrompt.run(args);
	}
	
	@Test
	public void testExample1() throws Exception{
		//-i BondFeature.gb -a
		String[] args = new String[7];
		args[0] = "-i";
		args[1] = "./src/test/resources/BondFeature.gb";;
		args[2] = "-a";
		args[3] = "-o";
		args[4] = "./src/test/resources/testgb.tsv";
		args[5] = "-f";
		args[6] = "tsv";
		CommandPrompt.run(args);
	}
	
	@Test
	public void testExample1WithCSV() throws Exception{
		//-i BondFeature.gb -a
		String[] args = new String[7];
		args[0] = "-i";
		args[1] = "./src/test/resources/BondFeature.gb";;
		args[2] = "-a";
		args[3] = "-o";
		args[4] = "./src/test/resources/testgb.csv";
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
