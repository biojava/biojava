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
package org.biojava.nbio.aaproperties.xml;

//import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GenerateJavaCodesFromTextTest {

	private final static Logger logger = LoggerFactory.getLogger(GenerateJavaCodesFromTextTest.class);

	/*
	 * Generate java codes from two text files; Symbol2Name.txt and Symbol2Weight.txt
	 */
	// commenting out since there's no assertions here - JD 2015-07-28
	//@Test
	public void generateCodes() throws IOException{
		BufferedReader input = new BufferedReader(new FileReader("./src/test/resources/Symbol2Name.txt"));
		Map<String, String> symbol2Name = new HashMap<String, String>();
		String line = input.readLine();//Header line is not required
		while((line = input.readLine()) != null){
			String[] sA = line.split("\t");
			symbol2Name.put(sA[1].trim(), sA[2].trim());
		}
		input.close();

		input = new BufferedReader(new FileReader("./src/test/resources/Symbol2Weight.txt"));
		List<String> elementNameList = new ArrayList<String>();
		String elementName = null;
		String elementNameLower = null;
		Double elementMass = null;
		while((line = input.readLine()) != null){
			String[] sA = line.split("\t");
			if(sA[0].length() > 0){
				//Elements
				logger.info("{}.setIsotopes(iList);", elementNameLower);
				//int decimalPlace = getDecimalPlace(elementMass + "");
				//System.out.println("assertEquals(" + elementMass + ", Utils.roundToDecimals(" + elementNameLower +
					//	".getMass(), " + decimalPlace + "));");

				String symbol = sA[1].trim();
				elementName = symbol2Name.get(symbol);
				elementNameLower = elementName.toLowerCase() ;
				int protonNumber = Integer.parseInt(sA[0].trim());
				elementMass = cleanNumber(sA[5]);
				if(protonNumber > 82) break;
				elementNameList.add(elementNameLower);

				logger.info("iList = new ArrayList<Isotope>();");
				logger.info("Element {} = new Element(\"{}\", \"{}\", {}, null, {});",
						elementNameLower, elementName, symbol, protonNumber, elementMass);
			}
			int neutronNumber = Integer.parseInt(sA[2]);
			double weight = cleanNumber(sA[3]);
			//if(sA.length > 4 && sA[4].length() > 0) abundance = cleanNumber(sA[4]);
			logger.info("iList.add(new Isotope(\"{}-{}\", {}, {}));",
					elementName, neutronNumber, neutronNumber, weight);
		}
		input.close();

		logger.info("List<Element> eList = new ArrayList<Element>();");
		for(String e:elementNameList) logger.info("eList.add({});", e);
	}

	private double cleanNumber(String s){
		int index = s.indexOf("(");
		if(index != -1) s = s.substring(0, index);
		return Double.parseDouble(s.replace(" ", "").replace("[", "").replace("]", ""));
	}

	/*private int getDecimalPlace(String s){
		int i = s.indexOf(".");
		return s.length() - (i + 1);
	}*/
}
