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
package org.biojava.nbio.aaproperties.profeat;

import org.biojava.nbio.aaproperties.profeat.IProfeatProperties.ATTRIBUTE;
import org.biojava.nbio.aaproperties.profeat.IProfeatProperties.DISTRIBUTION;
import org.biojava.nbio.aaproperties.profeat.IProfeatProperties.GROUPING;
import org.biojava.nbio.aaproperties.profeat.IProfeatProperties.TRANSITION;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Map;

//import org.junit.Test;

public class CookBookTest {

	private final static Logger logger = LoggerFactory.getLogger(CookBookTest.class);

	// TODO there's no assertions here, i.e. this is not a test! must fix! For the moment removed test tags - JD 2016-03-08

	public void shortExample1() throws Exception{
		/*
		 * Composition
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE, Map<GROUPING, Double>> attribute2Grouping2Double = ProfeatProperties.getComposition(sequence);
		for(Map.Entry<ATTRIBUTE, Map<GROUPING, Double>> entry : attribute2Grouping2Double.entrySet()){
			logger.info("======={}=======", entry.getKey());
			logger.info("GROUP1 = {}", entry.getValue().get(GROUPING.GROUP1));
			logger.info("GROUP2 = {}", entry.getValue().get(GROUPING.GROUP2));
			logger.info("GROUP3 = {}", entry.getValue().get(GROUPING.GROUP3));
		}
	}


	public void shortExample2() throws Exception{
		/*
		 * Transition
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE, Map<TRANSITION, Double>> attribute2Transition2Double = ProfeatProperties.getTransition(sequence);
		for(Map.Entry<ATTRIBUTE, Map<TRANSITION, Double>> entry : attribute2Transition2Double.entrySet()){
			logger.info("======={}=======", entry.getKey());
			logger.info("1<=>1 = {}", entry.getValue().get(TRANSITION.BETWEEN_11));
			logger.info("2<=>2 = {}", entry.getValue().get(TRANSITION.BETWEEN_22));
			logger.info("3<=>3 = {}", entry.getValue().get(TRANSITION.BETWEEN_33));
			logger.info("1<=>2 = {}", entry.getValue().get(TRANSITION.BETWEEN_12));
			logger.info("1<=>3 = {}", entry.getValue().get(TRANSITION.BETWEEN_13));
			logger.info("2<=>3 = {}", entry.getValue().get(TRANSITION.BETWEEN_23));
		}
	}


	public void shortExample3() throws Exception{
		/*
		 * Distribution
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE , Map<GROUPING, Map<DISTRIBUTION, Double>>> attribute2Grouping2Distribution2Double = ProfeatProperties.getDistributionPosition(sequence);
		for(Map.Entry<ATTRIBUTE, Map<GROUPING, Map<DISTRIBUTION, Double>>> entry : attribute2Grouping2Distribution2Double.entrySet()){
			logger.info("======={}=======", entry.getKey());
			logger.info("GROUP1 = {}", entry.getValue().get(GROUPING.GROUP1));
			logger.info("GROUP2 = {}", entry.getValue().get(GROUPING.GROUP2));
			logger.info("GROUP3 = {}", entry.getValue().get(GROUPING.GROUP3));
		}
	}
}
