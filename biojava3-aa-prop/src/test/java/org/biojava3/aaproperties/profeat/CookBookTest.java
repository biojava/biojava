package org.biojava3.aaproperties.profeat;

import java.util.Map;

import org.biojava3.aaproperties.profeat.IProfeatProperties.ATTRIBUTE;
import org.biojava3.aaproperties.profeat.IProfeatProperties.DISTRIBUTION;
import org.biojava3.aaproperties.profeat.IProfeatProperties.GROUPING;
import org.biojava3.aaproperties.profeat.IProfeatProperties.TRANSITION;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CookBookTest {
	
	private final static Logger logger = LoggerFactory.getLogger(CookBookTest.class);

	@Test
	public void shortExample1() throws Exception{
		/*
		 * Composition
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE, Map<GROUPING, Double>> attribute2Grouping2Double = ProfeatProperties.getComposition(sequence); 
		for(ATTRIBUTE a:attribute2Grouping2Double.keySet()){
			logger.info("======={}=======", a);
			logger.info("GROUP1 = {}", attribute2Grouping2Double.get(a).get(GROUPING.GROUP1));
			logger.info("GROUP2 = {}", attribute2Grouping2Double.get(a).get(GROUPING.GROUP2));
			logger.info("GROUP3 = {}", attribute2Grouping2Double.get(a).get(GROUPING.GROUP3));
		}
	}
	
	@Test
	public void shortExample2() throws Exception{
		/*
		 * Transition 
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE, Map<TRANSITION, Double>> attribute2Transition2Double = ProfeatProperties.getTransition(sequence); 
		for(ATTRIBUTE a:attribute2Transition2Double.keySet()){
			logger.info("======={}=======", a);
			logger.info("1<=>1 = {}", attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_11));
			logger.info("2<=>2 = {}", attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_22));
			logger.info("3<=>3 = {}", attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_33));
			logger.info("1<=>2 = {}", attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_12));
			logger.info("1<=>3 = {}", attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_13));
			logger.info("2<=>3 = {}", attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_23));
		}
	}
	
	@Test
	public void shortExample3() throws Exception{
		/*
		 * Distribution
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE , Map<GROUPING, Map<DISTRIBUTION, Double>>> attribute2Grouping2Distribution2Double = ProfeatProperties.getDistributionPosition(sequence); 
		for(ATTRIBUTE a:attribute2Grouping2Distribution2Double.keySet()){
			logger.info("======={}=======", a);
			logger.info("GROUP1 = {}", attribute2Grouping2Distribution2Double.get(a).get(GROUPING.GROUP1));
			logger.info("GROUP2 = {}", attribute2Grouping2Distribution2Double.get(a).get(GROUPING.GROUP2));
			logger.info("GROUP3 = {}", attribute2Grouping2Distribution2Double.get(a).get(GROUPING.GROUP3));
		}
	}
}
