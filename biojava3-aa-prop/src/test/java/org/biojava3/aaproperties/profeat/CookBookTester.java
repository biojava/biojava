package org.biojava3.aaproperties.profeat;

import java.util.Map;

import org.biojava3.aaproperties.profeat.IProfeatProperties.ATTRIBUTE;
import org.biojava3.aaproperties.profeat.IProfeatProperties.DISTRIBUTION;
import org.biojava3.aaproperties.profeat.IProfeatProperties.GROUPING;
import org.biojava3.aaproperties.profeat.IProfeatProperties.TRANSITION;
import org.junit.Test;

public class CookBookTester {
	@Test
	public void shortExample1() throws Exception{
		/*
		 * Composition
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE, Map<GROUPING, Double>> attribute2Grouping2Double = ProfeatProperties.getComposition(sequence); 
		for(ATTRIBUTE a:attribute2Grouping2Double.keySet()){
			System.out.println("=======" + a + "=======");
			System.out.println("GROUP1 = " + attribute2Grouping2Double.get(a).get(GROUPING.GROUP1));
			System.out.println("GROUP2 = " + attribute2Grouping2Double.get(a).get(GROUPING.GROUP2));
			System.out.println("GROUP3 = " + attribute2Grouping2Double.get(a).get(GROUPING.GROUP3));
			System.out.println();
		}
		System.out.println();
	}
	
	@Test
	public void shortExample2() throws Exception{
		/*
		 * Transition 
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE, Map<TRANSITION, Double>> attribute2Transition2Double = ProfeatProperties.getTransition(sequence); 
		for(ATTRIBUTE a:attribute2Transition2Double.keySet()){
			System.out.println("=======" + a + "=======");
			System.out.println("1<=>1 = " + attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_11));
			System.out.println("2<=>2 = " + attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_22));
			System.out.println("3<=>3 = " + attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_33));
			System.out.println("1<=>2 = " + attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_12));
			System.out.println("1<=>3 = " + attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_13));
			System.out.println("2<=>3 = " + attribute2Transition2Double.get(a).get(TRANSITION.BETWEEN_23));
			System.out.println();
		}
		System.out.println();
	}
	
	@Test
	public void shortExample3() throws Exception{
		/*
		 * Distribution
		 */
		String sequence = "QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE";
		Map<ATTRIBUTE , Map<GROUPING, Map<DISTRIBUTION, Double>>> attribute2Grouping2Distribution2Double = ProfeatProperties.getDistributionPosition(sequence); 
		for(ATTRIBUTE a:attribute2Grouping2Distribution2Double.keySet()){
			System.out.println("=======" + a + "=======");
			System.out.println("GROUP1 = " + attribute2Grouping2Distribution2Double.get(a).get(GROUPING.GROUP1));
			System.out.println("GROUP2 = " + attribute2Grouping2Distribution2Double.get(a).get(GROUPING.GROUP2));
			System.out.println("GROUP3 = " + attribute2Grouping2Distribution2Double.get(a).get(GROUPING.GROUP3));
		}
	}
}
