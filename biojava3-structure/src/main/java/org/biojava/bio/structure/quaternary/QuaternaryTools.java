package org.biojava.bio.structure.quaternary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



/**
 * @author Peter Rose
 * 
 *
 */
public class QuaternaryTools {
	
	
	/**
	 * Checks if the passed in expression is a unary operator expression
	 * Example: (1,2,3) or (1-60) are unary operator expressions
	 *          (1-60)(61-88) is a binary operator expression, representing
	 *          a cartesian product of the two parenthesised lists
	 *          
	 * @param expression
	 * @return true if expression is a unary operator expression
	 */
	public static boolean isUnaryExpression(String expression) {
		int first = expression.indexOf("(");
		int last = expression.lastIndexOf("(");
		if (first < 0 || last < 0) {
			return true;
		}
		return ! (first == 0 && last > first);
	}
	
	public static List<String> parseUnaryOperatorExpression(String operatorExpression) throws IllegalArgumentException {
		return parseSubExpression(operatorExpression);
	}
	
	private static List<String> parseSubExpression(String expression) throws IllegalArgumentException {
		// remove parenthesis, if any
		String tmp = expression.replace("(", "");
		tmp = tmp.replace(")", "");

		// separate the operators
		List<String> components = null;
		try {
			components = Arrays.asList(tmp.split(","));
		} catch (Exception e) {
			throw new IllegalArgumentException("Invalid oper_expression: " + expression);
		}

		// expand ranges if present, i.e. 1-60 -> 1, 2, 3, ..., 60
		List<String> operators = new ArrayList<String>();
		for (String component : components) {
			if (component.contains("-")) {
				operators.addAll(expandRange(component));
			} else {
				operators.add(component);
			}
		}
		return operators;
	}
	
	/**
	 * Expands a range expression, i.e. (1-6) to a list 1,2,3,4,5,6
	 * @param expression the expression to be expanded
	 * @return list of items in range
	 * @throws IllegalArgumentException
	 */
	private static List<String> expandRange(String expression) throws IllegalArgumentException {
		int first = 0;
		int last = 0;
		try {
			String[] range = expression.split("-");
			first = Integer.parseInt(range[0]);
			last = Integer.parseInt(range[1]);
		} catch (Exception e) {
			throw new IllegalArgumentException("Invalid range specification in oper_expression: " + expression);
		}
		
		List<String> expandedExpression = new ArrayList<String>(last-first+1);
		for (int i = first; i <= last; i++) {
			expandedExpression.add(String.valueOf(i));
		}
		return expandedExpression;
	}
	
	public static List<OrderedPair<String>> parseBinaryOperatorExpression(String expression) 
			throws IllegalArgumentException {
				// split operator expression, i.e. (1,2,3)(4,5) into two subexpressions
				String[] subExpressions = null;
				try {
					subExpressions = expression.split("\\)\\(");
				} catch (Exception e) {
					throw new IllegalArgumentException("Invalid oper_expression: " + expression);
				}
				if (subExpressions.length != 2) {
					throw new IllegalArgumentException("Invalid oper_expression: " + expression);
				}
				List<String> leftSide = parseSubExpression(subExpressions[0]);
				List<String> rightSide = parseSubExpression(subExpressions[1]);

				// form the cartesian product of the two lists
				CartesianProduct<String> product = new CartesianProduct<String>(leftSide, rightSide);
				return product.getOrderedPairs();
			}
}
