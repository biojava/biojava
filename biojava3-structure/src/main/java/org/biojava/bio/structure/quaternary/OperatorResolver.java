package org.biojava.bio.structure.quaternary;

import java.util.Collections;
import java.util.List;


/** A class to resolve the operators for transformations
 * 
 * @author Peter Rose
 *
 */
public class OperatorResolver {
	
	/**
	 * Unary operator expressions are parsed stored unary operations.
	 * For example the operator expression "(1,2,3,4)" is stored as a list 1,2,3,4
	 */
	private List<String> unaryOperators = Collections.emptyList();
	/**
	 * Binary Operator expressions are parsed and stored as ordered pairs of
	 * binary operators. For example the operator expression "(1-60)(61-88)"
	 * is saved as a list of pairs {1,61}, {1,62}, .., {1,88}, ... {60,88}.
	 */
	private List<OrderedPair<String>> binaryOperators = Collections.emptyList();

	
	/**
	 * Parses the operator expression and save the operators as a list
	 * of unary or binary operators (i.e. matrix multiplication, see below).
	 * Operation expressions are given in a compact notation and specify
	 * matrices from the operations list.
	 * An operation expression can be a comma-separated list 1, 5, 9,
	 * a dash-delimited range 1-60 or a matrix multiplication involving two
	 * or more lists or ranges. For instance, (X0)(1-20) specifies the 
	 * portion of the X174 procapsid crystal asymmetric unit belonging to 
	 * the first independent virus particle and corresponds
	 * to the 20 transformations [X0][1], [X0][2], ... , [X0][20].
	 * See C. Lawson, Acta Cryst., D64, 874-882, 2008.
	 *   
	 * @param operatorExpression the operator expression to be parsed
	 */
	public  void parseOperatorExpressionString(String operatorExpression) throws IllegalArgumentException {
		String expression = operatorExpression.trim();
		
		// remove single quotes, i.e. '(1-49)' in 1CGM
		expression = expression.replaceAll("'", "");

		if (BioAssemblyTools.isUnaryExpression(expression)) {
			unaryOperators = BioAssemblyTools.parseUnaryOperatorExpression(expression);
		} else {
			binaryOperators = BioAssemblyTools.parseBinaryOperatorExpression(expression);
		}
		
		//System.out.println("OperatorResolver: unary: " + unaryOperators + " | binary: " + binaryOperators);
	}


	


	public void setUnaryOperators(List<String> unaryOperators) {
		this.unaryOperators = unaryOperators;
	}


	


	public void setBinaryOperators(List<OrderedPair<String>> binaryOperators) {
		this.binaryOperators = binaryOperators;
	}
	
	/**
	 * Returns a list of operators for this assembly. The operators
	 * refer to the transformations that should be applied to
	 * the asym ids to generate this macromolecular assembly.
	 * @return the unary operators for this assembly
	 */
	public List<String> getUnaryOperators() {
		return unaryOperators;
	}

	/**
	 * Returns a list of operators for this assembly. The operators
	 * refer to the transformations that should be applied to
	 * the asym ids to generate this macromolecular assembly. 
	 * Each ordered pair refers to the multiplication 
	 * of the two transformation matrices in the
	 * pdbx_structure_oper_list category.
	 * @return the binary operators for this assembly
	 */
	public List<OrderedPair<String>> getBinaryOperators() {
		return binaryOperators;
	}
	
	
}
