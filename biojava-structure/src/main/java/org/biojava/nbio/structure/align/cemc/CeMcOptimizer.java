package org.biojava.nbio.structure.align.cemc;

import java.util.concurrent.Callable;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;

/**
 * This class takes a MultipleAlignment seed previously generated and runs a Monte Carlo optimization.
 * It implements Callable in order to be run in parallel along other optimization instances.
 * 
 * @author Aleix Lafita
 *
 */
public class CeMcOptimizer implements Callable<MultipleAlignment> {

	@Override
	public MultipleAlignment call() throws StructureException, StructureAlignmentException{
		// TODO Auto-generated method stub
		return null;
	}

}
