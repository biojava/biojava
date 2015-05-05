package org.biojava.nbio.structure.align.cemc;

import java.util.concurrent.Callable;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.CallableStructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * Simple Callable Class that calculates an alignment.
 * It is designed to be run in a different thread, so that multiple pairwise alignments can be run in parallel.
 * Similar to {@link CallableStructureAlignment}, but more general and simpler (the other contains DB search details). 
 * 
 * @author Aleix Lafita
 */
public class ParallelAlignment implements Callable<AFPChain> {

	private StructureAlignment algorithm;
	private Atom[] ca1;
	private Atom[] ca2;
	
	/**
	 * Constructor.
	 * @param ca1 Atoms to align of the first structure
	 * @param ca2 Atoms to align of the second structure
	 * @param algorithm the alignment algorithm to use
	 */
	public ParallelAlignment(Atom[] ca1, Atom[] ca2, StructureAlignment algorithm){
		this.ca1 = ca1;
		this.ca2 = ca2;
		this.algorithm = algorithm;
	}

	@Override
	public AFPChain call() throws StructureException {
		//Perform and return the alignment
		return algorithm.align(ca1, ca2);
	}

}