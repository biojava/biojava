/*
 * BioJava development code This code may be freely distributed and modified
 * under the terms of the GNU Lesser General Public Licence. This should be
 * distributed with the code. If you do not have a copy, see:
 * http://www.gnu.org/copyleft/lesser.html Copyright for this code is held
 * jointly by the individual authors. These should be listed in @author doc
 * comments. For more information on the BioJava project and its aims, or to
 * join the biojava-l mailing list, visit the home page at:
 * http://www.biojava.org/
 */

package org.biojavax.ga.impl;

import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.Organism;

/**
 * Abstract implementation of Organism. Most implementations would want to
 * inherit from here.
 *
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public abstract class AbstractOrganism extends AbstractChangeable implements
    Organism {

	private double[]	fitness;

	protected SymbolList[]	   chromosomes;

	String	         name;

	protected AbstractOrganism() {
		chromosomes = new SymbolList[0];
		name = "";
	}

	protected AbstractOrganism(Organism org, String name) {
		chromosomes = org.getChromosomes();
		this.name = name;
	}

	public final SymbolList[] getChromosomes() {
		return chromosomes;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see org.biojavax.ga.Organism#getFitness()
	 */
	public final double[] getFitness() {
		return fitness;
	}

	public String getName() {
		return name;
	}

	public abstract boolean isHaploid();

	protected abstract void setChromImpl(SymbolList[] chromosomes);

	public final void setChromosomes(SymbolList[] chromosomes)
	    throws ChangeVetoException {
		if (!hasListeners()) {
			setChromImpl(chromosomes);
		} else {
			ChangeEvent ce = new ChangeEvent(this, Organism.CHROMOSOMES, chromosomes,
			    this.chromosomes);
			ChangeSupport changeSupport = super
			    .getChangeSupport(Organism.CHROMOSOMES);
			synchronized (changeSupport) {
				changeSupport.firePreChangeEvent(ce);
				this.chromosomes = chromosomes;
				changeSupport.firePostChangeEvent(ce);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see org.biojavax.ga.Organism#setFitness(double[])
	 */
	public final void setFitness(double[] fitness) {
		this.fitness = fitness.clone();
	}

	public final void setName(String name) throws ChangeVetoException {
		if (!hasListeners()) {
			this.name = name;
		} else {
			ChangeEvent ce = new ChangeEvent(this, Organism.NAME, chromosomes,
			    this.chromosomes);
			ChangeSupport changeSupport = super.getChangeSupport(Organism.NAME);
			synchronized (changeSupport) {
				changeSupport.firePreChangeEvent(ce);
				this.name = name;
				changeSupport.firePostChangeEvent(ce);
			}
		}
	}

}
