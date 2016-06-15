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
package org.biojava.nbio.structure.align;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;

public abstract class AbstractStructureAlignment implements StructureAlignment{

	public static String newline = System.getProperty("line.separator");

	@Override
	abstract public  AFPChain align(Atom[] ca1, Atom[] ca2)
			throws StructureException;

	@Override
	abstract public AFPChain align(Atom[] ca1, Atom[] ca2, Object params)
			throws StructureException;

	@Override
	abstract public String getAlgorithmName() ;

	@Override
	abstract public ConfigStrucAligParams getParameters() ;

	@Override
	abstract public String getVersion() ;

	@Override
	abstract public void setParameters(ConfigStrucAligParams parameters);

}
