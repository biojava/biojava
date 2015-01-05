package org.biojava.bio.structure.align;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.model.AFPChain;

public abstract class AbstractStructureAlignment implements StructureAlignment {

	public static String newline = System.getProperty("line.separator");

	@Override
	abstract public  AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException ;

	@Override
	abstract public AFPChain align(Atom[] ca1, Atom[] ca2, Object params) throws StructureException;

	@Override
	abstract public String getAlgorithmName() ;

	@Override
	abstract public ConfigStrucAligParams getParameters() ;

	@Override
	abstract public String getVersion() ;

	@Override
	abstract public void setParameters(ConfigStrucAligParams parameters);


}
