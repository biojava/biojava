package org.biojava.bio.structure.align;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.model.AFPChain;

public abstract class AbstractStructureAlignment implements StructureAlignment {

	public static String newline = System.getProperty("line.separator");

	abstract public  AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException ;

	abstract public AFPChain align(Atom[] ca1, Atom[] ca2, Object params) throws StructureException;

	abstract public String getAlgorithmName() ;

	abstract public ConfigStrucAligParams getParameters() ;

	abstract public String getVersion() ;

	abstract public void setParameters(ConfigStrucAligParams parameters);


}
