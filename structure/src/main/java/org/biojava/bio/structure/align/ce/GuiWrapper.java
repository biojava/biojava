package org.biojava.bio.structure.align.ce;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.align.model.AFPChain;

/** A class to wrap some of the strucutre.gui classes using Reflection
 *  
 * @author Andreas Prlic
 *
 */
public class GuiWrapper {

	static final String guiPackage = "org.biojava.bio.structure.gui";
	
	static final String strucAlignmentDisplay = "org.biojava.bio.structure.align.gui.StructureAlignmentDisplay";
	static final String displayAFP   = "org.biojava.bio.structure.align.gui.DisplayAFP" ;
	static final String alignmentGUI = "org.biojava.bio.structure.align.gui.AlignmentGui";
	static final String strucAligJmol = "org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol";
	
	public static boolean isGuiModuleInstalled(){
		String className = alignmentGUI;
		try {
			Class c = Class.forName(className);
		} catch (ClassNotFoundException ex){
			return false;
		}
		return true;
	}
	
	@SuppressWarnings("unchecked")
	public static Object display(AFPChain afpChain, Atom[] ca1, Atom[] ca2,
			List<Group> hetatms1, List<Group> nucs1, List<Group> hetatms2,
			List<Group> nucs2) 
	throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException{
	
		Class c = Class.forName(strucAlignmentDisplay);
				
		Method display = c.getMethod("display", new Class[]{AFPChain.class, Atom[].class, 
			Atom[].class, List.class, List.class, List.class,List.class});
			
		Object structureAlignmentJmol = display.invoke(null, afpChain,ca1,ca2, hetatms1, nucs1, hetatms2, nucs2);
		
		return structureAlignmentJmol;
	}

	@SuppressWarnings("unchecked")
	public static void showAlignmentImage(AFPChain afpChain, Atom[] ca1,
			Atom[] ca2, Object jmol)
	throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException{
		
	   Class structureAlignmentJmol = Class.forName(strucAligJmol);
	   
		Class c = Class.forName(displayAFP);		
		Method show = c.getMethod("showAlignmentImage", new Class[] {AFPChain.class, Atom[].class, Atom[].class, structureAlignmentJmol});
		
		show.invoke(null,afpChain, ca1, ca2, jmol);
	}

	@SuppressWarnings("unchecked")
	public static void showAlignmentGUI()
	throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
		// proxy for AlignmentGui.getInstance();
		
		
		Class c = Class.forName(alignmentGUI);
		Method m = c.getMethod("getInstance", null);
		m.invoke(c,null);
	}
}
