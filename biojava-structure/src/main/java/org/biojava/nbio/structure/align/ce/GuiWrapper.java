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
package org.biojava.nbio.structure.align.ce;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.jama.Matrix;

import javax.swing.*;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.List;

/** A class to wrap some of the strucutre.gui classes using Reflection
 *  
 * @author Andreas Prlic
 *
 */

public class GuiWrapper {

	static final String guiPackage = "org.biojava.nbio.structure.gui";

	static final String strucAlignmentDisplay = "org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay";
	static final String displayAFP   = "org.biojava.nbio.structure.align.gui.DisplayAFP" ;
	static final String alignmentGUI = "org.biojava.nbio.structure.align.gui.AlignmentGui";
	static final String strucAligJmol = "org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol";

	static final String scaleMatrixPanel = "org.biojava.nbio.structure.gui.ScaleableMatrixPanel";

	@SuppressWarnings("rawtypes")
	public static boolean isGuiModuleInstalled(){
		String className = displayAFP;
		try {
			@SuppressWarnings("unused")
			Class c = Class.forName(className);
		} catch (ClassNotFoundException ex){
			return false;
		}
		return true;
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static Object display(AFPChain afpChain, Atom[] ca1, Atom[] ca2) 
			throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException{

		Class c = Class.forName(strucAlignmentDisplay);

		Method display = c.getMethod("display", new Class[]{AFPChain.class, Atom[].class, 
				Atom[].class});

		Object structureAlignmentJmol = display.invoke(null, afpChain,ca1,ca2);

		return structureAlignmentJmol;
	}



	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void showAlignmentImage(AFPChain afpChain, Atom[] ca1,
			Atom[] ca2, Object jmol)
					throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException{

		Class structureAlignmentJmol = Class.forName(strucAligJmol);

		Class c = Class.forName(displayAFP);		
		Method show = c.getMethod("showAlignmentImage", new Class[] {AFPChain.class, Atom[].class, Atom[].class, structureAlignmentJmol});

		show.invoke(null,afpChain, ca1, ca2, jmol);
	}


	/** Shows a structure in Jmol
	 * @since 3.0.5
	 */
	public static void showStructure(Structure structure)
			throws ClassNotFoundException, NoSuchMethodException, 
			InvocationTargetException, IllegalAccessException, InstantiationException{

		Class<?> structureAlignmentJmol = Class.forName(strucAligJmol);

		Object strucAligJ = structureAlignmentJmol.newInstance();		

		Method setS = structureAlignmentJmol.getMethod("setStructure", new Class[] {Structure.class});

		setS.invoke(strucAligJ,structure);
	}


	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void showAlignmentGUI()
			throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
		// proxy for AlignmentGui.getInstance();


		Class c = Class.forName(alignmentGUI);
		Method m = c.getMethod("getInstance", (Class[])null);
		m.invoke(c,(Object[])null);
	}

	@SuppressWarnings({"unused" })
	public static Structure getAlignedStructure(Atom[] ca1, Atom[] ca2)
			throws ClassNotFoundException, NoSuchMethodException,
			InvocationTargetException, IllegalAccessException{

		Class<?> structureAlignmentJmol = Class.forName(strucAligJmol);

		Class<?> c = Class.forName(displayAFP);		
		Method show = c.getMethod("getAlignedStructure", new Class[] { Atom[].class, Atom[].class});

		Structure s = (Structure) show.invoke(null, ca1, ca2);

		return s;

	}

	public static JPanel getScaleableMatrixPanel(Matrix m)
			throws ClassNotFoundException, NoSuchMethodException,
			InvocationTargetException, IllegalAccessException, InstantiationException{

		Class<?> scaleMatrixPanelC = Class.forName(scaleMatrixPanel);

		Method setMatrix = scaleMatrixPanelC.getMethod("setMatrix", new Class[] { Matrix.class});

		JPanel panel = (JPanel) scaleMatrixPanelC.newInstance();

		setMatrix.invoke(panel, m);

		return panel;

	}

	/**
	 * 
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return
	 * @throws ClassNotFoundException If an error occurs when invoking jmol
	 * @throws NoSuchMethodException If an error occurs when invoking jmol
	 * @throws InvocationTargetException If an error occurs when invoking jmol
	 * @throws IllegalAccessException If an error occurs when invoking jmol
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static Group[] prepareGroupsForDisplay(AFPChain afpChain, Atom[] ca1,
			Atom[] ca2)
					throws ClassNotFoundException, NoSuchMethodException,
					InvocationTargetException, IllegalAccessException{
		Class c = Class.forName(strucAlignmentDisplay);

		Method display = c.getMethod("prepareGroupsForDisplay", new Class[]{AFPChain.class, Atom[].class, 
				Atom[].class});

		Object groups = display.invoke(null, afpChain,ca1,ca2);

		return (Group[]) groups;
	}

	@SuppressWarnings({ "rawtypes", "unchecked", "unused" })
	public static Atom[] getAtomArray(Atom[] ca, List<Group> hetatoms, List<Group> nucs)
			throws ClassNotFoundException, NoSuchMethodException,
			InvocationTargetException, IllegalAccessException{

		Class structureAlignmentJmol = Class.forName(strucAligJmol);

		Class c = Class.forName(displayAFP);		
		Method show = c.getMethod("getAtomArray", new Class[] { Atom[].class, List.class, List.class});

		Atom[] atoms = (Atom[]) show.invoke(null, ca, hetatoms, nucs);

		return atoms;

	}
	
	/**
	 * @since 3.0.5
	 */
	public static void showDBResults(StartupParameters params) {
		//System.err.println("not implemented full yet");
		
		// We want to do this, but because we don't know if structure-gui.jar is in the classpath we use reflection to hide the calls
		
		UserConfiguration config = UserConfiguration.fromStartupParams(params);
		
		String tableClass = "org.biojava.nbio.structure.align.gui.DBResultTable";
		
		try {
			Class<?> c = Class.forName(tableClass);
			Object table = c.newInstance();
			
			Method show = c.getMethod("show", new Class[]{File.class, UserConfiguration.class });
		
			show.invoke(table, new File(params.getShowDBresult()),config);
			
		} catch (Exception e){
			e.printStackTrace();
			
			System.err.println("Probably structure-gui.jar is not in the classpath, can't show results...");
		}
		
		//DBResultTable table = new DBResultTable();
		
		//table.show(new File(params.getShowDBresult()),config);
		
	}

}
