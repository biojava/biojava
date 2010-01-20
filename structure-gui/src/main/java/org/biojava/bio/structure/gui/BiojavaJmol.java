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
 * Created on 24.05.2004
 * @author Andreas Prlic
 *
 */

package org.biojava.bio.structure.gui;


import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JTextField;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.gui.util.MenuCreator;
import org.biojava.bio.structure.io.PDBFileReader;


/** A class that provides a simple GUI for Jmol
 * 
 * @author Andreas Prlic
 * @since 1.6
 *
 *
 *
 */
public class BiojavaJmol {

	public static final String viewer       = "org.jmol.api.JmolSimpleViewer";
	public static final String adapter      = "org.jmol.api.JmolAdapter";
	public static final String smartAdapter = "org.jmol.adapter.smarter.SmarterJmolAdapter";

	Structure structure; 

	JmolPanel jmolPanel;
	JFrame frame ;
	
	
	public static void main(String[] args){
		try {

			PDBFileReader pdbr = new PDBFileReader();   
			//pdbr.setAutoFetch(true);
			pdbr.setPath("/Users/andreas/WORK/PDB/");

			String pdbCode = "5pti";

			Structure struc = pdbr.getStructureById(pdbCode);

			BiojavaJmol jmolPanel = new BiojavaJmol();

			jmolPanel.setStructure(struc);

			// send some RASMOL style commands to Jmol
			jmolPanel.evalString("select * ; color chain;");
			jmolPanel.evalString("select *; spacefill off; wireframe off; backbone 0.4;  ");

		} catch (Exception e){
			e.printStackTrace();
		}
	}

	
	

	public BiojavaJmol() {		
		
		frame = new JFrame();
		
		JMenuBar menu = MenuCreator.initMenu();
		

		frame.setJMenuBar(menu);
		
		frame.addWindowListener( new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				frame.dispose();
				//System.exit(0);
			}
		});
		
		Container contentPane = frame.getContentPane();
				
		Box vBox = Box.createVerticalBox();

		try {
			
			jmolPanel = new JmolPanel();
			
		} catch (ClassNotFoundException e){
			e.printStackTrace();
			System.err.println("Could not find Jmol in classpath, please install first. http://www.jmol.org");
			return;
		}
		jmolPanel.setPreferredSize(new Dimension(500,500));
		vBox.add(jmolPanel);


		JTextField field = new JTextField();

		field.setMaximumSize(new Dimension(Short.MAX_VALUE,30));   
		field.setText("enter RASMOL like command...");
		RasmolCommandListener listener = new RasmolCommandListener(jmolPanel,field) ;

		field.addActionListener(listener);
		field.addMouseListener(listener);
		field.addKeyListener(listener);
		vBox.add(field);

		contentPane.add(vBox);



		frame.pack();
		frame.setVisible(true); 

	}

	
	
	/** returns true if Jmol can be found in the classpath, otherwise false.
	 * 
	 * @return true/false depending if Jmol can be found
	 */
	public static boolean jmolInClassPath(){
		try {
			Class.forName(viewer);		
		} catch (ClassNotFoundException e){
			e.printStackTrace();			
			return false;
		}
		return true;
	}

	public void evalString(String rasmolScript){
		if ( jmolPanel == null ){
			System.err.println("please install Jmol first");
			return;
		}
		jmolPanel.evalString(rasmolScript);
	}

	public void setStructure(Structure s) {

		if ( jmolPanel == null ){
			System.err.println("please install Jmol first");
			return;
		}
		
		setTitle(s.getPDBCode());
		
		// actually this is very simple
		// just convert the structure to a PDB file

		String pdb = s.toPDB();	
		//System.out.println(s.isNmr());

		//System.out.println(pdb);
		// Jmol could also read the file directly from your file system
		//viewer.openFile("/Path/To/PDB/1tim.pdb");

		//System.out.println(pdb);
		jmolPanel.openStringInline(pdb);

		// send the PDB file to Jmol.
		// there are also other ways to interact with Jmol, e.g make it directly
		// access the biojava structure object, but they require more
		// code. See the SPICE code repository for how to do this.




	}

	public void setTitle(String label){
		frame.setTitle(label);
		frame.repaint();
	}




	

	static class JmolPanel extends JPanel {
		/**
		 * 
		 */
		private static final long serialVersionUID = -3661941083797644242L;

		Class viewerC;
		Class adapterC;
		Class smartAdapterC;

		Object viewerO;
		Object adapterO;

		Method evalString;
		Method renderScreenImage;
		Method openStringInline;

		//JmolSimpleViewer viewer;
		//JmolAdapter adapter;
		JmolPanel() throws ClassNotFoundException {

			try {
				viewerC = Class.forName(viewer);

				adapterC = Class.forName(adapter);
				smartAdapterC = Class.forName(smartAdapter);

				Method m = viewerC.getMethod("allocateSimpleViewer", new Class[]{Component.class,adapterC});

				Constructor constructor = smartAdapterC.getConstructor(new Class[]{});
				adapterO = constructor.newInstance(new Object[]{});

				//viewerC = JmolSimpleViewer.allocateSimpleViewer(this, adapter);
				viewerO = m.invoke(viewerC, this, adapterO);

				evalString = viewerC.getMethod("evalString",String.class);

				renderScreenImage = viewerC.getMethod("renderScreenImage",
						new Class[]{Graphics.class,Dimension.class, Rectangle.class});

				openStringInline = viewerC.getMethod("openStringInline", new Class[]{String.class});

			} catch (InstantiationException e){
				e.printStackTrace();
			} catch (NoSuchMethodException e){
				e.printStackTrace();        		
			} catch ( InvocationTargetException e){
				e.printStackTrace();
			} catch ( IllegalAccessException e){
				e.printStackTrace();
			}
			
			evalString("set scriptQueue on;");

		}

		public Class getViewer() {
			return viewerC;
		}

		public void evalString(String rasmolScript){
			try {
				evalString.invoke(viewerO, rasmolScript);
			} catch (Exception e){
				e.printStackTrace();
			}
		}

		public void openStringInline(String pdbFile){
			try {
				openStringInline.invoke(viewerO, pdbFile);
			} catch (Exception e){
				e.printStackTrace();
			}
		}

		public void executeCmd(String rasmolScript){
			try {
				evalString.invoke(viewerO, rasmolScript);
			} catch (Exception e){
				e.printStackTrace();
			}
		}


		final Dimension currentSize = new Dimension();
		final Rectangle rectClip = new Rectangle();

		public void paint(Graphics g) {
			getSize(currentSize);
			g.getClipBounds(rectClip);
			//viewer.renderScreenImage(g, currentSize, rectClip);

			try {
				renderScreenImage.invoke(viewerO,g,currentSize,rectClip);
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}
	
	
}