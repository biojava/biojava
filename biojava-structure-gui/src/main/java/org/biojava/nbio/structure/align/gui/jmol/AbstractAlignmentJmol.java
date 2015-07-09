package org.biojava.nbio.structure.align.gui.jmol;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JTextField;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.align.util.ResourceManager;
import org.biojava.nbio.structure.jama.Matrix;
import org.jcolorbrewer.ColorBrewer;
import org.jmol.api.JmolViewer;

/**
 * An Abstract Class to generalize the visualization of AFP and 
 * MultipleAlignment structure alignments in Jmol.
 * 
 * @author Aleix Lafita
 *
 */
public abstract class AbstractAlignmentJmol 
implements MouseMotionListener, MouseListener, WindowListener, ActionListener {

	protected Structure structure; 
	protected ColorBrewer colorPalette = ColorBrewer.Spectral;

	protected JmolPanel jmolPanel;
	protected JFrame frame;
	protected JTextField text ;
	protected JTextField status;

	protected static final String COMMAND_LINE_HELP = 
			"enter Jmol scripting command...";

	protected static final int DEFAULT_HEIGHT = 500;
	protected static final int DEFAULT_WIDTH = 500;
	protected static final String DEFAULT_SCRIPT = 
			ResourceManager.getResourceManager("ce").
			getString("default.alignment.jmol.script");

	protected static int nrOpenWindows = 0;

	/**
	 * Display the structures after the variable initialization in the
	 * constructor.
	 */
	protected abstract void initCoords();

	/**
	 * Set all the member variables to null.
	 */
	public void destroy(){
		System.err.println("cleaning up AlignmentJmol window");
		jmolPanel.removeMouseListener(this);
		jmolPanel.removeMouseMotionListener(this);
		jmolPanel.destroy();
	}

	/**
	 * Return to the initial state of the alignment visualization.
	 */
	public abstract void resetDisplay();

	/**
	 * Create and set a new structure from a given atom array.
	 * @param atoms
	 */
	public void setAtoms(Atom[] atoms){
		Structure s = new StructureImpl();
		Chain c = new ChainImpl();
		c.setChainID("A");
		for (Atom a: atoms){
			c.addGroup(a.getGroup());
		}
		s.addChain(c);
		setStructure(s);
	}

	/**
	 * Return the color pattelete of the AlignmentJmol instance.
	 */
	public ColorBrewer getColorPalette() {
		return colorPalette;
	}

	/**
	 * Return the jmolPanel instance of the AlignmentJmol.
	 */
	public JmolPanel getJmolPanel() {
		return jmolPanel;
	}

	/**
	 * Set the jmolPanel of the AlignmentJmol instance.
	 * @param jmolPanel
	 */
	public void setJmolPanel(JmolPanel jmolPanel) {
		this.jmolPanel = jmolPanel;
	}

	/**
	 * Execute a command String in the current Jmol panel.
	 * @param rasmolScript
	 */
	public void evalString(String rasmolScript){
		if ( jmolPanel == null ){
			System.err.println("please install Jmol first");
			return;
		}
		jmolPanel.evalString(rasmolScript);
	}

	/**
	 * Set a new Structure to visualize in the AlignmentJmol window.
	 * @param s
	 */
	public void setStructure(Structure s) {

		if (jmolPanel == null){
			System.err.println("please install Jmol first");
			return;
		}
		setTitle(s.getPDBCode());
		jmolPanel.setStructure(s);

		// actually this is very simple
		// just convert the structure to a PDB file

		//String pdb = s.toPDB();	
		//System.out.println(s.isNmr());

		//System.out.println(pdb);
		// Jmol could also read the file directly from your file system
		//viewer.openFile("/Path/To/PDB/1tim.pdb");

		//System.out.println(pdb);
		//jmolPanel.openStringInline(pdb);

		// send the PDB file to Jmol.
		// there are also other ways to interact with Jmol, 
		// e.g make it directly
		// access the biojava structure object, but they require more
		// code. See the SPICE code repository for how to do this.

		structure = s;
	}

	/**
	 * Return the current Structure in the AlignmentJmol instance.
	 */
	public Structure getStructure(){
		return structure;
	}

	/**
	 * Returns a List of internal Distance Matrices, 
	 * one for each structure in the alignment.
	 * Returns null if no alignment is being displayed.
	 */
	public abstract List<Matrix> getDistanceMatrices();

	/**
	 * Set the title of the AlignmentJmol window.
	 * @param label
	 */
	public void setTitle(String title){
		frame.setTitle(title);
		frame.repaint();
	}

	/**
	 * Return the title of the AlignmentJmol window.
	 */
	public String getTitle(){
		return frame.getTitle();
	}

	@Override
	public void mouseDragged(MouseEvent e) {}

	@Override
	public void mouseMoved(MouseEvent e) {

		JmolViewer viewer = jmolPanel.getViewer();
		int pos = viewer.findNearestAtomIndex( e.getX(), e.getY() );
		if ( pos == -1 ) { return ; }

		String atomInfo = viewer.getAtomInfo(pos);
		text.setText(atomInfo);

	}

	@Override
	public void mouseClicked(MouseEvent e) {}

	@Override
	public void mouseEntered(MouseEvent e) {}

	@Override
	public void mouseExited(MouseEvent e) {}

	@Override
	public void mousePressed(MouseEvent e) {}

	@Override
	public void mouseReleased(MouseEvent e) {

		JmolViewer viewer = jmolPanel.getViewer();
		int pos = viewer.findNearestAtomIndex(e.getX(), e.getY());
		if (pos == -1) return;

		String atomInfo = viewer.getAtomInfo(pos);
		status.setText("clicked: " + atomInfo);
		AtomInfo ai = AtomInfoParser.parse(atomInfo);

		String cmd = "select " + ai.getResidueNumber()+":" 
				+ai.getChainId()+"/"+ai.getModelNumber() 
				+ "; set display selected;";
		evalString(cmd);
	}

	@Override
	public void windowActivated(WindowEvent e) {}

	@Override
	public void windowClosed(WindowEvent e) {}

	@Override
	public void windowClosing(WindowEvent e) {
		destroy();
	}

	@Override
	public void windowDeactivated(WindowEvent e) {}

	@Override
	public void windowDeiconified(WindowEvent e) {}

	@Override
	public void windowIconified(WindowEvent e) {}

	@Override
	public void windowOpened(WindowEvent e) {}

	@Override
	public abstract void actionPerformed(ActionEvent e);

}
