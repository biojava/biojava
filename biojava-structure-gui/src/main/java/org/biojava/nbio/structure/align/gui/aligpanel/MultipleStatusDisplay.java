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
package org.biojava.nbio.structure.align.gui.aligpanel;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

import javax.swing.JTextField;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentTools;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;

/**
 * This class provides information of the selected positions in the
 * {@link MultipleAligPanel}.
 * <p>
 * It has to be linked to a {@link MultipleAligPanel} in order to obtain
 * the raw information and convert the mouse position to a String.
 *
 * @author Aleix Lafita
 *
 */
public class MultipleStatusDisplay extends JTextField
implements AlignmentPositionListener, WindowListener {

	private static final long serialVersionUID = 6939947266417830429L;
	private MultipleAligPanel panel;

	private MultipleStatusDisplay(){
		super();
		this.setBackground(Color.white);
		this.setEditable(false);
		this.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
	}

	/**
	 * Constructor using an alignment panel as a parent, where the
	 * information will be displayed.
	 *
	 * @param panel
	 */
	public MultipleStatusDisplay(MultipleAligPanel panel) {
		this();
		this.panel = panel;
	}

	public void destroy(){
		panel = null;
	}

	@Override
	public void mouseOverPosition(AlignedPosition p) {

		if (panel == null) return;

		try {
			String msg = "alig pos";
			for (int str=0; str<panel.size; str++) {

				String alnseq  = panel.getAlnSequences().get(str);
				char c = alnseq.charAt(p.getPos1());

				Atom a = MultipleAlignmentTools.getAtomForSequencePosition(
						panel.getMultipleAlignment(),
						panel.getMapSeqToStruct(),
						str, p.getPos1());

				String pdbInfo = JmolTools.getPdbInfo(a);
				msg += ": "+pdbInfo + " ("+c+") ";
			}
			this.setText(msg);
			this.repaint();

		} catch (Exception e){
			e.printStackTrace();
		}
	}

	@Override
	public void positionSelected(AlignedPosition p) {
		mouseOverPosition(p);
	}

	@Override
	public void toggleSelection(AlignedPosition p) {

		if (panel == null) return;

		try {
			String msg = "Clicked pos";
			for (int str=0; str<panel.size; str++) {

				String alnseq  = panel.getAlnSequences().get(str);
				char c = alnseq.charAt(p.getPos1());

				Atom a = MultipleAlignmentTools.getAtomForSequencePosition(
						panel.getMultipleAlignment(),
						panel.getMapSeqToStruct(),
						str, p.getPos1());

				String pdbInfo = JmolTools.getPdbInfo(a);

				msg += ": "+pdbInfo + " ("+c+") ";
			}
			this.setText(msg);

		} catch (Exception e){
			e.printStackTrace();
		}
	}

	@Override
	public void rangeSelected(AlignedPosition start, AlignedPosition end) {

		try {
			String msg =  "Selected:";
			for (int str=0; str<panel.size; str++) {

				String alnseq  = panel.getAlnSequences().get(str);
				char c1 = alnseq.charAt(start.getPos1());
				char c2 = alnseq.charAt(end.getPos1());

				Atom a1 = MultipleAlignmentTools.getAtomForSequencePosition(
						panel.getMultipleAlignment(),
						panel.getMapSeqToStruct(),
						str, start.getPos1());

				Atom a2 = MultipleAlignmentTools.getAtomForSequencePosition(
						panel.getMultipleAlignment(),
						panel.getMapSeqToStruct(),
						str, end.getPos1());

				String pdbInfo1 = JmolTools.getPdbInfo(a1);
				String pdbInfo2 = JmolTools.getPdbInfo(a2);

				msg +=  " range"+str+": " + pdbInfo1 +
						" ("+c1+") - " + pdbInfo2 + " ("+c2+")";
			}
			this.setText(msg);

		} catch (Exception e){
			e.printStackTrace();
		}
	}

	@Override
	public void windowClosing(WindowEvent e) {
		destroy();
	}
	@Override
	public void selectionLocked() {}
	@Override
	public void selectionUnlocked() {}
	@Override
	public void windowActivated(WindowEvent e) {}
	@Override
	public void windowClosed(WindowEvent e) {}
	@Override
	public void windowDeactivated(WindowEvent e) {}
	@Override
	public void windowDeiconified(WindowEvent e) {}
	@Override
	public void windowIconified(WindowEvent e) {}
	@Override
	public void windowOpened(WindowEvent e) {}
}
