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

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;

import javax.swing.*;

import java.awt.*;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

public class StatusDisplay extends JTextField implements AlignmentPositionListener, WindowListener {

	private static final long serialVersionUID = 6939947266417830429L;
	MultipleAlignment multAln;

	public StatusDisplay(){
		super();
		this.setBackground(Color.white);
		this.setEditable(false);
		this.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
	}
	
	public StatusDisplay(MultipleAlignment multAln){
		this();
		this.multAln = multAln;
	}
	
	public void destroy(){
		multAln = null;
	}

	@Override
	public void mouseOverPosition(AlignedPosition p) {

		if (multAln == null) return;
				
		try {
			String msg = "alig pos";
			for (int str=0; str<multAln.size(); str++) {
			
				String alnseq  = multAln.getAlnSequences().get(str);
				char c = alnseq.charAt(p.getPos1());
		
				Atom a = DisplayAFP.getAtomForAligPos(multAln, str, p.getPos1());
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
		
		if (multAln == null) return;
		
		try {
			String msg = "Clicked pos";
			for (int str=0; str<multAln.size(); str++) {
			
				String alnseq  = multAln.getAlnSequences().get(str);
				char c = alnseq.charAt(p.getPos1());
		
				Atom a = DisplayAFP.getAtomForAligPos(multAln, str, p.getPos1());
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
			for (int str=0; str<multAln.size(); str++) {
				
				String alnseq  = multAln.getAlnSequences().get(str);
				char c1 = alnseq.charAt(start.getPos1());
				char c2 = alnseq.charAt(end.getPos1());
		
				Atom a1 = DisplayAFP.getAtomForAligPos(multAln, str, start.getPos1());
				Atom a2 = DisplayAFP.getAtomForAligPos(multAln, str, end.getPos1());
				
				String pdbInfo1 = JmolTools.getPdbInfo(a1);
				String pdbInfo2 = JmolTools.getPdbInfo(a2);

				msg +=  " range"+str+": " + pdbInfo1 + " ("+c1+") - " + pdbInfo2 + " ("+c2+")";
			}
			this.setText(msg);

		} catch (Exception e){
			e.printStackTrace();
		}
	}

	@Override
	public void selectionLocked() {
		// TODO Auto-generated method stub

	}

	@Override
	public void selectionUnlocked() {
		// TODO Auto-generated method stub

	}

	@Override
	public void windowActivated(WindowEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void windowClosed(WindowEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void windowClosing(WindowEvent e) {
		destroy();
		
	}
	@Override
	public void windowDeactivated(WindowEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void windowDeiconified(WindowEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void windowIconified(WindowEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void windowOpened(WindowEvent e) {
		// TODO Auto-generated method stub
		
	}
}
