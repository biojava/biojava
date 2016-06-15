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
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

public class StatusDisplay extends JTextField implements AlignmentPositionListener, WindowListener  {

	/**
	 *
	 */
	private static final long serialVersionUID = 6939947266417830429L;

	AFPChain afpChain;
	Atom[] ca1;
	Atom[] ca2;



	public StatusDisplay(){
		super();

		this.setBackground(Color.white);
		this.setEditable(false);
		this.setMaximumSize(new Dimension(Short.MAX_VALUE,30));

	}
	public void destroy(){
		afpChain = null;
		ca1= null;
		ca2 = null;

	}

	@Override
	public void mouseOverPosition(AlignedPosition p) {

		if ( afpChain == null)
			return;

		char[] aligs1  = afpChain.getAlnseq1();
		char[] aligs2  = afpChain.getAlnseq2();

		char c1 = aligs1[p.getPos1()];
		char c2 = aligs2[p.getPos2()];

		try {
			Atom a1 = DisplayAFP.getAtomForAligPos(afpChain, 0, p.getPos1(), ca1,false);
			Atom a2 = DisplayAFP.getAtomForAligPos(afpChain, 1, p.getPos2(), ca2,true);

			String pdbInfo1 = JmolTools.getPdbInfo(a1);
			String pdbInfo2 = JmolTools.getPdbInfo(a2);

			String msg = "alig pos:" + p.getPos1()+ " " +  pdbInfo1 + " ("+c1+") : " + pdbInfo2 + " ("+c2+")";

			this.setText(msg);



		} catch (StructureException e){
			e.printStackTrace();
		}

		this.repaint();

	}



	@Override
	public void positionSelected(AlignedPosition p) {
		mouseOverPosition(p);

	}
	@Override
	public void toggleSelection(AlignedPosition p) {
		if ( afpChain == null)
			return;

		char[] aligs1  = afpChain.getAlnseq1();
		char[] aligs2  = afpChain.getAlnseq2();

		char c1 = aligs1[p.getPos1()];
		char c2 = aligs2[p.getPos2()];

		try {
			Atom a1 = DisplayAFP.getAtomForAligPos(afpChain, 0, p.getPos1(), ca1,false);
			Atom a2 = DisplayAFP.getAtomForAligPos(afpChain, 1, p.getPos2(), ca2,true);

			String pdbInfo1 = JmolTools.getPdbInfo(a1);
			String pdbInfo2 = JmolTools.getPdbInfo(a2);

			String msg = "Clicked pos:" + p.getPos1()+ " " + pdbInfo1 + " ("+c1+") : " + pdbInfo2 + " ("+c2+")";

			this.setText(msg);
		} catch (StructureException e){
			e.printStackTrace();
		}

	}


	@Override
	public void rangeSelected(AlignedPosition start, AlignedPosition end) {
		char[] aligs1  = afpChain.getAlnseq1();
		char[] aligs2  = afpChain.getAlnseq2();

		char c1 = aligs1[start.getPos1()];
		char c3 = aligs1[end.getPos1()];

		char c2 = aligs2[start.getPos2()];
		char c4 = aligs2[end.getPos2()];

		try {
			Atom a1 = DisplayAFP.getAtomForAligPos(afpChain, 0, start.getPos1(), ca1,false);
			Atom a2 = DisplayAFP.getAtomForAligPos(afpChain, 1, start.getPos2(), ca2,true);

			Atom a3 = DisplayAFP.getAtomForAligPos(afpChain, 0, end.getPos1(), ca1,false);
			Atom a4 = DisplayAFP.getAtomForAligPos(afpChain, 1, end.getPos2(), ca2,true);

			String pdbInfo1 = JmolTools.getPdbInfo(a1);
			String pdbInfo2 = JmolTools.getPdbInfo(a2);

			String pdbInfo3 = JmolTools.getPdbInfo(a3);
			String pdbInfo4 = JmolTools.getPdbInfo(a4);

			String msg =  "Selected range1: " + pdbInfo1 + " ("+c1+") - " + pdbInfo3 + " ("+c3+")";
			msg       +=  " range2: "         + pdbInfo2 + " ("+c2+") - " + pdbInfo4 + " ("+c4+")";


			this.setText(msg);
		} catch (StructureException e){
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

	public AFPChain getAfpChain() {
		return afpChain;
	}

	public void setAfpChain(AFPChain afpChain) {
		this.afpChain = afpChain;
	}

	public Atom[] getCa1() {
		return ca1;
	}

	public void setCa1(Atom[] ca1) {
		this.ca1 = ca1;
	}

	public Atom[] getCa2() {
		return ca2;
	}

	public void setCa2(Atom[] ca2) {
		this.ca2 = ca2;
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
