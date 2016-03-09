/*
 *                    PDB web development code
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
 *
 * Created on Jul 28, 2009
 * Created by ap3
 *
 */

package org.biojava.nbio.structure.align.gui;

import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.gui.ScaleableMatrixPanel;
import org.biojava.nbio.structure.jama.Matrix;

import javax.swing.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

/**
 * Shows the interatomic Distance Matrices of all the Structures aligned in different Frames.
 */
public class MyDistMaxListener implements ActionListener{

	AbstractAlignmentJmol parent;

	public MyDistMaxListener(AbstractAlignmentJmol parent){
		this.parent = parent;
	}

	@Override
	public void actionPerformed(ActionEvent a) {

		System.out.println("Show interatomic Distance Matrices");

		if (parent.getDistanceMatrices() == null) {
			System.err.println("Not displaying any alignment currently!");
			return;
		}
		for (int i=0; i<parent.getDistanceMatrices().size(); i++){
			if (parent.getDistanceMatrices().get(i)!=null)
				showMatrix(parent.getDistanceMatrices().get(i), "Internal Distances for Structure "+(i+1));
		}
	}

	private void showMatrix(Matrix m, String title){
		ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
		JFrame frame = new JFrame(title);
		frame.addWindowListener(new WindowAdapter(){
			@Override
			public void windowClosing(WindowEvent e){
	            JFrame f = (JFrame) e.getSource();
	            f.setVisible(false);
	            f.dispose();
			}
		});

		smp.setMatrix(m);
		frame.getContentPane().add(smp);
		frame.pack();
		frame.setVisible(true);
	}
}
