/*
 *                  BioJava development code
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
 * Created on Aug 3, 2007
 */
package org.biojava.bio.structure.gui.util.color;

import java.awt.BasicStroke;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JPanel;

public class GradientPanel extends JPanel {
	private static final long serialVersionUID = -6387922432121206731L;
	private ContinuousColorMapper mapper;
	private double min, max;
	
	
	public GradientPanel(ContinuousColorMapper mapper, double min, double max) {
		this.min = min; 
		this.max = max;
		this.mapper = mapper;
		this.setPreferredSize(new Dimension(100,20));
	}
	
	public void paintComponent(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;
		int w = getWidth();
		int h = getHeight();
		
		g2.setStroke(new BasicStroke(1.0f));
		for(int i=0;i<w;i++) {
			double val = (max-min)*i/w+min;
			g2.setColor(mapper.getColor(val));
			g2.drawLine(i, 0, i, h);
		}
	}
	
}