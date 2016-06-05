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
package org.biojava.nbio.survival.kaplanmeier.figure;

import org.biojava.nbio.survival.cox.StrataInfo;
import org.biojava.nbio.survival.cox.SurvFitInfo;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class NumbersAtRiskPanel extends JPanel {

	private static final long serialVersionUID = 1L;

	KaplanMeierFigure kmf = null;
	Double timePercentage = .20;

	public NumbersAtRiskPanel() {
		this.setPreferredSize(new Dimension(400,100));
		this.setSize(400, 100);
	}

	/**
	 * Pick up needed info and details from the KM Figure
	 * @param kmf
	 */
	public void setKaplanMeierFigure(KaplanMeierFigure kmf) {
		this.kmf = kmf;

		int numRows = kmf.getSurvivalFitInfo().getStrataInfoHashMap().size();
		int height = (numRows + 1) * getFontMetrics(getFont()).getHeight();
		int width = kmf.getWidth();
		setPreferredSize(new Dimension(width,height));
		this.setSize(width, height);

	}

	private void paintTable(Graphics g) {

		if(kmf == null)
			return;
		KMFigureInfo kmfi = kmf.getKMFigureInfo();
		Graphics2D g2 = (Graphics2D) g;
		g2.setStroke(kmfi.kmStroke);
		SurvFitInfo sfi = kmf.getSurvivalFitInfo();

		LinkedHashMap<String, StrataInfo> sfiHashMap = new LinkedHashMap<String, StrataInfo>();
		if(sfi.isWeighted()){
			sfiHashMap = sfi.getUnweightedStrataInfoHashMap();
		}else{
			sfiHashMap = sfi.getStrataInfoHashMap();
		}

		if(sfiHashMap.size() == 0)
			return;
		//int height = this.getHeight();

		int row = 0;
		int left = kmf.getLeft();
		//int right = kmf.getRight();
		//int width = right - left;
		Font f = g2.getFont();
		Font nf = new Font(f.getName(), Font.BOLD, f.getSize());
		g2.setFont(nf);
		FontMetrics fm = getFontMetrics(nf);
		int index = 0;
		int fontHeight = getFontMetrics(getFont()).getHeight();
		int increment = fontHeight;
		ArrayList<Double> xaxisTimeValues = kmf.getxAxisTimeValues();
		ArrayList<Integer> xAxisTimeCoordinates = kmf.getxAxisTimeCoordinates();

		ArrayList<String> labels = new ArrayList<String>(sfiHashMap.keySet());
		Collections.sort(labels);

		for (String group : labels) {
			row = row + increment;
			g2.setColor(kmfi.getColor(index));
			index++;
			g2.drawLine(15, row - fontHeight/2, left - 5, row - fontHeight/2);
			g2.setColor(Color.BLACK);
			StrataInfo si = sfiHashMap.get(group);
			if(kmf.title.toString().equals("[APOBEC1 Transhera Observation Arm]")){
				//int dummy = 1;
			}
//           System.out.println(kmf.title + " Group " +  group);
//           System.out.println(si);
			for(int i = 0; i < xaxisTimeValues.size(); i++){
				Double time = xaxisTimeValues.get(i);
				int xvalue = xAxisTimeCoordinates.get(i);
				Double value = si.getNearestAtRisk(time);
				String nrisk = "";
				if(value == null){
					nrisk = "";
				}else{
					nrisk = String.valueOf(value.intValue());
				}
				if(time == 0.0){
					 g2.drawString(nrisk , xvalue, row);
				}else{
					int w = fm.stringWidth(nrisk );
					g2.drawString(nrisk , xvalue - w/2, row);
				}
			}

		}
	}

	@Override
	protected void paintComponent(Graphics g) {

		super.paintComponent(g); //To change body of generated methods, choose Tools | Templates.
		g.setColor(Color.white);
		g.fillRect(0, 0, this.getWidth(), this.getHeight());

		this.paintTable(g);
	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
	}
}
