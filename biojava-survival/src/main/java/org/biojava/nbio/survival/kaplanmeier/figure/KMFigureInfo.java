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

import java.awt.*;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Properties;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class KMFigureInfo {

	/**
	 *
	 */
	public final int titleHeight = 40;
	/**
	 *
	 */
	public final int padding = 20;
	/**
	 *
	 */
	public final Integer width = 600;
	/**
	 *
	 */
	public final Integer height = 400;
	/**
	 *
	 */
	public final double timeScale = 1.0; //multiplier to change time from days to months etc
	/**
	 *
	 */
	public final double yaxisPercentIncrement = .2;
	/**
	 *
	 */
	public final double xaxisPercentIncrement = .25;
	/**
	 *
	 */
	public double legendUpperPercentX = .95;
	/**
	 *
	 */
	public double legendUpperPercentY = .95;
	/**
	 *
	 */
	public final double figureLineInfoLowerPercentX = .01;
	/**
	 *
	 */
	public final double figureLineInfoLowerPercentY = .01;
	/**
	 *
	 */
	public final BasicStroke axisStroke = new BasicStroke(2);
	/**
	 *
	 */
	public final BasicStroke kmStroke = new BasicStroke(3);
	/**
	 *
	 */
	public final Color[] legendColor = {Color.RED, Color.BLUE, Color.GREEN, Color.CYAN, Color.ORANGE, Color.YELLOW, Color.MAGENTA, Color.PINK};
	public final ArrayList<Double> xAxisLabels = new ArrayList<>();//new ArrayList<Double>(Arrays.asList(0.0, 5.0, 10.0, 15.0, 20.0));
	public String xAxisLegend = "";
	public String yAxisLegend = "";
	public Color getColor(int index) {
		return legendColor[index];
	}

	/**
	 *
	 * @param propertyFileName
	 * @throws Exception
	 */
	public void init(String propertyFileName) throws Exception {
		Properties properties = new Properties();
		properties.load(new FileInputStream(propertyFileName));

		if (properties.containsKey("legendUpperPercentX")) {
			legendUpperPercentX = Double.parseDouble(properties.getProperty("legendUpperPercentX"));
		}

		if (properties.containsKey("legendUpperPercentY")) {
			legendUpperPercentY = Double.parseDouble(properties.getProperty("legendUpperPercentY"));
		}

		if (properties.containsKey("xAxisLabels")) {
			String values = properties.getProperty("xAxisLabels").trim();
			if (values.startsWith("[") || values.startsWith("(") || values.startsWith("{")) {
				values = values.substring(1).trim();
			}
			if (values.endsWith("]") || values.endsWith(")") || values.endsWith("}")) {
				values = values.substring(0, values.length() - 1).trim();
			}
			xAxisLabels.clear();
			String[] data = values.split(",");
			for (String d : data) {
				try {
					Double v = Double.parseDouble(d.trim());
					xAxisLabels.add(v);
				} catch (Exception e) {
				}
			}
		}

		if (properties.containsKey("xAxisLegend")) {
			xAxisLegend = properties.getProperty("xAxisLegend");
		}

				if (properties.containsKey("yAxisLegend")) {
			yAxisLegend = properties.getProperty("yAxisLegend");
		}
	}
}
