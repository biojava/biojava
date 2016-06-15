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

import org.biojava.nbio.survival.cox.*;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class KaplanMeierFigure extends JPanel {

	private static final long serialVersionUID = 1L;

	ArrayList<String> title = new ArrayList<String>();
	/**
	 *
	 */
	private int top;
	/**
	 *
	 */
	private int bottom;
	/**
	 *
	 */
	private int left;
	private int yaxisLabel = 20;
	/**
	 *
	 */
	private int right;
	int titleHeight;
	int xAxisLabelHeight;
	int labelWidth;
	double minTime = 0.0;
	double maxTime = 10.0;
	double minPercentage = 0.0;
	double maxPercentage = 1.0;
	FontMetrics fm;
	KMFigureInfo kmfi = new KMFigureInfo();
	LinkedHashMap<String, ArrayList<CensorStatus>> survivalData = new LinkedHashMap<String, ArrayList<CensorStatus>>();
	ArrayList<String> lineInfoList = new ArrayList<String>();
	SurvFitInfo sfi = new SurvFitInfo();
	private String fileName = "";
	private ArrayList<Double> xAxisTimeValues = new ArrayList<Double>();
	private ArrayList<Integer> xAxisTimeCoordinates = new ArrayList<Integer>();

	/**
	 *
	 */
	public KaplanMeierFigure() {
		super();
		setSize(500, 400);
		setBackground(Color.WHITE);
	}

	/**
	 * Get the name of the groups that are being plotted in the figure
	 *
	 * @return
	 */
	public ArrayList<String> getGroups() {
		return new ArrayList<String>(survivalData.keySet());
	}

	/**
	 * To get the median percentile for a particular group pass the value of
	 * .50.
	 *
	 * @param group
	 * @param percentile
	 * @return
	 */
	public Double getSurvivalTimePercentile(String group, double percentile) {

		StrataInfo si = sfi.getStrataInfoHashMap().get(group);
		ArrayList<Double> percentage = si.getSurv();
		Integer percentileIndex = null;
		for (int i = 0; i < percentage.size(); i++) {
			if (percentage.get(i) == percentile) {
				if (i + 1 < percentage.size()) {
					percentileIndex = i + 1;
				}
				break;
			} else if (percentage.get(i) < percentile) {
				percentileIndex = i;
				break;
			}
		}
		if (percentileIndex != null) {
			return si.getTime().get(percentileIndex);
		} else {
			return null;
		}
	}

	/**
	 *
	 * @param kmfi
	 */
	public void setKMFigureInfo(KMFigureInfo kmfi) {
		this.kmfi = kmfi;
		if (kmfi.width != null && kmfi.height != null) {
			this.setSize(kmfi.width, kmfi.height);
		}
	}

	public KMFigureInfo getKMFigureInfo() {
		return kmfi;
	}

	/**
	 *
	 * @param lineInfoList
	 */
	public void setFigureLineInfo(ArrayList<String> lineInfoList) {
		this.lineInfoList = lineInfoList;
		this.repaint();
	}

	/**
	 *
	 * @param title Title of figures
	 * @param ci
	 * @param strataVariable The column that based on value will do a figure
	 * line
	 * @param legendMap Map the value in the column to something readable
	 * @param useWeighted
	 * @throws Exception
	 */
	public void setCoxInfo(ArrayList<String> title, CoxInfo ci, String strataVariable, LinkedHashMap<String, String> legendMap, Boolean useWeighted) throws Exception {
		LinkedHashMap<String, ArrayList<CensorStatus>> survivalData = new LinkedHashMap<String, ArrayList<CensorStatus>>();
		ArrayList<SurvivalInfo> siList = ci.getSurvivalInfoList();
		int n = 0;
		int event = 0;
		for (SurvivalInfo si : siList) {
			String strata = si.getOriginalMetaData(strataVariable);
			String legend = legendMap.get(strata);
			if (legend == null) {

				legend = strata;
			}
			ArrayList<CensorStatus> censorStatusList = survivalData.get(legend);
			if (censorStatusList == null) {
				censorStatusList = new ArrayList<CensorStatus>();
				survivalData.put(legend, censorStatusList);
			}
			CensorStatus cs = new CensorStatus(strata, si.getTime(), si.getStatus() + "");
			cs.weight = si.getWeight();
			censorStatusList.add(cs);
			n++;
			if (si.getStatus() == 1) {
				event++;
			}
		}

		setSurvivalData(title, survivalData, useWeighted);
		CoxCoefficient cc = ci.getCoefficient(strataVariable);
		//DecimalFormat df = new DecimalFormat("#.##");
		String line1 = "HR=" + fmt(cc.getHazardRatio(), 2, 0) + " (CI:" + fmt(cc.getHazardRatioLoCI(), 2, 0) + "-" + fmt(cc.getHazardRatioHiCI(), 2, 0) + ")";
		String line2 = "p=" + fmt(cc.getPvalue(), 3, 0);
	   // String line2 = "logrank P=" + fmt(ci.getScoreLogrankTestpvalue(), 3, 0);
		String line3 = "n=" + n + " events=" + event;
//        System.out.println("setCoxInfo=" + cc.pvalue + " " + title);


		ArrayList<String> lines = new ArrayList<String>();
		lines.add(line1);
		lines.add(line2);
		lines.add(line3);
		setFigureLineInfo(lines);
	}

	/**
	 *
	 * @param d
	 * @param precision
	 * @param pad
	 * @return
	 */
	public static String fmt(Double d, int precision, int pad) {
		String value = "";
		DecimalFormat dfe = new DecimalFormat("0.00E0");
		String dpad = "0.";
		double p = 1.0;
		for (int i = 0; i < (precision); i++) {
			dpad = dpad + "0";
			p = p / 10.0;
		}
		DecimalFormat df = new DecimalFormat(dpad);
		if (Math.abs(d) >= p) {
			value = df.format(d);
		} else {
			value = dfe.format(d);
		}
		int length = value.length();
		int extra = pad - length;
		if (extra > 0) {
			for (int i = 0; i < extra; i++) {
				value = " " + value;
			}
		}
		return value;
	}

	/**
	 *
	 * @return
	 */
	public SurvFitInfo getSurvivalFitInfo() {
		return sfi;
	}

	/**
	 * Allow setting of points in the figure where weighted correction has been
	 * done and percentage has already been calculated.
	 *
	 * @param title
	 * @param sfi
	 * @param userSetMaxTime
	 */
	public void setSurvivalData(ArrayList<String> title, SurvFitInfo sfi, Double userSetMaxTime) {
		this.title = title;
		LinkedHashMap<String, StrataInfo> strataInfoHashMap = sfi.getStrataInfoHashMap();
		Double mTime = null;
		for (StrataInfo si : strataInfoHashMap.values()) {
			for (double t : si.getTime()) {
				if (mTime == null || t > mTime) {
					mTime = t;
				}
			}
		}

		int evenCheck = Math.round(mTime.floatValue());
		if (evenCheck % 2 == 1) {
			evenCheck = evenCheck + 1;
		}
		this.maxTime = evenCheck;

		if (userSetMaxTime != null && userSetMaxTime > maxTime) {
			this.maxTime = userSetMaxTime;
		}
		this.sfi = sfi;
		if (sfi.getStrataInfoHashMap().size() == 1) {
			return;
		}
		this.repaint();
	}

	/**
	 * The data will set the max time which will result in off time points for
	 * tick marks
	 *
	 * @param title
	 * @param survivalData
	 * @param useWeighted
	 * @throws Exception
	 */
	public void setSurvivalData(ArrayList<String> title, LinkedHashMap<String, ArrayList<CensorStatus>> survivalData, Boolean useWeighted) throws Exception {
		this.setSurvivalData(title, survivalData, null, useWeighted);
	}

	/**
	 *
	 * @param title
	 * @param survivalData
	 * @param userSetMaxTime
	 * @param useWeighted
	 * @throws Exception
	 */
	public void setSurvivalData(ArrayList<String> title, LinkedHashMap<String, ArrayList<CensorStatus>> survivalData, Double userSetMaxTime, Boolean useWeighted) throws Exception {
		this.title = title;
		this.survivalData = survivalData;
		Double mTime = null;
		ArrayList<String> labels = new ArrayList<String>(survivalData.keySet());
		Collections.sort(labels);
		for (String legend : labels) {
			ArrayList<CensorStatus> censorStatusList = survivalData.get(legend);
			for (CensorStatus cs : censorStatusList) {

				if (mTime == null || cs.time > mTime) {
					mTime = cs.time;
				}
			}
		}

		int evenCheck = Math.round(mTime.floatValue());
		if (evenCheck % 2 == 1) {
			evenCheck = evenCheck + 1;
		}
		this.maxTime = evenCheck;

		if (userSetMaxTime != null && userSetMaxTime > maxTime) {
			this.maxTime = userSetMaxTime;
		}

		//calculate percentages
		SurvFitKM survFitKM = new SurvFitKM();
		sfi = survFitKM.process(survivalData, useWeighted);
		this.repaint();
	}

	/**
	 * Save data from survival curve to text file
	 *
	 * @param fileName
	 * @throws Exception
	 */
	public void saveSurvivalData(String fileName) throws Exception {
		FileWriter fw = new FileWriter(fileName);
		fw.write("index\tTIME\tSTATUS\tGROUP\r\n");
		int index = 0;
		for (String group : survivalData.keySet()) {
			ArrayList<CensorStatus> sd = survivalData.get(group);
			for (CensorStatus cs : sd) {
				String line = index + "\t" + cs.time + "\t" + cs.censored + "\t" + cs.group + "\r\n";
				index++;
				fw.write(line);
			}
		}
		fw.close();
	}
	DecimalFormat df = new DecimalFormat("#.#");

	@Override
	public void paintComponent(Graphics g) // draw graphics in the panel
	{
		int width = getWidth();             // width of window in pixels
		int height = getHeight();           // height of window in pixels
		setFigureDimensions();
		g.setColor(Color.white);
		g.clearRect(0, 0, width, height);

		super.paintComponent(g);            // call superclass to make panel display correctly

		drawLegend(g);
		drawSurvivalCurves(g);
		drawFigureLineInfo(g);
		// Drawing code goes here
	}

	private void drawFigureLineInfo(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;
		setRenderingHints(g2);
		g2.setColor(Color.BLACK);
		fm = getFontMetrics(getFont());
		int yoffset = fm.getHeight() * lineInfoList.size();

		int x = getTimeX(kmfi.figureLineInfoLowerPercentX * maxTime);
		int y = getPercentageY(kmfi.figureLineInfoLowerPercentY) - yoffset;

		for (String line : lineInfoList) {
			g2.drawString(line, x, y);
			y = y + fm.getHeight();
		}

	}

	private void drawSurvivalCurves(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;
		setRenderingHints(g2);
		g2.setStroke(kmfi.kmStroke);


		int colorIndex = 0;
		ArrayList<String> labels = new ArrayList<String>(sfi.getStrataInfoHashMap().keySet());
		Collections.sort(labels);

		LinkedHashMap<String, StrataInfo> strataInfoHashMap = sfi.getStrataInfoHashMap();

		for (String legend : labels) {
			StrataInfo si = strataInfoHashMap.get(legend);
			g2.setColor(kmfi.legendColor[colorIndex]);
			colorIndex++;

			for (int i = 0; i < si.getSurv().size() - 1; i++) {
				double p0time = si.getTime().get(i);
				double p1time = si.getTime().get(i + 1);
				double p0percentage = si.getSurv().get(i);
				double p1percentage = si.getSurv().get(i + 1);
				if (i == 0) {
					g2.drawLine(getTimeX(0), getPercentageY(1), getTimeX(p0time), getPercentageY(1));
					g2.drawLine(getTimeX(p0time), getPercentageY(1), getTimeX(p0time), getPercentageY(p0percentage));
				}
				g2.drawLine(getTimeX(p0time), getPercentageY(p0percentage), getTimeX(p1time), getPercentageY(p0percentage));

				g2.drawLine(getTimeX(p1time), getPercentageY(p0percentage), getTimeX(p1time), getPercentageY(p1percentage));
				// if (si.getStatus().get(i) == 0) {
				if (i > 0 && si.getNcens().get(i) > 0) {
					g2.drawLine(getTimeX(p0time), getPercentageY(p0percentage) - 4, getTimeX(p0time), getPercentageY(p0percentage) + 4);
					g2.drawLine(getTimeX(p0time) - 4, getPercentageY(p0percentage), getTimeX(p0time) + 4, getPercentageY(p0percentage));
				}
			}


		}

		String maxString = "";
		for (String legend : labels) {
			if (legend.length() > maxString.length()) {
				maxString = legend;
			}
		}

		int offset = fm.stringWidth(maxString);
		int x = getTimeX(kmfi.legendUpperPercentX * maxTime) - offset;
		int y = getPercentageY(kmfi.legendUpperPercentY);

		colorIndex = 0;
		for (String legend : labels) {
			g2.setColor(kmfi.legendColor[colorIndex]);
			colorIndex++;
			g2.drawLine(x - 20, y - (fm.getHeight() / 3), x - 5, y - (fm.getHeight() / 3));
			g2.drawString(legend, x, y);
			y = y + fm.getHeight();
		}


	}

	/**
	 * Get the X coordinate based on a time value
	 *
	 * @param value
	 * @return
	 */
	private int getTimeX(double value) {
		double d = left + (((right - left) * value) / (maxTime - minTime));
		return (int) d;
	}

	/**
	 * Get the Y coordinate based on percent value 0.0-1.0
	 *
	 * @param value
	 * @return
	 */
	private int getPercentageY(double value) {
		value = 1.0 - value;
		double d = top + (((bottom - top) * value) / (maxPercentage - minPercentage));
		return (int) d;
	}

	/**
	 * @return the fileName
	 */
	public String getFileName() {
		return fileName;
	}

	/**
	 * @return the top
	 */
	public int getTop() {
		return top;
	}

	/**
	 * @return the bottom
	 */
	public int getBottom() {
		return bottom;
	}

	/**
	 * @return the left
	 */
	public int getLeft() {
		return left;
	}

	/**
	 * @return the right
	 */
	public int getRight() {
		return right;
	}

	/**
	 * @return the xAxisTimeValues
	 */
	public ArrayList<Double> getxAxisTimeValues() {
		return xAxisTimeValues;
	}

	/**
	 * @return the xAxisTimeValues
	 */
	public ArrayList<Integer> getxAxisTimeCoordinates() {
		return xAxisTimeCoordinates;
	}

	class PlotInfo {

		double time;
		double atRisk;
		double censored;
		double events;
		double percentage;

		@Override
		public String toString() {
			return time + "\t" + atRisk + "\t" + censored + "\t" + events + "\t" + (atRisk - events) + "\t" + percentage;
		}
	}

	/**
	 * Do higher quality rendering options
	 *
	 * @param g
	 */
	private void setRenderingHints(Graphics2D g) {
		RenderingHints rh = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		rh.put(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_ENABLE);
		rh.put(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_GASP);
		rh.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);

		g.setRenderingHints(rh);

	}

	/**
	 * Setup the axis, labels etc
	 *
	 * @param g
	 */
	private void drawLegend(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;
		setRenderingHints(g2);
		g2.setColor(Color.BLACK);
		Font font = g2.getFont();
		Font f = new Font(font.getFontName(), Font.BOLD, font.getSize());
		g2.setFont(f);
		fm = getFontMetrics(f);
		int fontHeight = fm.getHeight();
		for (int i = 0; i < title.size(); i++) {
			if (fm.stringWidth(title.get(i)) > .8 * this.getWidth()) {
				f = new Font(font.getFontName(), Font.BOLD, 10);
				g2.setFont(f);
				fm = getFontMetrics(f);
			}
			g2.drawString(title.get(i), (getSize().width - fm.stringWidth(title.get(i))) / 2, ((i + 1) * fontHeight));
			// g2.setFont(font);
		}
		// draw the maxPercentage and minPercentage values
		String label = df.format(minPercentage);
		g2.drawString(label, left - 5 - (fm.stringWidth(label)), bottom + titleHeight / 6);
		g2.drawLine(left - 5, bottom, left, bottom);
		double d = minPercentage + kmfi.yaxisPercentIncrement;
		//double graphHeight = top - bottom;

		while (d < maxPercentage) {
			int yvalue = bottom - (int) (d * (bottom - top));
			label = df.format(d * 100);
			g2.drawString(label, left - 5 - (fm.stringWidth(label)), yvalue + titleHeight / 6); //

			g2.drawLine(left - 5, yvalue, left, yvalue);
			d = d + kmfi.yaxisPercentIncrement;
		}

		label = df.format(maxPercentage * 100);
		g2.drawString(label, left - 5 - (fm.stringWidth(label)), top + (titleHeight) / 6);
		g2.drawLine(left - 5, top, left, top);

		// Create a rotation transformation for the font.
		AffineTransform fontAT = new AffineTransform();


		// Derive a new font using a rotatation transform
		fontAT.rotate(270 * java.lang.Math.PI / 180);
		Font theDerivedFont = f.deriveFont(fontAT);

		// set the derived font in the Graphics2D context
		g2.setFont(theDerivedFont);

		// Render a string using the derived font
		int yaxisHeight = fm.stringWidth(kmfi.yAxisLegend);
		g2.drawString(kmfi.yAxisLegend, yaxisLabel, (bottom - (int) (.5 * (bottom - top))) + yaxisHeight / 2);

		// put the original font back
		g2.setFont(f);



		double timeDistance = maxTime - minTime;
		double timeIncrement = timeDistance * kmfi.xaxisPercentIncrement;
		double timeInt = (int) Math.floor(timeIncrement);
		if (timeInt < 1.0) {
			timeInt = 1.0;
		}
		adjustedPercentIncrement = timeInt / timeDistance;

		d = adjustedPercentIncrement; //kmfi.xaxisPercentIncrement;
		xAxisTimeValues.clear();
		xAxisTimeCoordinates.clear();

		//if we don't have time values then use percentage to set time. Not perfect but allows different tics
		if (kmfi.xAxisLabels.isEmpty()) {
			xAxisTimeValues.add(minTime);
			xAxisTimeCoordinates.add(left);
			while (d <= 1.0) {
				double xaxisTime = ((minTime * kmfi.timeScale) + d * ((maxTime - minTime) * kmfi.timeScale)); //
				xAxisTimeValues.add(xaxisTime);

				Integer coordinate = left + (int) (d * (right - left));
				xAxisTimeCoordinates.add(coordinate);
				//       System.out.println(d + " " + left + " " + right + " " + coordinate + " " + minTime + " " + maxTime);
				d = d + adjustedPercentIncrement; //kmfi.xaxisPercentIncrement;
			}
		} else {
			minTime = kmfi.xAxisLabels.get(0);
			maxTime = kmfi.xAxisLabels.get(kmfi.xAxisLabels.size() - 1);
			for (Double xaxisTime : kmfi.xAxisLabels) {
				xAxisTimeValues.add(xaxisTime);
				d = (xaxisTime - minTime) / (maxTime - minTime);
				Integer coordinate = left + (int) (d * (right - left));
				xAxisTimeCoordinates.add(coordinate);
			}
		}

		for (int i = 0; i < xAxisTimeValues.size(); i++) {
			Double xaxisTime = xAxisTimeValues.get(i);
			Integer xCoordinate = xAxisTimeCoordinates.get(i);
			label = df.format(xaxisTime);
			if (i == xAxisTimeValues.size() - 1) {
				g2.drawString(label, xCoordinate - (fm.stringWidth(label)), bottom + fm.getHeight() + 5);
			} else {
				g2.drawString(label, xCoordinate - (fm.stringWidth(label) / 2), bottom + fm.getHeight() + 5);
			}
			g2.drawLine(xCoordinate, bottom, xCoordinate, bottom + 5);
		}

		// draw the vertical and horizontal lines
		g2.setStroke(kmfi.axisStroke);
		g2.drawLine(left, top, left, bottom);
		g2.drawLine(left, bottom, right, bottom);

		// draw xAxis legend
		g2.drawString(kmfi.xAxisLegend, getSize().width / 2 - (fm.stringWidth(kmfi.xAxisLegend) / 2), bottom + 2 * fm.getHeight() + 10);
	}
	Double adjustedPercentIncrement = 0.0;

	/**
	 * Get the percentage increment for the time axis
	 *
	 * @return
	 */
	public Double getTimeAxisIncrementPercentage() {
		return adjustedPercentIncrement;
	}

	/**
	 * Reset the various bounds used to draw graph
	 */
	private void setFigureDimensions() {
		fm = getFontMetrics(getFont());
		titleHeight = kmfi.titleHeight;//fm.getHeight();
		xAxisLabelHeight = titleHeight;
		labelWidth = Math.max(fm.stringWidth(df.format(minPercentage)),
				fm.stringWidth(df.format(maxPercentage))) + 5;
		top = kmfi.padding + titleHeight;
		bottom = this.getHeight() - kmfi.padding - xAxisLabelHeight;
		left = kmfi.padding + labelWidth + yaxisLabel;
		right = this.getWidth() - kmfi.padding;

	}

	/**
	 * Combine the KM and Num risk into one image
	 *
	 * @param fileName
	 */
	public void savePNGKMNumRisk(String fileName) {
		if (fileName.startsWith("null") || fileName.startsWith("Null") || fileName.startsWith("NULL")) {
			return;
		}
		this.fileName = fileName;

		NumbersAtRiskPanel numbersAtRiskPanel = new NumbersAtRiskPanel();
		numbersAtRiskPanel.setKaplanMeierFigure(this);
		numbersAtRiskPanel.setSize(this.getWidth(), numbersAtRiskPanel.getHeight());
		BufferedImage imageKM = new BufferedImage(this.getWidth(), this.getHeight(), BufferedImage.TYPE_INT_RGB);
		Graphics2D graphics2D = imageKM.createGraphics();

		this.paint(graphics2D);

		BufferedImage imageNumRisk = new BufferedImage(numbersAtRiskPanel.getWidth(), numbersAtRiskPanel.getHeight(), BufferedImage.TYPE_INT_RGB);
		Graphics2D graphics2DNumRisk = imageNumRisk.createGraphics();
		numbersAtRiskPanel.paint(graphics2DNumRisk);


		BufferedImage image = new BufferedImage(numbersAtRiskPanel.getWidth(), numbersAtRiskPanel.getHeight() + this.getHeight(), BufferedImage.TYPE_INT_RGB);
		Graphics2D g = image.createGraphics();

		g.drawImage(imageKM, 0, 0, null);
		g.drawImage(imageNumRisk, 0, this.getHeight(), null);

		try {
			ImageIO.write(image, "png", new File(fileName));
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	/**
	 *
	 * @param fileName
	 */
	public void savePNG(String fileName) {
		if (fileName.startsWith("null") || fileName.startsWith("Null") || fileName.startsWith("NULL")) {
			return;
		}
		this.fileName = fileName;
		BufferedImage image = new BufferedImage(this.getWidth(), this.getHeight(), BufferedImage.TYPE_INT_RGB);
		Graphics2D graphics2D = image.createGraphics();

		this.paint(graphics2D);
		try {
			ImageIO.write(image, "png", new File(fileName));
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
		try {

			KaplanMeierFigure kaplanMeierFigure = new KaplanMeierFigure();
			LinkedHashMap<String, ArrayList<CensorStatus>> survivalDataHashMap = new LinkedHashMap<String, ArrayList<CensorStatus>>();

//            if (false) { //http://sph.bu.edu/otlt/MPH-Modules/BS/BS704_Survival/
//                ArrayList<CensorStatus> graph1 = new ArrayList<CensorStatus>();
//                graph1.add(new CensorStatus("A", 24.0, "0"));
//                graph1.add(new CensorStatus("A", 3.0, "1"));
//                graph1.add(new CensorStatus("A", 11.0, "0"));
//                graph1.add(new CensorStatus("A", 19.0, "0"));
//                graph1.add(new CensorStatus("A", 24.0, "0"));
//                graph1.add(new CensorStatus("A", 13.0, "0"));
//
//                graph1.add(new CensorStatus("A", 14.0, "1"));
//                graph1.add(new CensorStatus("A", 2.0, "0"));
//                graph1.add(new CensorStatus("A", 18.0, "0"));
//                graph1.add(new CensorStatus("A", 17.0, "0"));
//                graph1.add(new CensorStatus("A", 24.0, "0"));
//                graph1.add(new CensorStatus("A", 21.0, "0"));
//                graph1.add(new CensorStatus("A", 12.0, "0"));
//
//                graph1.add(new CensorStatus("A", 1.0, "1"));
//                graph1.add(new CensorStatus("A", 10.0, "0"));
//                graph1.add(new CensorStatus("A", 23.0, "1"));
//                graph1.add(new CensorStatus("A", 6.0, "0"));
//                graph1.add(new CensorStatus("A", 5.0, "1"));
//                graph1.add(new CensorStatus("A", 9.0, "0"));
//                graph1.add(new CensorStatus("A", 17.0, "1"));
//
//                survivalDataHashMap.put("Label 1", graph1);
//
//
//
//            }


			if (true) {



				ArrayList<CensorStatus> graph1 = new ArrayList<CensorStatus>();
				graph1.add(new CensorStatus("A", 1.0, "1"));
				graph1.add(new CensorStatus("A", 1.0, "1"));
				graph1.add(new CensorStatus("A", 1.0, "1"));
				graph1.add(new CensorStatus("A", 2.0, "1"));
				graph1.add(new CensorStatus("A", 2.0, "1"));
				graph1.add(new CensorStatus("A", 3.0, "1"));

				graph1.add(new CensorStatus("A", 4.0, "1"));
				graph1.add(new CensorStatus("A", 4.0, "1"));
				graph1.add(new CensorStatus("A", 4.0, "1"));
				graph1.add(new CensorStatus("A", 4.0, "1"));
				graph1.add(new CensorStatus("A", 4.0, "1"));
				graph1.add(new CensorStatus("A", 4.0, "1"));
				graph1.add(new CensorStatus("A", 4.0, "0"));

				graph1.add(new CensorStatus("A", 5.0, "1"));
				graph1.add(new CensorStatus("A", 5.0, "1"));

				graph1.add(new CensorStatus("A", 8.0, "0"));
				graph1.add(new CensorStatus("A", 8.0, "0"));
				graph1.add(new CensorStatus("A", 8.0, "0"));
				graph1.add(new CensorStatus("A", 8.0, "0"));
				graph1.add(new CensorStatus("A", 8.0, "0"));
				graph1.add(new CensorStatus("A", 8.0, "0"));
				graph1.add(new CensorStatus("A", 8.0, "1"));

				graph1.add(new CensorStatus("A", 9.0, "1"));
				graph1.add(new CensorStatus("A", 9.0, "1"));
				graph1.add(new CensorStatus("A", 9.0, "1"));
				graph1.add(new CensorStatus("A", 9.0, "1"));
				graph1.add(new CensorStatus("A", 9.0, "1"));


				graph1.add(new CensorStatus("A", 13.0, "0"));
				graph1.add(new CensorStatus("A", 13.0, "0"));
				graph1.add(new CensorStatus("A", 13.0, "1"));

				survivalDataHashMap.put("Label 1", graph1);

				ArrayList<CensorStatus> graph2 = new ArrayList<CensorStatus>();
				graph2.add(new CensorStatus("A", 1.0, "1"));
				graph2.add(new CensorStatus("A", 1.0, "1"));
				graph2.add(new CensorStatus("A", 1.0, "0"));
				graph2.add(new CensorStatus("A", 3.0, "0"));
				graph2.add(new CensorStatus("A", 3.0, "1"));
				graph2.add(new CensorStatus("A", 4.0, "1"));

				graph2.add(new CensorStatus("A", 4.0, "1"));
				graph2.add(new CensorStatus("A", 4.0, "1"));
				graph2.add(new CensorStatus("A", 4.0, "1"));
				graph2.add(new CensorStatus("A", 5.0, "1"));
				graph2.add(new CensorStatus("A", 5.0, "0"));
				graph2.add(new CensorStatus("A", 5.0, "0"));
				graph2.add(new CensorStatus("A", 5.0, "0"));

				graph2.add(new CensorStatus("A", 6.0, "1"));
				graph2.add(new CensorStatus("A", 6.0, "0"));

				graph2.add(new CensorStatus("A", 7.0, "0"));
				graph2.add(new CensorStatus("A", 7.0, "0"));
				graph2.add(new CensorStatus("A", 7.0, "0"));
				graph2.add(new CensorStatus("A", 7.0, "0"));
				graph2.add(new CensorStatus("A", 8.0, "1"));
				graph2.add(new CensorStatus("A", 8.0, "1"));
				graph2.add(new CensorStatus("A", 8.0, "1"));

				graph2.add(new CensorStatus("A", 8.0, "1"));
				graph2.add(new CensorStatus("A", 8.0, "1"));
				graph2.add(new CensorStatus("A", 8.0, "0"));
				graph2.add(new CensorStatus("A", 9.0, "0"));
				graph2.add(new CensorStatus("A", 9.0, "1"));


				graph2.add(new CensorStatus("A", 10.0, "0"));
				graph2.add(new CensorStatus("A", 10.0, "0"));
				graph2.add(new CensorStatus("A", 10.0, "0"));

				survivalDataHashMap.put("Label 2", graph2);
			}

			ArrayList<String> figureInfo = new ArrayList<String>();
			//DecimalFormat dfe = new DecimalFormat("0.00E0");
			//DecimalFormat df = new DecimalFormat("0.00");



			JFrame application = new JFrame();
			application.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			application.add(kaplanMeierFigure);
			kaplanMeierFigure.setSize(500, 400);

			application.setSize(500, 400);         // window is 500 pixels wide, 400 high
			application.setVisible(true);

			ArrayList<String> titles = new ArrayList<String>();
			titles.add("Line 1");
			titles.add("line 2");
			kaplanMeierFigure.setSurvivalData(titles, survivalDataHashMap, true);

			//   figureInfo.add("HR=2.1 95% CI(1.8-2.5)");
			//   figureInfo.add("p-value=.001");
			kaplanMeierFigure.setFigureLineInfo(figureInfo);

			kaplanMeierFigure.savePNGKMNumRisk("/Users/Scooter/Downloads/test.png");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
