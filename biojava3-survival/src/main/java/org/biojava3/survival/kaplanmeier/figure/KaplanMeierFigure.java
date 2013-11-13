/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.kaplanmeier.figure;


import org.biojava3.survival.cox.CoxCoefficient;
import org.biojava3.survival.cox.CoxInfo;
import org.biojava3.survival.cox.SurvFitInfo;

import org.biojava3.survival.cox.SurvivalInfo;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;
import org.biojava3.survival.cox.StrataInfo;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class KaplanMeierFigure extends JPanel {

    ArrayList<String> title = new ArrayList<String>();
    /**
     *
     */
    public int top;
    /**
     *
     */
    public int bottom;
    /**
     *
     */
    public int left;
    /**
     *
     */
    public int right;
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

    /**
     *
     */
    public KaplanMeierFigure() {
        super();
        setSize(500, 400);
        setBackground(Color.WHITE);
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
            if(legend == null){
                
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
        DecimalFormat df = new DecimalFormat("#.##");
        String line1 = "HR=" + fmt(cc.getHazardRatio(), 2, 0) + " (CI:" + fmt(cc.getHazardRatioLoCI(), 2, 0) + "-" + fmt(cc.getHazardRatioHiCI(), 2, 0) + ")";
        String line2 = "p=" + fmt(cc.getPvalue(), 3, 0);
        String line3 = "n=" + n + " e=" + event;
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
    public SurvFitInfo getSurvivalFitInfo(){
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

        int evenCheck = Math.round((float) mTime.floatValue());
        if (evenCheck % 2 == 1) {
            evenCheck = evenCheck + 1;
        }
        this.maxTime = (double) evenCheck;

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

        int evenCheck = Math.round((float) mTime.floatValue());
        if (evenCheck % 2 == 1) {
            evenCheck = evenCheck + 1;
        }
        this.maxTime = (double) evenCheck;

        if (userSetMaxTime != null && userSetMaxTime > maxTime) {
            this.maxTime = userSetMaxTime;
        }

        //calculate percentages
        SurvFitKM survFitKM = new SurvFitKM();
        sfi = survFitKM.process(survivalData, useWeighted);





        this.repaint();
    }
    DecimalFormat df = new DecimalFormat("#.#");

    public void paintComponent(Graphics g) // draw graphics in the panel
    {
        int width = getWidth();             // width of window in pixels
        int height = getHeight();           // height of window in pixels
        setFigureDimensions();
        super.paintComponent(g);            // call superclass to make panel display correctly

        drawLegend(g);
        drawSurvivalCurves(g);
        drawFigureLineInfo(g);
        // Drawing code goes here
    }

    private void drawFigureLineInfo(Graphics g) {
        g.setColor(Color.BLACK);
        fm = getFontMetrics(getFont());
        int yoffset = fm.getHeight() * lineInfoList.size();

        int x = getTimeX(kmfi.figureLineInfoLowerPercentX * maxTime);
        int y = getPercentageY(kmfi.figureLineInfoLowerPercentY) - yoffset;

        for (String line : lineInfoList) {
            g.drawString(line, x, y);
            y = y + fm.getHeight();
        }

    }

    private void drawSurvivalCurves(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
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

                g2.drawLine(getTimeX(p0time), getPercentageY(p0percentage), getTimeX(p1time), getPercentageY(p0percentage));

                g2.drawLine(getTimeX(p1time), getPercentageY(p0percentage), getTimeX(p1time), getPercentageY(p1percentage));
                if (si.getStatus().get(i) == 0) {
                    g2.drawLine(getTimeX(p0time), getPercentageY(p0percentage) - 4, getTimeX(p0time), getPercentageY(p0percentage) + 4);
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

    private void drawSurvivalCurvesOld(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
        g2.setStroke(kmfi.kmStroke);

        int colorIndex = 0;
        ArrayList<String> labels = new ArrayList<String>(survivalData.keySet());
        Collections.sort(labels);

        for (String legend : labels) {
            g2.setColor(kmfi.legendColor[colorIndex]);
            colorIndex++;
            ArrayList<CensorStatus> censorStatusList = survivalData.get(legend);
            Collections.sort(censorStatusList);
            ArrayList<PlotInfo> plotInfoList = new ArrayList<PlotInfo>();
            PlotInfo pi = new PlotInfo();
            pi.time = 0;
            pi.percentage = 1.0;
            pi.atRisk = censorStatusList.size();
            pi.censored = 0;
            pi.events = 0;
            plotInfoList.add(pi);
            pi = null;
            double lastTime = -1.0;

            double lastPercentage = 1.0;
            double atRisk = 0;
            boolean useExternalPercentage = false;
            for (CensorStatus cs : censorStatusList) {
                if (lastTime < 0 || lastTime != cs.time) {
                    if (pi != null) {
                        if (pi.atRisk == 0) {
                            pi.percentage = lastPercentage;
                        } else {
                            Double percentage = null;
                            if (cs.getPercentage() == null) {
                                percentage = ((pi.atRisk - pi.events) / pi.atRisk) * lastPercentage;
                            } else {
                                percentage = cs.getPercentage();
                                pi.censored = cs.ncens;
                                pi.events = cs.nevents;
                                pi.atRisk = cs.nrisk;
                                useExternalPercentage = true;
                            }
                            pi.percentage = percentage;
                            lastPercentage = percentage;
                        }
                        atRisk = pi.atRisk - pi.events - pi.censored;
                    }
                    pi = new PlotInfo();
                    pi.time = cs.time;
                    plotInfoList.add(pi);
                    if (lastTime < 0) {
                        pi.atRisk = censorStatusList.size();
                    } else {
                        pi.atRisk = atRisk;
                    }
                    lastTime = cs.time;
                }
                if (cs.censored.equals("0")) {
                    pi.censored++;
                    //pi.atRisk--;
                } else {
                    pi.events++;
                }
                if (cs.getPercentage() != null) {
                    useExternalPercentage = true;
                    pi.percentage = cs.getPercentage();
                    pi.censored = cs.ncens;
                    pi.events = cs.nevents;
                    pi.atRisk = cs.nrisk;

                }
            }

            if (useExternalPercentage == false) { //don't need to correct for last one if set externally
                if (pi.atRisk == 0) {
                    pi.percentage = lastPercentage;
                } else {
                    pi.percentage = ((pi.atRisk - pi.events) / pi.atRisk) * lastPercentage;
                }
            }

            if (true) {
                for (PlotInfo p : plotInfoList) {
                    //                   System.out.println(p);
                }
            }
            for (int i = 0; i < plotInfoList.size() - 1; i++) {
                PlotInfo p0 = plotInfoList.get(i);
                PlotInfo p1 = plotInfoList.get(i + 1);
                g2.drawLine(getTimeX(p0.time), getPercentageY(p0.percentage), getTimeX(p1.time), getPercentageY(p0.percentage));

                g2.drawLine(getTimeX(p1.time), getPercentageY(p0.percentage), getTimeX(p1.time), getPercentageY(p1.percentage));
            }
//            System.out.println();
            for (CensorStatus cs : censorStatusList) {
                if (cs.censored.equals("0")) {
                    PlotInfo p0 = null;
                    for (PlotInfo p : plotInfoList) {
                        if (p.time == cs.time) {
                            p0 = p;
                            break;
                        }
                    }
                    g2.drawLine(getTimeX(cs.time), getPercentageY(p0.percentage) - 4, getTimeX(cs.time), getPercentageY(p0.percentage) + 4);
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

    private int getTimeX(double value) {
        double d = left + (((right - left) * value) / (maxTime - minTime));
        return (int) d;
    }

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

    private void drawLegend(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
        Font font = g.getFont();
        fm = getFontMetrics(font);
        int fontHeight = fm.getHeight();
        for (int i = 0; i < title.size(); i++) {
            if (fm.stringWidth(title.get(i)) > .8 * this.getWidth()) {

                Font f = new Font(font.getFontName(), Font.PLAIN, 10);
                g.setFont(f);
                fm = getFontMetrics(f);
            }
            g.drawString(title.get(i), (size().width - fm.stringWidth(title.get(i))) / 2, ((i + 1) * fontHeight));
            g.setFont(font);
        }
        // draw the maxPercentage and minPercentage values
        g.drawString(df.format(minPercentage), kmfi.padding, bottom + titleHeight / 6);
        g.drawLine(left - 5, bottom, left, bottom);
        double d = minPercentage + kmfi.yaxisPercentIncrement;
        double graphHeight = top - bottom;
        String label = "";
        while (d < maxPercentage) {
            int yvalue = bottom - (int) (d * (bottom - top));
            label = df.format(d);
            g.drawString(label, kmfi.padding - (int) (.6 * fm.stringWidth(label)), yvalue + titleHeight / 6); //

            g.drawLine(left - 5, yvalue, left, yvalue);
            d = d + kmfi.yaxisPercentIncrement;
        }

        label = df.format(maxPercentage);
        g.drawString(label, kmfi.padding - (int) (.6 * fm.stringWidth(label)), top + (titleHeight) / 6);
        g.drawLine(left - 5, top, left, top);

        double timeDistance = maxTime - minTime;
        double timeIncrement = timeDistance * kmfi.xaxisPercentIncrement;
        double timeInt = (int) timeIncrement;
        if (timeInt < 1.0) {
            timeInt = 1.0;
        }
        double adjustedPercentIncrement = timeInt / timeDistance;

        d = adjustedPercentIncrement; //kmfi.xaxisPercentIncrement;
        while (d <= 1.0) {
            label = df.format((minTime * kmfi.timeScale) + d * ((maxTime - minTime) * kmfi.timeScale)); //
            if (d + adjustedPercentIncrement > 1.0) { //if this is the last one then adjust
                g.drawString(label, left + (int) (d * (right - left)) - (int) (.5 * fm.stringWidth(label)), bottom + fm.getHeight() + 5);
            } else {
                g.drawString(label, left + (int) (d * (right - left)) - (fm.stringWidth(label) / 2), bottom + fm.getHeight() + 5);
            }
            g.drawLine(left + (int) (d * (right - left)), bottom, left + (int) (d * (right - left)), bottom + 5);
            d = d + adjustedPercentIncrement; //kmfi.xaxisPercentIncrement;
        }


        // draw the vertical and horizontal lines
        g2.setStroke(kmfi.axisStroke);
        g.drawLine(left, top, left, bottom);
        g.drawLine(left, bottom, right, bottom);
    }

    private void setFigureDimensions() {
        fm = getFontMetrics(getFont());
        titleHeight = kmfi.titleHeight;//fm.getHeight();
        xAxisLabelHeight = titleHeight;
        labelWidth = Math.max(fm.stringWidth(df.format(minPercentage)),
                fm.stringWidth(df.format(maxPercentage))) + 5;
        top = kmfi.padding + titleHeight;
        bottom = this.getHeight() - kmfi.padding - xAxisLabelHeight;
        left = kmfi.padding + labelWidth;
        right = this.getWidth() - kmfi.padding;

    }

    /**
     *
     * @param fileName
     */
    public void savePNG(String fileName) {
        if(fileName.startsWith("null") || fileName.startsWith("Null") || fileName.startsWith("NULL"))
            return;
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

            if (false) { //http://sph.bu.edu/otlt/MPH-Modules/BS/BS704_Survival/
                ArrayList<CensorStatus> graph1 = new ArrayList<CensorStatus>();
                graph1.add(new CensorStatus("A", 24.0, "0"));
                graph1.add(new CensorStatus("A", 3.0, "1"));
                graph1.add(new CensorStatus("A", 11.0, "0"));
                graph1.add(new CensorStatus("A", 19.0, "0"));
                graph1.add(new CensorStatus("A", 24.0, "0"));
                graph1.add(new CensorStatus("A", 13.0, "0"));

                graph1.add(new CensorStatus("A", 14.0, "1"));
                graph1.add(new CensorStatus("A", 2.0, "0"));
                graph1.add(new CensorStatus("A", 18.0, "0"));
                graph1.add(new CensorStatus("A", 17.0, "0"));
                graph1.add(new CensorStatus("A", 24.0, "0"));
                graph1.add(new CensorStatus("A", 21.0, "0"));
                graph1.add(new CensorStatus("A", 12.0, "0"));

                graph1.add(new CensorStatus("A", 1.0, "1"));
                graph1.add(new CensorStatus("A", 10.0, "0"));
                graph1.add(new CensorStatus("A", 23.0, "1"));
                graph1.add(new CensorStatus("A", 6.0, "0"));
                graph1.add(new CensorStatus("A", 5.0, "1"));
                graph1.add(new CensorStatus("A", 9.0, "0"));
                graph1.add(new CensorStatus("A", 17.0, "1"));

                survivalDataHashMap.put("Label 1", graph1);



            }


            if (false) {



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
            DecimalFormat dfe = new DecimalFormat("0.00E0");
            DecimalFormat df = new DecimalFormat("0.00");



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

            kaplanMeierFigure.savePNG("/Users/Scooter/Downloads/test.png");

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
