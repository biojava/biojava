/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.kaplanmeier.figure;


import java.awt.BasicStroke;
import java.awt.Color;
import java.io.FileInputStream;
import java.util.Properties;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class KMFigureInfo {

    /**
     *
     */
    public int titleHeight = 40;
    /**
     *
     */
    public int padding = 20;
    /**
     *
     */
    public Integer width = 600;
    /**
     *
     */
    public Integer height = 400;
    /**
     *
     */
    public double timeScale = 1.0; //multiplier to change time from days to months etc
    /**
     *
     */
    public double yaxisPercentIncrement = .2;
    /**
     *
     */
    public double xaxisPercentIncrement = .1;
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
    public double figureLineInfoLowerPercentX = .01;
    /**
     *
     */
    public double figureLineInfoLowerPercentY = .01;
    /**
     *
     */
    public BasicStroke axisStroke = new BasicStroke(2);
    /**
     *
     */
    public BasicStroke kmStroke = new BasicStroke(3);
    /**
     *
     */
    public Color[] legendColor = {Color.RED, Color.BLUE, Color.GREEN, Color.CYAN, Color.ORANGE, Color.YELLOW, Color.MAGENTA, Color.PINK};

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

    }
}
